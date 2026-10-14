"""Validated import and descriptive analysis of 101-node bilateral CMJ exports."""
from dataclasses import dataclass
import csv
import hashlib
import io
import re
import numpy as np
import pandas as pd

PERCENT = np.arange(101)

@dataclass
class Export:
    name: str
    digest: str
    values: np.ndarray
    mapping: pd.DataFrame


def read_export(raw: bytes, name: str) -> Export:
    text = raw.decode('utf-8-sig')
    rows = list(csv.reader(io.StringIO(text), delimiter='\t'))
    rows = [r for r in rows if any(v.strip() for v in r)]
    if len(rows) != 103:
        raise ValueError('Expected two header rows and exactly 101 data rows; no silent resampling is performed.')
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError('Header/data column counts differ.')
    try:
        values = np.array(rows[2:], dtype=float)
    except ValueError as exc:
        raise ValueError('Non-numeric or missing sample found; correct the source export.') from exc
    if not np.isfinite(values).all():
        raise ValueError('NaN or infinite values are not supported.')
    if not (np.allclose(values[:, 0], np.arange(1,102)) or np.allclose(values[:,0], PERCENT)):
        raise ValueError('First column must be the ordered 1–101 or 0–100 node index.')
    counts = {}
    records = []
    for col, (label, axis) in enumerate(zip(rows[0][1:], rows[1][1:]), 1):
        label, axis = label.strip(), axis.strip().upper()
        match = re.match(r'^(Left |Right |L_|R_)(.*)$', label, re.I)
        limb = 'Left' if match and match[1].lower() in ('left ', 'l_') else 'Right' if match else 'Unknown'
        signal = match[2].strip() if match else label
        key = (limb, signal, axis)
        counts[key] = counts.get(key, 0) + 1
        occurrence = counts[key]
        records.append(dict(include=occurrence <= 3, column=col, source=label,
                            limb=limb, signal=signal, axis=axis, trial=occurrence))
    return Export(name, hashlib.sha256(raw).hexdigest(), values, pd.DataFrame(records))


def curves(export: Export, mapping: pd.DataFrame, trials: list[int]):
    selected = mapping[mapping['include'] & mapping.trial.isin(trials)].copy()
    if len(trials) < 2 or len(set(trials)) != len(trials):
        raise ValueError('Select at least two distinct trial IDs.')
    if selected.empty:
        raise ValueError('No selected channels.')
    if selected[['limb','signal','axis','trial','column']].isna().any().any():
        raise ValueError('Selected mapping fields cannot be blank.')
    if any(not isinstance(v,str) or not v.strip() for col in ['limb','signal','axis'] for v in selected[col]):
        raise ValueError('Selected limb, signal and axis labels must be nonempty text.')
    if any(float(t) != int(t) or int(t) < 1 for t in selected.trial):
        raise ValueError('Trial IDs must be positive integers.')
    if selected.duplicated(['limb','signal','axis','trial']).any():
        raise ValueError('Each limb/signal/axis must have one column per trial.')
    result = {}
    for key, group in selected.groupby(['limb','signal','axis'], sort=False):
        if set(group.trial) != set(trials):
            raise ValueError(f'{key} is missing a selected trial. Exclude this channel or fix the mapping.')
        cols = group.set_index('trial').loc[trials, 'column'].astype(int).to_numpy()
        result[key] = export.values[:, cols].T
    return result


def phase_candidates(data, trial_index, threshold_fraction=.05):
    """Conservative force-contact candidates. Never labels an initial flight as takeoff."""
    left = data.get(('Left','GRF','Z'))
    right = data.get(('Right','GRF','Z'))
    if left is None or right is None:
        return {}, 'Bilateral vertical GRF is missing; manual events required.'
    force = left[trial_index] + right[trial_index]
    peak = force.max()
    if peak <= 0:
        return {}, 'No positive vertical GRF; check axis/sign.'
    contact = force > threshold_fraction * peak
    if not contact[:3].all():
        return {}, 'Recording does not begin in sustained contact. Full CMJ phases cannot be inferred.'
    runs = [(m.start(), m.end()) for m in re.finditer('0{3,}', ''.join(contact.astype(int).astype(str)))]
    for start, end in runs:
        if start >= 5 and end <= 98 and contact[end:end+3].all():
            events = {'takeoff': start, 'landing': end}
            knees = [data.get((s,'Knee Angles','X')) for s in ['Left','Right']]
            if all(k is not None for k in knees):
                k = np.mean([v[trial_index] for v in knees], axis=0)
                bottom = int(np.argmax(k[:start]))
                if 0 < bottom < start:
                    events['bottom'] = bottom
            return events, 'Candidates only: force threshold is relative to peak; bottom assumes positive knee flexion. Review every trial.'
    return {}, 'No sustained flight bracketed by contact; manual events required.'


def validate_events(events):
    required = ['onset','bottom','takeoff','landing','end']
    if not all(k in events for k in required):
        raise ValueError('Supply onset, bottom, takeoff, landing and end.')
    vals = [float(events[k]) for k in required]
    if not all(np.isfinite(vals)) or not 0 <= vals[0] < vals[1] < vals[2] < vals[3] < vals[4] <= 100:
        raise ValueError('Events must satisfy 0 ≤ onset < bottom < takeoff < landing < end ≤ 100.')
    return dict(zip(required, vals))


def windows(events):
    out = {'whole_cycle': (0,100)}
    if events:
        e = validate_events(events)
        out.update(descent=(e['onset'],e['bottom']), ascent=(e['bottom'],e['takeoff']),
                   flight=(e['takeoff'],e['landing']), landing_absorption=(e['landing'],e['end']))
    return out


def metrics(subject, data, trial_ids, events_by_trial, units):
    rows = []
    for (limb, signal, axis), y in data.items():
        for i, trial in enumerate(trial_ids):
            for phase, (lo,hi) in windows(events_by_trial.get(trial)).items():
                x = np.unique(np.r_[lo, PERCENT[(PERCENT>lo)&(PERCENT<hi)], hi])
                v = np.interp(x,PERCENT,y[i])
                row = dict(subject=subject, trial=trial, limb=limb, signal=signal, axis=axis,
                           phase=phase, unit=units.get(signal,'source units'))
                row.update(mean=float(np.trapezoid(v,x)/(hi-lo)), minimum=float(v.min()), maximum=float(v.max()),
                           rom=float(np.ptp(v)), rms=float(np.sqrt(np.trapezoid(v*v,x)/(hi-lo))),
                           peak_abs=float(np.max(np.abs(v))), peak_abs_percent=float(x[np.argmax(np.abs(v))]),
                           area_unit_percent=float(np.trapezoid(v,x)))
                rows.append(row)
    trial_table = pd.DataFrame(rows)
    keys = ['subject','limb','signal','axis','phase','unit']
    values = ['mean','minimum','maximum','rom','rms','peak_abs','peak_abs_percent','area_unit_percent']
    summary = trial_table.groupby(keys)[values].agg(['mean','std','count'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    return trial_table, summary.reset_index()


def average_table(subject, data):
    frames=[]
    for (limb,signal,axis), y in data.items():
        frames.append(pd.DataFrame(dict(subject=subject,limb=limb,signal=signal,axis=axis,
                         percent=PERCENT,mean=y.mean(0),sd=y.std(0,ddof=1),n_trials=len(y))))
    return pd.concat(frames,ignore_index=True)


def differences(subject, data):
    rows=[]
    for limb,signal,axis in data:
        if limb != 'Left' or ('Right',signal,axis) not in data:
            continue
        l=data[(limb,signal,axis)].mean(0); r=data[('Right',signal,axis)].mean(0)
        denom=np.abs(r)+np.abs(l)
        si=np.divide(200*(r-l),denom,out=np.full(101,np.nan),where=denom>1e-12)
        rows.append(pd.DataFrame(dict(subject=subject,signal=signal,axis=axis,percent=PERCENT,
                    right_minus_left=r-l,signed_symmetry_percent=si,rmse=np.sqrt(np.mean((r-l)**2)))))
    return pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()


KNEE_PLANES = {'X': 'Flexion', 'Y': 'Varus', 'Z': 'Rotation'}


def channel_label(signal, axis):
    return f'{signal} {axis}' + (f' · {KNEE_PLANES[axis]}' if signal == 'Knee Angles' and axis in KNEE_PLANES else '')


def confirmed_three_block_layout(export):
    """User-confirmed mapping, applied only to the matching original 45-column layout."""
    expected = [
        ('Left','Gluteus Medius','X'), ('Right','Gluteus Maximus','X'),
        ('Right','Gluteus Medius','X'), ('Left','GRF','Z'), ('Left','Hip Angles','X'),
        ('Left','Knee Angles','X'), ('Left','Knee Angles','Y'), ('Left','Knee Angles','Z'),
        ('Right','GRF','Z'), ('Right','Hip Angles','X'), ('Right','Knee Angles','X'),
        ('Right','Knee Angles','Y'), ('Right','Knee Angles','Z'),
        ('Right','Knee Power','X'), ('Left','Knee Power','X')]
    actual=list(export.mapping[['limb','signal','axis']].head(45).itertuples(index=False,name=None))
    return actual == expected * 3
