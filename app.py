"""Run: python -m streamlit run app.py"""
from pathlib import Path
import io
import hashlib
import json
import re
import zipfile
import numpy as np
import pandas as pd
import streamlit as st
from cmj.core import read_export, curves, phase_candidates, validate_events, metrics, average_table, differences, confirmed_three_block_layout, KNEE_PLANES, channel_label
from cmj.statistics import paired_spm, exact_max_t
from cmj.plots import comparison, spm_figure, phase_review_figure
import matplotlib.pyplot as plt

st.set_page_config(page_title='CMJ bilateral analysis',layout='wide')
st.title('Countermovement jump · bilateral analysis')
st.caption('Import normalized 101-point exports • review trial mapping • compare limbs • export results')
st.info('Samples are normalized movement percentages, not seconds. Configured for your confirmed first three trial blocks, synchronized limbs, GRF in N, and knee X/Y/Z = flexion/varus/rotation. Verify settings for other export formats.')


def safe(value):
    return re.sub(r'[^A-Za-z0-9_.-]+','_',value).strip('.') or 'subject'


def csv_bytes(table):
    return table.to_csv(index=False).encode('utf-8')


def figure_bytes(fig,fmt):
    buf=io.BytesIO(); fig.savefig(buf,format=fmt,dpi=300,bbox_inches='tight'); return buf.getvalue()


def add_stats(label,right,left,plane,alpha,family,exploratory,assets,figures,records):
    base=dict(subject=label,plane=plane,contrast='Right - Left')
    if len(right)<=12:
        try:
            result,p=exact_max_t(right,left,alpha,family)
            records.append(base|result)
            assets[f'{safe(label)}/exact_{plane}_node_p.csv']=csv_bytes(pd.DataFrame({'percent':np.arange(101),'p_max_t':p}))
            st.write(f'Knee {plane} exact permutation: global p = {result["p_global"]:.4g}; minimum attainable two-sided p = {result["minimum_two_sided_p"]:.4g}.')
        except ValueError as exc:
            records.append(base|dict(method='exact',status=str(exc)))
            st.warning(str(exc))
    if exploratory:
        try:
            result,inference=paired_spm(right,left,alpha,family)
            records.append(base|result)
            st.write(f'Knee {plane} paired SPM: threshold {result["threshold"]:.3g}; significant clusters: {len(result["clusters"])}.')
            fig=spm_figure(f'{label} · knee {plane} ({KNEE_PLANES[plane]}) · paired SPM',inference)
            figures[f'{safe(label)}/SPM_knee_{plane}']=fig
            st.pyplot(fig)
            assets[f'{safe(label)}/spm_{plane}_field.csv']=csv_bytes(pd.DataFrame({'percent':np.arange(101),'t':inference.z,'threshold':result['threshold']}))
        except (ValueError,AssertionError,RuntimeError,np.linalg.LinAlgError) as exc:
            records.append(base|dict(method='paired parametric SPM',status=str(exc)))
            st.warning(f'Knee {plane}: {exc}')

with st.sidebar:
    st.header('Analysis settings')
    alpha=st.number_input('Family alpha',.001,.10,.05,.001,format='%.3f')
    show_trials=st.checkbox('Show individual trials',True)
    parametric=st.checkbox('Enable exploratory parametric SPM',False)
    st.caption('Two/three trials cannot establish robust subject-level significance. Exact sign flips cannot reach p < .05 with n = 2 or 3. Parametric SPM needs ≥3 pairs and strong assumptions; small-n results are exploratory.')
    family_scope=st.selectbox('Multiplicity family',['Three planes per comparison','All subjects plus group × three planes'])
    angle_unit=st.text_input('Angle units','degrees (confirm)')
    force_unit=st.text_input('GRF units','N')
    power_unit=st.text_input('Power units','source units (unconfirmed)')
    emg_unit=st.selectbox('EMG output scale',['MVC-normalized (scale unconfirmed)','MVC ratio','%MVC','mV (unnormalized)'])
    st.caption('Raw EMG was measured in mV. After division by MVC it is a ratio, or %MVC if multiplied by 100. Selecting a label does not rescale the values.')
    threshold=st.slider('Candidate contact threshold (% of peak total GRF)',1,15,5)/100
    st.caption('Relative GRF threshold proposes events only. Review against the acquisition protocol; it is not a validated force threshold in newtons.')

uploads=st.file_uploader('Select subject text files',type=['txt','tsv'],accept_multiple_files=True)
if not uploads:
    st.write('Select one or more files. Each file should contain two tab-delimited header rows and 101 numeric rows. The first column is ITEM (1–101 or 0–100).')
    st.stop()

assets={}; figures={}; records=[]; trial_tables=[]; summary_tables=[]; average_tables=[]; diff_tables=[]
subjects=[]; seen={}; names=set(); metadata=[]
family=3 if family_scope.startswith('Three') else 3*(len(uploads)+1)
for file_index,upload in enumerate(uploads):
    with st.expander(upload.name,expanded=True):
        try: exported=read_export(upload.getvalue(),upload.name)
        except (ValueError,UnicodeError) as exc:
            st.error(f'Import failed: {exc}'); continue
        known_layout=confirmed_three_block_layout(exported)
        if known_layout:
            exported.mapping['include']=exported.mapping.column.between(1,45)
        key=f'{file_index}_{exported.digest[:10]}'
        subject=st.text_input('Subject ID',Path(upload.name).stem,key=key+'id').strip()
        if not subject or safe(subject) in names:
            st.error('Use a nonempty, unique subject ID (also unique after filename sanitization).'); continue
        names.add(safe(subject))
        duplicated=exported.digest in seen
        if duplicated:
            st.warning(f'Identical file content to {seen[exported.digest]}. Descriptive exports are available, but this copy is excluded from group inference.')
        else: seen[exported.digest]=subject
        surgical=st.selectbox('Operated limb',['Unknown','Left','Right','Bilateral'],key=key+'surgery')
        st.caption('Your confirmed layout uses columns 1–45: three complete trial blocks. Other layouts default to the first three occurrences per channel and require review. Extra columns remain excluded. Trial numbers must identify the same jump across all channels and both limbs.')
        with st.expander('Review / edit all column mappings'):
            mapping=st.data_editor(exported.mapping,hide_index=True,disabled=['column','source'],key=key+'mapping',
                column_config={'limb':st.column_config.SelectboxColumn(options=['Left','Right','Unknown']),
                               'trial':st.column_config.NumberColumn(min_value=1,step=1)})
        try:
            mapped_trials=pd.to_numeric(mapping.loc[mapping.include,'trial'],errors='raise')
            if mapped_trials.isna().any() or (mapped_trials<1).any() or (mapped_trials%1 != 0).any():
                raise ValueError('Trial IDs must be positive integers.')
            trial_options=sorted(set(mapped_trials.astype(int)))
        except (ValueError,TypeError) as exc:
            st.error(f'Invalid trial mapping: {exc}'); continue
        trials=st.multiselect('Trials to include',trial_options,default=trial_options[:3],key=key+'trials')
        mapped=st.checkbox('I verified the selected columns and paired trial IDs',value=known_layout,key=key+'mapped')
        if not mapped:
            st.info('Review the mapping and confirm it to analyze this subject.'); continue
        try: data=curves(exported,mapping,trials)
        except (ValueError,KeyError) as exc:
            st.error(str(exc)); continue
        # Axis sign changes are explicit, recorded and applied equally to all selected right trials.
        flips=st.multiselect('Right knee axes to multiply by −1 (only if required by your coordinate convention)',
                            ['X','Y','Z'],key=key+'flips')
        for axis in flips:
            k=('Right','Knee Angles',axis)
            if k in data: data[k]=-data[k]
        synchronized=st.checkbox('Both limbs share the same movement normalization (confirmed for these exports)',value=True,key=key+'sync')
        st.caption('Knee X = flexion; Y = varus; Z = rotation. Verify that left/right axis signs represent comparable anatomical directions before interpreting differences.')
        if not synchronized:
            st.warning('Plots and differences are descriptive only until alignment/conventions are confirmed; inference is withheld.')
        units={signal:angle_unit if 'Angles' in signal else force_unit if signal=='GRF' else power_unit if 'Power' in signal else emg_unit
               for _,signal,_ in data}
        channels=sorted(set((signal,axis) for _,signal,axis in data))
        events={}; candidate_notes={}
        phase_key=key+hashlib.sha256((mapping.to_csv(index=False)+str(trials)+str(flips)+str(threshold)).encode()).hexdigest()[:12]
        with st.expander('Phase identification: review events separately for each trial'):
            st.write('Descent: onset → deepest position; ascent: deepest position → takeoff; flight: takeoff → landing; landing absorption: landing → end. These are movement phases, not force-derived unweighting/braking phases.')
            st.caption('Automatic candidates require initial force contact and a sustained flight interval. Knee X maximum is only a proxy for the deepest position. If the recording omits a phase, leave it unconfirmed; whole-cycle metrics remain available.')
            phase_rows=[]
            for i,trial in enumerate(trials):
                candidate,note=phase_candidates(data,i,threshold); candidate_notes[trial]=note
                st.write(f'Trial {trial}: {note}')
                phase_rows.append(dict(trial=trial,confirmed=False,onset=None,bottom=candidate.get('bottom'),
                                       takeoff=candidate.get('takeoff'),landing=candidate.get('landing'),end=None))
            phase_df=st.data_editor(pd.DataFrame(phase_rows),hide_index=True,disabled=['trial'],key=phase_key+'events',
                column_config={k:st.column_config.NumberColumn(min_value=0,max_value=100,step=1) for k in ['onset','bottom','takeoff','landing','end']})
            for row in phase_df.to_dict('records'):
                if row['confirmed']:
                    try: events[int(row['trial'])]=validate_events(row)
                    except (ValueError,TypeError) as exc: st.error(f'Trial {row["trial"]}: {exc}')
            review_trial=st.selectbox('Trial to inspect for phase review',trials,key=key+'review_trial')
            review_candidate,_=phase_candidates(data,trials.index(review_trial),threshold)
            review_fig=phase_review_figure(subject,data,trials.index(review_trial),review_trial,threshold,events.get(review_trial,review_candidate))
            st.pyplot(review_fig)
            figures[f'{safe(subject)}/phase_review_trial_{review_trial}']=review_fig
            st.caption('Only valid, confirmed rows contribute phase metrics. Counts in the summary show how many trials contributed. Markers on average plots are mean event locations, not events detected from the mean waveform.')
        mean_events={k:float(np.mean([e[k] for e in events.values()])) for k in ['onset','bottom','takeoff','landing','end']} if events else None
        preview=st.multiselect('Graphs to preview',channels,default=[c for c in channels if c[0]=='Knee Angles'],
                              format_func=lambda c:channel_label(*c),key=key+'preview')
        for signal,axis in channels:
            fig=comparison(subject,signal,axis,data,units[signal],mean_events,show_trials)
            figures[f'{safe(subject)}/{safe(signal)}_{safe(axis)}']=fig
            if (signal,axis) in preview: st.pyplot(fig)
        trial_table,summary=metrics(subject,data,trials,events,units)
        for table in [trial_table,summary]:
            table['operated_limb']=surgical
            table['limb_role']=table.limb.map(lambda limb: 'operated' if surgical==limb else 'contralateral' if surgical in ('Left','Right') and limb in ('Left','Right') else 'unspecified')
        assets[f'{safe(subject)}/trial_metrics.csv']=csv_bytes(trial_table)
        assets[f'{safe(subject)}/summary_metrics.csv']=csv_bytes(summary)
        trial_tables.append(trial_table); summary_tables.append(summary)
        average_tables.append(average_table(subject,data)); diff_tables.append(differences(subject,data))
        raw_frames=[]
        for (limb,signal,axis),arr in data.items():
            for trial,values in zip(trials,arr):
                raw_frames.append(pd.DataFrame(dict(subject=subject,trial=trial,limb=limb,signal=signal,axis=axis,percent=np.arange(101),value=values)))
        assets[f'{safe(subject)}/selected_trial_waveforms.csv']=csv_bytes(pd.concat(raw_frames,ignore_index=True))
        st.write('Subject metrics: mean and SD of trial-level metrics')
        st.dataframe(summary,hide_index=True)
        st.caption('Normalized areas have units × percentage points, not impulse/work/iEMG in seconds. EMG RMS describes the supplied waveform; no raw-EMG filtering or MVC normalization is inferred.')
        assets[f'{safe(subject)}/column_mapping.csv']=csv_bytes(mapping)
        assets[f'{safe(subject)}/phase_review.csv']=csv_bytes(phase_df)
        meta=dict(subject=subject,source=upload.name,sha256=exported.digest,operated_limb=surgical,
                  trials=trials,right_knee_sign_flips=flips,aligned=synchronized,units=units,knee_axis_meanings=KNEE_PLANES,confirmed_original_layout=known_layout,
                  phase_events=events,candidate_notes=candidate_notes,duplicate=duplicated)
        metadata.append(meta)
        if synchronized:
            with st.expander('Subject-level knee inference',expanded=parametric):
                for axis in ['X','Y','Z']:
                    l=data.get(('Left','Knee Angles',axis)); r=data.get(('Right','Knee Angles',axis))
                    if l is not None and r is not None:
                        add_stats(subject,r,l,axis,alpha,family,parametric,assets,figures,records)
                    else: st.warning(f'Knee {axis}: missing paired channels.')
        if not duplicated and synchronized:
            subjects.append((subject,data,exported.digest,surgical))

if not summary_tables:
    st.stop()

st.subheader('Group knee comparison')
st.caption('Each independent subject contributes one mean curve per limb. Trials are never pooled as independent subjects. Contrast is Right − Left; mixed surgery sides are not an operated-versus-unoperated contrast.')
run_group=st.checkbox('Run group paired SPM (confirm the included IDs are independent subjects)')
if run_group:
    if len(subjects)<3:
        st.warning('At least three unique, aligned subjects are required for group parametric SPM. Identical uploads are excluded.')
    else:
        with st.expander('Group results',expanded=True):
            for axis in ['X','Y','Z']:
                eligible=[(s,d) for s,d,_,_ in subjects if all((side,'Knee Angles',axis) in d for side in ['Left','Right'])]
                if len(eligible)<3:
                    st.warning(f'Knee {axis}: fewer than three complete subjects.'); continue
                l=np.array([d[('Left','Knee Angles',axis)].mean(0) for _,d in eligible])
                r=np.array([d[('Right','Knee Angles',axis)].mean(0) for _,d in eligible])
                st.write(f'Knee {axis}: included IDs: {", ".join(s for s,_ in eligible)}')
                add_stats('GROUP',r,l,axis,alpha,family,True,assets,figures,records)
                for side,arr in [('Left',l),('Right',r)]:
                    table=pd.DataFrame(arr,columns=[f'p{i}' for i in range(101)])
                    table.insert(0,'subject',[s for s,_ in eligible])
                    assets[f'GROUP/knee_{axis}_{side}_input.csv']=csv_bytes(table)

st.subheader('Export results and selected figures')
selection=st.multiselect('Figures to export',list(figures),default=list(figures))
formats=st.multiselect('Image formats',['png','svg','pdf'],default=['png'])
assets['all_subjects_trial_metrics.csv']=csv_bytes(pd.concat(trial_tables,ignore_index=True))
assets['all_subjects_summary_metrics.csv']=csv_bytes(pd.concat(summary_tables,ignore_index=True))
assets['all_subjects_average_waveforms.csv']=csv_bytes(pd.concat(average_tables,ignore_index=True))
assets['all_subjects_limb_differences.csv']=csv_bytes(pd.concat(diff_tables,ignore_index=True))
assets['inference_results.csv']=csv_bytes(pd.DataFrame([{k:v for k,v in record.items() if k!='clusters'} | {'cluster_count':len(record.get('clusters',[])) if 'clusters' in record else None} for record in records]))
cluster_rows=[dict(subject=record['subject'],plane=record['plane'],alpha_plane=record['alpha_plane'],**cluster) for record in records for cluster in record.get('clusters',[])]
assets['spm_clusters.csv']=csv_bytes(pd.DataFrame(cluster_rows,columns=['subject','plane','alpha_plane','start_percent','end_percent','p_cluster']))
assets['analysis_manifest.json']=json.dumps(dict(subjects=metadata,alpha=alpha,family_size=family,
    family_scope=family_scope,spm_enabled=parametric,group_enabled=run_group,
    force_threshold_fraction=threshold,software='CMJ bilateral analysis 1.0',inference=records),indent=2).encode()
assets['READ_ME.txt']=b'''Curves: pointwise arithmetic mean and sample SD across selected trials.
Summary metrics: mean, sample SD and count of trial-level metrics, not extrema of the average curve.
Events: only reviewed valid trial events are used. All time coordinates are normalized percentages.
SPM: paired, two-sided; Bonferroni across three planes (or chosen larger family). Cluster p-values
are random-field corrected within the plane at the recorded alpha_plane. No cluster is not proof of equivalence.
Exact comparison: paired sign-flip max-|t| correction across nodes plus plane Bonferroni.
With 2 or 3 trials, minimum two-sided exact p is 0.5 or 0.25. Parametric small-n SPM is exploratory.
Symmetry = 200*(Right-Left)/(|Right|+|Left|); undefined at zero denominator, unstable near zero.
RMSE is between the mean curves. This is not a return-to-sport clearance score.
Missing contralateral channels are not imputed. Source EMG processing/normalization must be verified.
No seconds-based impulse, work, jump height or phase duration is inferred from 101-point files.
'''
if st.button('Prepare download',type='primary'):
    with st.spinner('Rendering selected figures…'):
        out=io.BytesIO()
        with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
            for name,content in assets.items(): z.writestr(name,content)
            for name in selection:
                for fmt in formats: z.writestr(f'{name}.{fmt}',figure_bytes(figures[name],fmt))
        st.download_button('Download CMJ results ZIP',out.getvalue(),'cmj_results.zip','application/zip')
for fig in figures.values(): plt.close(fig)
