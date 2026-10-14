import io
from pathlib import Path
import numpy as np
import pytest
from cmj.core import read_export, curves, phase_candidates, metrics, validate_events, differences, confirmed_three_block_layout
from cmj.statistics import exact_max_t, paired_spm
from cmj.plots import comparison, spm_figure


def sample_bytes():
    headers=['Left Knee Angles','Right Knee Angles','Left GRF','Right GRF']*3
    axes=['X','X','Z','Z']*3
    rows=['\t'+'\t'.join(headers),'ITEM\t'+'\t'.join(axes)]
    for i in range(101):
        vals=[]
        for trial in range(3): vals += [i+trial,i+trial+2,0 if 50<=i<70 else 1,0 if 50<=i<70 else 1]
        rows.append('\t'.join(map(str,[i+1]+vals)))
    return '\n'.join(rows).encode()


def test_import_average_and_pairing():
    e=read_export(sample_bytes(),'example')
    d=curves(e,e.mapping,[1,2,3])
    np.testing.assert_allclose(d[('Left','Knee Angles','X')].mean(0),np.arange(101)+1)
    assert np.all(differences('S',d).query("signal == 'Knee Angles'").right_minus_left==2)
    e.mapping.loc[0,'trial']=2
    with pytest.raises(ValueError,match='one column'): curves(e,e.mapping,[1,2,3])


def test_invalid_rows_and_missing_trial():
    with pytest.raises(ValueError,match='101'): read_export(b'bad','x')
    e=read_export(sample_bytes(),'x'); e.mapping.loc[0,'include']=False
    with pytest.raises(ValueError,match='missing'): curves(e,e.mapping,[1,2,3])


def test_phase_detection_rejects_initial_flight():
    e=read_export(sample_bytes(),'x'); d=curves(e,e.mapping,[1,2,3])
    events,_=phase_candidates(d,0)
    assert events['takeoff']==50 and events['landing']==70
    d[('Left','GRF','Z')][0,0:4]=0; d[('Right','GRF','Z')][0,0:4]=0
    events,note=phase_candidates(d,0)
    assert not events and 'does not begin' in note


def test_metrics_use_trial_peaks_not_peak_of_mean():
    y=np.zeros((2,101)); y[0,20]=10; y[1,80]=10
    trial,summary=metrics('S',{('Left','EMG','X'):y},[1,2],{}, {})
    assert summary.maximum_mean.iloc[0]==10
    assert y.mean(0).max()==5
    with pytest.raises(ValueError): validate_events(dict(onset=0,bottom=40,takeoff=30,landing=70,end=100))


def test_fractional_phase_integration_and_counts():
    y=np.ones((3,101))*2
    events={1:dict(onset=0,bottom=20.5,takeoff=50,landing=70,end=100)}
    trial,summary=metrics('S',{('Left','EMG','X'):y},[1,2,3],events,{})
    descent=trial.query("phase=='descent'").iloc[0]
    assert descent.area_unit_percent==41
    assert summary.query("phase=='descent'").maximum_count.iloc[0]==1


def test_small_n_exact_and_parametric_guards():
    rng=np.random.default_rng(42); l=rng.normal(size=(3,101)); r=l+rng.normal(size=(3,101))+.8
    result,p=exact_max_t(r,l)
    assert result['minimum_two_sided_p']==.25 and p.min()>=.25
    assert not result['reject']
    with pytest.raises(ValueError,match='three pairs'): paired_spm(r[:2],l[:2])
    with pytest.raises(ValueError,match='Zero'): paired_spm(l,l)


def test_paired_spm_and_figure_exports():
    from scipy.ndimage import gaussian_filter1d
    rng=np.random.default_rng(3)
    l=gaussian_filter1d(rng.normal(size=(12,101)),4,axis=1)
    r=l+gaussian_filter1d(rng.normal(size=(12,101)),4,axis=1)+.6
    result,inf=paired_spm(r,l)
    assert result['n_pairs']==12 and result['df']==11 and result['reject']
    for fig in [comparison('test','Knee Angles','X',{('Left','Knee Angles','X'):l,('Right','Knee Angles','X'):r},'degrees'),spm_figure('test',inf)]:
        buf=io.BytesIO(); fig.savefig(buf,format='png'); assert buf.getvalue().startswith(b'\x89PNG')


def test_uploaded_examples_when_available():
    files=list(Path('/workspace/attachments').rglob('SUBID* ACL.txt'))
    if not files: pytest.skip('Private sample uploads are intentionally not bundled.')
    parsed=[read_export(p.read_bytes(),p.name) for p in files]
    assert len(parsed[0].mapping)==98
    assert confirmed_three_block_layout(parsed[0])
    assert not confirmed_three_block_layout(read_export(sample_bytes(),'different-layout'))
    assert parsed[0].digest==parsed[1].digest
    d=curves(parsed[0],parsed[0].mapping,[1,2,3])
    assert len(d)==15 and ('Left','Gluteus Maximus','X') not in d
    assert d[('Left','Knee Angles','X')][0,0]==pytest.approx(4.031282)
    for axis in 'XYZ':
        result,inf=paired_spm(d[('Right','Knee Angles',axis)],d[('Left','Knee Angles',axis)])
        assert result['n_pairs']==3 and len(inf.z)==101
