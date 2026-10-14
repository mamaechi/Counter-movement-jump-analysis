from pathlib import Path
from unittest.mock import patch
import io
import zipfile
import json
import streamlit
import pytest
from streamlit.testing.v1 import AppTest


@pytest.mark.parametrize('use_uploaded',[False,True])
def test_app_upload_review_analyze_and_export(use_uploaded):
    paths=list(Path('/workspace/attachments').rglob('SUBID* ACL.txt'))
    if use_uploaded and not paths:
        pytest.skip('Private sample uploads are not bundled.')
    if not use_uploaded:
        from test_analysis import sample_bytes
        contents=sample_bytes()
    else: contents=paths[0].read_bytes()
    uploads=[]
    for name in ['subject1.txt','subject2.txt']:
        f=io.BytesIO(contents); f.name=name; uploads.append(f)
    with patch('streamlit.file_uploader',return_value=uploads):
        app=AppTest.from_file(Path(__file__).resolve().parents[1] / 'app.py',default_timeout=60).run()
        assert not app.exception
        for widget in app.checkbox:
            if 'verified the selected columns' in widget.label: widget.set_value(True)
        app.run()
        assert not app.exception
        for widget in app.checkbox:
            if 'Both limbs share' in widget.label or 'Enable exploratory' in widget.label: widget.set_value(True)
        app.run()
        assert not app.exception
        assert any('Identical file content' in w.value for w in app.warning)
        for widget in app.checkbox:
            if widget.label.startswith('Run group'): widget.set_value(True)
        app.run()
        assert not app.exception
        assert any('three unique' in w.value for w in app.warning)
        choices=next(w for w in app.multiselect if w.label=='Figures to export')
        selected=choices.options[0]
        choices.set_value([selected])
        next(w for w in app.multiselect if w.label=='Image formats').set_value(['png','svg','pdf'])
        next(w for w in app.button if w.label=='Prepare download').click()
        with patch('streamlit.download_button',wraps=streamlit.download_button) as download:
            app.run()
            assert not app.exception
            archive=zipfile.ZipFile(io.BytesIO(download.call_args.args[1]))
            assert archive.read(selected+'.png').startswith(b'\x89PNG')
            assert archive.read(selected+'.pdf').startswith(b'%PDF')
            assert b'<svg' in archive.read(selected+'.svg')
            manifest=json.loads(archive.read('analysis_manifest.json'))
            assert len(manifest['subjects'])==2
            assert manifest['subjects'][1]['duplicate']
            assert 'all_subjects_summary_metrics.csv' in archive.namelist()
        assert len(app.get('download_button'))==1
