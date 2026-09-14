"""Opt-in audit on three existing, read-only local P. aeruginosa assemblies.

These are a workflow-validation cohort, not an epidemiologically defined outbreak.
Run with WMLSTUDIO_REALDATA=1; reports go only to the test output directory.
"""

import hashlib
import json
import os
from pathlib import Path

import pytest

from wmlstudio.app import MainWindow
from wmlstudio.export import write_json, write_review_report
from wmlstudio.investigation import InvestigationStore
from wmlstudio.typing import call_assembly, load_scheme


@pytest.mark.realdata
def test_real_assemblies_weekly_review_and_selected_report(qtbot, tmp_path):
    if os.environ.get('WMLSTUDIO_REALDATA') != '1':
        pytest.skip('Opt in to known local sequence validation with WMLSTUDIO_REALDATA=1.')
    repository = Path(__file__).resolve().parents[1]
    root = Path('/home/iowa/Desktop/tesseract/work/eval/paeruginosa_final')
    known = [('GCF059737285v1', '155'), ('GCF059544455v1', '2952'), ('GCF049679795v1', '1858')]
    paths = [root / name / 'spades' / 'contigs.fasta' for name, _ in known]
    if not all(path.is_file() for path in paths):
        pytest.skip('Known local validation assemblies are not available on this machine.')
    before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    scheme = load_scheme(repository / 'src/wmlstudio/resources/schemes/paeruginosa')
    output = Path(os.environ.get('WMLSTUDIO_REALDATA_OUTPUT', str(tmp_path / 'evidence'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    window = MainWindow(storage_root=tmp_path / 'app')
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(40)
    identifiers = []
    store = InvestigationStore(window.project)
    try:
        for path, (name, expected_st) in zip(paths, known):
            result = call_assembly(path, scheme)
            assert result['st'] == expected_st and result['status'] == 'complete'
            sid = window.project.add_sample(path, name=name)
            window.project.set_result(sid, result)
            window.project.set_metadata(sid, {'organism': {'genus': 'Pseudomonas', 'species': 'aeruginosa'},
                'annotations': {'validation_scope': 'Independent local assemblies; no outbreak ground truth'}})
            identifiers.append(sid)
        plan = store.save('Real-data workflow audit · classical MLST', identifiers[:2], scheme=scheme.name,
                          scheme_digest=scheme.digest, threshold=1, min_overlap=1.0,
                          protocol='Software validation only. Seven-locus MLST is not a cgMLST outbreak investigation.', include_new=False)
        window.select_investigation(plan['id'])
        first = store.snapshot(plan['id'])
        assert first['reuse']['computed_pairs'] == 1
        store.save(plan['name'], identifiers, investigation_id=plan['id'], scheme=scheme.name,
                   scheme_digest=scheme.digest, threshold=1, min_overlap=1.0, protocol=plan['protocol'])
        window.select_investigation(plan['id'])
        latest = store.snapshot(plan['id'])
        assert latest['reuse']['reused_pairs'] == 1 and latest['reuse']['computed_pairs'] == 2
        assert all(pair['shared_loci'] == pair['total_loci'] == 7 for pair in latest['pairs'])
        window.report_isolate_proximity(identifiers[0])
        window.write_pdf_report(output / 'real-isolate-proximity.pdf', selected_ids={identifiers[0]})
        write_review_report(window.report_records(), output / 'real-isolate-proximity.html',
                            selected_ids={identifiers[0]}, investigation=latest,
                            options=window.report_options(), graph_png=window.report_graph_png())
        write_json(window.report_records(), output / 'real-investigation-evidence.json',
                   selected_ids=identifiers, investigation=latest)
        document = json.loads((output / 'real-investigation-evidence.json').read_text())
        assert [row['st'] for row in document['samples']] == [st for _, st in known]
        assert len(document['investigation']['proximity']) == 3
        assert [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths] == before
        assert window.project.get_sample(identifiers[0])['input_path'] == str(paths[0])
    finally:
        window.close()
