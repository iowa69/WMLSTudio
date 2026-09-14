import sqlite3

import pytest

from wmlstudio.project import Project


def test_summary_queries_do_not_decode_large_locus_arrays(tmp_path, monkeypatch):
    with Project(tmp_path / 'summary.wmlstudio') as project:
        sid = project.add_profile('A', {'scheme': 'Primary', 'scheme_digest': 'primary', 'alleles': {'a': '1'}, 'st': '20'})
        project.set_analysis(sid, {'scheme': 'Core', 'scheme_digest': 'core', 'input_sha256': 'a' * 64,
                                  'status': 'complete', 'analysis_kind': 'cgmlst', 'scheme_path': '/reference',
                                  'alleles': {f'locus_{i}': '1' for i in range(3000)}, 'calls': [{'large': 'A' * 100000}]})
        monkeypatch.setattr('wmlstudio.project.json.loads', lambda *args: pytest.fail('Large JSON was decoded in Python for a GUI label'))
        summaries = project.analysis_summaries(sid)
        core = next(row for row in summaries if row['scheme_digest'] == 'core')
        assert core['locus_count'] == 3000 and core['scheme'] == 'Core'
        assert core['input_sha256'] == 'a' * 64 and core['analysis_kind'] == 'cgmlst'
        assert 'alleles' not in core and 'calls' not in core
        with pytest.raises(KeyError):
            project.analysis_summaries('missing')


def test_jsonless_sqlite_fallback_has_identical_summaries(tmp_path):
    with Project(tmp_path / 'fallback.wmlstudio') as project:
        sid = project.add_profile('A', {'scheme': 'One', 'scheme_digest': 'd', 'alleles': {'a': '1', 'b': None}})
        expected = project.analysis_summaries(sid)
        connection = project._connection
        class JsonlessConnection:
            def execute(self, sql, parameters=()):
                if 'json_extract' in sql:
                    raise sqlite3.OperationalError('no such function: json_extract')
                return connection.execute(sql, parameters)
        project._connection = JsonlessConnection()
        try:
            assert project.analysis_summaries(sid) == expected
        finally:
            project._connection = connection
