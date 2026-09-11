# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Proof that no attacker-controlled string can execute in the HTML report.

Covers docs/ARCHITECTURE.md section 11.1 rules 3-8. A FASTA arrives from a
sequencing provider, the report opens from ``file://`` with a full local
origin, and it gets emailed onward; contig ids, filenames, labels, locus names
and scheme names are therefore treated as hostile input.

Runs under ``python3 -m pytest tests/test_report_escaping.py`` and as a plain
script (``python3 tests/test_report_escaping.py``).
"""

from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_report_html import walk
from tests.test_tabular import mk_result, mk_sample
from wmlst import report
from wmlst.engine import (
    AlleleCall,
    Hit,
    NovelAllele,
    RunConfig,
    SampleResult,
    SchemeScore,
)

DBDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "db")

#: Every payload is injected into every untrusted field.
PAYLOADS = [
    "<script>alert(1)</script>",
    "><script>alert(1)</script>",
    "\"><img src=x onerror=alert(1)>",
    "'><svg onload=alert(1)>",
    "</script><script>alert(1)</script>",
    "</title></head><body onload=alert(1)>",
    "javascript:alert(1)",
    "` onmouseover=alert(1) x=",
    "</style><style>body{display:none}</style>",
    "\" style=\"position:fixed;inset:0\" data-x=\"",
    "</pre></details><script>alert(1)</script>",
    "]]><script>alert(1)</script>",
    " alert(1) ",
    "a\x00b\x07c\x1bd",
    "&lt;already escaped&gt;",
]


def poisoned_sample(payload: str) -> SampleResult:
    """One sample with the payload in EVERY untrusted position."""
    hit = Hit(
        sseqid=payload, scheme=payload, locus=payload, allele=payload,
        slen=465, length=465, nident=465, qseqid=payload, qstart=1, qend=465,
        qseq="ACGT", sstrand="plus", pct_identity=100.0, pct_coverage=100.0,
        call_kind="exact", index=1,
    )
    call = AlleleCall(locus=payload, code=payload, symbol="exact", best=hit, hits=(hit,))
    novel = NovelAllele(scheme=payload, locus=payload, md5="0" * 32, seq="ACGTACGT",
                        source_label=payload, fasta_id=payload, nearest_allele=payload,
                        length_bp=8)
    return SampleResult(
        path=payload, label=payload, scheme=payload, st=payload,
        signature=payload, score=100, status="PERFECT",
        alleles=(call,),
        candidates=(SchemeScore(payload, payload, payload, 100, 1),
                    SchemeScore(payload + "2", payload, payload, 100, 1)),
        novel=(novel,),
        warnings=("WARNING: " + payload,),
        n_contigs=1, total_bp=100, hits_seen=1, hits_kept=1, elapsed_s=0.1,
    )


def poisoned_html(payload: str) -> str:
    cfg = RunConfig(dbdir=DBDIR, files=(payload,), label=payload,
                    scheme=payload, exclude=frozenset({payload}))
    result = mk_result([poisoned_sample(payload), poisoned_sample(payload + "-b")],
                       cfg, novel=poisoned_sample(payload).novel)
    return report.render_html(result, report.HtmlOptions(title=payload))


def all_payload_html():
    return [(p, poisoned_html(p)) for p in PAYLOADS]


# ---------------------------------------------------------------------------
# h() / attr() / js_json() -- the interpolators themselves
# ---------------------------------------------------------------------------

def test_h_escapes_the_seven_dangerous_characters():
    assert report.h("&") == "&amp;"
    assert report.h("<") == "&lt;"
    assert report.h(">") == "&gt;"
    assert report.h('"') == "&quot;"
    assert report.h("'") == "&#x27;"
    assert report.h("`") == "&#x60;"
    assert report.h("=") == "&#x3d;"


def test_h_neutralises_controls_but_keeps_tab_and_newline():
    assert report.h("a\x00b") == "a�b"
    assert report.h("a\x1bb") == "a�b"
    assert report.h("a\x9bb") == "a�b"
    assert report.h("a b") == "a�b"
    assert report.h("a b") == "a�b"
    assert report.h("a\rb") == "a�b"
    assert report.h("a\tb\nc") == "a\tb\nc"


def test_h_handles_none_and_non_strings():
    assert report.h(None) == ""
    assert report.h(17) == "17"


def test_attr_is_always_double_quoted():
    for payload in PAYLOADS:
        value = report.attr(payload)
        assert value.startswith('"') and value.endswith('"')
        assert '"' not in value[1:-1]


def test_js_json_cannot_close_a_script_element():
    for payload in PAYLOADS:
        blob = report.js_json({"x": payload})
        assert "</script" not in blob.lower()
        assert "<" not in blob and ">" not in blob and "&" not in blob
        assert json.loads(blob)["x"] == payload


# ---------------------------------------------------------------------------
# the rendered document
# ---------------------------------------------------------------------------

def test_no_payload_survives_as_live_markup():
    # An escaped payload legitimately still contains its own words as inert
    # TEXT ("javascript:", "onerror=") -- what must never appear is the markup.
    for payload, html in all_payload_html():
        for needle in ("<script>alert", "<img src=x", "<svg onload",
                       "<body onload", "</style><style>"):
            assert needle not in html, (payload, needle)
        for _tag, name, value in walk(html).attrs:
            if name in ("href", "src", "action", "formaction", "xlink:href"):
                assert not (value or "").strip().lower().startswith("javascript:"), payload


def test_the_only_scripts_are_ours_and_they_contain_no_payload():
    for payload, html in all_payload_html():
        parser = walk(html)
        assert parser.tags.count("script") == 2, payload
        island = [s for s in parser.scripts if s.strip().startswith("[")]
        code = [s for s in parser.scripts if s not in island]
        assert len(island) == 1 and len(code) == 1, payload
        # the executable script is verbatim JS_TEXT: no payload can reach it
        assert code[0].strip() == report.JS_TEXT.strip(), payload
        # the island is inert data -- it may quote the payload, but only in a
        # form that cannot close the element or be read as markup
        assert "</script" not in island[0].lower(), payload
        assert "<" not in island[0] and ">" not in island[0], payload
        json.loads(island[0])


def test_no_event_handler_attribute_is_ever_produced():
    for payload, html in all_payload_html():
        for tag, name, _value in walk(html).attrs:
            assert not name.startswith("on"), (payload, tag, name)


def test_no_payload_reaches_a_style_attribute_or_a_style_element():
    for payload, html in all_payload_html():
        parser = walk(html)
        assert parser.tags.count("style") == 1, payload
        for _tag, name, value in parser.attrs:
            if name == "style":
                assert re.fullmatch(r"width:\d{1,3}%", value), (payload, value)
        block = html[html.index("<style>"):html.index("</style>")]
        assert "alert" not in block and payload not in block, payload


def test_document_still_parses_and_balances_under_every_payload():
    for payload, html in all_payload_html():
        parser = walk(html)
        assert parser.errors == [], (payload, parser.errors)
        assert parser.stack == [], (payload, parser.stack)


def test_contig_id_specifically():
    payload = "<script>alert(1)</script>"
    html = poisoned_html(payload)
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    parser = walk(html)
    # the contig appears both as cell text and inside a title attribute
    titles = [v for t, n, v in parser.attrs if n == "title"]
    assert payload in titles           # decoded by the parser == inert text
    assert "<script>alert(1)</script>" not in html.replace(
        "&lt;script&gt;alert(1)&lt;/script&gt;", "")


def test_filename_and_label_specifically():
    payload = "\"><img src=x onerror=alert(1)>"
    html = poisoned_html(payload)
    assert "&quot;&gt;&lt;img" in html
    assert "<img src=x" not in html
    for _tag, name, _value in walk(html).attrs:
        assert name != "onerror"


def test_locus_name_specifically():
    payload = "</td></tr><script>alert(1)</script>"
    html = poisoned_html(payload)
    parser = walk(html)
    assert parser.errors == [], parser.errors
    assert parser.tags.count("script") == 2
    assert "&lt;/td&gt;&lt;/tr&gt;" in html


def test_scheme_name_specifically():
    payload = "</style><style>body{display:none}</style>"
    html = poisoned_html(payload)
    assert walk(html).tags.count("style") == 1
    css = html[html.index("<style>"):html.index("</style>")]
    assert "body{display:none}" not in css
    assert "</style><style>" not in html


def test_novel_header_and_sequence_are_escaped():
    payload = "><script>alert(1)</script>"
    html = poisoned_html(payload)
    assert "&gt;&lt;script&gt;" in html
    bad = NovelAllele("s", "g", "0" * 32, "ACGT<script>alert(1)</script>", payload,
                      "s.g-" + "0" * 32, "1", 29)
    html2 = report.render_html(
        mk_result([mk_sample("x.fa", "-", "-", "NONE", 0)], RunConfig(dbdir=DBDIR),
                  novel=[bad]))
    assert "not shown here" in html2        # not DNA -> replaced by a notice
    assert "ACGT&lt;script&gt;" not in html2
    assert '<pre class="seq mono" translate="no">' not in html2


def test_warning_text_is_escaped():
    payload = "<b>bold</b>"
    html = poisoned_html(payload)
    assert "<b>bold</b>" not in html
    assert "&lt;b&gt;bold&lt;/b&gt;" in html


def test_batch_search_attribute_cannot_break_out():
    payload = "\" onmouseover=\"alert(1)"
    html = poisoned_html(payload)
    for _tag, name, _value in walk(html).attrs:
        assert not name.startswith("on"), name
    assert "onmouseover=\"alert" not in html


def test_title_element_cannot_be_closed_early():
    payload = "</title></head><body onload=alert(1)>"
    html = poisoned_html(payload)
    title = re.search(r"<title>(.*?)</title>", html, re.S)
    assert title
    assert "<" not in title.group(1) and ">" not in title.group(1)
    assert html.count("<title>") == 1
    assert html.count("</title>") == 1


def test_unicode_line_separators_cannot_split_a_js_string():
    html = poisoned_html(" alert(1) ")
    assert " " not in html and " " not in html


def test_compat_outputs_are_untouched_by_escaping():
    # h() must never leak into a byte-identity surface.
    import io
    payload = "<script>alert(1)</script>"
    result = mk_result([mk_sample(payload, "saureus", "1", "PERFECT", 100, "arcC(1)")],
                       RunConfig(full=True))
    buf = io.StringIO()
    report.write_tsv(result, buf)
    assert buf.getvalue().startswith("FILE\t")
    assert payload in buf.getvalue()
    assert "&lt;" not in buf.getvalue()


# ---------------------------------------------------------------------------

def _main():
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("ok   " + name)
            except AssertionError as exc:
                failures += 1
                print("FAIL " + name + ": " + str(exc))
            except Exception as exc:
                failures += 1
                print("ERROR " + name + ": " + repr(exc))
    print("{0} failure(s)".format(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(_main())
