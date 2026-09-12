# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Structural, accessibility and self-containment tests for the HTML report.

Covers docs/ARCHITECTURE.md section 11. Fixtures are built by hand from the
frozen dataclasses in ``wmlst.engine``.

Runs under ``python3 -m pytest tests/test_report_html.py`` and as a plain
script (``python3 tests/test_report_html.py``).
"""

from __future__ import annotations

import os
import re
import sys
import tempfile
from html.parser import HTMLParser

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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

HERE = os.path.dirname(os.path.abspath(__file__))
DBDIR = os.path.join(os.path.dirname(HERE), "db")

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link",
        "meta", "param", "source", "track", "wbr"}


class Walker(HTMLParser):
    """Tracks tag balance and collects attributes and script bodies."""

    def __init__(self):
        HTMLParser.__init__(self, convert_charrefs=True)
        self.stack = []
        self.errors = []
        self.attrs = []          # (tag, name, value)
        self.scripts = []
        self.tags = []
        self._script = None

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)
        for name, value in attrs:
            self.attrs.append((tag, name, value))
        if tag == "script":
            self._script = []
        if tag not in VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        self.tags.append(tag)
        for name, value in attrs:
            self.attrs.append((tag, name, value))

    def handle_endtag(self, tag):
        if tag == "script" and self._script is not None:
            self.scripts.append("".join(self._script))
            self._script = None
        if tag in VOID:
            return
        if not self.stack:
            self.errors.append("closing </{0}> with an empty stack".format(tag))
            return
        if self.stack[-1] != tag:
            self.errors.append("</{0}> closes <{1}>".format(tag, self.stack[-1]))
            return
        self.stack.pop()

    def handle_data(self, data):
        if self._script is not None:
            self._script.append(data)


def walk(html: str) -> Walker:
    parser = Walker()
    parser.feed(html)
    parser.close()
    return parser


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def mk_hit(**kw):
    base = {"sseqid": "sepidermidis.arcC_16", "scheme": "sepidermidis", "locus": "arcC",
                "allele": "16", "slen": 465, "length": 465, "nident": 465, "qseqid": "contig_00001",
                "qstart": 1200, "qend": 1664, "qseq": "ACGTACGT", "sstrand": "plus",
                "pct_identity": 100.0, "pct_coverage": 100.0, "call_kind": "exact", "index": 1}
    base.update(kw)
    return Hit(**base)


def rich_sample():
    """One sample exercising every symbol, a tie, a warning and evidence."""
    hit_exact = mk_hit()
    hit_novel = mk_hit(allele="80", locus="infB", nident=460, call_kind="novel",
                       pct_identity=98.92473118279571, pct_coverage=98.92473118279571,
                       sstrand="minus", qstart=900, qend=440, index=2,
                       qseqid="NODE_2_length_5000_cov_3.1")
    hit_part = mk_hit(allele="114", locus="6pgd", length=300, nident=300,
                      call_kind="partial", pct_identity=100.0,
                      pct_coverage=64.51612903225806, index=3)
    hit_a = mk_hit(allele="1", locus="MLST_atpA", index=4)
    hit_b = mk_hit(allele="2", locus="MLST_atpA", index=5, qseqid="contig_00002")
    calls = (
        AlleleCall("arcC", "16", "exact", hit_exact, (hit_exact,)),
        AlleleCall("infB", "~80", "novel", hit_novel, (hit_novel,)),
        AlleleCall("6pgd", "114?", "partial", hit_part, (hit_part,)),
        AlleleCall("MLST_atpA", "1,2", "multiple", hit_a, (hit_a, hit_b)),
        AlleleCall("g3pd", "-", "missing", None, ()),
        AlleleCall("frdB", "0", "null", None, ()),
    )
    candidates = (
        SchemeScore("sepidermidis", "184", "16/~80/114?/1,2/-/0", 100, 6,
                    n_exact=2, n_novel=1, n_partial=1, n_missing=1, n_null=1),
        SchemeScore("saureus", "-", "-/-/-/-/-/-", 100, 6, n_missing=6),
        SchemeScore("ecoli", "-", "-/-/-/-/-/-", 20, 6, below_minscore=True),
        SchemeScore("vcholerae_2", "-", "-/-/-", 0, 3, excluded=True),
    )
    novel = (NovelAllele("sepidermidis", "infB", "d974fcafa3e20f56dd36e6671f968ff5",
                         "ACGT" * 40, "sample one.fna",
                         "sepidermidis.infB-d974fcafa3e20f56dd36e6671f968ff5",
                         "80", 160),)
    return SampleResult(
        path="/data/sample one.fna", label="sample one.fna", scheme="sepidermidis",
        st="184", signature="16/~80/114?/1,2/-/0", score=100, status="PERFECT",
        alleles=calls, candidates=candidates, tied=candidates[:2], novel=novel,
        warnings=("WARNING: sepidermidis(184)==saureus(-) score=100 /data/sample one.fna",),
        n_contigs=42, total_bp=2_812_345, hits_seen=1423, hits_kept=88,
        elapsed_s=1.5,
    )


def one_result():
    cfg = RunConfig(full=True, dbdir=DBDIR, html_path="r.html")
    result = mk_result([rich_sample()], cfg)
    return result


def batch_result():
    cfg = RunConfig(full=True, dbdir=DBDIR)
    samples = [
        rich_sample(),
        mk_sample("messy.fa", "hparasuis", "-", "MISSING", 51,
                  "atpD(~66);g3pd(-)"),
        mk_sample("none.fa", "-", "-", "NONE", 0, ""),
        mk_sample("null.fa", "-", "-", "NONE", 0, "", failed=True,
                  error_text="The input appears to be empty"),
        mk_sample("dup.fa", "hparasuis", "-", "MISSING", 51, "atpD(~66);g3pd(-)"),
    ]
    return mk_result(samples, cfg)


# ---------------------------------------------------------------------------
# 11.1  self-containment
# ---------------------------------------------------------------------------

def test_document_shape():
    html = report.render_html(one_result())
    assert html.startswith("<!DOCTYPE html>\n<html lang=\"en\">")
    assert html.rstrip().endswith("</html>")
    assert "<title>" in html


def test_html_parses_and_tags_balance():
    for result in (one_result(), batch_result(), mk_result([], RunConfig(dbdir=DBDIR))):
        parser = walk(report.render_html(result))
        assert parser.errors == [], parser.errors
        assert parser.stack == [], parser.stack


def test_no_external_references():
    html = report.render_html(batch_result())
    for needle in ("http://", "https://", "//cdn", "@import",
                   '<link rel="stylesheet"', "<script src", "srcset="):
        assert needle not in html, needle
    # the only url()/data: use is the base64 vendor mark
    for match in re.findall(r"url\(\s*([^)]*)", html):
        assert match.strip().strip("\"'").startswith("data:"), match
    for _tag, name, value in walk(html).attrs:
        if name == "src":
            assert value.startswith("data:image/svg+xml;base64,"), value


def test_csp_and_referrer_locks():
    html = report.render_html(one_result())
    assert "Content-Security-Policy" in html
    assert "default-src 'none'" in html
    assert 'name="referrer" content="no-referrer"' in html


def test_branding_is_present():
    html = report.render_html(one_result())
    assert "IOWA-Tech" in html
    assert "Giovanni Lorenzin" in html
    assert "WMLST" in html
    assert "Torsten Seemann" in html
    assert "PubMLST" in html


def test_render_is_deterministic():
    result = batch_result()
    assert report.render_html(result) == report.render_html(result)


# ---------------------------------------------------------------------------
# 11.1 rule 4 / 11.6  the JavaScript contract
# ---------------------------------------------------------------------------

def test_no_forbidden_javascript_apis():
    html = report.render_html(batch_result())
    for needle in ("innerHTML", "outerHTML", "insertAdjacentHTML", "document.write",
                   "eval(", "new Function", "XMLHttpRequest", "WebSocket",
                   "fetch(", "javascript:"):
        assert needle not in html, needle
    assert not re.search(r"setTimeout\(\s*['\"]", html)
    assert not re.search(r"\bimport\s*\(", html)


def test_no_inline_event_handlers():
    for _tag, name, _value in walk(report.render_html(batch_result())).attrs:
        assert not name.startswith("on"), name


def test_js_text_budget_and_shape():
    lines = report.JS_TEXT.splitlines()
    assert len(lines) <= 120, len(lines)
    assert lines[0].startswith("(function ()")
    assert "'use strict';" in report.JS_TEXT
    assert report.JS_TEXT.rstrip().endswith("}());")


def test_exactly_two_script_elements_and_the_island_is_json():
    import json
    html = report.render_html(one_result())
    parser = walk(html)
    assert parser.tags.count("script") == 2
    island = [s for s in parser.scripts if s.strip().startswith("[")]
    assert island, "the JSON island is missing"
    payload = json.loads(island[0])
    assert payload[0]["scheme"] == "sepidermidis"


def test_progressive_enhancement_controls_start_hidden():
    html = report.render_html(batch_result())
    for control in ("theme-toggle", "expand-all", "reset-order"):
        match = re.search(r'id="{0}"[^>]*>'.format(control), html)
        assert match, control
        assert "hidden" in match.group(0), control


# ---------------------------------------------------------------------------
# 11.3  status badges
# ---------------------------------------------------------------------------

def test_every_status_has_glyph_word_class_and_data_attribute():
    html = report.render_html(one_result())
    for status, (glyph, cls, style, why) in report.STATUS_INFO.items():
        assert 'data-status="{0}"'.format(status) in html, status
        assert glyph in html, status
        assert why in html or report.h(why) in html, status
        assert cls in report.CSS_TEXT, cls
        assert style in ("solid", "double", "dashed", "thick")


def test_glyphs_are_aria_hidden_so_the_word_is_heard_once():
    html = report.render_html(one_result())
    for glyph, _cls, _style, _why in report.STATUS_INFO.values():
        for match in re.finditer(re.escape(glyph), html):
            start = html.rfind("<span", 0, match.start())
            assert 'aria-hidden="true"' in html[start:match.start()], glyph


def test_unknown_status_falls_back_without_a_colour_lookup():
    sample = mk_sample("x.fa", "saureus", "1", "SOMETHING_NEW", 42, "arcC(1)")
    html = report.render_html(mk_result([sample], RunConfig(dbdir=DBDIR)))
    assert "s-unknown" in html
    assert "SOMETHING_NEW" in html
    assert "s-unknown-bg" in report.CSS_TEXT


def test_card_border_style_differs_per_status_not_only_colour():
    styles = {info[2] for info in report.STATUS_INFO.values()}
    assert styles == {"solid", "double", "dashed", "thick"}
    for style in styles:
        assert ".card.b-{0}".format(style) in report.CSS_TEXT


# ---------------------------------------------------------------------------
# 11.4  locus table and evidence
# ---------------------------------------------------------------------------

def test_locus_table_is_a_real_table_with_scoped_headers():
    html = report.render_html(one_result())
    assert "<caption" in html
    assert '<th scope="col">Locus</th>' in html
    assert re.search(r'<th scope="row"[^>]*>arcC</th>', html)
    assert 'role="table"' not in html


def test_evidence_columns_show_the_right_arithmetic():
    html = report.render_html(one_result())
    assert "contig_00001" in html
    assert "1,200–1,664" in html          # en-dash, thousands separators
    assert "465 / 465" in html                  # alignment length vs allele length
    assert "98.9 %" in html                     # floored, never 99.0
    assert "64.5 %" in html                     # coverage = nident/slen
    assert "the quantity the --mincov threshold is applied to" in html


def test_minus_strand_renders_a_word_not_only_a_sign():
    html = report.render_html(one_result())
    assert "reverse" in html and "forward" in html
    assert "minus strand," in html


def test_multiple_code_gets_one_evidence_row_per_component():
    html = report.render_html(one_result())
    card = html[html.index('id="sample-1"'):]
    assert 'rowspan="2"' in card
    assert "contig_00002" in card


def test_codes_are_rendered_verbatim():
    html = report.render_html(one_result())
    for code in ("16", "~80", "114?", "1,2", "-", "0"):
        assert ">{0}<".format(report.h(code)) in html, code


def test_legend_always_lists_all_six_symbols():
    for result in (one_result(), mk_result([], RunConfig(dbdir=DBDIR))):
        html = report.render_html(result)
        for code, _symbol, text in report.SYMBOL_INFO:
            assert report.h(code) in html, code
            assert report.h(text)[:40] in html, code


def test_runner_up_table_explains_why_the_scheme_won():
    html = report.render_html(one_result())
    assert "Why this scheme won" in html
    assert "below --minscore 50" in html
    assert "Dropped by --exclude before scoring" in html
    assert "vcholerae_2" in html


def test_tie_is_surfaced_outside_a_fold():
    html = report.render_html(one_result())
    index = html.index("Tie:")
    card_start = html.index('id="sample-1"')
    fold_start = html.index("<details", card_start)
    assert card_start < index < fold_start


def test_tie_names_both_schemes_with_their_sequence_types():
    """5.14a -- the ambiguity must be legible, not just flagged."""
    html = report.render_html(one_result())
    assert "sepidermidis ST 184" in html
    assert "saureus (no ST)" in html
    assert "fit this assembly equally well at score 100" in html
    assert "confirm the species by another method" in html.lower()


def test_tied_rows_are_flagged_in_the_runner_up_table():
    html = report.render_html(one_result())
    assert "ties with the reported scheme" in html


def test_no_tie_block_when_the_winner_stands_alone():
    sample = mk_sample("solo.fa", "sepidermidis", "184", "PERFECT", 100, "arcC(16)")
    html = report.render_html(mk_result([sample], RunConfig(dbdir=DBDIR)))
    assert "fit this assembly equally well" not in html


def test_species_line_and_caveat():
    html = report.render_html(one_result())
    assert "Staphylococcus epidermidis" in html
    assert "not an independent identification" in html


def test_unmapped_scheme_says_so_instead_of_guessing():
    sample = mk_sample("x.fa", "abetadefensin", "1", "PERFECT", 100, "a(1)")
    html = report.render_html(mk_result([sample], RunConfig(dbdir=DBDIR)))
    assert ("not mapped to a species" in html) or ("Homo" in html)


def test_species_label_rules():
    assert report.species_label_for(DBDIR, "sepidermidis") == "Staphylococcus epidermidis"
    assert report.species_label_for(DBDIR, "no_such_scheme_at_all") is None
    label = report.species_label_for(DBDIR, "ecoli")
    assert label and label.startswith("Escherichia")


# ---------------------------------------------------------------------------
# 11.2  batch section
# ---------------------------------------------------------------------------

def test_batch_section_only_when_more_than_one_sample():
    assert 'id="batch-table"' not in report.render_html(one_result())
    html = report.render_html(batch_result())
    assert 'id="batch-table"' in html
    assert 'id="filter-count"' in html and 'role="status"' in html


def test_batch_rows_are_in_input_order_with_a_stable_index():
    html = report.render_html(batch_result())
    indices = re.findall(r'<tr data-index="(\d+)"', html)
    assert indices == ["1", "2", "3", "4", "5"]


def test_batch_search_attribute_is_prelowercased():
    html = report.render_html(batch_result())
    for value in re.findall(r'data-search="([^"]*)"', html):
        assert value == value.lower()


def test_rollup_groups_and_carries_the_outbreak_caveat():
    html = report.render_html(batch_result())
    assert "Sequence-type roll-up" in html
    assert "does not prove an outbreak link" in html


def test_failed_sample_gets_a_card_but_no_call():
    html = report.render_html(batch_result())
    assert "Not analysed" in html
    assert "The input appears to be empty" in html


def test_failed_card_is_styled_as_failed_not_as_a_real_status():
    html = report.render_html(batch_result())
    card = re.search(r'<article class="card ([^"]+)" id="sample-4" data-status="([^"]+)"', html)
    assert card, "the failed sample has no card"
    assert card.group(1) == "s-bad b-thick"
    assert card.group(2) == "FAILED"
    assert "FAILED" not in report.STATUS_INFO   # not part of the mlst vocabulary


def test_empty_run_renders_an_empty_state():
    html = report.render_html(mk_result([], RunConfig(dbdir=DBDIR)))
    assert "No input files were analysed" in html
    assert walk(html).errors == []


# ---------------------------------------------------------------------------
# novel alleles
# ---------------------------------------------------------------------------

def test_novel_section_wraps_at_60_and_warns_about_revcom():
    result = one_result()
    result = mk_result(list(result.samples), result.meta.config,
                       novel=result.samples[0].novel)
    html = report.render_html(result)
    assert "already reverse-complemented" in html
    assert "&gt;sepidermidis.infB-d974fcafa3e20f56dd36e6671f968ff5" in html
    body = html[html.index("<pre class=\"seq mono\" translate=\"no\">"):]
    seq_lines = [ln for ln in body.split("</pre>")[0].split("\n")[1:] if ln]
    assert all(len(ln) <= 60 for ln in seq_lines)


def test_non_dna_sequence_is_not_inlined():
    bad = NovelAllele("s", "g", "0" * 32, "ACGT<script>", "lbl", "s.g-" + "0" * 32,
                      "1", 12)
    html = report.render_html(mk_result([mk_sample("x.fa", "-", "-", "NONE", 0)],
                                        RunConfig(dbdir=DBDIR), novel=[bad]))
    assert "not shown here" in html
    assert "<script>ACGT" not in html


def test_oversize_novel_sequence_shows_a_notice():
    big = NovelAllele("s", "g", "1" * 32, "A" * 500, "lbl", "s.g-" + "1" * 32, "1", 500)
    html = report.render_html(
        mk_result([mk_sample("x.fa", "-", "-", "NONE", 0)], RunConfig(dbdir=DBDIR),
                  novel=[big]),
        report.HtmlOptions(max_novel_bp=100))
    assert "not shown (longer than" in html
    assert "A" * 200 not in html


def test_hint_when_novel_was_not_requested():
    assert "--novel FILE.fa" in report.render_html(one_result())
    # ...but not when --novel WAS given and simply found nothing
    cfg = RunConfig(full=True, dbdir=DBDIR, novel_path="n.fa")
    assert "--novel FILE.fa" not in report.render_html(mk_result([rich_sample()], cfg))


# ---------------------------------------------------------------------------
# real BLAST evidence (tests/golden/example.blast.tsv, 1423 real HSPs)
# ---------------------------------------------------------------------------

def hits_from_golden_blast(locus: str):
    """Build Hit objects from the real blastn output of tests/data/example.fna."""
    hits = []
    path = os.path.join(HERE, "golden", "example.blast.tsv")
    with open(path, encoding="utf-8") as fh:
        for i, line in enumerate(fh, start=1):
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != 9:
                continue
            sseqid, slen, length, nident, qseqid, qstart, qend, qseq, sstrand = cols
            scheme, _, rest = sseqid.partition(".")
            gene, _, allele = rest.rpartition("_")
            if scheme != "sepidermidis" or gene != locus:
                continue
            slen, length, nident = int(slen), int(length), int(nident)
            qstart, qend = int(qstart), int(qend)
            kind = ("exact" if slen == length and nident == slen else
                    "novel" if slen == length else "partial")
            hits.append(Hit(sseqid=sseqid, scheme=scheme, locus=gene, allele=allele,
                            slen=slen, length=length, nident=nident, qseqid=qseqid,
                            qstart=qstart, qend=qend, qseq=qseq, sstrand=sstrand,
                            pct_identity=100.0 * nident / length,
                            pct_coverage=100.0 * nident / slen,
                            call_kind=kind, index=i))
    return hits


def test_real_blast_evidence_renders_and_the_arithmetic_adds_up():
    import math
    hits = hits_from_golden_blast("arcC")
    assert len(hits) > 5, len(hits)
    exact = [hh for hh in hits if hh.call_kind == "exact"]
    assert exact, "no exact arcC hit in the golden BLAST output"
    best = exact[0]
    call = AlleleCall("arcC", best.allele, "exact", best, tuple(hits[:6]))
    sample = SampleResult(path="example.fna", label="example.fna",
                          scheme="sepidermidis", st="184", signature=best.allele,
                          score=100, status="PERFECT", alleles=(call,),
                          n_contigs=99, total_bp=2_499_279, hits_seen=1423,
                          hits_kept=len(hits))
    cfg = RunConfig(full=True, dbdir=DBDIR)
    html = report.render_html(mk_result([sample], cfg),
                              report.HtmlOptions(evidence="all"))
    assert walk(html).errors == []
    for hh in hits[:6]:
        assert report.h(hh.qseqid) in html
        assert "{0} / {1}".format(hh.length, hh.slen) in html
        if hh.call_kind == "exact":
            expect_id = expect_cov = "100.0"
        else:
            expect_id = "{0:.1f}".format(math.floor(hh.pct_identity * 10) / 10)
            expect_cov = "{0:.1f}".format(math.floor(hh.pct_coverage * 10) / 10)
        assert "{0} %".format(expect_id) in html, hh.sseqid
        assert "{0} %".format(expect_cov) in html, hh.sseqid
    # the minus-strand hits in this file must be labelled as reverse
    if any(hh.sstrand == "minus" for hh in hits[:6]):
        assert "reverse" in html


# ---------------------------------------------------------------------------
# 11.5  CSS contract
# ---------------------------------------------------------------------------

def _lin(channel: float) -> float:
    c = channel / 255.0
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def luminance(hexcolor: str) -> float:
    hexcolor = hexcolor.lstrip("#")
    r, g, b = (int(hexcolor[i:i + 2], 16) for i in (0, 2, 4))
    return 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)


def contrast(fg: str, bg: str) -> float:
    a, b = luminance(fg), luminance(bg)
    lo, hi = min(a, b), max(a, b)
    return (hi + 0.05) / (lo + 0.05)


def test_every_declared_contrast_pair_meets_wcag_aa():
    for theme, tokens in (("light", report.TOKENS_LIGHT), ("dark", report.TOKENS_DARK)):
        for fg, bg in report.CONTRAST_PAIRS:
            ratio = contrast(tokens[fg], tokens[bg])
            assert ratio >= 4.5, "{0}: {1} on {2} = {3:.2f}".format(theme, fg, bg, ratio)


def test_both_themes_declare_the_same_tokens():
    assert set(report.TOKENS_LIGHT) == set(report.TOKENS_DARK)


def test_no_colour_is_defined_only_inside_a_media_or_theme_block():
    root = report.CSS_TEXT.split("}", 1)[0]
    for name in report.TOKENS_LIGHT:
        assert "--{0}:".format(name) in root, name


def test_theme_override_wins_in_both_directions():
    for needle in ('@media (prefers-color-scheme: dark)',
                   ':root:not([data-theme="light"])',
                   ':root[data-theme="dark"]',
                   ':root[data-theme="light"]'):
        assert needle in report.CSS_TEXT, needle


def test_required_css_blocks_exist():
    for needle in ("@page { size: A4 portrait", "@media print",
                   "@media (prefers-reduced-motion: reduce)",
                   "@media (prefers-contrast: more)",
                   "@media (forced-colors: active)",
                   "@media (max-width: 760px)", "@media (max-width: 420px)",
                   "break-inside: avoid", "display: table-header-group",
                   "repeating-linear-gradient", ".table-wrap { overflow-x: auto",
                   "position: sticky"):
        assert needle in report.CSS_TEXT, needle


def test_print_expands_every_details_block():
    print_block = report.CSS_TEXT[report.CSS_TEXT.index("@media print"):]
    assert "details > *:not(summary) { display: block !important; }" in print_block
    assert ".no-print { display: none !important; }" in print_block


def test_body_has_an_explicit_token_background():
    assert "background: var(--bg);" in report.CSS_TEXT


def test_no_user_string_can_reach_a_style_attribute():
    html = report.render_html(batch_result())
    for _tag, name, value in walk(html).attrs:
        if name == "style":
            assert re.fullmatch(r"width:\d{1,3}%", value), value


# ---------------------------------------------------------------------------
# write_html
# ---------------------------------------------------------------------------

def test_write_html_is_atomic_utf8_lf():
    result = one_result()
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "report.html")
        assert report.write_html(result, path) == path
        assert not os.path.exists(path + ".tmp")
        with open(path, "rb") as fh:
            raw = fh.read()
    assert b"\r\n" not in raw
    assert not raw.startswith(b"\xef\xbb\xbf")
    assert raw.decode("utf-8").startswith("<!DOCTYPE html>")


def test_write_html_keeps_the_previous_report_when_the_rename_fails():
    # Finding 26: a failed write must not destroy the file the user already
    # had, and must not leave the staging file behind.
    result = one_result()
    real_replace = os.replace

    def boom(src, dst):
        raise OSError(27, "File too large")

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "report.html")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("PREVIOUS-HTML")
        report.os.replace = boom
        try:
            report.write_html(result, path)
        except OSError:
            pass
        else:
            raise AssertionError("write_html swallowed the failure")
        finally:
            report.os.replace = real_replace
        assert open(path, encoding="utf-8").read() == "PREVIOUS-HTML"
        assert not os.path.exists(path + ".tmp"), "stray .tmp left behind"


def test_evidence_mode_none_drops_the_evidence_columns():
    html = report.render_html(one_result(), report.HtmlOptions(evidence="none"))
    assert "contig_00001" not in html
    assert ">Locus</th>" in html


def test_evidence_mode_all_shows_every_hit():
    html = report.render_html(one_result(), report.HtmlOptions(evidence="all"))
    assert "contig_00002" in html


def test_run_meta_always_prints_the_settings():
    html = report.render_html(one_result())
    for needle in ("Minimum identity", "Minimum coverage", "Minimum score",
                   "Excluded schemes", "BLAST threads", "Scheme", "Command line"):
        assert needle in html, needle
    assert "auto-detect" in html


def test_forced_scheme_shows_the_effective_minscore():
    cfg = RunConfig(dbdir=DBDIR, scheme="saureus", minscore=0.0, exclude=frozenset())
    html = report.render_html(mk_result([mk_sample("x.fa", "saureus", "-", "NONE", 0)], cfg))
    assert "0 (forced by --scheme)" in html


def test_repair_locus_ids_is_labelled_loudly():
    cfg = RunConfig(dbdir=DBDIR, repair_locus_ids=True)
    html = report.render_html(mk_result([mk_sample("x.fa", "saureus", "1", "PERFECT", 100,
                                                   "arcC(1)")], cfg))
    assert "Results will NOT match tseemann/mlst" in html


def test_title_option_is_honoured_and_escaped():
    html = report.render_html(one_result(), report.HtmlOptions(title="Ward <6> run"))
    assert "<title>Ward &lt;6&gt; run</title>" in html


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
