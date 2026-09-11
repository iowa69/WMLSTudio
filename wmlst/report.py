# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""WMLST rendering — every byte the user reads.

Pure rendering. No biology, no recomputation: every number in this module comes
straight off the frozen dataclasses in :mod:`wmlst.engine`.

Implements docs/ARCHITECTURE.md section 4.7 (public API), section 11 (the
self-contained HTML report) and section 12 (the byte-identity text surfaces).

Ported from bin/mlst:96-119 (--list/--longlist/--info), bin/mlst:167-250
(the row, JSON and novel-FASTA writers) and bin/mlst:417-428 (print_row).

Two families of output live here and they must never contaminate each other:

* the **compat emitters** (:func:`print_row`, :func:`write_tsv`,
  :func:`write_json`, :func:`write_novel_fasta`, :func:`write_list`,
  :func:`write_longlist`, :func:`write_info`) reproduce ``tseemann/mlst``
  2.35.0 byte for byte, quirks included;
* the **WMLST-only surfaces** (:func:`write_evidence_tsv`,
  :func:`render_html`, :func:`write_html`) are additive supersets.
"""

from __future__ import annotations

import base64
import io
import json
import math
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from . import branding as _branding
from .engine import RunResult, out_sep
from .version import __version__

__all__ = [
    "CONTRAST_PAIRS",
    # constants
    "CSS_TEXT",
    "JS_TEXT",
    "LOGO_SVG_DATA_URI",
    "STATUS_INFO",
    "SYMBOL_INFO",
    "TOKENS_DARK",
    "TOKENS_LIGHT",
    # section 11 -- the HTML report
    "HtmlOptions",
    "attr",
    # the only sanctioned interpolators
    "h",
    "js_json",
    # section 12 -- byte-identity surfaces
    "print_row",
    "render_html",
    # helpers the GUI shares
    "species_label_for",
    "tsv_text",
    "write_evidence_tsv",
    "write_html",
    "write_info",
    "write_json",
    "write_list",
    "write_longlist",
    "write_novel_fasta",
    "write_tsv",
]


# =============================================================================
# Section 12 -- the byte-identity text surfaces
# =============================================================================

def print_row(cols, sep: str, out) -> None:
    """Write one row, reproducing ``bin/mlst:417-428`` exactly (section 12.1).

    Quote a field iff it contains a double quote OR the *current* separator;
    then double the quotes, THEN double the backslashes, then wrap. The
    backslash doubling is upstream's non-standard escape and happens only
    inside quoted fields and only after quote doubling, so ``a"b\\c`` becomes
    ``"a""b\\\\c"``. A comma is not quoted in TSV mode and a tab is not quoted
    in CSV mode, because the trigger is the current ``$OUTSEP``. The ``csv``
    module cannot reproduce this and must not be used.
    """
    fields = []
    for col in cols:
        c = "" if col is None else str(col)
        if '"' in c or sep in c:
            c = c.replace('"', '""')
            c = c.replace("\\", "\\\\")
            c = '"' + c + '"'
        fields.append(c)
    out.write(sep.join(fields) + "\n")


def _profile_header_genes(datadir: str, scheme: str) -> Tuple[str, ...]:
    """Gene names from ``<datadir>/<scheme>/<scheme>.txt``'s first line.

    Last-resort fallback for the ``--legacy`` header when the run produced no
    usable sample to take the gene order from (section 12.2). Mirrors the
    ``MLST::Scheme`` filter regex; it reads a header line, never a profile.
    """
    path = os.path.join(datadir, scheme, scheme + ".txt")
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
    except OSError:
        return ()
    cols = first.rstrip("\r\n").split("\t")
    drop = re.compile(r"^(ST|mlst_clade|clonal_complex|species|CC|Lineage)$", re.ASCII)
    return tuple(c for c in cols if not drop.match(c))


def _legacy_genes(result: RunResult, catalog=None) -> Tuple[str, ...]:
    """Locus names for the ``--legacy`` header row, in profile-header order."""
    for sample in result.samples:
        if not sample.failed and sample.alleles:
            return tuple(call.locus for call in sample.alleles)
    scheme = result.meta.config.scheme
    if not scheme:
        return ()
    if catalog is not None:
        try:
            return tuple(catalog[scheme].genes)
        except Exception:
            pass
    try:  # the catalog is owner A's; degrade to reading the header directly
        from .schemes import SchemeCatalog

        return tuple(SchemeCatalog(result.meta.config.datadir)[scheme].genes)
    except Exception:
        return _profile_header_genes(result.meta.config.datadir, scheme)


def _compat_samples(result: RunResult):
    """The samples that reach a compat output.

    ``failed=True`` rows never appear in TSV/CSV/JSON: upstream dies on a bad
    file and emits no row for it, so neither do we (section 3.5, D10).
    """
    return [s for s in result.samples if not s.failed]


def write_tsv(result: RunResult, out, catalog=None) -> None:
    """Write the default / ``--full`` / ``--legacy`` layout (sections 12.2, 12.3).

    The separator comes from ``cfg.csv``. ``--legacy`` wins over ``--full``
    because upstream's header chain is an ``elsif`` and the row printer tests
    legacy first. Default mode emits no header at all.
    """
    cfg = result.meta.config
    sep = out_sep(cfg)
    legacy = bool(cfg.scheme) and cfg.legacy

    if legacy:
        print_row(("FILE", "SCHEME", "ST", *_legacy_genes(result, catalog)), sep, out)
    elif cfg.full:
        print_row(("FILE", "SCHEME", "ST", "STATUS", "SCORE", "ALLELES"), sep, out)

    for sample in _compat_samples(result):
        if legacy:
            print_row(
                [sample.label, sample.scheme, sample.st]
                + [call.code for call in sample.alleles],
                sep, out,
            )
        elif cfg.full:
            # scheme == "-" yields an EMPTY ALLELES field and therefore a
            # trailing separator. Do not substitute "-".
            print_row(
                [
                    sample.label, sample.scheme, sample.st, sample.status,
                    str(int(sample.score)),
                    ";".join("{0}({1})".format(c.locus, c.code) for c in sample.alleles),
                ],
                sep, out,
            )
        else:
            print_row(
                [sample.label, sample.scheme, sample.st]
                + ["{0}({1})".format(c.locus, c.code) for c in sample.alleles],
                sep, out,
            )


def tsv_text(result: RunResult, catalog=None) -> str:
    """The exact bytes :func:`write_tsv` would emit, as a ``str``."""
    buf = io.StringIO()
    write_tsv(result, buf, catalog)
    return buf.getvalue()


def _json_records(result: RunResult) -> List[Dict[str, Any]]:
    """The ``--json`` payload: five keys per file, every value a string (12.4)."""
    records = []
    for sample in _compat_samples(result):
        alleles = None
        if sample.scheme != "-":
            alleles = {}
            for call in sample.alleles:
                alleles[call.locus] = call.code
        records.append({
            "id": sample.label,
            "filename": sample.path,
            "scheme": sample.scheme,
            "sequence_type": sample.st,
            "alleles": alleles,
        })
    return records


def json_text(result: RunResult) -> str:
    """The ``--json`` document as a ``str`` (section 12.4).

    Ends with a single LF: Perl's ``to_json(..., {pretty=>1})`` emits one, and
    the golden ``example.json`` confirms it. Spec C5 claimed otherwise on the
    grounds of byte-identity, which key order (D2) makes unreachable anyway.
    """
    return json.dumps(
        _json_records(result), indent=3, separators=(",", " : "), ensure_ascii=True
    ) + "\n"


def write_json(result: RunResult, path: str) -> None:
    """Write the ``--json`` file (section 12.4).

    ``indent=3`` with a ``" : "`` key separator, ``ensure_ascii=True``, UTF-8,
    LF, no BOM and a single trailing LF, as Perl emits. Key order is fixed by WMLST
    to id, filename, scheme, sequence_type, alleles (D2); upstream's is a
    randomised Perl hash. ``alleles`` is JSON ``null`` when no scheme matched.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json_text(result))


def write_novel_fasta(result: RunResult, path: str) -> None:
    """Write the ``--novel`` FASTA (section 12.6).

    ``">{fasta_id} {source_label}\\n{seq}\\n"`` per allele, one unwrapped
    sequence line, description is the LABEL and not the path, deduplicated by
    the full id with first occurrence winning. An empty file is still created
    when there are no novel alleles.
    """
    seen = set()
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for allele in result.novel:
            if allele.fasta_id in seen:
                continue
            seen.add(allele.fasta_id)
            fh.write(">{0} {1}\n{2}\n".format(
                allele.fasta_id, allele.source_label, allele.seq))


EVIDENCE_COLUMNS = (
    "LABEL", "SCHEME", "LOCUS", "ALLELE", "CALL_KIND", "USED", "CONTIG",
    "QSTART", "QEND", "STRAND", "NIDENT", "ALIGN_LEN", "ALLELE_LEN",
    "PCT_IDENTITY", "PCT_COVERAGE",
)


def write_evidence_tsv(result: RunResult, path: str) -> None:
    """Write the WMLST-only long-format hit dump (section 12.7).

    One row per surviving BLAST hit of every locus of the reported scheme.
    Unlike the compat surfaces this includes ``failed=True`` samples, which
    contribute a single row recording the failure. Never affects a compat
    output.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        print_row(EVIDENCE_COLUMNS, "\t", fh)
        for sample in result.samples:
            if sample.failed:
                print_row(
                    [sample.label, "-", "-", "-", "failed", "-",
                     sample.error_text or "failed", "", "", "", "", "", "", "", ""],
                    "\t", fh,
                )
                continue
            for call in sample.alleles:
                if not call.hits:
                    print_row(
                        [sample.label, sample.scheme, call.locus, call.code,
                         call.symbol, "-", "", "", "", "", "", "", "", "", ""],
                        "\t", fh,
                    )
                    continue
                for hit in call.hits:
                    print_row(
                        [
                            sample.label, sample.scheme, call.locus, hit.allele,
                            hit.call_kind,
                            "yes" if _is_best(call, hit) else "no",
                            hit.qseqid, hit.qstart, hit.qend, hit.sstrand,
                            hit.nident, hit.length, hit.slen,
                            _pct_text(hit.pct_identity, hit.call_kind == "exact"),
                            _pct_text(hit.pct_coverage, hit.call_kind == "exact"),
                        ],
                        "\t", fh,
                    )


def write_list(catalog, out) -> None:
    """``--list``: space-joined scheme names plus a newline (section 12.5).

    Goes to stdout directly, not through ``print_row`` and not to ``$OUTFH``.
    """
    out.write(" ".join(catalog.names()) + "\n")


def write_longlist(catalog, out, sep: str = "\t") -> None:
    """``--longlist``: one ``print_row(name, *genes)`` per scheme (section 12.5)."""
    for name in catalog.names():
        print_row([name, *list(catalog[name].genes)], sep, out)


def write_info(catalog, out, sep: str = "\t") -> None:
    """``--info``: six columns, upstream's ``LOCII`` misspelling kept (12.5).

    The ``,,`` at ``bin/mlst:112`` is a Perl no-op, so the header and the rows
    both carry exactly six fields. ``LOCII_NAMES`` joins gene names with a
    single space; ``DATE`` is ``Unknown`` for all 162 shipped schemes.
    """
    print_row(("SCHEME", "LOCII", "TYPES", "ALLELES", "DATE", "LOCII_NAMES"), sep, out)
    for name in sorted(catalog.names()):
        scheme = catalog[name]
        print_row(
            [name, scheme.num_genes, scheme.num_genotypes, scheme.num_alleles,
             scheme.last_updated, " ".join(scheme.genes)],
            sep, out,
        )


# =============================================================================
# Section 11.1 -- the only sanctioned interpolators
# =============================================================================

#: C0 and C1 controls except tab and newline, plus the two line separators that
#: terminate a JavaScript string literal.
_CONTROL_RE = re.compile(
    "[\x00-\x08\x0b-\x1f\x7f-\x9f\u2028\u2029]"
)

_ESCAPE_MAP = {
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#x27;",
    "`": "&#x60;",
    "=": "&#x3d;",
}


def h(value) -> str:
    """HTML-escape for element text AND quoted attribute values (section 11.1).

    Escapes ``& < > " ' ` =`` and replaces C0/C1 controls (except tab and
    newline) plus U+2028/U+2029 with U+FFFD. ``html.escape()`` is NOT
    sufficient — it leaves the backtick and the equals sign, which are enough
    to break out of an unquoted attribute in several parsers — and must not be
    used directly. Every value that originated in a FASTA header, a filename,
    a label or a scheme name goes through this function.
    """
    text = "" if value is None else str(value)
    text = _CONTROL_RE.sub("\ufffd", text)
    return "".join(_ESCAPE_MAP.get(ch, ch) for ch in text)


def attr(value) -> str:
    """``'"' + h(value) + '"'``. Attribute values are ALWAYS double-quoted."""
    return '"' + h(value) + '"'


def js_json(obj) -> str:
    """Serialise ``obj`` for embedding inside a ``<script>`` element (4.7).

    ``<``, ``>``, ``&``, U+2028 and U+2029 are ``\\u``-escaped, so the byte
    sequence ``</script`` can never appear inside the embedded JSON island.
    """
    text = json.dumps(obj, ensure_ascii=True, separators=(",", ":"))
    return (text.replace("<", "\\u003c").replace(">", "\\u003e")
                .replace("&", "\\u0026")
                .replace("\u2028", "\\u2028").replace("\u2029", "\\u2029"))


# =============================================================================
# Section 11.3 / 11.4 -- the vocabulary the report renders
# =============================================================================

#: status -> (glyph, css class, left-border style, plain-English explanation).
#: The seven values are the complete vocabulary (section 5.15 step 52).
STATUS_INFO: Dict[str, Tuple[str, str, str, str]] = {
    "PERFECT": (
        "\u2714", "s-perfect", "solid",
        "Every locus matched a known allele exactly, and the combination is a "
        "recognised sequence type.",
    ),
    "NOVEL": (
        "\u271a", "s-novel", "solid",
        "Every locus matched a known allele exactly, but this combination is "
        "not yet a named sequence type \u2014 it may be a new ST worth "
        "submitting to PubMLST.",
    ),
    "MIXED": (
        "\u29c9", "s-mixed", "double",
        "At least one locus matched two or more different alleles equally "
        "well. This usually means the assembly contains more than one strain, "
        "or a duplicated gene \u2014 check culture purity before reporting.",
    ),
    "MISSING": (
        "\u25cc", "s-missing", "dashed",
        "At least one locus could not be found in this assembly. It may be "
        "incomplete, or the locus may genuinely be absent \u2014 an ST cannot "
        "be assigned.",
    ),
    "BAD": (
        "\u26a0", "s-bad", "thick",
        "The best-matching scheme scored below 70 out of 100. Treat this "
        "result as unreliable: wrong organism, heavily fragmented, or "
        "contaminated.",
    ),
    "OK": (
        "\u25cf", "s-ok", "solid",
        "The scheme matched reasonably well, but at least one locus is only "
        "an approximate or partial match, so no exact sequence type could be "
        "assigned.",
    ),
    "NONE": (
        "\u25cb", "s-none", "dashed",
        "No MLST scheme matched this file at all. Check that it really "
        "contains assembled contigs for a species covered by one of the "
        "bundled schemes.",
    ),
}

#: An unrecognised status never reaches a colour lookup (section 11.1 rule 5).
_STATUS_UNKNOWN = (
    "?", "s-unknown", "dashed",
    "WMLST does not recognise this status value. The result is shown "
    "unaltered; treat it as unverified.",
)

#: The legend, in the order section 11.4 requires it to always be rendered.
SYMBOL_INFO: Tuple[Tuple[str, str, str], ...] = (
    ("16", "exact",
     "Exact match: the locus was found at full length and is identical to "
     "known allele 16."),
    ("~16", "novel",
     "Inexact match: the locus was found at full length and allele 16 is the "
     "closest known one, but the sequence is not identical \u2014 possibly a "
     "new allele."),
    ("16?", "partial",
     "Partial match: only part of the locus was found, usually because the "
     "assembly breaks inside the gene. The allele number is a best guess."),
    ("-", "missing",
     "Missing: the locus was not found in this assembly at all."),
    ("1,2", "multiple",
     "Multiple matches: two or more different alleles matched equally well. "
     "Suspect a mixed culture or a duplicated gene."),
    ("0", "null",
     "Null: the locus was not found, but the sequence type that was matched "
     "records no allele at this position, so the gap is expected by the "
     "scheme rather than a failure of the assembly."),
)


#: Not an mlst status: WMLST's own marker for a file that could not be read
#: (D10). It never reaches a compat output.
_FAILED_INFO = (
    "\u26a0", "s-bad", "thick",
    "This file could not be analysed, so no typing was attempted. The rest of "
    "the run continued without it.",
)


def _status_info(status: str) -> Tuple[str, str, str, str]:
    if status == "FAILED":
        return _FAILED_INFO
    return STATUS_INFO.get(status, _STATUS_UNKNOWN)


def _is_best(call, hit) -> bool:
    """Whether ``hit`` is the hit that set the call's code."""
    if call.best is None:
        return False
    return hit is call.best or hit == call.best


# =============================================================================
# Section 11.5 -- the design tokens (single source of truth for CSS + tests)
# =============================================================================

TOKENS_LIGHT: Dict[str, str] = {
    "bg": "#ffffff",
    "surface": "#f4f6f9",
    "surface-2": "#e7ecf2",
    "text": "#11161c",
    "muted": "#4d5763",
    "border": "#c2ccd6",
    "border-strong": "#7b8794",
    "accent": "#0a4c99",
    "accent-fg": "#ffffff",
    "link": "#0a4c99",
    "meter-track": "#dde3ea",
    "meter-high": "#1f7a45",
    "meter-mid": "#8a5a06",
    "meter-low": "#a12222",
    "warn-bg": "#fdf3df",
    "warn-fg": "#5e3705",
    "warn-line": "#a97a14",
    "code-bg": "#eef1f5",
    "s-perfect-bg": "#dcf1e3", "s-perfect-fg": "#0e4a29", "s-perfect-line": "#1f7a45",
    "s-novel-bg": "#dfe9fb", "s-novel-fg": "#0d3b7a", "s-novel-line": "#2a5fae",
    "s-mixed-bg": "#fbeddb", "s-mixed-fg": "#5e3705", "s-mixed-line": "#a3701a",
    "s-missing-bg": "#fbe6de", "s-missing-fg": "#6d2a10", "s-missing-line": "#b4562a",
    "s-bad-bg": "#fbe1e1", "s-bad-fg": "#7a1414", "s-bad-line": "#b02020",
    "s-ok-bg": "#e3ecf2", "s-ok-fg": "#1c4054", "s-ok-line": "#3a6f8c",
    "s-none-bg": "#e9ebee", "s-none-fg": "#31383f", "s-none-line": "#6b747d",
    "s-unknown-bg": "#e9ebee", "s-unknown-fg": "#31383f", "s-unknown-line": "#6b747d",
}

TOKENS_DARK: Dict[str, str] = {
    "bg": "#10141a",
    "surface": "#171d25",
    "surface-2": "#202834",
    "text": "#e9eef5",
    "muted": "#a8b4c2",
    "border": "#323c49",
    "border-strong": "#54616f",
    "accent": "#79b3f5",
    "accent-fg": "#0b1017",
    "link": "#8cc0ff",
    "meter-track": "#2a3341",
    "meter-high": "#5fd08a",
    "meter-mid": "#e8b552",
    "meter-low": "#ff8d8d",
    "warn-bg": "#2e2410",
    "warn-fg": "#f3d69a",
    "warn-line": "#a97a14",
    "code-bg": "#1d242e",
    "s-perfect-bg": "#102c1d", "s-perfect-fg": "#83e0a8", "s-perfect-line": "#3fae69",
    "s-novel-bg": "#121f38", "s-novel-fg": "#a6c8ff", "s-novel-line": "#4a7fd0",
    "s-mixed-bg": "#2e2410", "s-mixed-fg": "#f0c979", "s-mixed-line": "#b4881f",
    "s-missing-bg": "#301c12", "s-missing-fg": "#f5b394", "s-missing-line": "#c06a3c",
    "s-bad-bg": "#331414", "s-bad-fg": "#ff9f9f", "s-bad-line": "#c94040",
    "s-ok-bg": "#152531", "s-ok-fg": "#a3cce4", "s-ok-line": "#4d87a8",
    "s-none-bg": "#1c222a", "s-none-fg": "#c2cad3", "s-none-line": "#78828d",
    "s-unknown-bg": "#1c222a", "s-unknown-fg": "#c2cad3", "s-unknown-line": "#78828d",
}

#: (foreground token, background token) pairs whose WCAG contrast is a checked-in
#: contract. tests/test_report_html.py recomputes every one of these in both
#: themes and fails the build on anything below 4.5:1.
CONTRAST_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("text", "bg"), ("text", "surface"), ("text", "surface-2"),
    ("muted", "bg"), ("muted", "surface"), ("muted", "surface-2"),
    ("link", "bg"), ("link", "surface"),
    ("accent-fg", "accent"),
    ("warn-fg", "warn-bg"),
    ("text", "code-bg"),
    ("s-perfect-fg", "s-perfect-bg"),
    ("s-novel-fg", "s-novel-bg"),
    ("s-mixed-fg", "s-mixed-bg"),
    ("s-missing-fg", "s-missing-bg"),
    ("s-bad-fg", "s-bad-bg"),
    ("s-ok-fg", "s-ok-bg"),
    ("s-none-fg", "s-none-bg"),
    ("s-unknown-fg", "s-unknown-bg"),
)


def _token_block(tokens: Dict[str, str], indent: str = "  ") -> str:
    return "\n".join(
        "{0}--{1}: {2};".format(indent, name, tokens[name]) for name in sorted(tokens)
    )


_CSS_BODY = r"""
/* ---- reset / base ------------------------------------------------------- */
*, *::before, *::after { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0;
  padding: 0;
  background: var(--bg);
  color: var(--text);
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue",
               Arial, "Noto Sans", sans-serif;
  font-size: 16px;
  line-height: 1.5;
}
.mono, code, pre, .code {
  font-family: ui-monospace, "Cascadia Mono", "Consolas", "SFMono-Regular",
               "Liberation Mono", "Menlo", monospace;
}
h1, h2, h3, h4 { line-height: 1.2; margin: 0 0 .4em; font-weight: 650; }
h1 { font-size: 1.5rem; }
h2 { font-size: 1.2rem; }
h3 { font-size: 1.05rem; }
p { margin: 0 0 .75em; }
a { color: var(--link); }
.wrap { max-width: 1180px; margin: 0 auto; padding: 0 20px; }
.muted { color: var(--muted); }
.small { font-size: .82rem; }
.visually-hidden {
  position: absolute !important; width: 1px; height: 1px; padding: 0;
  margin: -1px; overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap;
  border: 0;
}
.skip {
  position: absolute; left: -9999px; top: 0; background: var(--accent);
  color: var(--accent-fg); padding: 10px 14px; z-index: 20;
}
.skip:focus { left: 0; }
:focus-visible { outline: 3px solid var(--accent); outline-offset: 2px; }

/* ---- masthead ----------------------------------------------------------- */
.masthead {
  background: var(--surface);
  border-bottom: 3px solid var(--accent);
  padding: 16px 0;
}
.masthead .wrap {
  display: flex; flex-wrap: wrap; gap: 14px; align-items: center;
}
.masthead img { width: 44px; height: 44px; flex: none; }
.brand { flex: 1 1 320px; min-width: 0; }
.brand .product { font-size: 1.35rem; font-weight: 700; letter-spacing: .01em; }
.brand .vendor { font-size: .9rem; color: var(--muted); }
.masthead .tools { display: flex; gap: 8px; flex-wrap: wrap; }
button {
  font: inherit; font-size: .85rem; cursor: pointer; color: var(--text);
  background: var(--bg); border: 1px solid var(--border-strong);
  border-radius: 6px; padding: 6px 11px;
}
button:hover { background: var(--surface-2); }

/* ---- run metadata ------------------------------------------------------- */
.run-meta { margin: 22px 0; }
.run-meta dl {
  display: grid; grid-template-columns: repeat(auto-fill, minmax(230px, 1fr));
  gap: 2px 18px; margin: 0;
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 14px 16px;
}
.run-meta dt {
  font-size: .74rem; text-transform: uppercase; letter-spacing: .06em;
  color: var(--muted); margin-top: 8px;
}
.run-meta dd { margin: 0 0 2px; font-size: .9rem; word-break: break-word; }
.dagger { color: var(--warn-line); font-weight: 700; }
.notice {
  background: var(--warn-bg); color: var(--warn-fg);
  border: 1px solid var(--warn-line); border-left-width: 6px;
  border-radius: 6px; padding: 10px 14px; margin: 14px 0;
}

/* ---- cards -------------------------------------------------------------- */
.card {
  background: var(--surface); border: 1px solid var(--border);
  border-left: 8px solid var(--border-strong);
  border-radius: 10px; padding: 18px 20px; margin: 18px 0;
}
.card.b-solid  { border-left-style: solid; }
.card.b-double { border-left-style: double; border-left-width: 10px; }
.card.b-dashed { border-left-style: dashed; }
.card.b-thick  { border-left-style: solid; border-left-width: 14px; }
.card.s-perfect { border-left-color: var(--s-perfect-line); }
.card.s-novel   { border-left-color: var(--s-novel-line); }
.card.s-mixed   { border-left-color: var(--s-mixed-line); }
.card.s-missing { border-left-color: var(--s-missing-line); }
.card.s-bad     { border-left-color: var(--s-bad-line); }
.card.s-ok      { border-left-color: var(--s-ok-line); }
.card.s-none    { border-left-color: var(--s-none-line); }
.card.s-unknown { border-left-color: var(--s-unknown-line); }
.card-head { display: flex; flex-wrap: wrap; gap: 14px; align-items: flex-start; }
.card-head .who { flex: 1 1 320px; min-width: 0; }
.filename { font-size: .85rem; color: var(--muted); word-break: break-all; }
.headline { display: flex; flex-wrap: wrap; align-items: baseline; gap: 10px; margin: 4px 0; }
.headline .scheme { font-size: 1.5rem; font-weight: 700; word-break: break-word; }
.headline .st { font-size: 2.3rem; font-weight: 800; letter-spacing: -.01em; }
.headline .st-label {
  font-size: .74rem; text-transform: uppercase; letter-spacing: .07em;
  color: var(--muted);
}
.species { font-size: .95rem; }
.stats { font-size: .82rem; color: var(--muted); margin-top: 6px; }
.stats span + span::before { content: " \00b7 "; }

/* ---- badges ------------------------------------------------------------- */
.badge {
  display: inline-flex; align-items: center; gap: 7px;
  border: 2px solid currentColor; border-radius: 999px;
  padding: 5px 13px; font-weight: 700; font-size: .84rem;
  letter-spacing: .07em; white-space: nowrap;
}
.badge .glyph { font-size: 1rem; line-height: 1; }
.s-perfect .badge, .badge.s-perfect { background: var(--s-perfect-bg); color: var(--s-perfect-fg); }
.s-novel   .badge, .badge.s-novel   { background: var(--s-novel-bg);   color: var(--s-novel-fg); }
.s-mixed   .badge, .badge.s-mixed   { background: var(--s-mixed-bg);   color: var(--s-mixed-fg); }
.s-missing .badge, .badge.s-missing { background: var(--s-missing-bg); color: var(--s-missing-fg); }
.s-bad     .badge, .badge.s-bad     { background: var(--s-bad-bg);     color: var(--s-bad-fg); }
.s-ok      .badge, .badge.s-ok      { background: var(--s-ok-bg);      color: var(--s-ok-fg); }
.s-none    .badge, .badge.s-none    { background: var(--s-none-bg);    color: var(--s-none-fg); }
.s-unknown .badge, .badge.s-unknown { background: var(--s-unknown-bg); color: var(--s-unknown-fg); }
.status-box { flex: 0 1 380px; }
.status-why { font-size: .86rem; margin: 8px 0 0; }

/* ---- score meter -------------------------------------------------------- */
.meter { margin: 14px 0 6px; }
.meter-label { display: flex; justify-content: space-between; font-size: .8rem; color: var(--muted); }
.meter-track {
  height: 16px; border-radius: 8px; background: var(--meter-track);
  border: 1px solid var(--border-strong); overflow: hidden;
}
.meter-fill { height: 100%; }
.meter-fill.high { background: var(--meter-high); }
/* the hatch keeps the low/mid bands distinguishable in greyscale print */
.meter-fill.mid {
  background: repeating-linear-gradient(45deg,
    var(--meter-mid) 0 6px, var(--meter-track) 6px 10px);
}
.meter-fill.low {
  background: repeating-linear-gradient(-45deg,
    var(--meter-low) 0 4px, var(--meter-track) 4px 8px);
}

/* ---- tables ------------------------------------------------------------- */
.table-wrap { overflow-x: auto; margin: 12px 0; border: 1px solid var(--border); border-radius: 8px; }
table { border-collapse: collapse; width: 100%; font-size: .87rem; background: var(--bg); }
caption { text-align: left; padding: 10px 12px; font-weight: 650; background: var(--surface-2); }
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--border); vertical-align: top; }
thead th {
  position: sticky; top: 0; background: var(--surface-2); z-index: 1;
  border-bottom: 2px solid var(--border-strong); white-space: nowrap;
}
tbody th { font-weight: 650; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
tbody tr:last-child td, tbody tr:last-child th { border-bottom: 0; }
th[data-sort] { cursor: pointer; }
th[aria-sort="ascending"]::after { content: " \25b2"; }
th[aria-sort="descending"]::after { content: " \25bc"; }
.code { font-weight: 650; white-space: nowrap; }
.sym { display: block; font-size: .76rem; color: var(--muted); font-weight: 400; }
.contig { word-break: break-all; max-width: 22ch; }

/* ---- details / folds ---------------------------------------------------- */
details {
  border: 1px solid var(--border); border-radius: 8px;
  margin: 12px 0; background: var(--bg);
}
summary { cursor: pointer; padding: 9px 12px; font-weight: 650; font-size: .9rem; }
summary::marker { color: var(--muted); }
details > *:not(summary) { padding: 0 12px 12px; }
details .table-wrap { padding: 0; margin: 0 12px 12px; }
.pill {
  display: inline-block; font-size: .72rem; padding: 1px 8px; border-radius: 999px;
  border: 1px solid var(--border-strong); color: var(--muted);
}
tr.dropped th, tr.dropped td { color: var(--muted); }
tr.winner th, tr.winner td { font-weight: 700; }

/* ---- batch -------------------------------------------------------------- */
.chips { display: flex; flex-wrap: wrap; gap: 8px; margin: 10px 0; }
.chip {
  display: inline-flex; align-items: center; gap: 6px; border-radius: 999px;
  border: 2px solid currentColor; padding: 3px 11px; font-size: .8rem; font-weight: 650;
}
.filters { display: flex; flex-wrap: wrap; gap: 10px; align-items: end; margin: 12px 0; }
.filters label { display: block; font-size: .78rem; color: var(--muted); }
.filters input, .filters select {
  font: inherit; font-size: .88rem; padding: 6px 9px; border-radius: 6px;
  border: 1px solid var(--border-strong); background: var(--bg); color: var(--text);
  min-width: 180px;
}

/* ---- legend / novel / empty -------------------------------------------- */
.legend { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px 20px; }
.legend dl { display: grid; grid-template-columns: max-content 1fr; gap: 6px 16px; margin: 0; }
.legend dt { font-weight: 700; }
.legend dd { margin: 0; }
pre.seq {
  background: var(--code-bg); border: 1px solid var(--border); border-radius: 6px;
  padding: 10px 12px; overflow-x: auto; font-size: .8rem; line-height: 1.45; margin: 6px 0;
}
.seq-head { color: var(--muted); }
.empty {
  border: 2px dashed var(--border-strong); border-radius: 10px;
  padding: 26px; text-align: center; color: var(--muted);
}
.colophon {
  border-top: 1px solid var(--border); margin-top: 30px; padding: 18px 0 34px;
  font-size: .82rem; color: var(--muted);
}
.colophon ul { margin: .4em 0; padding-left: 1.2em; }

/* ---- responsive --------------------------------------------------------- */
@media (max-width: 760px) {
  body { font-size: 15px; }
  .wrap { padding: 0 14px; }
  .headline .st { font-size: 1.9rem; }
  .run-meta dl { grid-template-columns: 1fr 1fr; }
  .legend dl { grid-template-columns: 1fr; }
  .legend dt { margin-top: 8px; }
}
@media (max-width: 420px) {
  .wrap { padding: 0 10px; }
  .card { padding: 14px; border-left-width: 6px; }
  .run-meta dl { grid-template-columns: 1fr; }
  .headline .scheme { font-size: 1.2rem; }
  .headline .st { font-size: 1.6rem; }
  .filters input, .filters select { min-width: 0; width: 100%; }
}

/* ---- user preferences --------------------------------------------------- */
@media (prefers-reduced-motion: reduce) {
  *, *::before, *::after { animation: none !important; transition: none !important; }
}
@media (prefers-contrast: more) {
  .card, .legend, .table-wrap, details, .run-meta dl { border-color: var(--text); }
  .badge, .chip { border-width: 3px; }
  .muted, .small, .stats, .sym { color: var(--text); }
}
@media (forced-colors: active) {
  .badge, .chip, .card, .table-wrap, details, .meter-track { border-color: CanvasText; }
  .meter-fill { background: Highlight !important; forced-color-adjust: none; }
  .card { border-left-color: CanvasText; }
}

/* ---- print -------------------------------------------------------------- */
@page { size: A4 portrait; margin: 14mm 12mm 16mm; }
@media print {
  :root {
    color-scheme: light;
COLOR_TOKENS_PRINT
  }
  body { background: #ffffff; color: #11161c; font-size: 10.5pt; }
  .no-print { display: none !important; }
  .wrap { max-width: none; padding: 0; }
  .masthead { border-bottom: 2pt solid #11161c; }
  .card, tr, .legend, .notice { break-inside: avoid; page-break-inside: avoid; }
  thead { display: table-header-group; }
  tfoot { display: table-footer-group; }
  details { break-inside: auto; }
  details > *:not(summary) { display: block !important; }
  .table-wrap { overflow: visible; }
  thead th { position: static; }
  a { text-decoration: none; color: inherit; }
}
"""


def _build_css() -> str:
    light = _token_block(TOKENS_LIGHT)
    dark = _token_block(TOKENS_DARK)
    css = (
        ":root {\n  color-scheme: light dark;\n" + light + "\n}\n"
        "@media (prefers-color-scheme: dark) {\n"
        "  :root:not([data-theme=\"light\"]) {\n" + _token_block(TOKENS_DARK, "    ")
        + "\n  }\n}\n"
        ":root[data-theme=\"dark\"] {\n" + dark + "\n}\n"
        ":root[data-theme=\"light\"] {\n" + light + "\n}\n"
    )
    # The print block forces the light palette so a dark-mode screen never
    # prints white-on-black onto a laboratory record.
    return css + _CSS_BODY.replace(
        "COLOR_TOKENS_PRINT", _token_block(TOKENS_LIGHT, "    ")
    )


CSS_TEXT = _build_css()


JS_TEXT = r"""(function () {
  'use strict';
  var root = document.documentElement, KEY = 'wmlst-theme';
  function show(el) { if (el) { el.hidden = false; } }
  function each(list, fn) { Array.prototype.forEach.call(list, fn); }
  function saved() { try { return localStorage.getItem(KEY); } catch (e) { return null; } }
  function remember(v) { try { localStorage.setItem(KEY, v); } catch (e) { /* file:// */ } }
  var pref = saved();
  if (pref === 'dark' || pref === 'light') { root.setAttribute('data-theme', pref); }
  var toggle = document.getElementById('theme-toggle');
  if (toggle) {
    show(toggle);
    toggle.addEventListener('click', function () {
      var dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
      var cur = root.getAttribute('data-theme') || (dark ? 'dark' : 'light');
      var next = cur === 'dark' ? 'light' : 'dark';
      root.setAttribute('data-theme', next);
      remember(next);
      toggle.setAttribute('aria-pressed', next === 'dark' ? 'true' : 'false');
    });
  }
  var expand = document.getElementById('expand-all');
  if (expand) {
    show(expand);
    expand.addEventListener('click', function () {
      var open = expand.getAttribute('data-open') !== 'yes', all = document.querySelectorAll('details'), i;
      for (i = 0; i < all.length; i++) { all[i].open = open; }
      expand.setAttribute('data-open', open ? 'yes' : 'no');
      expand.textContent = open ? 'Collapse all sections' : 'Expand all sections';
    });
  }
  var body = document.getElementById('batch-body');
  var text = document.getElementById('filter-text');
  var status = document.getElementById('filter-status');
  var count = document.getElementById('filter-count');
  function applyFilter() {
    if (!body) { return; }
    var term = text ? text.value.trim().toLowerCase() : '', want = status ? status.value : '';
    var rows = body.rows, shown = 0, i, row, okText, okStatus;
    for (i = 0; i < rows.length; i++) {
      row = rows[i];
      okText = !term || (row.getAttribute('data-search') || '').indexOf(term) >= 0;
      okStatus = !want || row.getAttribute('data-status') === want;
      row.hidden = !(okText && okStatus);
      if (okText && okStatus) { shown++; }
    }
    if (count) { count.textContent = shown + ' of ' + rows.length + ' samples shown'; }
  }
  if (text) { show(text.parentNode); text.addEventListener('input', applyFilter); }
  if (status) { show(status.parentNode); status.addEventListener('change', applyFilter); }
  function value(row, idx, numeric) {
    var cell = row.cells[idx];
    if (!cell) { return numeric ? -1 : ''; }
    var raw = cell.getAttribute('data-value');
    if (raw === null) { raw = cell.textContent; }
    if (!numeric) { return raw.trim().toLowerCase(); }
    return parseFloat(raw.replace(/[^0-9.eE+-]/g, '')) || -1;
  }
  function byIndex(a, b) { return (+a.getAttribute('data-index')) - (+b.getAttribute('data-index')); }
  function place(rows) { for (var i = 0; i < rows.length; i++) { body.appendChild(rows[i]); } }
  var table = document.getElementById('batch-table');
  if (table && body) {
    var heads = table.querySelectorAll('th[data-sort]');
    var clearSort = function () { each(heads, function (o) { o.setAttribute('aria-sort', 'none'); }); };
    each(heads, function (th) {
      th.tabIndex = 0;
      function run() {
        var dir = th.getAttribute('aria-sort') === 'ascending' ? -1 : 1;
        clearSort();
        th.setAttribute('aria-sort', dir === 1 ? 'ascending' : 'descending');
        var idx = th.cellIndex, numeric = th.getAttribute('data-sort') === 'number';
        var rows = Array.prototype.slice.call(body.rows);
        rows.sort(function (a, b) {
          var av = value(a, idx, numeric), bv = value(b, idx, numeric);
          if (av < bv) { return -dir; }
          if (av > bv) { return dir; }
          return byIndex(a, b);
        });
        place(rows);
      }
      th.addEventListener('click', run);
      th.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); run(); }
      });
    });
    var reset = document.getElementById('reset-order');
    if (reset) {
      show(reset);
      reset.addEventListener('click', function () {
        clearSort();
        place(Array.prototype.slice.call(body.rows).sort(byIndex));
      });
    }
  }
  function legacyCopy(str) {
    var area = document.createElement('textarea'), ok = false;
    area.value = str;
    document.body.appendChild(area);
    area.select();
    try { ok = document.execCommand('copy'); } catch (e) { ok = false; }
    document.body.removeChild(area);
    return ok;
  }
  each(document.querySelectorAll('button[data-copy]'), function (btn) {
    var target = document.getElementById(btn.getAttribute('data-copy'));
    if (!target) { return; }
    show(btn);
    btn.addEventListener('click', function () {
      var out = document.getElementById('copy-status'), str = target.textContent;
      function done(ok) {
        if (!out) { return; }
        out.textContent = ok ? 'Copied to the clipboard.' : 'Copy failed — press Ctrl+C.';
      }
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(str).then(function () { done(true); },
                                               function () { done(legacyCopy(str)); });
      } else { done(legacyCopy(str)); }
    });
  });
}());"""


_LOGO_SVG = (
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" '
    'role="img" aria-label="%s %s">' % (_branding.VENDOR, _branding.APP_NAME) +
    '<rect width="64" height="64" rx="12" fill="#0a4c99"/>'
    '<path d="M20 12c0 12 24 12 24 24S20 48 20 52" fill="none" stroke="#ffffff" '
    'stroke-width="4" stroke-linecap="round"/>'
    '<path d="M44 12c0 12-24 12-24 24s24 12 24 16" fill="none" stroke="#8cc0ff" '
    'stroke-width="4" stroke-linecap="round"/>'
    '<g stroke="#ffffff" stroke-width="2.5" stroke-linecap="round">'
    '<path d="M23 20h18"/><path d="M22 28h20"/><path d="M22 36h20"/><path d="M23 44h18"/>'
    '</g></svg>'
)

#: The only permitted non-inline URI in the whole document (section 11.1).
LOGO_SVG_DATA_URI = "data:image/svg+xml;base64," + base64.b64encode(
    _LOGO_SVG.encode("utf-8")
).decode("ascii")


# =============================================================================
# Section 11 -- small rendering helpers
# =============================================================================

_DNA_RE = re.compile(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn-]*")


def _pct_text(value: float, exact: bool) -> str:
    """Format a percentage to 1 dp (section 11.4).

    Exact calls print ``100.0`` unconditionally; everything else FLOORS, so an
    inexact or partial match can never be rounded up into a misleading
    ``100.0 %`` (C14).
    """
    if exact:
        return "100.0"
    try:
        return "{0:.1f}".format(math.floor(float(value) * 10) / 10.0)
    except (TypeError, ValueError):
        return "-"


def _int_text(value) -> str:
    try:
        return "{0:,}".format(int(value))
    except (TypeError, ValueError):
        return str(value)


def _num(value) -> str:
    """Compact numeric rendering for the settings block (95.0 -> '95')."""
    try:
        return "{0:g}".format(float(value))
    except (TypeError, ValueError):
        return str(value)


_URL_RE = re.compile(r"\s*https?://\S+")
_PUBMED_RE = re.compile(r"https?://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")


def _no_urls(text: str) -> str:
    """Strip external URLs from a citation string.

    The report must contain zero references to an external origin (section
    11.1 rule 1) and CI greps the rendered file for ``http://``/``https://``,
    so a PubMed link becomes a PMID and every other URL is dropped.
    """
    text = _PUBMED_RE.sub(r"PMID \1", text)
    return _URL_RE.sub("", text).strip()


def _wrap60(seq: str) -> str:
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def _default_branding() -> Dict[str, str]:
    return {
        "app_name": _branding.APP_NAME,
        "tagline": _branding.APP_TAGLINE,
        "vendor": _branding.VENDOR,
        "author": _branding.AUTHOR,
        "attribution": _branding.ATTRIBUTION,
        "copyright": _branding.COPYRIGHT,
        "version": __version__,
        "credit_html": _branding.CREDIT_HTML,
        "citations": tuple(_no_urls(c) for c in _branding.CITATIONS),
    }


# ---- species resolution (db/scheme_species_map.tab, section 6.4) ------------

_SPECIES_CACHE: Dict[str, Dict[str, Tuple[Tuple[str, str], ...]]] = {}


def _species_map(dbdir: str) -> Dict[str, Tuple[Tuple[str, str], ...]]:
    """Parse ``<dbdir>/scheme_species_map.tab`` (section 6.4).

    Several GENUS values carry a trailing space and three schemes appear on
    more than one row, so every field is stripped and rows accumulate. A
    missing or unreadable file yields an empty map — the report says
    "not mapped" rather than guessing.
    """
    key = os.path.abspath(dbdir or "")
    cached = _SPECIES_CACHE.get(key)
    if cached is not None:
        return cached
    out: Dict[str, Tuple[Tuple[str, str], ...]] = {}
    path = os.path.join(dbdir or "", "scheme_species_map.tab")
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                name = parts[0].strip()
                genus = parts[1].strip() if len(parts) > 1 else ""
                species = parts[2].strip() if len(parts) > 2 else ""
                if not name:
                    continue
                out[name] = (*out.get(name, ()), (genus, species))
    except OSError:
        out = {}
    _SPECIES_CACHE[key] = out
    return out


def species_label_for(dbdir: str, scheme: str) -> Optional[str]:
    """Human species label for a scheme, or ``None`` when unmapped (4.5).

    0 rows -> ``None``; 1 row with a species -> ``Genus species``; 1 row
    without -> ``Genus sp.``; several rows sharing a genus -> ``Genus a/b``;
    otherwise ``Genus a / Other b``. Never fabricates.
    """
    rows = _species_map(dbdir).get(scheme, ())
    if not rows:
        return None
    labels = []
    genera = {genus for genus, _ in rows}
    if len(rows) > 1 and len(genera) == 1:
        genus = rows[0][0]
        species = [sp for _, sp in rows if sp]
        if not species:
            return "{0} sp.".format(genus)
        return "{0} {1}".format(genus, "/".join(sorted(set(species))))
    for genus, species in rows:
        labels.append("{0} {1}".format(genus, species) if species else "{0} sp.".format(genus))
    seen: List[str] = []
    for label in labels:
        if label not in seen:
            seen.append(label)
    return " / ".join(seen)


# =============================================================================
# Section 11 -- HtmlOptions and the renderer
# =============================================================================

@dataclass(frozen=True)
class HtmlOptions:
    """Knobs for :func:`render_html` (section 4.7)."""

    title: Optional[str] = None
    evidence: str = "best"          # "none" | "best" | "all"
    max_novel_bp: int = 100_000
    open_after: bool = False


_SPECIES_CAVEAT = (
    "The species shown is what this typing scheme is designed for \u2014 it is "
    "not an independent identification of your isolate."
)

_OUTBREAK_CAVEAT = (
    "Samples sharing a scheme and sequence type carry the same alleles. That "
    "is consistent with a common source, but MLST alone does not prove an "
    "outbreak link."
)

_COVERAGE_FOOTNOTE = (
    "Coverage is identical bases \u00f7 allele length \u00d7 100 \u2014 the "
    "same quantity the --mincov threshold is applied to. Identity is "
    "identical bases \u00f7 alignment length \u00d7 100. Percentages for "
    "inexact and partial calls are rounded DOWN, so an imperfect match can "
    "never be shown as 100.0 %."
)


def _evidence_rows(call, mode: str):
    """(component-code, hit-or-None) pairs for one locus row group (11.4)."""
    if mode == "none":
        return []
    if mode == "all":
        return [(None, hit) for hit in call.hits] or [(None, None)]
    if call.symbol == "multiple":
        rows = []
        for part in call.code.split(","):
            chosen = None
            for hit in call.hits:
                if hit.allele == part:
                    chosen = hit
                    break
            rows.append((part, chosen))
        return rows
    return [(None, call.best)]


def _hit_cells(hit) -> str:
    """The six evidence cells of one locus row (section 11.4 table)."""
    if hit is None:
        return (
            '<td class="muted">\u2014</td><td class="muted num">\u2014</td>'
            '<td class="muted">\u2014</td><td class="muted num">\u2014</td>'
            '<td class="muted num">\u2014</td><td class="muted num">\u2014</td>'
        )
    lo, hi = min(hit.qstart, hit.qend), max(hit.qstart, hit.qend)
    if hit.sstrand == "minus":
        strand = ('<span aria-hidden="true">\u2212</span>'
                  '<span class="visually-hidden">minus strand, </span>reverse')
    elif hit.sstrand == "plus":
        strand = ('<span aria-hidden="true">+</span>'
                  '<span class="visually-hidden">plus strand, </span>forward')
    else:
        strand = h(hit.sstrand)
    exact = hit.call_kind == "exact"
    return (
        '<td class="contig mono" translate="no" title={title}>{contig}</td>'
        '<td class="num mono">{lo}\u2013{hi}</td>'
        '<td>{strand}</td>'
        '<td class="num">{ident} %</td>'
        '<td class="num">{alen} / {slen}</td>'
        '<td class="num">{cov} %</td>'
    ).format(
        title=attr(hit.qseqid),
        contig=h(hit.qseqid),
        lo=_int_text(lo), hi=_int_text(hi),
        strand=strand,
        ident=h(_pct_text(hit.pct_identity, exact)),
        alen=_int_text(hit.length), slen=_int_text(hit.slen),
        cov=h(_pct_text(hit.pct_coverage, exact)),
    )


def _locus_table(sample, mode: str, index: int) -> str:
    out = io.StringIO()
    show_evidence = mode != "none" and any(call.hits for call in sample.alleles)
    out.write('<div class="table-wrap"><table>')
    out.write('<caption id="loci-{0}">Allele calls for every locus of {1}</caption>'.format(
        index, h(sample.scheme)))
    out.write('<thead><tr>')
    out.write('<th scope="col">Locus</th><th scope="col">Call</th>')
    if show_evidence:
        out.write(
            '<th scope="col" title="the BLAST qseqid of the matching contig">Contig</th>'
            '<th scope="col" class="num" title="qstart/qend, 1-based inclusive">Position</th>'
            '<th scope="col" title="BLAST sstrand">Strand</th>'
            '<th scope="col" class="num" title="identical bases / alignment length">Identity</th>'
            '<th scope="col" class="num" title="alignment length vs allele length">Aligned / allele</th>'
            '<th scope="col" class="num" title="the quantity the --mincov threshold is applied to">Coverage</th>'
        )
    out.write('</tr></thead><tbody>')
    for call in sample.alleles:
        rows = _evidence_rows(call, mode) if show_evidence else []
        span = max(1, len(rows))
        first = True
        if not rows:
            rows = [(None, None)]
        for part, hit in rows:
            out.write('<tr>')
            if first:
                out.write(
                    '<th scope="row" rowspan="{0}" class="mono" translate="no">{1}</th>'.format(
                        span, h(call.locus))
                )
                out.write(
                    '<td rowspan="{0}"><span class="code mono" translate="no">{1}</span>'
                    '<span class="sym">{2}</span></td>'.format(
                        span, h(call.code), h(_symbol_word(call.symbol)))
                )
                first = False
            if show_evidence:
                if part is not None and hit is None:
                    out.write('<td class="mono" colspan="6">allele {0}: no surviving hit</td>'.format(
                        h(part)))
                else:
                    out.write(_hit_cells(hit))
            out.write('</tr>')
    out.write('</tbody></table></div>')
    return out.getvalue()


_SYMBOL_WORDS = {
    "exact": "exact match",
    "novel": "inexact, full length",
    "partial": "partial match",
    "missing": "not found",
    "null": "null allele",
    "multiple": "multiple matches",
}


def _symbol_word(symbol: str) -> str:
    return _SYMBOL_WORDS.get(symbol, symbol)


def _meter(score: int) -> str:
    value = max(0, min(100, int(score)))
    band = "high" if value >= 90 else ("mid" if value >= 70 else "low")
    return (
        '<div class="meter">'
        '<div class="meter-label"><span>Match score</span>'
        '<span><strong>{0}</strong> / 100</span></div>'
        '<div class="meter-track" role="meter" aria-valuemin="0" aria-valuemax="100" '
        'aria-valuenow="{0}" aria-label="Match score {0} out of 100">'
        '<div class="meter-fill {1}" style="width:{0}%"></div></div>'
        '</div>'
    ).format(value, band)


def _badge(status: str) -> str:
    glyph, cls, _style, _why = _status_info(status)
    return (
        '<span class="badge {cls}" data-status={ds}>'
        '<span class="glyph" aria-hidden="true">{glyph}</span>{word}</span>'
    ).format(cls=cls, ds=attr(status), glyph=h(glyph), word=h(status))


def _ties(sample) -> List[str]:
    if not sample.candidates:
        return []
    top = sample.candidates[0]
    return [c.scheme for c in sample.candidates[1:] if c.score == top.score]


def _runner_up_table(sample, cfg) -> str:
    if not sample.candidates:
        return ""
    out = io.StringIO()
    out.write('<details><summary>Why this scheme won \u2014 all {0} scored candidates</summary>'.format(
        len(sample.candidates)))
    out.write('<div class="table-wrap"><table>')
    out.write('<caption>Every scheme scored for this sample, best first</caption>')
    out.write('<thead><tr><th scope="col">Scheme</th><th scope="col">ST</th>'
              '<th scope="col" class="num">Score</th><th scope="col" class="num">Loci</th>'
              '<th scope="col">Allele signature</th></tr></thead><tbody>')
    for i, cand in enumerate(sample.candidates):
        classes = []
        if i == 0:
            classes.append("winner")
        if cand.below_minscore or cand.excluded:
            classes.append("dropped")
        cls = ' class="{0}"'.format(" ".join(classes)) if classes else ""
        pills = ""
        if i == 0:
            pills += ' <span class="pill">reported</span>'
        if cand.below_minscore:
            pills += ' <span class="pill">below --minscore {0}</span>'.format(
                h(_num(cfg.minscore)))
        if cand.excluded:
            pills += ' <span class="pill">excluded</span>'
        out.write(
            '<tr{cls}><th scope="row" class="mono" translate="no">{name}{pills}</th>'
            '<td class="mono">{st}</td><td class="num">{score}</td>'
            '<td class="num">{loci}</td><td class="mono">{sig}</td></tr>'.format(
                cls=cls, name=h(cand.scheme), pills=pills, st=h(cand.st),
                score=h(int(cand.score)), loci=h(cand.num_loci), sig=h(cand.signature))
        )
    out.write('</tbody></table></div>')
    dropped = sorted({c.scheme for c in sample.candidates if c.excluded})
    if dropped:
        out.write('<p class="small muted">Dropped by --exclude before scoring: '
                  '<span class="mono">{0}</span>.</p>'.format(h(", ".join(dropped))))
    out.write('</details>')
    return out.getvalue()


def _novel_block(alleles, opts: HtmlOptions, summary: str) -> str:
    if not alleles:
        return ""
    out = io.StringIO()
    out.write('<details class="novel"><summary>{0}</summary>'.format(h(summary)))
    out.write('<div>')
    out.write('<p class="small">These sequences are the assembly\u2019s version of a locus '
              'where no known allele matched exactly. A capture from the minus strand is shown '
              'already reverse-complemented, so it will read in the opposite direction to the '
              'coordinates in the locus table above. The header line is byte-identical to what '
              '<span class="mono">--novel</span> writes.</p>')
    for allele in alleles:
        out.write('<h4 class="mono" translate="no">{0}</h4>'.format(h(allele.fasta_id)))
        out.write(
            '<p class="small muted">scheme {0}, locus {1}, {2} bp, nearest known allele {3}, '
            'from {4}</p>'.format(
                h(allele.scheme), h(allele.locus), _int_text(allele.length_bp),
                h(allele.nearest_allele), h(allele.source_label))
        )
        if allele.length_bp > int(opts.max_novel_bp):
            out.write('<p class="notice">Sequence of {0} bp not shown (longer than the '
                      '{1} bp report limit). Use --novel to write it to a FASTA file.</p>'.format(
                          _int_text(allele.length_bp), _int_text(opts.max_novel_bp)))
            continue
        if _DNA_RE.fullmatch(allele.seq or "") is None:
            out.write('<p class="notice">This sequence contains characters that are not '
                      'nucleotides, so it is not shown here.</p>')
            continue
        out.write(
            '<pre class="seq mono" translate="no"><span class="seq-head">&gt;{0} {1}</span>\n'
            '{2}</pre>'.format(h(allele.fasta_id), h(allele.source_label),
                               h(_wrap60(allele.seq)))
        )
    out.write('</div></details>')
    return out.getvalue()


def _sample_card(sample, index: int, opts: HtmlOptions, dbdir: str, cfg) -> str:
    status = "FAILED" if sample.failed else sample.status
    _glyph, cls, style, why = _status_info(status)
    out = io.StringIO()
    out.write('<article class="card {cls} b-{style}" id="sample-{i}" data-status={ds}>'.format(
        cls=cls, style=style, i=index, ds=attr(status)))
    out.write('<div class="card-head"><div class="who">')
    out.write('<div class="filename mono" translate="no">{0}</div>'.format(h(sample.label)))
    if sample.path != sample.label:
        out.write('<div class="filename mono small" translate="no">file: {0}</div>'.format(
            h(sample.path)))
    if sample.failed:
        out.write('<div class="headline"><span class="scheme">Not analysed</span></div>')
        out.write('<p class="notice">{0}</p>'.format(
            h(sample.error_text or "This file could not be analysed.")))
        out.write('</div><div class="status-box">{0}<p class="status-why">{1}</p></div>'.format(
            _badge(status), h(why)))
        out.write('</div></article>')
        return out.getvalue()

    out.write('<div class="headline">')
    out.write('<span class="scheme mono" translate="no">{0}</span>'.format(h(sample.scheme)))
    out.write('<span class="st-label">sequence type</span>')
    out.write('<span class="st mono">{0}</span>'.format(h(sample.st)))
    out.write('</div>')
    species = species_label_for(dbdir, sample.scheme) if sample.scheme != "-" else None
    if sample.scheme == "-":
        out.write('<p class="species muted">No scheme matched, so there is no species to show.</p>')
    elif species:
        out.write('<p class="species">Scheme designed for <strong>{0}</strong>. '
                  '<span class="small muted">{1}</span></p>'.format(
                      h(species), h(_SPECIES_CAVEAT)))
    else:
        out.write('<p class="species muted">This scheme is not mapped to a species in '
                  'scheme_species_map.tab. <span class="small">{0}</span></p>'.format(
                      h(_SPECIES_CAVEAT)))
    out.write('<p class="stats small">')
    out.write('<span>{0} contigs</span><span>{1} bp</span>'.format(
        _int_text(sample.n_contigs), _int_text(sample.total_bp)))
    out.write('<span>{0} BLAST hits parsed, {1} kept</span>'.format(
        _int_text(sample.hits_seen), _int_text(sample.hits_kept)))
    out.write('<span>{0:.2f} s</span>'.format(float(sample.elapsed_s)))
    out.write('</p></div>')

    out.write('<div class="status-box">')
    out.write(_badge(sample.status))
    out.write('<p class="status-why">{0}</p>'.format(h(why)))
    out.write(_meter(sample.score))
    out.write('</div></div>')

    ties = _ties(sample)
    if ties:
        out.write(
            '<p class="notice"><strong><span aria-hidden="true">\u26a0</span> Tie:</strong> '
            '{names} also scored {score}. {app} reports the alphabetically first scheme. '
            'Confirm the species by another method before using this ST.</p>'.format(
                names=h(", ".join(ties)), score=h(int(sample.score)),
                app=h(_branding.APP_NAME))
        )
    for warning in sample.warnings:
        out.write('<p class="notice small mono">{0}</p>'.format(h(warning)))

    if sample.alleles:
        out.write(_locus_table(sample, opts.evidence, index))
        out.write('<p class="small muted">{0}</p>'.format(h(_COVERAGE_FOOTNOTE)))
    else:
        out.write('<div class="empty">No loci were called for this sample.</div>')

    out.write(_runner_up_table(sample, cfg))
    if sample.novel:
        out.write(_novel_block(
            sample.novel, opts,
            "Novel allele sequences from this sample ({0})".format(len(sample.novel))))
    out.write('</article>')
    return out.getvalue()


def _batch_section(result: RunResult, dbdir: str) -> str:
    samples = result.samples
    out = io.StringIO()
    out.write('<section class="batch" aria-labelledby="batch-h">')
    out.write('<h2 id="batch-h">All {0} samples</h2>'.format(len(samples)))

    counts: Dict[str, int] = {}
    for sample in samples:
        key = "FAILED" if sample.failed else sample.status
        counts[key] = counts.get(key, 0) + 1
    out.write('<div class="chips">')
    for status in [*list(STATUS_INFO), "FAILED"]:
        if status not in counts:
            continue
        glyph, cls, _s, _w = _status_info(status)
        out.write('<span class="chip {cls}"><span aria-hidden="true">{g}</span>'
                  '{word} {n}</span>'.format(cls=cls, g=h(glyph), word=h(status),
                                             n=counts[status]))
    out.write('</div>')

    out.write('<div class="filters no-print">')
    out.write('<div hidden><label for="filter-text">Filter by text</label>'
              '<input type="search" id="filter-text" placeholder="sample, scheme, ST\u2026"></div>')
    out.write('<div hidden><label for="filter-status">Filter by status</label>'
              '<select id="filter-status"><option value="">Any status</option>')
    for status in [*list(STATUS_INFO), "FAILED"]:
        if status in counts:
            out.write('<option value={0}>{1}</option>'.format(attr(status), h(status)))
    out.write('</select></div>')
    out.write('<button type="button" id="reset-order" hidden>Reset order</button>')
    out.write('</div>')
    out.write('<p class="small muted" id="filter-count" role="status">{0} samples</p>'.format(
        len(samples)))

    out.write('<div class="table-wrap"><table id="batch-table">')
    out.write('<caption>One row per input file, in the order they were given</caption>')
    out.write('<thead><tr>'
              '<th scope="col" data-sort="text" aria-sort="none">Sample</th>'
              '<th scope="col" data-sort="text" aria-sort="none">Scheme</th>'
              '<th scope="col" data-sort="text" aria-sort="none">Species</th>'
              '<th scope="col" data-sort="text" aria-sort="none">ST</th>'
              '<th scope="col" data-sort="text" aria-sort="none">Status</th>'
              '<th scope="col" class="num" data-sort="number" aria-sort="none">Score</th>'
              '<th scope="col">Alleles</th>'
              '</tr></thead><tbody id="batch-body">')
    for i, sample in enumerate(samples, start=1):
        status = "FAILED" if sample.failed else sample.status
        species = species_label_for(dbdir, sample.scheme) or ""
        alleles = ";".join("{0}({1})".format(c.locus, c.code) for c in sample.alleles)
        search = " ".join([sample.label, sample.path, sample.scheme, sample.st,
                           status, species, alleles]).lower()
        out.write('<tr data-index="{i}" data-status={ds} data-search={search}>'.format(
            i=i, ds=attr(status), search=attr(search)))
        out.write('<td class="mono" translate="no"><a href="#sample-{i}">{label}</a></td>'.format(
            i=i, label=h(sample.label)))
        out.write('<td class="mono" translate="no">{0}</td>'.format(h(sample.scheme)))
        out.write('<td>{0}</td>'.format(h(species) if species else '<span class="muted">\u2014</span>'))
        out.write('<td class="mono">{0}</td>'.format(h(sample.st)))
        out.write('<td>{0}</td>'.format(_badge(status)))
        out.write('<td class="num" data-value={0}>{1}</td>'.format(
            attr(int(sample.score)), h(int(sample.score))))
        out.write('<td class="mono small">{0}</td>'.format(
            h(alleles) if alleles else '<span class="muted">\u2014</span>'))
        out.write('</tr>')
    out.write('</tbody></table></div>')

    groups: Dict[Tuple[str, str], int] = {}
    for sample in samples:
        if sample.failed:
            continue
        key = (sample.scheme, sample.st)
        groups[key] = groups.get(key, 0) + 1
    if groups:
        out.write('<details><summary>Sequence-type roll-up</summary>')
        out.write('<div class="table-wrap"><table>')
        out.write('<caption>How many samples share each scheme and sequence type</caption>')
        out.write('<thead><tr><th scope="col">Scheme</th><th scope="col">ST</th>'
                  '<th scope="col" class="num">Samples</th></tr></thead><tbody>')
        for (scheme, st) in sorted(groups):
            out.write('<tr><th scope="row" class="mono" translate="no">{0}</th>'
                      '<td class="mono">{1}</td><td class="num">{2}</td></tr>'.format(
                          h(scheme), h(st), groups[(scheme, st)]))
        out.write('</tbody></table></div>')
        out.write('<p class="small">{0}</p></details>'.format(h(_OUTBREAK_CAVEAT)))
    out.write('</section>')
    return out.getvalue()


def _run_meta_section(result: RunResult, brand: Dict[str, Any]) -> str:
    meta = result.meta
    cfg = meta.config
    out = io.StringIO()

    def row(term: str, value: str, flagged: bool = False) -> None:
        mark = ' <span class="dagger" title="not the default">\u2020</span>' if flagged else ""
        out.write('<dt>{0}</dt><dd>{1}{2}</dd>'.format(h(term), value, mark))

    out.write('<section class="run-meta" aria-labelledby="meta-h">')
    out.write('<h2 id="meta-h" class="visually-hidden">Run details</h2><dl>')
    row("Started", h(meta.started_local))
    row("Started (UTC)", h(meta.started_utc))
    row("Duration", h("{0:.2f} s".format(float(meta.duration_s))))
    row("Samples", h(len(result.samples)))
    row("{0} version".format(brand["app_name"]), h(meta.wmlst_version))
    row("mlst compatibility", h(meta.mlst_compat))
    row("Allele database", h("{0} ({1} schemes)".format(meta.db_version, meta.db_scheme_count)))
    row("BLAST+", h(meta.blast_version))
    row("Computer", h(meta.hostname))
    row("Platform", h(meta.platform))
    row("Python", h(meta.python_version))
    row("Minimum identity", h("{0} %".format(_num(cfg.minid))), float(cfg.minid) != 95.0)
    row("Minimum coverage", h("{0} %".format(_num(cfg.mincov))), float(cfg.mincov) != 50.0)
    if cfg.scheme:
        row("Minimum score", h("0 (forced by --scheme)"), True)
    else:
        row("Minimum score", h(_num(cfg.minscore)), float(cfg.minscore) != 50.0)
    row("Scheme", h(cfg.scheme) if cfg.scheme else h("auto-detect"), bool(cfg.scheme))
    excluded = ", ".join(sorted(cfg.exclude)) if cfg.exclude else "none"
    row("Excluded schemes", h(excluded),
        frozenset(cfg.exclude) != frozenset(
            {"ecoli", "abaumannii", "vcholerae_2", "senterica_achtman_2"}))
    row("BLAST threads", h(cfg.threads), int(cfg.threads) != 1)
    row("Parallel files", h(cfg.jobs), int(cfg.jobs) != 1)
    row("Database directory", '<span class="mono">{0}</span>'.format(h(meta.dbdir)))
    row("Command line", '<span class="mono">{0}</span>'.format(h(" ".join(meta.argv))))
    out.write('</dl>')
    out.write('<p class="small muted">\u2020 marks a setting that differs from the default.</p>')
    if getattr(cfg, "repair_locus_ids", False):
        out.write('<p class="notice"><strong>--repair-locus-ids was used.</strong> '
                  'Results will NOT match tseemann/mlst. This option recovers alleles whose '
                  'identifiers the upstream parser discards; do not compare these calls with '
                  'output from the original tool.</p>')
    out.write('</section>')
    return out.getvalue()


def _legend_section() -> str:
    out = io.StringIO()
    out.write('<section class="legend" aria-labelledby="legend-h">')
    out.write('<h2 id="legend-h">What the allele codes mean</h2><dl>')
    for code, symbol, text in SYMBOL_INFO:
        out.write('<dt class="mono" translate="no">{0}</dt><dd>{1} <span class="muted small">'
                  '({2})</span></dd>'.format(h(code), h(text), h(symbol)))
    out.write('</dl>')
    out.write('<h3>What the status words mean</h3><dl>')
    for status in STATUS_INFO:
        _glyph, _cls, _s, why = STATUS_INFO[status]
        out.write('<dt>{0}</dt><dd>{1}</dd>'.format(_badge(status), h(why)))
    out.write('</dl></section>')
    return out.getvalue()


def _method_section(result: RunResult, brand: Dict[str, Any]) -> str:
    cfg = result.meta.config
    return (
        '<section class="method" aria-labelledby="method-h">'
        '<h2 id="method-h">How these results were produced</h2>'
        '<p class="small">Each input file was converted to FASTA and searched with '
        'BLAST+ {blast} against a database of {n} PubMLST typing schemes '
        '(snapshot {db}). Hits were kept when they covered at least {mincov} % of a '
        'reference allele at {minid} % identity or better. Every scheme with a '
        'surviving hit was scored out of 100 \u2014 90 points for locus quality, plus '
        '10 when the resulting allele combination is a named sequence type \u2014 and '
        'the highest-scoring scheme is the one reported. {app} {ver} reproduces the '
        'scoring arithmetic of mlst {compat} exactly.</p></section>'
    ).format(
        blast=h(result.meta.blast_version), n=h(result.meta.db_scheme_count),
        db=h(result.meta.db_version), mincov=h(_num(cfg.mincov)), minid=h(_num(cfg.minid)),
        app=h(brand["app_name"]), ver=h(brand["version"]), compat=h(result.meta.mlst_compat),
    )


def _colophon(brand: Dict[str, Any]) -> str:
    out = io.StringIO()
    out.write('<footer class="colophon"><div class="wrap">')
    out.write('<p>{0}</p>'.format(brand["credit_html"]))
    out.write('<p><strong>Please cite:</strong></p><ul>')
    for citation in brand["citations"]:
        out.write('<li>{0}</li>'.format(h(citation)))
    out.write('</ul>')
    out.write('<p class="small">{0}. This report is a record of a computation, not a '
              'clinical diagnosis. Confirm any result that will inform patient care or an '
              'outbreak investigation by an independent method.</p>'.format(
                  h(brand["copyright"])))
    out.write('</div></footer>')
    return out.getvalue()


def render_html(result: RunResult, opts: Optional[HtmlOptions] = None,
                branding: Optional[Dict[str, Any]] = None) -> str:
    """Render the complete, self-contained HTML report (section 11).

    Returns one document as a ``str``: no external CSS, JS, fonts or images,
    the single permitted non-inline URI being the base64 vendor mark. Every
    value that came from a FASTA header, a filename, a label or a scheme name
    passes through :func:`h`. Deterministic: two renders of the same
    :class:`~wmlst.engine.RunResult` are byte-identical.
    """
    opts = opts or HtmlOptions()
    brand = _default_branding()
    if branding:
        brand.update(branding)
    meta = result.meta
    dbdir = meta.dbdir
    samples = result.samples
    title = opts.title or "{0} report \u2014 {1} sample{2}".format(
        brand["app_name"], len(samples), "" if len(samples) == 1 else "s")

    out = io.StringIO()
    out.write('<!DOCTYPE html>\n<html lang="en">\n<head>\n')
    out.write('<meta charset="utf-8">\n')
    out.write('<meta name="viewport" content="width=device-width, initial-scale=1">\n')
    # Second lock behind the escaping: file:// user agents honour CSP
    # inconsistently, so this never replaces h(), it only backs it up.
    out.write('<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; '
              'style-src \'unsafe-inline\'; script-src \'unsafe-inline\'; img-src data:; '
              'font-src \'none\'; connect-src \'none\'; form-action \'none\'; '
              'base-uri \'none\'; frame-ancestors \'none\'">\n')
    out.write('<meta name="referrer" content="no-referrer">\n')
    out.write('<meta name="generator" content={0}>\n'.format(
        attr("{0} {1} \u2014 {2}".format(brand["app_name"], brand["version"],
                                         brand["attribution"]))))
    out.write('<title>{0}</title>\n'.format(h(title)))
    out.write('<style>\n{0}\n</style>\n'.format(CSS_TEXT))
    out.write('</head>\n<body>\n')
    out.write('<a class="skip" href="#main">Skip to the results</a>\n')

    out.write('<header class="masthead"><div class="wrap">')
    out.write('<img src={0} alt="" width="44" height="44">'.format(attr(LOGO_SVG_DATA_URI)))
    out.write('<div class="brand"><div class="product">{0} <span class="muted">{1}</span></div>'
              '<div class="vendor">{2}</div></div>'.format(
                  h(brand["app_name"]), h("v" + str(brand["version"])),
                  h(brand["attribution"])))
    out.write('<div class="tools no-print">')
    out.write('<button type="button" id="theme-toggle" aria-pressed="false" hidden>'
              'Light / dark</button>')
    out.write('<button type="button" id="expand-all" data-open="no" hidden>'
              'Expand all sections</button>')
    out.write('<button type="button" data-copy="tsv-block" hidden>Copy results as TSV</button>')
    out.write('<button type="button" data-copy="wmlst-data" hidden>Copy results as JSON</button>')
    out.write('</div></div>')
    out.write('<div class="wrap"><h1>{0}</h1>'
              '<p class="small muted" id="copy-status" role="status"></p></div>'.format(h(title)))
    out.write('</header>\n')

    out.write('<main id="main" class="wrap">\n')
    out.write(_run_meta_section(result, brand))
    if not samples:
        out.write('<div class="empty">No input files were analysed in this run.</div>')
    if len(samples) > 1:
        out.write(_batch_section(result, dbdir))
    out.write('<section aria-labelledby="samples-h"><h2 id="samples-h">Results by sample</h2>')
    for i, sample in enumerate(samples, start=1):
        out.write(_sample_card(sample, i, opts, dbdir, meta.config))
    out.write('</section>')

    if result.novel:
        out.write('<section aria-labelledby="novel-h"><h2 id="novel-h">Novel alleles</h2>')
        out.write(_novel_block(
            result.novel, opts,
            "All {0} novel allele sequences from this run".format(len(result.novel))))
        out.write('</section>')
    elif not meta.config.novel_path and any(
            c.symbol == "novel" for s in samples for c in s.alleles):
        out.write('<p class="small muted">Inexact allele matches were found but their sequences '
                  'were not captured. Re-run with <span class="mono">--novel FILE.fa</span> to '
                  'write them out (they are also included in this report when that option is '
                  'used).</p>')

    out.write(_legend_section())
    out.write(_method_section(result, brand))
    out.write('<details class="no-print"><summary>Plain-text results (the same rows the '
              'command line prints)</summary><pre class="seq mono" id="tsv-block">{0}</pre>'
              '</details>'.format(h(tsv_text(result))))
    out.write('\n</main>\n')

    out.write(_colophon(brand))
    out.write('\n<script type="application/json" id="wmlst-data">{0}</script>\n'.format(
        js_json(_json_records(result))))
    out.write('<script>\n{0}\n</script>\n'.format(JS_TEXT))
    out.write('</body>\n</html>\n')
    return out.getvalue()


def write_html(result: RunResult, path: str, opts: Optional[HtmlOptions] = None,
               branding: Optional[Dict[str, Any]] = None) -> str:
    """Render and write the HTML report atomically; return the final path (4.7).

    Writes ``path + ".tmp"`` in the same directory and then :func:`os.replace`,
    so a half-written report can never be opened. UTF-8, LF, no BOM.
    """
    text = render_html(result, opts, branding)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    os.replace(tmp, path)
    return path
