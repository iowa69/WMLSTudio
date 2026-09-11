# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Every vendor-facing string, defined exactly once.

No other module may hard-code the product, vendor or author name; import from
here so the GUI title bar, the HTML report, the CLI banner and the Windows
installer can never drift apart.
"""

from .version import UPSTREAM_MLST_VERSION, __version__

#: Product name.
APP_NAME = "WMLST"

#: What the product does, in one line.
APP_TAGLINE = "MLST typing for Windows"

#: Vendor / organisation.
VENDOR = "IOWA-Tech"

#: Author.
AUTHOR = "Giovanni Lorenzin"

#: Combined attribution, the canonical visible form.
ATTRIBUTION = f"{VENDOR} — {AUTHOR}"

#: Project home.
HOMEPAGE = "https://github.com/iowa69/WMLST"

#: Upstream tool this is a port of, and its author.
UPSTREAM_NAME = "mlst"
UPSTREAM_AUTHOR = "Torsten Seemann"
UPSTREAM_URL = "https://github.com/tseemann/mlst"

#: Window / document title.
WINDOW_TITLE = f"{APP_NAME} — {APP_TAGLINE} — {VENDOR}"

#: Copyright line.
COPYRIGHT = f"© {VENDOR} — {AUTHOR}"

#: One-line CLI banner (stderr, suppressed by --quiet).
BANNER = (
    f"This is {APP_NAME} {__version__} by {ATTRIBUTION} "
    f"(compatible with {UPSTREAM_NAME} {UPSTREAM_MLST_VERSION})"
)

#: Citations users should include in publications.
CITATIONS = (
    "Jolley et al (2018) 'Open-access bacterial population genomics: BIGSdb "
    "software, the PubMLST.org website and their applications' "
    "https://pubmed.ncbi.nlm.nih.gov/30345391",
    f"Seemann T, {UPSTREAM_NAME}, {UPSTREAM_URL}",
    f"{APP_NAME} {__version__}, {ATTRIBUTION}, {HOMEPAGE}",
)

#: Short credit shown in the GUI About box and the HTML report footer.
CREDIT_HTML = (
    f"{APP_NAME} v{__version__} &middot; {VENDOR} &mdash; {AUTHOR}. "
    f"A Windows port of <em>{UPSTREAM_NAME}</em> by {UPSTREAM_AUTHOR}. "
    "Allele data from PubMLST (Jolley et al. 2018)."
)


#: The PubMLST citation, as a single sentence (section 4.2).
PUBMLST_CITATION = (
    "Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: "
    "BIGSdb software, the PubMLST.org website and their applications. "
    "Wellcome Open Res 2018;3:124. PMID 30345391."
)

#: The BLAST+ citation, as a single sentence (section 4.2).
BLAST_CITATION = (
    "Camacho C et al. BLAST+: architecture and applications. "
    "BMC Bioinformatics 2009;10:421."
)

#: Alias kept because section 4.2 names the product constant ``APP``.
APP = APP_NAME


def banner_lines(db_version: str) -> list:
    """The three stderr banner lines printed at startup unless --quiet (4.2).

    Line 0 reproduces ``bin/mlst:55`` exactly; lines 1-2 carry the vendor and
    the upstream attribution.  ``db_version`` is accepted for callers that want
    to log it alongside; it is not part of the three lines, which must stay
    byte-stable against upstream's first line.
    """
    import sys as _sys

    return [
        "This is %s %s running on %s with Python %d.%d.%d"
        % (APP_NAME.lower(), __version__, _sys.platform,
           _sys.version_info[0], _sys.version_info[1], _sys.version_info[2]),
        BANNER,
        "A port of %s by %s (%s), packaged by %s"
        % (UPSTREAM_NAME, UPSTREAM_AUTHOR, UPSTREAM_URL, ATTRIBUTION),
    ]


def window_title(running=None) -> str:
    """The GUI title bar (section 4.2).

    Without ``running``: ``'WMLST 1.0.0 - MLST typing - IOWA-Tech - Giovanni
    Lorenzin'``.  With it, the middle segment becomes the progress text.
    """
    middle = running if running else APP_TAGLINE
    return "%s %s \u2014 %s \u2014 %s" % (APP_NAME, __version__, middle, ATTRIBUTION)


def html_branding() -> dict:
    """The mapping handed to :func:`wmlst.report.render_html` (section 4.2)."""
    return {
        "vendor": VENDOR,
        "author": AUTHOR,
        "app": APP_NAME,
        "version": __version__,
        "compat": UPSTREAM_MLST_VERSION,
        "homepage": HOMEPAGE,
        "upstream_url": UPSTREAM_URL,
        "upstream_author": UPSTREAM_AUTHOR,
        "pubmlst_citation": PUBMLST_CITATION,
        "blast_citation": BLAST_CITATION,
    }


def about_text(meta=None) -> str:
    """The About-dialog body (section 4.2).

    ``meta`` is any RunMeta-like object; its ``db_version``, ``blast_version``
    and ``platform`` attributes are used when present.  It may be ``None``.
    """
    def _get(name, default="unknown"):
        value = getattr(meta, name, None) if meta is not None else None
        return default if value in (None, "") else str(value)

    lines = [
        "%s %s" % (APP_NAME, __version__),
        APP_TAGLINE,
        COPYRIGHT,
        "",
        "A port of %s %s by %s." % (UPSTREAM_NAME, UPSTREAM_MLST_VERSION,
                                    UPSTREAM_AUTHOR),
        UPSTREAM_URL,
        "",
        "Allele database: %s" % _get("db_version"),
        "BLAST+: %s" % _get("blast_version"),
        "Platform: %s" % _get("platform"),
        "",
        "Please cite:",
    ]
    lines.extend("  " + c for c in CITATIONS)
    return "\n".join(lines)


__all__ = [
    "APP",
    "APP_NAME",
    "APP_TAGLINE",
    "ATTRIBUTION",
    "AUTHOR",
    "BANNER",
    "BLAST_CITATION",
    "CITATIONS",
    "COPYRIGHT",
    "CREDIT_HTML",
    "HOMEPAGE",
    "PUBMLST_CITATION",
    "UPSTREAM_AUTHOR",
    "UPSTREAM_NAME",
    "UPSTREAM_URL",
    "VENDOR",
    "WINDOW_TITLE",
    "about_text",
    "banner_lines",
    "html_branding",
    "window_title",
]
