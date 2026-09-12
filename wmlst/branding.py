# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
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
VENDOR = "IOWA-BioTech"

#: Author.
AUTHOR = "Giovanni Lorenzin"

#: Combined attribution, the canonical visible form.
ATTRIBUTION = f"{VENDOR} — {AUTHOR}"

#: Project home.
HOMEPAGE = "https://github.com/iowa69/WMLST"


#: Window / document title.
WINDOW_TITLE = f"{APP_NAME} — {APP_TAGLINE} — {VENDOR}"

#: Copyright line.
COPYRIGHT = f"© {VENDOR} — {AUTHOR}"

#: One-line CLI banner (stderr, suppressed by --quiet).
BANNER = f"This is {APP_NAME} {__version__} by {ATTRIBUTION}"

#: Citations users should include in publications.
CITATIONS = ()

#: Short credit shown in the GUI About box and the HTML report footer.
CREDIT_HTML = (
    f"{APP_NAME} v{__version__} &middot; {VENDOR} &mdash; {AUTHOR}."
)


#: The PubMLST citation, as a single sentence (section 4.2).
PUBMLST_CITATION = ""

#: The BLAST+ citation, as a single sentence (section 4.2).
BLAST_CITATION = ""

#: Alias kept because section 4.2 names the product constant ``APP``.
APP = APP_NAME


def banner_lines(db_version: str) -> list:
    """The stderr banner printed at start-up unless --quiet (4.2).

    ``db_version`` is accepted for callers that want to log it alongside; it is
    not part of the lines returned here.
    """
    import sys as _sys

    return [
        "This is %s %s running on %s with Python %d.%d.%d"
        % (APP_NAME.lower(), __version__, _sys.platform,
           _sys.version_info[0], _sys.version_info[1], _sys.version_info[2]),
        BANNER,
    ]


def window_title(running=None) -> str:
    """The GUI title bar (section 4.2).

    Without ``running``: ``'WMLST 1.0.0 - MLST typing - IOWA-BioTech - Giovanni
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
        "Allele database: %s" % _get("db_version"),
        "BLAST+: %s" % _get("blast_version"),
        "Platform: %s" % _get("platform"),
    ]
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
    "VENDOR",
    "WINDOW_TITLE",
    "about_text",
    "banner_lines",
    "html_branding",
    "window_title",
]
