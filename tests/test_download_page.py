# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""The GitHub Pages download page, and the workflow that publishes it.

``docs/index.html`` has exactly one job: put the *current* portable zip one
click away from somebody who has never opened a terminal. Everything that could
quietly break that job is asserted here.

The failure this file really exists to prevent is the silent one. The asset
filename carries the version, so it cannot live in the HTML; the deploy
workflow substitutes it in from the GitHub API by matching a handful of
sentinel strings. Edit the HTML carelessly, a sentinel stops matching, the
deploy fails (loudly, by design) or - worse, if the assertion were ever
relaxed - the page keeps shipping a link to last year's build. So the sentinels
are extracted from the workflow's own stamping script and checked against the
page, which means the two files can never drift apart unnoticed.

Runnable two ways::

    python -m pytest tests/test_download_page.py
    python tests/test_download_page.py
"""

from __future__ import annotations

import ast
import os
import re
import sys
import textwrap
from html.parser import HTMLParser

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

PAGE = os.path.join(REPO_ROOT, "docs", "index.html")
WORKFLOW = os.path.join(REPO_ROOT, ".github", "workflows", "pages.yml")
RELEASE_WORKFLOW = os.path.join(REPO_ROOT, ".github", "workflows", "release.yml")
INNO = os.path.join(REPO_ROOT, "packaging", "wmlst.iss")

#: Hosts the page is allowed to name. Everything here is a plain hyperlink or
#: the GitHub REST endpoint the refresh script calls; none of them is a
#: sub-resource, so the page still renders with the network unplugged.
ALLOWED_HOSTS = frozenset({"github.com", "api.github.com", "pubmlst.org"})

#: Hosts that would mean somebody reached for a CDN or a web font.
FORBIDDEN_HOST_SUBSTRINGS = (
    "cdnjs", "jsdelivr", "unpkg", "cdn.", "bootstrapcdn", "googleapis",
    "gstatic", "typekit", "fontawesome", "jquery.com", "polyfill.io",
    "google-analytics", "googletagmanager",
)

#: ``1.0.2``, ``v1.0`` - anything that would go stale the moment we release.
RE_VERSIONISH = re.compile(r"\b\d+\.\d+\.\d+\b|\bv\d+\.\d+\b")

#: The two release assets the page may offer, as built by release.yml and
#: packaging/wmlst.iss. Lower-cased, because the stamper matches case-blind.
PORTABLE_SUFFIX = "-win64-portable.zip"
INSTALLER_SUFFIX = "-win64-setup.exe"


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


HTML = read(PAGE)


# --------------------------------------------------------------------------
# A minimal DOM, because the page is deliberately small enough not to need one.
# --------------------------------------------------------------------------
class _Collector(HTMLParser):
    """Record every start tag with its attributes, plus the <style> bodies."""

    def __init__(self):
        HTMLParser.__init__(self, convert_charrefs=True)
        self.tags = []          # list of (tagname, {attr: value})
        self.styles = []
        self.scripts = []
        self._in = None

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, {k: (v or "") for k, v in attrs}))
        if tag in ("style", "script"):
            self._in = tag

    def handle_endtag(self, tag):
        if tag == self._in:
            self._in = None

    def handle_data(self, data):
        if self._in == "style":
            self.styles.append(data)
        elif self._in == "script":
            self.scripts.append(data)


def parse(html=HTML):
    collector = _Collector()
    collector.feed(html)
    collector.close()
    return collector


DOM = parse()


def tags(name):
    return [attrs for tag, attrs in DOM.tags if tag == name]


def classes(attrs):
    return attrs.get("class", "").split()


# --------------------------------------------------------------------------
# 1. Exactly one primary download control.
# --------------------------------------------------------------------------
def test_page_exists_and_is_html():
    assert HTML.lstrip().lower().startswith("<!doctype html>")
    assert "<html lang=\"en\">" in HTML


def test_exactly_one_primary_download_control():
    """One button. Not two, not a button plus a second equally loud link."""
    buttons = [a for a in tags("a") if "btn" in classes(a)]
    assert len(buttons) == 1, (
        "the page must offer exactly ONE primary download control, found %d: %r"
        % (len(buttons), [b.get("id") or b.get("href") for b in buttons])
    )
    button = buttons[0]
    assert button.get("id") == "get-portable"
    assert "portable" in button["href"]

    # No <button>/<input type=submit> sneaking in as a second call to action.
    assert not tags("button")
    assert not [i for i in tags("input") if i.get("type") in ("submit", "button")]


def test_the_one_button_says_what_it_gives_you():
    text = re.search(r'id="get-portable"[^>]*>(.*?)</a>', HTML, re.S)
    assert text is not None
    label = re.sub(r"\s+", " ", text.group(1)).strip()
    assert label == "Download WMLST for Windows (portable)"


def test_installer_is_a_small_secondary_link_not_a_button():
    installer = [a for a in tags("a") if a.get("id") == "get-installer"]
    assert len(installer) == 1
    assert "btn" not in classes(installer[0])
    # It lives inside the <small class="secondary"> block, below the button.
    assert re.search(
        r'<small class="secondary">.*?id="get-installer".*?</small>', HTML, re.S
    ), "the installer link must sit inside the secondary block"
    assert HTML.index('id="get-portable"') < HTML.index('id="get-installer"')


def test_nothing_else_competes_with_the_portable_zip():
    """No wheel, no sdist, no source tarball, no PyPI - those are distractions."""
    hrefs = " ".join(a.get("href", "") for a in tags("a")).lower()
    for distraction in (".whl", ".tar.gz", "pypi.org", "/archive/", "pip install"):
        assert distraction not in hrefs, "%r must not be offered here" % distraction
    assert "pip install" not in HTML.lower()


# --------------------------------------------------------------------------
# 2. No hard-coded version anywhere in the page.
# --------------------------------------------------------------------------
def test_no_hard_coded_version_string():
    hits = RE_VERSIONISH.findall(HTML)
    assert hits == [], (
        "docs/index.html must contain no version number - the link is resolved "
        "at deploy time and in the browser. Found: %r" % hits
    )


def test_no_hard_coded_asset_filename():
    lowered = HTML.lower()
    for suffix in (PORTABLE_SUFFIX, INSTALLER_SUFFIX):
        # The bare suffix is fine (the browser script matches on it); a full
        # versioned filename is not.
        for match in re.finditer(re.escape(suffix), lowered):
            before = lowered[max(0, match.start() - 40):match.start()]
            assert not re.search(r"wmlst-\d", before), (
                "a versioned asset filename is baked into the page near %r" % before
            )


def test_the_fallback_href_is_the_version_less_latest_redirect():
    """With no stamping and no JS, the button must still reach the download."""
    for attrs in tags("a"):
        if attrs.get("id") in ("get-portable", "get-installer"):
            href = attrs["href"]
            assert href.startswith("https://github.com/iowa69/WMLST/releases/latest"), href
            assert not RE_VERSIONISH.search(href)


# --------------------------------------------------------------------------
# 3. No external resources. The page must render offline, from one file.
# --------------------------------------------------------------------------
def _hosts(text):
    return {m.group(1).lower() for m in re.finditer(r"https?://([^/\s\"')]+)", text)}


def test_no_external_sub_resources():
    for attrs in tags("script"):
        assert "src" not in attrs, "the script must be inline, not fetched: %r" % attrs
    for attrs in tags("link"):
        rel = attrs.get("rel", "").lower()
        assert rel not in ("stylesheet", "preload", "prefetch", "preconnect",
                           "dns-prefetch", "modulepreload"), (
            "no external stylesheet or preconnect: %r" % attrs
        )
    for tag in ("img", "iframe", "video", "audio", "source", "object", "embed"):
        for attrs in tags(tag):
            src = attrs.get("src", "") or attrs.get("data", "")
            assert not src.lower().startswith(("http://", "https://", "//")), (
                "<%s> loads an off-site resource: %r" % (tag, src)
            )


def test_css_is_self_contained():
    css = "\n".join(DOM.styles)
    assert css.strip(), "the page must carry its own stylesheet"
    assert "@import" not in css
    assert not re.search(r"url\(\s*['\"]?(?:https?:)?//", css), "no remote url() in CSS"
    assert "@font-face" not in css, "no web fonts; use the system font stack"


def test_only_allowlisted_hosts_are_named_at_all():
    unexpected = {h for h in _hosts(HTML) if h not in ALLOWED_HOSTS}
    assert not unexpected, "unexpected host(s) in the page: %r" % sorted(unexpected)


def test_no_cdn_or_tracker_anywhere():
    lowered = HTML.lower()
    for needle in FORBIDDEN_HOST_SUBSTRINGS:
        assert needle not in lowered, "%r must not appear in the page" % needle


# --------------------------------------------------------------------------
# 4. The copy a non-technical user actually needs.
# --------------------------------------------------------------------------
def test_smartscreen_note_is_in_plain_language():
    lowered = re.sub(r"\s+", " ", HTML.lower())
    assert "protected your pc" in lowered
    assert "more info" in lowered
    assert "run anyway" in lowered
    # And it must not tell people to turn their defences off.
    assert not re.search(r"(disable|turn off|switch off) (smartscreen|defender)", lowered)


def test_the_three_steps_are_present_and_in_order():
    steps = re.findall(r"<li>(.*?)</li>", HTML, re.S)
    assert len(steps) >= 3
    unzip, run, drop = (re.sub(r"<[^>]+>|\s+", " ", s).lower() for s in steps[:3])
    assert "unzip" in unzip
    assert "wmlst.exe" in run and "double-click" in run
    assert "fasta" in drop and ("drag" in drop or "drop" in drop)


def test_metadata_fields_have_sensible_unstamped_defaults():
    for cls, default in (("f-version", "the latest release"),
                         ("f-size", "ZIP archive"),
                         ("f-date", "published on GitHub")):
        assert '<span class="%s">%s</span>' % (cls, default) in HTML


def test_branding_and_attribution():
    assert "IOWA-Tech" in HTML
    assert "Giovanni Lorenzin" in HTML
    assert "Torsten Seemann" in HTML
    assert "GPL-2.0-only" in HTML
    assert "PubMLST" in HTML
    assert "Not a medical device" in HTML


def test_responsive_and_both_colour_schemes():
    viewport = [m for m in tags("meta") if m.get("name") == "viewport"]
    assert viewport and "width=device-width" in viewport[0]["content"]
    css = "\n".join(DOM.styles)
    assert "prefers-color-scheme: dark" in css
    assert re.search(r"@media\s*\(max-width", css), "no phone-width breakpoint"
    assert 'name="color-scheme"' in HTML or "color-scheme" in css
    assert "focus-visible" in css, "keyboard focus must stay visible"


# --------------------------------------------------------------------------
# 5. The deploy workflow, and its coupling to the page.
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def workflow():
    yaml = pytest.importorskip("yaml")
    return yaml.safe_load(read(WORKFLOW))


def test_workflow_yaml_parses(workflow):
    assert isinstance(workflow, dict)
    assert workflow["name"] == "Pages"
    # PyYAML resolves the bare key `on` to the boolean True (YAML 1.1).
    triggers = workflow.get("on", workflow.get(True))
    assert isinstance(triggers, dict)
    assert triggers["push"]["branches"] == ["main"]
    assert "docs/**" in triggers["push"]["paths"]
    # A release published after the last docs push must re-stamp the page.
    assert triggers["release"]["types"] == ["published"]
    assert "workflow_dispatch" in triggers


def test_workflow_permissions_are_exactly_what_pages_needs(workflow):
    assert workflow["permissions"] == {
        "contents": "read", "pages": "write", "id-token": "write",
    }
    assert workflow["concurrency"]["group"] == "pages"
    assert workflow["concurrency"]["cancel-in-progress"] is False


def test_workflow_uses_the_official_pages_actions(workflow):
    uses = [
        step.get("uses", "")
        for job in workflow["jobs"].values()
        for step in job["steps"]
    ]
    joined = " ".join(uses)
    for action in ("actions/configure-pages@", "actions/upload-pages-artifact@",
                   "actions/deploy-pages@"):
        assert action in joined, "%s is missing from pages.yml" % action
    deploy = workflow["jobs"]["deploy"]
    assert deploy["needs"] == "build" or "build" in deploy["needs"]
    assert deploy["environment"]["name"] == "github-pages"


def stamper_source(workflow):
    """The inline Python the workflow uses to rewrite the page."""
    for job in workflow["jobs"].values():
        for step in job["steps"]:
            run = step.get("run", "")
            if "<<'PY'" in run:
                body = run.split("<<'PY'\n", 1)[1].rsplit("PY\n", 1)[0]
                return textwrap.dedent(body)
    raise AssertionError("pages.yml no longer contains an inline stamping script")


def test_stamper_is_valid_stdlib_only_python(workflow):
    source = stamper_source(workflow)
    tree = ast.parse(source)          # raises SyntaxError if it ever breaks
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert imported <= {"json", "os", "sys", "urllib", "datetime"}, imported


def test_every_sentinel_the_stamper_replaces_occurs_exactly_once(workflow):
    """The one test that keeps docs/index.html and pages.yml in lock-step."""
    tree = ast.parse(stamper_source(workflow))
    sentinels = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        value = node.value
        if "%s" in value or "%d" in value:
            continue          # a replacement template, not a sentinel
        if value.startswith("https://github.com/iowa69/WMLST/releases/latest#") or (
            value.startswith('<span class="f-') and value.endswith("</span>")
        ):
            sentinels.add(value)

    assert len(sentinels) == 5, (
        "expected 5 sentinels (portable href, installer href, version, size, "
        "date), found %r" % sorted(sentinels)
    )
    for sentinel in sorted(sentinels):
        assert HTML.count(sentinel) == 1, (
            "pages.yml substitutes %r, but docs/index.html contains it %d time(s); "
            "the deploy would fail or, worse, stamp the wrong element"
            % (sentinel, HTML.count(sentinel))
        )


def test_stamper_matches_the_asset_names_release_yml_actually_builds(workflow):
    source = stamper_source(workflow)
    assert PORTABLE_SUFFIX in source
    assert INSTALLER_SUFFIX in source
    # release.yml builds WMLST-$v-win64-portable.zip ...
    assert "win64-portable.zip" in read(RELEASE_WORKFLOW)
    # ... and packaging/wmlst.iss emits WMLST-{version}-win64-setup.exe.
    assert "win64-setup" in read(INNO)


def test_browser_refresh_script_matches_the_same_asset_names():
    script = "\n".join(DOM.scripts)
    assert "-win64-portable" in script and "-win64-setup" in script
    assert "api.github.com/repos/iowa69/WMLST/releases/latest" in script
    # Failure must be a no-op, never a broken link.
    assert ".catch(" in script


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
