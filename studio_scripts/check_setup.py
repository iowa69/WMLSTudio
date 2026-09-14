"""Print what this installation has, what it is missing, and what to do about it.

Run this to diagnose a downloaded build without opening the desktop. It reads the
installation; it never downloads, installs or repairs anything. The exit status is
1 when something required is not ready, so a support script can act on it.

    python studio_scripts/check_setup.py
    python studio_scripts/check_setup.py --json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from wmlstudio.provisioning import report

# Fixed width so the states line up in a terminal and in a pasted support message.
BADGES = {"ready": "[ READY   ]", "partial": "[ PARTIAL ]",
          "missing": "[ MISSING ]", "unusable": "[ BROKEN  ]"}
INDENT = " " * 12


def wrap(text, width=88, indent=INDENT):
    """Fold one paragraph onto indented lines without importing a formatter."""
    words, lines, current = str(text).split(), [], indent
    for word in words:
        candidate = word if current.strip() == "" else current + " " + word
        if len(candidate) > width and current.strip():
            lines.append(current)
            current = indent + word
        else:
            current = candidate if current.strip() else indent + word
    if current.strip():
        lines.append(current)
    return lines


def describe(item) -> list[str]:
    """One requirement as a person reads it: state, what breaks, why, and the fix."""
    scope = "required" if item["required"] else "optional"
    lines = [f"{BADGES[item['state']]} {item['title']} ({scope})"]
    if not item["ready"]:
        lines += wrap(item["consequence"])
    lines += wrap("Why: " + item["reason"])
    action = item["action"]
    if action and not item["ready"]:
        owner = "the application can do this" if action["automatic"] else "you need to do this"
        lines += wrap(f"Fix ({owner}): {action['label']}")
        if action["ui_route"]:
            lines += wrap("In the app: " + action["ui_route"], indent=INDENT + "  ")
        if action["detail"]:
            lines += wrap(action["detail"], indent=INDENT + "  ")
        if action["command"]:
            lines += wrap("Terminal: " + action["command"], indent=INDENT + "  ")
    return lines


def render(data) -> list[str]:
    lines = [f"WMLSTudio {data['application_version']} installation check",
             f"Data folder:    {data['data_root']}",
             f"Scheme folders: {len(data['scheme_paths'])} searched",
             f"Checked:        {data['generated_utc']}"
             + ("  (full verification)" if data["verified"] else "  (quick check, no hashing)"),
             ""]
    for item in data["items"]:
        lines += describe(item)
        lines.append("")
    lines.append("Can you run an investigation right now?")
    lines += wrap(data["summary"], indent="  ")
    lines.append("")
    lines.append("What each menu can do right now:")
    for entry in data["capabilities"].values():
        lines.append(f"  {'ready  ' if entry['ready'] else 'blocked'}  {entry['title']}")
        if not entry["ready"]:
            lines += wrap("waiting on: " + ", ".join(entry["blocked_by"]), indent=" " * 12)
    if data["next_action"]:
        action = data["next_action"]
        lines += ["", "Next step: " + action["label"]]
        if action["ui_route"]:
            lines.append("           " + action["ui_route"])
    return lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path,
                        help="Inspect this data folder instead of the one the desktop uses")
    parser.add_argument("--scheme-root", type=Path, action="append", dest="scheme_roots",
                        help="A scheme folder to inspect; repeatable. Defaults to the "
                             "application's own scheme locations")
    parser.add_argument("--hydra-database", type=Path,
                        help="The AMR database store a project has selected")
    parser.add_argument("--quick", action="store_true",
                        help="Skip content hashing and tool execution; states say which ran")
    parser.add_argument("--json", action="store_true", dest="as_json",
                        help="Print the machine-readable report instead of the text one")
    args = parser.parse_args(argv)
    scheme_paths = None
    if args.scheme_roots:
        scheme_paths = [path for root in args.scheme_roots
                        for path in (sorted(p for p in root.iterdir() if p.is_dir())
                                     if root.is_dir() else [])]
    data = report(data_root=args.data_root, scheme_paths=scheme_paths,
                  hydra_database_root=args.hydra_database, verify=not args.quick)
    if args.as_json:
        print(json.dumps(data, indent=2))
    else:
        # A Windows console is not always UTF-8; the report is ASCII, but a path
        # or an upstream message inside it may not be.
        for line in render(data):
            print(line.encode(sys.stdout.encoding or "utf-8", "replace")
                  .decode(sys.stdout.encoding or "utf-8", "replace"))
    return 0 if data["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
