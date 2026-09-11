# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, GPL-2.0-only)
#
# Direct translation of tseemann/mlst 2.35.0 `bin/mlst` lines 40-155, 417-513
# and of `perl5/MLST/Logger.pm`.
"""WMLST command line interface.

Implements docs/ARCHITECTURE.md sections 4.8 (public API), 5.2-5.3 (startup and
per-file preflight), 12 (output surfaces) and 14 (exit codes).

The option table mirrors `mlst` 2.35.0 exactly, including the order and the
wording of the six help sections; the WMLST-only switches live in a seventh
section after them so that a user diffing `--help` against upstream sees the
familiar text first.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter

from . import branding
from .engine import (
    RunConfig,
    WmlstError,
)
from .version import BUNDLED_DB_VERSION, __version__

__all__ = [
    "EXE",
    "Logger",
    "WmlstArgumentParser",
    "build_parser",
    "config_from_args",
    "main",
    "main_make_blast_db",
    "main_show_seqs",
    "main_update_db",
    "reconfigure_streams",
    "usage_text",
]

#: Name the program calls itself in messages and in `--version` (C15).
EXE = "wmlst"

#: Upstream's default --exclude string (bin/mlst:26), kept as a string because
#: that is what the user sees in --help and what --exclude replaces wholesale.
EXCLUDE_DEFAULT = "ecoli,abaumannii,vcholerae_2,senterica_achtman_2"

_DB_HINT = "WMLST_DBDIR"


# ---------------------------------------------------------------------------
# Logger  --  port of perl5/MLST/Logger.pm (section 4.8)
# ---------------------------------------------------------------------------
class Logger:
    """Port of ``MLST::Logger``. Module-global state, set once from RunConfig.

    Implements docs/ARCHITECTURE.md section 4.8. Every method writes to
    ``sys.stderr`` as it is bound *at call time* so that a test harness which
    swaps the stream still captures the output.
    """

    _quiet = False
    _debug = False

    @classmethod
    def configure(cls, quiet: bool = False, debug: bool = False) -> None:
        """Set the global quiet/debug flags (bin/mlst:52-53)."""
        cls._quiet = bool(quiet)
        cls._debug = bool(debug)

    @classmethod
    def quiet(cls) -> bool:
        """True when ``msg`` output is suppressed."""
        return cls._quiet

    @classmethod
    def debug(cls) -> bool:
        """True when ``dbg`` output is enabled."""
        return cls._debug

    @staticmethod
    def _emit(text: str) -> None:
        stream = sys.stderr
        if stream is None:  # pythonw.exe before the runtime hook installs a stub
            return
        stream.write(text)
        try:
            stream.flush()
        except ValueError:  # pragma: no cover - closed stream at interpreter exit
            pass

    @classmethod
    def msg(cls, *args) -> None:
        """Informational line on stderr. SUPPRESSED by --quiet."""
        if cls._quiet:
            return
        cls._emit(" ".join(str(a) for a in args) + "\n")

    @classmethod
    def notice(cls, *args) -> None:
        """An output-location notice on stderr. NOT suppressed by --quiet (12.7)."""
        cls._emit(" ".join(str(a) for a in args) + "\n")

    @classmethod
    def wrn(cls, *args) -> None:
        """``WARNING: ...`` on stderr. NOT suppressed by --quiet."""
        cls._emit("WARNING: " + " ".join(str(a) for a in args) + "\n")

    @classmethod
    def err(cls, *args) -> None:
        """``ERROR: ...`` on stderr, then exit 1. Not suppressed (divergence D3)."""
        cls._emit("ERROR: " + " ".join(str(a) for a in args) + "\n")
        sys.exit(1)

    @classmethod
    def dbg(cls, *args) -> None:
        """Debug block on stderr, only with --debug. test.sh:103-106 greps ``=== DEBUG``."""
        if not cls._debug:
            return
        payload = " ".join(str(a) for a in args)
        cls._emit(
            "=" * 10 + " DEBUG START " + "=" * 10 + "\n"
            + payload + "\n"
            + "=" * 10 + " DEBUG END " + "=" * 10 + "\n"
        )


def reconfigure_streams() -> None:
    """Force UTF-8 and LF on the standard streams (section 4.8).

    Without this a non-cp1252 filename in the FILE column kills the run at the
    final print, after every expensive step has already succeeded.
    """
    # stderr gets LF too, not the platform default. It is a byte-identity
    # surface: the tie WARNING and the msg() lines are part of the golden corpus,
    # and leaving stderr on Windows' default put a CR on every one of them while
    # stdout stayed clean -- the two streams of one run disagreeing.
    for stream in (sys.stdout, sys.stderr):
        if stream is None or not hasattr(stream, "reconfigure"):
            continue
        try:
            stream.reconfigure(encoding="utf-8", errors="replace", newline="\n")
        except (ValueError, AttributeError, OSError):  # pragma: no cover
            # A redirected or already-detached stream; nothing to reconfigure.
            pass


# ---------------------------------------------------------------------------
# The option table  --  bin/mlst:440-473
# ---------------------------------------------------------------------------
class _Opt:
    """One row of the upstream ``@Options`` table.

    ``spec`` uses Getopt::Long syntax (``name!``/``name=s``/``name=i``/``name=f``)
    so that the help renderer can reproduce ``bin/mlst:496-505`` character for
    character.
    """

    __slots__ = ("aliases", "default", "desc", "dest", "help_default", "spec")

    def __init__(self, spec, dest=None, default=None, desc="", aliases=(),
                 help_default=None):
        self.spec = spec
        self.dest = dest
        self.default = default
        self.desc = desc
        self.aliases = tuple(aliases)
        self.help_default = help_default

    @property
    def name(self):
        """The bare option name, e.g. ``threads``."""
        return self.spec.rstrip("!").split("=")[0]

    @property
    def kind(self):
        """``"bool"``, ``"str"``, ``"int"``, ``"num"`` or ``"action"``."""
        if self.spec.endswith("!"):
            return "bool"
        if self.spec.endswith("=s"):
            return "str"
        if self.spec.endswith("=i"):
            return "int"
        if self.spec.endswith("=f"):
            return "num"
        return "action"


def _opts():
    """Build the option table: headings (str) interleaved with :class:`_Opt`."""
    return [
        "GENERAL",
        _Opt("help", desc="This help"),
        _Opt("version", desc="Print version and exit"),
        _Opt("check!", "check", False, "Just check dependencies and exit"),
        _Opt("skipcheck!", "skipcheck", False, "Skip dependency check at runtime"),
        _Opt("quiet!", "quiet", False, "Quiet - no stderr output", aliases=("-q",)),
        _Opt("threads=i", "threads", 1,
             "Number of BLAST threads (suggest GNU Parallel instead)"),
        _Opt("debug!", "debug", False, "Verbose debug output to stderr"),
        _Opt("fofn=s", "fofn", "", "File of input filenames to use"),
        "SCHEME",
        _Opt("scheme=s", "scheme", "", "Do not auto-detect, force this scheme"),
        _Opt("info!", "info", False, "More information about schemes"),
        _Opt("list!", "list", False, "List available MLST scheme names"),
        _Opt("longlist!", "longlist", False, "List allelles for all MLST schemes"),
        _Opt("exclude=s", "exclude", EXCLUDE_DEFAULT, "Ignore these schemes"),
        "OUTPUT FORMAT",
        _Opt("full!", "full", False, "More detailed output format (recommended)"),
        _Opt("legacy!", "legacy", False,
             "Use old legacy output with allele header row (requires --scheme)"),
        _Opt("csv!", "csv", False, "Output CSV instead of TSV"),
        _Opt("label=s", "label", "", "Replace FILE with this name instead"),
        _Opt("nopath!", "nopath", False, "Strip filename paths from FILE column"),
        "OUTPUT FILES",
        _Opt("outfile=s", "outfile", "", "Save output to this file [STDOUT]"),
        _Opt("novel=s", "novel_path", "", "Save novel alleles to this FASTA file"),
        _Opt("json=s", "json_path", "",
             "Also write results to this file in JSON format"),
        "SCORING",
        _Opt("minid=f", "minid", 95.0,
             "DNA %identity of full allelle to consider similar/'~''"),
        _Opt("mincov=f", "mincov", 50.0,
             "DNA %cov to report partial/'?' allele"),
        _Opt("minscore=f", "minscore", 50.0,
             "Minumum score out of 100 to match a scheme (when auto --scheme)"),
        "DATABASE",
        _Opt("blastdb=s", "blastdb", "",
             "BLAST database prefix", help_default="$%s/blast/mlst.fa" % _DB_HINT),
        _Opt("datadir=s", "datadir", "",
             "PubMLST data folders", help_default="$%s/pubmlst" % _DB_HINT),
        "WMLST EXTRAS",
        _Opt("html=s", "html_path", "", "Write a self-contained HTML report here"),
        _Opt("evidence-tsv=s", "evidence_tsv_path", "",
             "Write the per-hit evidence table here"),
        _Opt("jobs=i", "jobs", 1, "Number of input files to analyse at once"),
        _Opt("gui!", "gui", False, "Launch the graphical interface"),
        _Opt("update-db!", "update_db", False,
             "Refresh the allele database from PubMLST and Institut Pasteur"),
        _Opt("check-only!", "check_only", False,
             "With --update-db, only report what is out of date"),
        _Opt("rollback-db!", "rollback_db", False,
             "Restore the database snapshot the last update kept"),
        _Opt("make-blast-db!", "make_blast_db", False,
             "Rebuild the derived BLAST index from db/pubmlst"),
        _Opt("bootstrap-blast!", "bootstrap_blast", False,
             "Download and verify NCBI BLAST+ into your user profile"),
        _Opt("blast-timeout=f", "blast_timeout_s", 900.0,
             "Seconds before a stuck blastn is killed"),
        _Opt("repair-locus-ids!", "repair_locus_ids", False,
             "Recover alleles whose id the sseqid regex drops"
             " - results will NOT match tseemann/mlst"),
        _Opt("html-evidence=s", "html_evidence", "best",
             "How much BLAST evidence the HTML report shows: none|best|all"),
    ]


OPTIONS = _opts()


# ---------------------------------------------------------------------------
# --help  --  bin/mlst:485-513
# ---------------------------------------------------------------------------
def usage_text() -> str:
    """Render the help screen exactly like ``bin/mlst:485-513``.

    Implements docs/ARCHITECTURE.md section 13.4: the six upstream section
    headings appear in order with the upstream misspellings intact, and the
    output is snapshot-compared against ``tests/golden/help.txt``.
    """
    out = []
    out.append("SYNOPSIS\n  Automatic MLST calling from assembled contigs\n")
    out.append("USAGE\n")
    out.append("  %% %s --longlist                                        "
               "# list known schemes\n" % EXE)
    out.append("  %% %s [options] <contigs.{fasta,gbk,embl}[.gz]          "
               "# auto-detect scheme\n" % EXE)
    out.append("  %% %s --scheme <scheme> <contigs.{fasta,gbk,embl}[.gz]> "
               "# force a scheme\n" % EXE)
    for entry in OPTIONS:
        if isinstance(entry, str):
            out.append(entry + "\n")
            continue
        kind = entry.kind
        if kind == "bool":
            # Perl renders a defined-but-false DEFAULT as " (default '0')", which
            # is truthy, so every negatable flag prints "(default OFF)".
            suffix = " (default OFF)"
        elif entry.help_default is not None:
            suffix = " (default '%s')" % entry.help_default
        elif entry.default is None:
            suffix = ""
        else:
            suffix = " (default '%s')" % _perl_str(entry.default)
        opt = entry.name
        if kind == "str":
            opt += " STR"
        elif kind == "int":
            opt += " INT"
        elif kind == "num":
            opt += " NUM"
        out.append("  --%-15s %s%s\n" % (opt, entry.desc, suffix))
    out.append("HOMEPAGE\n  %s - %s\n" % (branding.HOMEPAGE, branding.ATTRIBUTION))
    out.append("UPSTREAM\n  %s - %s\n"
               % (branding.UPSTREAM_URL, branding.UPSTREAM_AUTHOR))
    return "".join(out)


def _perl_str(value) -> str:
    """Render a default the way Perl interpolates it (95.0 -> 95, '' -> '')."""
    if isinstance(value, float) and value == int(value):
        return str(int(value))
    return str(value)


class _HelpAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        sys.stdout.write(usage_text())
        parser.exit(0)


class _VersionAction(argparse.Action):
    def __init__(self, option_strings, dest=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        # C15: the bare tool name and the WMLST version; the upstream-compat
        # version lives on the banner and in --check, never here.
        sys.stdout.write("%s %s\n" % (EXE, __version__))
        parser.exit(0)


class WmlstArgumentParser(argparse.ArgumentParser):
    """argparse with the two behaviours upstream's Getopt::Long has (section 4.8).

    ``error()`` turns argparse's "unrecognized arguments" into
    ``Unknown option: <name>`` with NO usage block (test.sh:32-36), and
    ``print_usage()`` is suppressed entirely.
    """

    def error(self, message):
        text = str(message)
        unknown = None
        if text.startswith("unrecognized arguments:"):
            unknown = text.split(":", 1)[1].strip().split()[0]
        elif text.startswith("no such option:"):
            unknown = text.split(":", 1)[1].strip().split()[0]
        if unknown is not None:
            # Getopt::Long names the OPTION, never the value glued to it:
            # `mlst --bogus=1` prints "Unknown option: bogus" (finding 25).
            name = unknown.lstrip("-").partition("=")[0]
            sys.stderr.write("Unknown option: %s\n" % name)
        else:
            sys.stderr.write("ERROR: %s\n" % text)
        self.exit(1)

    def print_usage(self, file=None):
        return

    def print_help(self, file=None):
        (file or sys.stdout).write(usage_text())


def build_parser() -> WmlstArgumentParser:
    """Build the argument parser from :data:`OPTIONS` (section 4.8).

    Every negatable boolean uses :class:`argparse.BooleanOptionalAction` so that
    ``--no-quiet``, ``--no-full``, ... behave like Getopt::Long's ``name!``
    (test.sh:63,96).
    """
    # allow_abbrev=False: prefix resolution is done by _expand_abbrevs()
    # instead, because argparse has no way to make an upstream option win
    # a prefix it shares with one of the WMLST extras (finding 15).
    parser = WmlstArgumentParser(prog=EXE, add_help=False, allow_abbrev=False)
    parser.add_argument("--help", "-h", action=_HelpAction)
    parser.add_argument("--version", action=_VersionAction)
    for entry in OPTIONS:
        if isinstance(entry, str) or entry.kind == "action":
            continue
        flags = ["--" + entry.name, *list(entry.aliases)]
        if entry.kind == "bool":
            parser.add_argument(*flags, dest=entry.dest, default=False,
                                action=argparse.BooleanOptionalAction)
        elif entry.kind == "int":
            parser.add_argument(*flags, dest=entry.dest, default=entry.default,
                                type=int, metavar="INT")
        elif entry.kind == "num":
            parser.add_argument(*flags, dest=entry.dest, default=entry.default,
                                type=float, metavar="NUM")
        else:
            parser.add_argument(*flags, dest=entry.dest, default=entry.default,
                                metavar="STR")
    # Accepted spelling from the divergence register (D14); the documented flag
    # is --evidence-tsv.
    parser.add_argument("--tsv-evidence", dest="evidence_tsv_path",
                        default=argparse.SUPPRESS, metavar="STR",
                        help=argparse.SUPPRESS)
    parser.add_argument("files", nargs="*", metavar="FILE")
    return parser


def _upstream_long_names():
    """Every ``--name`` upstream's own @Options table defines (bin/mlst:440-473).

    Includes the ``--no-<flag>`` spellings Getopt::Long's ``name!`` implies, so
    that a prefix such as ``--no-q`` still resolves the upstream way.
    """
    names = {"help", "version"}
    for entry in OPTIONS:
        if isinstance(entry, str):
            if entry == "WMLST EXTRAS":
                break      # everything past this heading is WMLST-only
            continue
        names.add(entry.name)
        if entry.kind == "bool":
            names.add("no-" + entry.name)
    return frozenset(names)


#: Option names that exist in `mlst` 2.35.0 as well as in WMLST.
UPSTREAM_LONG_NAMES = _upstream_long_names()


def _ambiguous(name, hits) -> None:
    """Getopt::Long's wording for an unresolvable prefix, then exit 1."""
    sys.stderr.write("Option %s is ambiguous (%s)\n"
                     % (name, ", ".join(sorted(hits))))
    raise SystemExit(1)


def _option_namespace(parser):
    """Every option name Getopt::Long would answer to -> its argparse spelling.

    ``bin/mlst:438`` is a bare ``use Getopt::Long;`` with no
    ``Getopt::Long::Configure`` call anywhere in the file, so every default is
    live.  Two of them widen the name space (finding 25):

    * a ``name!`` option answers to the hyphen-free negation too, so ``csv!``
      takes ``--nocsv`` as well as ``--no-csv``;
    * ``ignore_case`` is on, so ``--CSV`` is ``--csv``.

    The value of each entry is the canonical argparse flag name, which is what
    gets substituted back into argv.
    """
    space = {}
    for action in parser._actions:
        for flag in action.option_strings:
            if flag.startswith("--") and len(flag) > 2:
                space[flag[2:]] = flag[2:]
    # `no-csv` -> also reachable as `nocsv`; never overwrite a real option
    # (`nopath` and `novel` are options in their own right).
    for name in list(space):
        if name.startswith("no-"):
            space.setdefault(name.replace("-", "", 1), name)
    return space


def _resolve_long_name(name, space):
    """Match one option name the way Getopt::Long does. -> canonical or None.

    Exact first, then case-insensitively, then as a unique prefix.  A prefix
    that is ambiguous exits 1 with upstream's wording; ``None`` means nothing
    matched, which leaves the token alone so the ``Unknown option:`` path still
    reports it.
    """
    if name in space:
        return space[name]
    lowered = name.lower()
    for spelling, canonical in space.items():
        if spelling.lower() == lowered:
            return canonical
    hits = sorted({canonical for spelling, canonical in space.items()
                   if spelling.lower().startswith(lowered)})
    # Upstream has no WMLST extras, so `--t`/`--j`/`--e`/`--b`/`--h` are
    # unambiguous there and must stay unambiguous here: when a prefix matches
    # exactly one upstream option, that option wins even though a WMLST-only
    # option shares the prefix (`--threads` over `--tsv-evidence`, `--help`
    # over `--html`).  A prefix that is ambiguous upstream too (`--l`) still
    # fails, with upstream's message.
    upstream = [n for n in hits if n in UPSTREAM_LONG_NAMES]
    if len(upstream) == 1:
        hits = upstream
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        _ambiguous(name, hits)
    return None


def _expand_abbrevs(argv, parser):
    """Rewrite argv into the spellings argparse knows (section 4.8, finding 25).

    Covers the three Getopt::Long leniencies argparse has no equivalent for --
    single-dash long options (``-csv``), case-insensitive names (``--CSV``) and
    hyphen-free negation (``--nocsv``) -- plus ``auto_abbrev``.  Parsing stops
    at ``--``, ``-`` (stdin) is never touched, and the registered short flags
    ``-q``/``-h`` keep their own meaning.
    """
    space = _option_namespace(parser)
    shorts = frozenset(flag for action in parser._actions
                       for flag in action.option_strings
                       if not flag.startswith("--"))
    out = []
    for i, tok in enumerate(argv):
        if tok == "--":
            out.extend(argv[i:])
            return out
        if tok.startswith("--") and len(tok) > 2:
            lead = "--"
        elif (tok.startswith("-") and len(tok) > 2
              and not tok.startswith("--") and tok not in shorts):
            # Getopt::Long is not bundling by default, so a multi-character
            # single-dash token IS a long option: `-csv`, `-scheme=bogus`.
            lead = "-"
        else:
            out.append(tok)
            continue
        name, sep, value = tok[len(lead):].partition("=")
        if name:
            canonical = _resolve_long_name(name, space)
            if canonical is not None:
                tok = "--" + canonical + (sep + value if sep else "")
        out.append(tok)
    return out


# ---------------------------------------------------------------------------
# argv -> RunConfig  --  bin/mlst:70-144, section 3.6.1
# ---------------------------------------------------------------------------
def _read_fofn(path: str) -> tuple:
    """Read a file of filenames (step 7, divergences D8 and D9)."""
    try:
        with open(path, encoding="utf-8", errors="replace", newline=None) as fh:
            raw = fh.read().splitlines()
    except OSError as exc:
        Logger.err("Unable to read from '%s': %s" % (path, exc.strerror or exc))
    base = os.path.dirname(os.path.abspath(path))
    files = []
    for line in raw:
        name = line.strip()
        if not name:
            continue  # D9: upstream keeps the blank and then dies on -r ''
        if name != "-" and not os.path.isabs(name) and not os.path.exists(name):
            # D8: a relative entry the CWD cannot resolve is retried against the
            # FOFN's own directory. Entries that DO resolve against the CWD are
            # left verbatim, so the FILE column still matches upstream's.
            here = os.path.join(base, name)
            if os.path.exists(here):
                name = here
        files.append(name)
    return tuple(files)


def _resolve_db(ns):
    """Resolve dbdir/datadir/blastdb (steps 9, 11; sections 4.5, 3.6.1)."""
    from . import schemes  # local: keeps --help working without the whole stack

    # bin/mlst:471-472 defaults --blastdb and --datadir INDEPENDENTLY off
    # $MLST_DBDIR, so supplying one of them never moves the other (nor the
    # root itself, which bin/mlst:78 tests unconditionally). Deriving a dbdir
    # "hint" from either flag made `--blastdb /elsewhere/blast/mlst.fa` look
    # for pubmlst/ next to it and abort.
    try:
        dbdir = schemes.resolve_dbdir(None)
    except WmlstError as exc:
        Logger.err(str(exc))
        raise  # pragma: no cover - err() never returns
    datadir = (os.path.abspath(ns.datadir) if ns.datadir
               else os.path.join(dbdir, "pubmlst"))
    blastdb = (os.path.abspath(ns.blastdb) if ns.blastdb
               else schemes.resolve_blastdb(dbdir))
    if not os.path.isdir(datadir):
        Logger.err("Database directory does not exist: %s" % datadir)
    return dbdir, datadir, blastdb


#: One SchemeCatalog per datadir; --scheme validation and --list would otherwise
#: each walk the 162 scheme directories.
_CATALOGS = {}


def _catalog_for(datadir: str):
    """Return the (memoised) :class:`~wmlst.schemes.SchemeCatalog` for ``datadir``."""
    from . import schemes

    catalog = _CATALOGS.get(datadir)
    if catalog is None:
        catalog = _CATALOGS[datadir] = schemes.SchemeCatalog(datadir)
    return catalog


def config_from_args(ns, argv=None) -> RunConfig:
    """Turn parsed options into a :class:`~wmlst.engine.RunConfig`.

    Applies every validation invariant of docs/ARCHITECTURE.md section 3.6.1 and
    the ``--scheme`` side effects, in the order upstream applies them
    (bin/mlst:70-144). Errors exit the process through :meth:`Logger.err`.
    """
    files = tuple(ns.files)
    if ns.fofn:
        Logger.msg("Using '%s' as source of input files" % ns.fofn)
        files = _read_fofn(ns.fofn)
        Logger.msg("Read", len(files), "filenames")

    # step 8 - checked BEFORE the database check (bin/mlst:76)
    if ns.label and len(files) > 1:
        Logger.err("Using --label when scanning multiple files does not make sense")

    dbdir, datadir, blastdb = _resolve_db(ns)

    # step 10 (bin/mlst:81)
    if ns.legacy and not ns.scheme:
        Logger.err("Must specify a --scheme for --legacy output mode")

    catalog = _catalog_for(datadir)

    scheme = ns.scheme or None
    if scheme is not None and scheme not in catalog:
        # Case-sensitive set membership, never a path test (step 12).
        Logger.err("Invalid --scheme '%s'. Check using --list" % scheme)

    minscore = float(ns.minscore)
    if scheme is not None and minscore != 0:
        Logger.msg("Setting --minscore=0 because user chose --scheme")
        minscore = 0.0

    listing = bool(ns.list or ns.longlist or ns.info)
    if not files and not listing and not ns.gui:
        # step 15 (bin/mlst:121)
        Logger.err("Please provide some FASTA/Genbank files to genotype (can be .gz)")

    exclude = frozenset(s for s in str(ns.exclude).split(",") if s)
    if scheme is not None:
        exclude = frozenset()  # bin/mlst:125
    if exclude and not listing:
        Logger.msg("Excluding", len(exclude), "schemes:", *sorted(exclude))

    threads, jobs = _resolve_concurrency(ns)

    return RunConfig(
        files=files,
        minid=float(ns.minid),
        mincov=float(ns.mincov),
        minscore=minscore,
        scheme=scheme,
        exclude=exclude,
        full=bool(ns.full),
        legacy=bool(ns.legacy),
        csv=bool(ns.csv),
        label=(ns.label or None),
        nopath=bool(ns.nopath),
        outfile=(ns.outfile or None),
        json_path=(ns.json_path or None),
        novel_path=(ns.novel_path or None),
        html_path=(ns.html_path or None),
        evidence_tsv_path=(getattr(ns, "evidence_tsv_path", "") or None),
        threads=threads,
        jobs=jobs,
        quiet=bool(ns.quiet),
        debug=bool(ns.debug),
        blast_timeout_s=float(ns.blast_timeout_s),
        dbdir=dbdir,
        datadir=datadir,
        blastdb=blastdb,
        blastn=None,
        repair_locus_ids=bool(ns.repair_locus_ids),
        html_evidence=_html_evidence(ns.html_evidence),
    )


def _html_evidence(value: str) -> str:
    """Validate --html-evidence (section 3.6)."""
    text = str(value).lower()
    if text not in ("none", "best", "all"):
        Logger.err("--html-evidence must be one of: none, best, all")
    return text


def _resolve_concurrency(ns):
    """Invariants 4 and 5 of section 3.6.1: clamp --threads and --jobs."""
    cpus = os.cpu_count() or 1
    threads = int(ns.threads)
    jobs = int(ns.jobs)
    if threads < 1:
        Logger.err("--threads must be 1 or more")
    if jobs < 1:
        Logger.err("--jobs must be 1 or more")
    ceiling = min(cpus, 64)
    if threads > ceiling:
        Logger.wrn("--threads %d is more than this machine can use; using %d"
                   % (threads, ceiling))
        threads = ceiling
    if jobs * threads > cpus:
        allowed = max(1, cpus // jobs)
        if allowed < threads:
            Logger.wrn("--jobs %d x --threads %d exceeds %d CPUs; using --threads %d"
                       % (jobs, threads, cpus, allowed))
            threads = allowed
    return threads, jobs


# ---------------------------------------------------------------------------
# Startup helpers  --  bin/mlst:52-68
# ---------------------------------------------------------------------------
def _banner() -> None:
    """Step 5: the stderr banner, suppressed by --quiet (sections 4.2, 5.2)."""
    for line in branding.banner_lines(BUNDLED_DB_VERSION):
        Logger.msg(line)


def _dependency_check(ns) -> None:
    """Step 6 (bin/mlst:58-68). Exits 0 when --check succeeds."""
    if ns.skipcheck:
        # --skipcheck wins: upstream never reaches the --check branch either.
        Logger.msg("Skipping dependency check due to --skipcheck")
        return
    Logger.msg("Checking %s dependencies:" % EXE)
    from . import blastbin

    try:
        tools = blastbin.find_blast()
    except WmlstError as exc:
        Logger.err(str(exc))
        return  # pragma: no cover - err() never returns
    Logger.msg("Found blastn - %s (%s)" % (tools.blastn, tools.version))
    if tools.version != blastbin.BLAST_PINNED_VERSION:
        # Worth saying, but only through Logger.msg: as a logging warning from
        # blastbin it survived --quiet and broke stderr byte-parity on every
        # distro-packaged BLAST. Supported, simply not the pinned build.
        Logger.msg("  (WMLST is validated against %s; %s is supported)"
                   % (blastbin.BLAST_PINNED_VERSION, tools.version))
    Logger.msg("any2fasta is built in - no external copy needed")
    if ns.check:
        _check_report(ns, tools)
        Logger.msg("OK.")
        sys.exit(0)


def _check_report(ns, tools) -> None:
    """The extra lines --check prints before ``OK.`` (step 6)."""
    dbdir, datadir, blastdb = _resolve_db(ns)
    Logger.msg("Database directory:", dbdir)
    Logger.msg("PubMLST data:", datadir)
    Logger.msg("BLAST database:", blastdb)
    version_file = os.path.join(dbdir, "VERSION.txt")
    db_version = "Unknown"
    if os.path.isfile(version_file):
        with open(version_file, encoding="utf-8", errors="replace") as fh:
            db_version = fh.readline().strip() or "Unknown"
    Logger.msg("Database version:", db_version)
    Logger.msg(branding.ATTRIBUTION)


def _listing(ns, cfg) -> bool:
    """Steps 14 and section 12.5: --list / --longlist / --info. True when handled."""
    if not (ns.list or ns.longlist or ns.info):
        return False
    from . import report

    catalog = _catalog_for(cfg.datadir)
    # Upstream's listing block ends in exit(0) at bin/mlst:118, twenty lines
    # BEFORE `$OUTSEP = ',' if $csv` at bin/mlst:136 - so --csv never reaches
    # --longlist/--info and these two surfaces are unconditionally tab
    # separated (and print_row therefore never quotes a bare comma).
    if ns.list:
        report.write_list(catalog, sys.stdout)
    elif ns.longlist:
        report.write_longlist(catalog, sys.stdout, "\t")
    else:
        report.write_info(catalog, sys.stdout, "\t")
    sys.stdout.flush()
    return True


def _preflight(cfg) -> None:
    """Step 20: readability, then the directory test (bin/mlst:157-158)."""
    for path in cfg.files:
        if path == "-":
            continue
        if not os.access(path, os.R_OK):
            Logger.err("Unable to read from '%s'" % path)
        if os.path.isdir(path):
            Logger.err("'%s' seems to be a directory, not a file" % path)


#: Upstream emits this one through msg() (bin/mlst:344), so --quiet hides it.
_DUPLICATE_EXACT = "found additional exact allele match"


def _emit_warning(text) -> None:
    """Route one engine warning the way upstream routed it.

    ``bin/mlst:344`` announces a duplicate exact allele through ``msg()``, so
    ``--quiet`` hides it and ``tests/golden/mixedzip.err`` is empty; that line
    now reaches stderr through :func:`_msg_sink`, in hit order, and is dropped
    here. The equal-score tie at ``bin/mlst:405`` goes through ``wrn()``, which
    ``--quiet`` must NOT hide - ``tests/golden/equality.err`` holds that line
    even though the golden was generated with ``--quiet``. Both arrive here
    already carrying the ``WARNING: `` prefix, so the text itself is the only
    discriminator.
    """
    text = str(text)
    body = text[9:] if text.startswith("WARNING: ") else text
    if body.startswith(_DUPLICATE_EXACT):
        # Already printed by _msg_sink, in hit order, from walk.messages - the
        # engine appends this one line to BOTH messages and warnings, so
        # emitting it here as well would double it.
        return
    Logger.wrn(body)


def _make_warn_sink():
    """Live warning callback plus the multiset of what it has already emitted."""
    seen = Counter()

    def sink(text):
        seen[str(text)] += 1
        _emit_warning(text)

    return sink, seen


def _drain_warnings(result, seen) -> None:
    """Emit any SampleResult.warnings the live callback never saw."""
    for sample in result.samples:
        for text in sample.warnings:
            if seen[str(text)] > 0:
                seen[str(text)] -= 1
                continue
            _emit_warning(text)


#: The engine announces the novel-allele count too (engine.py:1305); upstream
#: emits it exactly once and AFTER the result table (bin/mlst:247), which is
#: where :func:`_run` emits it, so the engine's copy is dropped here.
_ENGINE_NOVEL_RE = re.compile(r"^Found \d+ novel alleles$", re.ASCII)


def _msg_sink(text) -> None:
    """Route one engine ``msg()`` line, upstream's way (bin/mlst:344-349).

    ``Found exact allele match ...`` and ``WARNING: found additional exact
    allele match ...`` both leave upstream through ``msg()``, in hit order, so
    ``--quiet`` hides them and they interleave exactly as Perl prints them.
    """
    text = str(text)
    if _ENGINE_NOVEL_RE.match(text):
        return
    Logger.msg(text)


def _progress_sink(event) -> None:
    """Render engine progress at debug level only (section 5.16, C17)."""
    Logger.dbg("%s %5.1f%% %s" % (event.phase, event.percent, event.text))


# ---------------------------------------------------------------------------
# The run  --  bin/mlst:132-250
# ---------------------------------------------------------------------------
def _run(cfg, argv) -> int:
    """Analyse every input file and write every requested output (section 14)."""
    from . import engine as engine_mod
    from . import report

    eng = engine_mod.Engine(cfg)
    warn_sink, seen = _make_warn_sink()
    try:
        result = eng.analyse(progress=_progress_sink, warn=warn_sink,
                             msg=_msg_sink, dbg=Logger.dbg)
    finally:
        close = getattr(eng, "close", None)
        if close is not None:
            close()
    _drain_warnings(result, seen)

    if cfg.outfile:
        with open(cfg.outfile, "w", encoding="utf-8", newline="\n") as fh:
            report.write_tsv(result, fh)
    else:
        report.write_tsv(result, sys.stdout)
        sys.stdout.flush()

    if cfg.json_path:
        Logger.msg("Writing JSON: %s" % cfg.json_path)
        report.write_json(result, cfg.json_path)

    if cfg.novel_path:
        Logger.msg("Found", len(result.novel), "novel alleles")
        report.write_novel_fasta(result, cfg.novel_path)

    if cfg.evidence_tsv_path:
        report.write_evidence_tsv(result, cfg.evidence_tsv_path)
        Logger.msg("Wrote evidence table: %s" % cfg.evidence_tsv_path)

    if cfg.html_path:
        opts = report.HtmlOptions(evidence=cfg.html_evidence,
                                  max_novel_bp=cfg.html_max_novel_bp)
        written = report.write_html(result, cfg.html_path, opts,
                                    branding.html_branding())
        # An output-location notice, not a witticism: it survives --quiet.
        Logger.notice("Wrote HTML report: %s" % written)

    failed = [s for s in result.samples if s.failed]
    for sample in failed:
        # Upstream dies on the first bad file; WMLST reports each one and keeps
        # the batch alive (D10). The wording is the engine's, verbatim.
        Logger.notice("ERROR: %s" % (sample.error_text or "analysis failed"))

    Logger.msg("Done.")

    # Any unreadable input is a non-zero exit. WMLST keeps the batch alive where
    # upstream dies on the first bad file (D10), but a pipeline that ran
    # `wmlst *.fna > calls.tsv` must still notice that a sample is missing from
    # the table rather than reading a silent 0 as "all samples typed".
    # Upstream also exits non-zero here, so this is the more compatible choice.
    if failed:
        return 1
    return 0


def _forward(ns, flag, value, extra=()):
    """Build an argv for one of the maintenance entry points (section 4.8)."""
    argv = []
    if value:
        argv += [flag, value]
    if getattr(ns, "quiet", False):
        argv.append("--quiet")
    argv += list(extra)
    return argv


def _bootstrap_blast(ns) -> int:
    """``wmlst --bootstrap-blast``: install NCBI BLAST+ for this user.

    Windows has no package manager we can rely on, so the download is part of
    the product rather than a documented prerequisite.
    """
    from . import blastbin

    def progress(done, total, phase):
        if total:
            Logger.msg("%s: %d%%" % (phase, round(100.0 * done / total)))
        else:
            Logger.msg(phase)

    try:
        tools = blastbin.bootstrap(progress=progress)
    except WmlstError as exc:
        Logger.notice("ERROR: %s" % (exc.user_message or exc))
        return 1
    Logger.notice("Installed blastn %s at %s" % (tools.version, tools.blastn))
    return 0


def _launch_gui(ns) -> int:
    """--gui: hand over to the Tkinter application (section 10)."""
    try:
        from . import gui
    except ImportError as exc:
        Logger.err("The graphical interface is not available: %s" % exc)
        return 1  # pragma: no cover - err() never returns
    return int(gui.main() or 0)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    """Entry point for ``wmlst`` and ``python -m wmlst`` (sections 4.8, 5.2, 14).

    Returns a process exit code and never raises to the interpreter.
    """
    import multiprocessing

    multiprocessing.freeze_support()
    reconfigure_streams()

    if argv is None:
        argv = sys.argv[1:]
    argv = list(argv)

    try:
        parser = build_parser()
        # Upstream calls GetOptions() with no Getopt::Long::Configure, i.e. in
        # default PERMUTE mode: `mlst a.fna --quiet b.fna` collects BOTH files.
        # argparse stops the positional run at the first option, so gather the
        # leftovers and put them back in argv order (finding 16).
        ns, extra = parser.parse_known_args(_expand_abbrevs(argv, parser))
        for token in extra:
            if token.startswith("-") and token != "-":
                # A real unknown option still dies the upstream way.
                parser.error("unrecognized arguments: %s" % token)
        ns.files = list(ns.files) + list(extra)

        Logger.configure(quiet=ns.quiet, debug=ns.debug)
        _banner()
        Logger.dbg("argv=%r" % (argv,))

        if ns.gui:
            return _launch_gui(ns)

        # These three build or fetch what the engine needs, so they run before
        # the dependency check that would otherwise reject a machine with no
        # BLAST+ installed -- which is exactly when --bootstrap-blast is used.
        if ns.bootstrap_blast:
            return _bootstrap_blast(ns)
        if ns.make_blast_db:
            return main_make_blast_db(_forward(ns, "--dbdir", ns.datadir and
                                               os.path.dirname(ns.datadir)))
        if ns.rollback_db:
            return main_update_db(_forward(
                ns, "--dbdir", ns.datadir and os.path.dirname(ns.datadir),
                extra=("--rollback",)))
        if ns.update_db:
            return main_update_db(_forward(
                ns, "--dbdir", ns.datadir and os.path.dirname(ns.datadir),
                extra=("--check",) if ns.check_only else ()))

        _dependency_check(ns)

        cfg = config_from_args(ns, argv)
        Logger.dbg("config=%r" % (cfg,))

        if _listing(ns, cfg):
            return 0

        _preflight(cfg)
        return _run(cfg, argv)
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        sys.stderr.write("%s\n" % code)
        return 1
    except KeyboardInterrupt:
        Logger.notice("ERROR: Cancelled by the user.")
        return 1
    except WmlstError as exc:
        Logger.notice("ERROR: %s" % (exc.user_message or exc))
        return 1
    except BrokenPipeError:  # pragma: no cover - `wmlst --list | head`
        try:
            sys.stdout.close()
        except Exception:
            pass
        return 0
    except OSError as exc:
        Logger.notice("ERROR: %s" % exc)
        return 1


# ---------------------------------------------------------------------------
# The other three console entry points (pyproject [project.scripts])
# ---------------------------------------------------------------------------
SHOW_SEQS_USAGE = """SYNOPSIS
  Print allele seqs for a sequence type in FASTA format
USAGE
  wmlst-show-seqs -s [scheme] -t [ST] > alleles.ffn
OPTIONS
  -v         Show version and exit
  -h         This help
  -s SCHEME  MLST scheme name
  -t ST      Sequence type number
  -d DBDIR   PubMLST folder [<datadir>]
END
"""


def main_show_seqs(argv=None) -> int:
    """``wmlst-show-seqs``: port of ``scripts/mlst-show_seqs`` (bats case 40).

    Writes the allele sequences of one ST to stdout as FASTA; every progress
    line goes to stderr, exactly as the Perl does.
    """
    import getopt

    reconfigure_streams()
    if argv is None:
        argv = sys.argv[1:]
    try:
        opts, _rest = getopt.getopt(list(argv), "vhs:t:d:")
    except getopt.GetoptError as exc:
        sys.stderr.write("Unknown option: %s\n" % exc.opt)
        return 1
    flags = dict(opts)
    if "-v" in flags:
        sys.stderr.write("wmlst-show-seqs %s\n" % __version__)
        return 0
    if "-h" in flags:
        sys.stdout.write(SHOW_SEQS_USAGE)
        return 0

    def fail(text):
        sys.stderr.write("ERROR: %s\n" % text)
        return 1

    scheme = flags.get("-s")
    st = flags.get("-t")
    if not scheme:
        return fail("Please provide -s scheme")
    if not st:
        return fail("Please provide -t ST")
    datadir = flags.get("-d")
    if not datadir:
        from . import schemes

        try:
            datadir = schemes.resolve_datadir(schemes.resolve_dbdir())
        except WmlstError as exc:
            return fail(str(exc))
    if not os.path.isdir(datadir):
        return fail("Bad database dir -d %s" % datadir)
    folder = os.path.join(datadir, scheme)
    if not os.path.isdir(folder):
        return fail("No such folder: %s" % folder)
    if not re.match(r"^\d+$", str(st), re.ASCII):
        return fail("-t '%s' should be an integer" % st)
    profile = os.path.join(folder, scheme + ".txt")
    if not os.access(profile, os.R_OK):
        return fail("No schema file: %s" % profile)
    sys.stderr.write("Schema: %s\n" % profile)
    sys.stderr.write("Target: %s ST%s\n" % (scheme, st))

    genes = sorted(n[:-4] for n in os.listdir(folder) if n.endswith(".tfa"))
    sys.stderr.write("Genes: %s\n" % " ".join(genes))

    header = []
    need = {}
    with open(profile, encoding="utf-8", errors="replace", newline=None) as fh:
        for line in fh:
            cols = line.rstrip("\r\n").split("\t")
            if cols[0] == "ST":
                header = cols[1:len(genes) + 1]
                continue
            if not header:
                continue
            try:
                same = int(cols[0]) == int(st)
            except ValueError:
                same = False
            if same:
                need = {header[i]: cols[i + 1]
                        for i in range(len(header)) if i + 1 < len(cols)}
                break
    if not need:
        return fail("Could not find ST%s in %s" % (st, scheme))

    for gene in sorted(need):
        ident = "%s_%s" % (gene, need[gene])
        sys.stderr.write("Extracting: %s\n" % ident)
        _fasta_extract(os.path.join(folder, gene + ".tfa"), ident)
    sys.stderr.write("Done.\n")
    return 0


def _fasta_extract(path: str, ident: str) -> None:
    """Copy the one record whose id is ``ident`` to stdout, verbatim."""
    try:
        fh = open(path, encoding="utf-8", errors="replace", newline=None)
    except OSError:
        sys.stderr.write("ERROR: Can't open FASTA: %s\n" % path)
        sys.exit(1)
    with fh:
        keep = False
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0:1] == [ident]
            if keep:
                sys.stdout.write(line if line.endswith("\n") else line + "\n")


def main_make_blast_db(argv=None) -> int:
    """``wmlst-make-blast-db``: rebuild ``db/blast/mlst.fa`` (sections 4.9, D14)."""
    reconfigure_streams()
    parser = WmlstArgumentParser(prog="wmlst-make-blast-db", add_help=False)
    parser.add_argument("--help", "-h", action="store_true", dest="show_help")
    parser.add_argument("--dbdir", default=None)
    parser.add_argument("--quiet", "-q", action="store_true")
    ns = parser.parse_args(sys.argv[1:] if argv is None else list(argv))
    if ns.show_help:
        sys.stdout.write(
            "SYNOPSIS\n  Rebuild the WMLST BLAST database from db/pubmlst\n"
            "USAGE\n  wmlst-make-blast-db [--dbdir DIR] [--quiet]\n")
        return 0
    Logger.configure(quiet=ns.quiet)
    from . import schemes, updatedb

    try:
        dbdir = schemes.resolve_dbdir(ns.dbdir)
        Logger.msg("Rebuilding the BLAST database in", os.path.join(dbdir, "blast"))
        path = updatedb.build_blast_db(dbdir, progress=_updatedb_progress)
    except WmlstError as exc:
        sys.stderr.write("ERROR: %s\n" % exc)
        return 1
    Logger.msg("Wrote", path)
    return 0


def main_update_db(argv=None) -> int:
    """``wmlst-update-db``: refresh the PubMLST/Pasteur snapshot (sections 4.9, 7)."""
    reconfigure_streams()
    parser = WmlstArgumentParser(prog="wmlst-update-db", add_help=False)
    parser.add_argument("--help", "-h", action="store_true", dest="show_help")
    parser.add_argument("--dbdir", default=None)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--yes", "-y", action="store_true")
    parser.add_argument("--rollback", action="store_true")
    parser.add_argument("--no-backup", action="store_true")
    parser.add_argument("--write-version-files", action="store_true")
    parser.add_argument("--quiet", "-q", action="store_true")
    ns = parser.parse_args(sys.argv[1:] if argv is None else list(argv))
    if ns.show_help:
        sys.stdout.write(
            "SYNOPSIS\n  Refresh the bundled PubMLST allele database\n"
            "USAGE\n  wmlst-update-db [--check] [--yes] [--dbdir DIR]\n"
            "OPTIONS\n"
            "  --check                Report what changed and exit; write nothing\n"
            "  --yes                  Apply every change without asking\n"
            "  --rollback             Restore the copies the last update kept\n"
            "  --no-backup            Do not keep rollback copies\n"
            "  --write-version-files  Stamp database_version.txt (changes --info DATE)\n"
            "  --quiet                No progress output\n")
        return 0
    Logger.configure(quiet=ns.quiet)
    from . import schemes, updatedb

    try:
        dbdir = schemes.resolve_dbdir(ns.dbdir)
        if ns.rollback:
            Logger.msg("Database rolled back to", updatedb.rollback(dbdir))
            return 0
        Logger.msg("Checking every scheme against the PubMLST/Pasteur APIs...")
        plan = updatedb.check(dbdir, progress=_updatedb_progress)
        chosen = [u.name for u in plan.updates if u.status == "changed"]
        for update in plan.updates:
            if update.status != "up-to-date":
                Logger.msg("  %-28s %-26s %s"
                           % (update.name, update.status, update.detail))
        Logger.msg("%d scheme(s) would change; %d new scheme(s) upstream."
                   % (len(chosen), len(plan.new_schemes)))
        if ns.check:
            return 0
        if not chosen:
            Logger.msg("Nothing to do.")
            return 0
        if not ns.yes:
            sys.stderr.write("Apply these updates? [y/N] ")
            sys.stderr.flush()
            answer = sys.stdin.readline().strip().lower()
            if answer not in ("y", "yes"):
                Logger.msg("Cancelled.")
                return 0
        version = updatedb.apply(dbdir, plan, chosen, backup=not ns.no_backup,
                                 write_version_files=ns.write_version_files,
                                 progress=_updatedb_progress)
    except WmlstError as exc:
        sys.stderr.write("ERROR: %s\n" % exc)
        return 1
    except KeyboardInterrupt:
        sys.stderr.write("ERROR: Cancelled by the user.\n")
        return 1
    Logger.msg("Database updated to", version)
    return 0


def _updatedb_progress(fraction, message=None) -> None:
    """Render updatedb progress as one stderr line (suppressed by --quiet)."""
    if message is None:
        Logger.msg(str(fraction))
    else:
        Logger.msg("[%3.0f%%] %s" % (float(fraction) * 100.0, message))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
