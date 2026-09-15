"""How far a long step has got, and roughly how much is left, in plain words.

Every long-running step in this application had the same defect in a different
dialog: it named the item it had just finished and nothing else, so ten seconds
in and ten minutes in looked identical, and a download that was working perfectly
was read as a frozen one. These are the words all of them use now, so "is this
still doing something" has the same answer wherever it is asked.

An estimate is offered only once enough of the work has been done for it to mean
anything, and it is always called an estimate.
"""

from __future__ import annotations

#: Below this, a percentage has not seen enough of the work to predict the rest.
ESTIMATE_PERCENT = 10
#: Below this many seconds, the clock has not run long enough to divide by.
ESTIMATE_SECONDS = 5
#: How many items must be finished before a per-item rate is worth quoting.
ESTIMATE_ITEMS = 12


def duration_words(seconds) -> str:
    """A number of seconds as somebody would say it out loud."""
    seconds = max(0, int(seconds))
    if seconds < 60:
        return f"{seconds} seconds"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes} min" if not seconds else f"{minutes} min {seconds} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


def describe_progress(current, total, text, elapsed) -> str:
    """Say where a counted step has got to, and roughly how much is left."""
    if total <= 0 or current <= 0:
        return f"{text} · {duration_words(elapsed)} so far" if elapsed >= ESTIMATE_SECONDS else text
    current = min(current, total)
    percent = int(current * 100 / total)
    detail = f"{current:,} of {total:,} ({percent}%)"
    if current >= ESTIMATE_ITEMS and elapsed >= 2 * ESTIMATE_SECONDS and current < total:
        remaining = elapsed / current * (total - current)
        return f"{text} · {detail} · about {duration_words(remaining)} left"
    return f"{text} · {detail}"


def remaining_words(percent, elapsed) -> str:
    """How much longer, for a step that reports only a percentage.

    Empty until the guess would mean something, because "about 4 h left" read off
    the first two percent of a two-minute download is worse than saying nothing.
    """
    try:
        percent, elapsed = float(percent), float(elapsed)
    except (TypeError, ValueError):
        return ""
    if percent < ESTIMATE_PERCENT or percent >= 100 or elapsed < ESTIMATE_SECONDS:
        return ""
    return f"about {duration_words(elapsed * (100 - percent) / percent)} left"
