# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-BioTech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Generate the WMLST window icon as a PNG, with nothing but the stdlib.

Run ``python -m wmlst.assets.make_icon`` (or execute this file) to regenerate
``wmlst.png`` after changing the palette. The artwork is drawn arithmetically so
it stays crisp at every size and needs no image library at build or run time.

Implements the asset side of docs/ARCHITECTURE.md section 10.1 (no bitmaps in the
UI itself; the window icon is the single exception, because Tk requires one).
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

#: Accent blue, matching gui.PALETTE["accent"].
ACCENT = (0x1F, 0x6F, 0xEB)
INK = (0xFF, 0xFF, 0xFF)
SIZES = (16, 32, 48, 64, 128, 256)


def _rounded(x: int, y: int, size: int, radius: float) -> bool:
    """True when the pixel is inside a rounded square of side ``size``."""
    cx = min(max(x + 0.5, radius), size - radius)
    cy = min(max(y + 0.5, radius), size - radius)
    return ((x + 0.5 - cx) ** 2 + (y + 0.5 - cy) ** 2) <= radius * radius


def render(size: int) -> bytes:
    """Return raw RGBA rows for one icon size.

    The mark is four allele bars of different lengths — a typing profile, which
    is what the product actually does — on a rounded accent tile.
    """
    radius = size * 0.22
    bar_left = int(round(size * 0.22))
    bar_h = max(1, int(round(size * 0.10)))
    gap = max(1, int(round(size * 0.07)))
    block = bar_h + gap
    top = int(round((size - (4 * block - gap)) / 2.0))
    lengths = (0.56, 0.38, 0.50, 0.30)

    rows = bytearray()
    for y in range(size):
        rows.append(0)  # PNG filter type 0 for this scanline
        for x in range(size):
            if not _rounded(x, y, size, radius):
                rows.extend((0, 0, 0, 0))
                continue
            pixel = ACCENT
            for i, frac in enumerate(lengths):
                y0 = top + i * block
                if y0 <= y < y0 + bar_h:
                    if bar_left <= x < bar_left + int(round(size * frac)):
                        pixel = INK
                    break
            rows.extend((pixel[0], pixel[1], pixel[2], 255))
    return bytes(rows)


def _chunk(tag: bytes, payload: bytes) -> bytes:
    return (struct.pack(">I", len(payload)) + tag + payload
            + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF))


def write_png(path: Path, size: int = 64) -> Path:
    """Write a single-size RGBA PNG and return the path."""
    header = struct.pack(">IIBBBBB", size, size, 8, 6, 0, 0, 0)
    data = zlib.compress(render(size), 9)
    png = (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", header)
           + _chunk(b"IDAT", data) + _chunk(b"IEND", b""))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(png)
    return path


def main() -> int:
    """Regenerate every bundled icon size."""
    here = Path(__file__).resolve().parent
    write_png(here / "wmlst.png", 64)
    write_png(here / "wmlst-256.png", 256)
    print("wrote {} and {}".format(here / "wmlst.png", here / "wmlst-256.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
