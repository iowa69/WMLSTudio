#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025-2026 IOWA-Tech - Giovanni Lorenzin
# Copyright (C) Torsten Seemann (upstream `mlst`, from which WMLST is ported)
"""Generate ``packaging/wmlst.ico`` with nothing but the standard library.

Implements docs/ARCHITECTURE.md section 13.1 (packaging assets) and the
zero-dependency rule of section 13.3: Pillow is an optional *build* extra at
most and is explicitly excluded from the frozen bundle, so the application icon
cannot be produced by it. Everything here is ``zlib`` + ``struct`` + arithmetic.

What it draws
-------------
A rounded square with an indigo-to-teal vertical gradient and a white ``W``
drawn as a thick round-capped polyline. The mark is deliberately a single heavy
glyph: at 16x16 in the Windows taskbar anything finer turns to mush.

How it renders
--------------
The artwork is rasterised once at 1024x1024 with a boolean inside/outside test
per pixel, then box-downsampled to each icon size in premultiplied alpha. Box
downsampling from a single master is both the anti-aliasing and the guarantee
that every size in the file is the same picture.

Container format
----------------
A real multi-size ``.ico``::

    ICONDIR       6 bytes    reserved=0, type=1, count=N
    ICONDIRENTRY  16 bytes   x N, sorted smallest-first
    payloads                 BMP/DIB for 16..128, PNG for 256

Sizes 16..128 are stored as 32-bit bottom-up DIBs with the (mandatory, even
when unused) 1-bit AND mask; 256x256 is stored PNG-compressed, which is the
Vista+ convention and keeps the file around a tenth of the size it would
otherwise be.

Usage::

    python3 packaging/make_icon.py [-o packaging/wmlst.ico] [--verify]
"""

from __future__ import annotations

import argparse
import os
import struct
import sys
import zlib

#: Sizes the icon must contain. Asserted by tests/test_packaging.py.
ICON_SIZES = (16, 24, 32, 48, 64, 128, 256)

#: Sizes at or above this are stored PNG-compressed inside the ICO.
PNG_FROM = 256

#: Master raster edge. 1536 = 2^9 * 3, the smallest comfortable common multiple
#: of every size above -- 24 and 48 are not powers of two, so the obvious 1024
#: cannot divide them and would silently chamfer those two entries.
MASTER = 1536

#: Gradient stops, top to bottom (R, G, B).
TOP_COLOR = (0x1B, 0x2A, 0x6B)      # deep indigo
BOTTOM_COLOR = (0x0D, 0x94, 0x88)   # teal

#: Glyph colour.
INK = (0xFF, 0xFF, 0xFF)

#: Corner radius and glyph geometry, as fractions of the edge.
CORNER_RADIUS = 0.185
STROKE = 0.118
GLYPH = (
    (0.215, 0.285),
    (0.350, 0.730),
    (0.500, 0.455),
    (0.650, 0.730),
    (0.785, 0.285),
)


def _rounded_rect_mask(n, radius_frac=CORNER_RADIUS, inset_frac=0.0):
    """Return a bytearray of n*n 0/1 flags: 1 where the rounded square covers.

    The test is the standard squared-distance-to-corner-centre one, so the
    corners are true circular arcs rather than the chamfers a naive
    implementation produces.
    """
    inset = inset_frac * n
    lo = inset
    hi = n - inset
    radius = radius_frac * n
    rx0 = lo + radius
    rx1 = hi - radius
    r2 = radius * radius
    mask = bytearray(n * n)
    for y in range(n):
        py = y + 0.5
        if py < lo or py > hi:
            continue
        row = y * n
        if py < rx0:
            cy = rx0
        elif py > rx1:
            cy = rx1
        else:
            cy = None
        for x in range(n):
            px = x + 0.5
            if px < lo or px > hi:
                continue
            if cy is None:
                mask[row + x] = 1
                continue
            if px < rx0:
                cx = rx0
            elif px > rx1:
                cx = rx1
            else:
                mask[row + x] = 1
                continue
            dx = px - cx
            dy = py - cy
            if dx * dx + dy * dy <= r2:
                mask[row + x] = 1
    return mask


def _polyline_mask(n, points, stroke_frac=STROKE):
    """Return a bytearray of n*n 0/1 flags: 1 within stroke/2 of the polyline.

    Distance-to-segment gives round caps and round joins for free, which is
    exactly what the glyph wants and what a scanline polygon filler would make
    hard work of.
    """
    half = stroke_frac * n / 2.0
    h2 = half * half
    segs = []
    for i in range(len(points) - 1):
        ax, ay = points[i][0] * n, points[i][1] * n
        bx, by = points[i + 1][0] * n, points[i + 1][1] * n
        dx, dy = bx - ax, by - ay
        length2 = dx * dx + dy * dy
        # Bounding box of the segment, fattened by the stroke radius, so the
        # inner loop only visits pixels that can possibly be covered.
        x0 = max(0, int(min(ax, bx) - half) - 1)
        x1 = min(n - 1, int(max(ax, bx) + half) + 1)
        y0 = max(0, int(min(ay, by) - half) - 1)
        y1 = min(n - 1, int(max(ay, by) + half) + 1)
        segs.append((ax, ay, dx, dy, length2 or 1.0, x0, x1, y0, y1))

    mask = bytearray(n * n)
    for ax, ay, dx, dy, length2, x0, x1, y0, y1 in segs:
        for y in range(y0, y1 + 1):
            py = y + 0.5
            row = y * n
            for x in range(x0, x1 + 1):
                if mask[row + x]:
                    continue
                px = x + 0.5
                t = ((px - ax) * dx + (py - ay) * dy) / length2
                if t < 0.0:
                    t = 0.0
                elif t > 1.0:
                    t = 1.0
                ex = px - (ax + t * dx)
                ey = py - (ay + t * dy)
                if ex * ex + ey * ey <= h2:
                    mask[row + x] = 1
    return mask


def render_master(n=MASTER):
    """Rasterise the artwork at n x n as premultiplied RGBA bytes.

    Returns a bytearray of length n*n*4 in R,G,B,A order. Premultiplied,
    because averaging straight-alpha pixels across a transparent edge drags
    black into the fringe.
    """
    plate = _rounded_rect_mask(n)
    glyph = _polyline_mask(n, GLYPH)
    buf = bytearray(n * n * 4)
    tr, tg, tb = TOP_COLOR
    br, bg, bb = BOTTOM_COLOR
    ir, ig, ib = INK
    for y in range(n):
        f = (y + 0.5) / n
        gr = int(tr + (br - tr) * f + 0.5)
        gg = int(tg + (bg - tg) * f + 0.5)
        gb = int(tb + (bb - tb) * f + 0.5)
        row = y * n
        for x in range(n):
            i = row + x
            if not plate[i]:
                continue
            o = i * 4
            if glyph[i]:
                buf[o] = ir
                buf[o + 1] = ig
                buf[o + 2] = ib
            else:
                buf[o] = gr
                buf[o + 1] = gg
                buf[o + 2] = gb
            buf[o + 3] = 255
    return buf


def downsample(src, src_n, dst_n):
    """Box-average a premultiplied RGBA buffer from src_n to dst_n.

    ``src_n`` must be an exact multiple of ``dst_n``; the caller guarantees it
    by choosing MASTER as a common multiple of every icon size.
    """
    if src_n % dst_n:
        raise ValueError("%d is not an integer multiple of %d" % (src_n, dst_n))
    block = src_n // dst_n
    area = block * block
    out = bytearray(dst_n * dst_n * 4)
    for oy in range(dst_n):
        base_y = oy * block
        orow = oy * dst_n * 4
        for ox in range(dst_n):
            base_x = ox * block
            r = g = b = a = 0
            for by in range(block):
                o = ((base_y + by) * src_n + base_x) * 4
                for _ in range(block):
                    r += src[o]
                    g += src[o + 1]
                    b += src[o + 2]
                    a += src[o + 3]
                    o += 4
            o = orow + ox * 4
            out[o] = (r + area // 2) // area
            out[o + 1] = (g + area // 2) // area
            out[o + 2] = (b + area // 2) // area
            out[o + 3] = (a + area // 2) // area
    return out


def _downsample_chain(art, master, sizes):
    """Box-downsample the master to every requested size, cheaply.

    A box average of box averages over exact integer blocks is the same box
    average over the larger block, so each size is derived from the smallest
    already-computed raster that is an exact multiple of it. This turns seven
    full-master passes into roughly one and a half.
    """
    have = {master: art}
    out = {}
    for size in sorted(sizes, reverse=True):
        best = min(
            (n for n in have if n % size == 0 and n != size),
            key=lambda n: n,
            default=None,
        )
        if best is None:
            raise ValueError("no raster is an integer multiple of %d" % size)
        out[size] = downsample(have[best], best, size)
        have[size] = out[size]
    return out


def unpremultiply(buf):
    """Convert premultiplied RGBA to straight RGBA in place, returning buf."""
    for o in range(0, len(buf), 4):
        a = buf[o + 3]
        if a == 0:
            buf[o] = buf[o + 1] = buf[o + 2] = 0
        elif a < 255:
            for k in range(3):
                v = (buf[o + k] * 255 + a // 2) // a
                buf[o + k] = 255 if v > 255 else v
    return buf


def _png_chunk(tag, payload):
    """Return one length-tag-payload-CRC PNG chunk."""
    body = tag + payload
    return struct.pack(">I", len(payload)) + body + struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF)


def encode_png(rgba, n):
    """Encode straight-alpha RGBA bytes as a PNG (colour type 6, 8 bit)."""
    raw = bytearray()
    stride = n * 4
    for y in range(n):
        raw.append(0)  # filter type 0 (None): the artwork is flat, so filtering buys little
        raw += rgba[y * stride:(y + 1) * stride]
    out = bytearray(b"\x89PNG\r\n\x1a\n")
    out += _png_chunk(b"IHDR", struct.pack(">IIBBBBB", n, n, 8, 6, 0, 0, 0))
    out += _png_chunk(b"IDAT", zlib.compress(bytes(raw), 9))
    out += _png_chunk(b"IEND", b"")
    return bytes(out)


def encode_dib(rgba, n):
    """Encode straight-alpha RGBA as the BMP/DIB payload an ICO entry expects.

    A BITMAPINFOHEADER whose biHeight is *twice* the image height, then the
    bottom-up 32-bit BGRA colour data, then the 1-bit AND mask. The mask is
    unused for 32-bit entries but is not optional: omit it and Windows renders
    a garbage bottom half on some shell surfaces.
    """
    header = struct.pack(
        "<IiiHHIIiiII",
        40,        # biSize
        n,         # biWidth
        n * 2,     # biHeight: colour data plus AND mask
        1,         # biPlanes
        32,        # biBitCount
        0,         # biCompression = BI_RGB
        n * n * 4, # biSizeImage
        0, 0, 0, 0,
    )
    colour = bytearray()
    for y in range(n - 1, -1, -1):
        o = y * n * 4
        for _ in range(n):
            colour.append(rgba[o + 2])  # B
            colour.append(rgba[o + 1])  # G
            colour.append(rgba[o])      # R
            colour.append(rgba[o + 3])  # A
            o += 4
    mask_stride = ((n + 31) // 32) * 4
    mask = bytearray(mask_stride * n)
    return bytes(header) + bytes(colour) + bytes(mask)


def build_ico(sizes=ICON_SIZES, master=MASTER):
    """Render every size and return the complete ``.ico`` file as bytes."""
    art = render_master(master)
    rendered = _downsample_chain(art, master, sorted(sizes))
    payloads = []
    for size in sorted(sizes):
        rgba = unpremultiply(bytearray(rendered[size]))
        if size >= PNG_FROM:
            payloads.append((size, encode_png(rgba, size)))
        else:
            payloads.append((size, encode_dib(rgba, size)))

    count = len(payloads)
    offset = 6 + 16 * count
    directory = bytearray(struct.pack("<HHH", 0, 1, count))
    for size, data in payloads:
        dim = 0 if size >= 256 else size  # 0 is how an ICO spells 256
        directory += struct.pack(
            "<BBBBHHII", dim, dim, 0, 0, 1, 32, len(data), offset
        )
        offset += len(data)
    return bytes(directory) + b"".join(data for _, data in payloads)


def verify_ico(path):
    """Parse an ``.ico`` and return the sorted tuple of sizes it contains.

    Raises ValueError on anything structurally wrong. Used by
    tests/test_packaging.py and by ``--verify``.
    """
    with open(path, "rb") as fh:
        blob = fh.read()
    if len(blob) < 6:
        raise ValueError("%s: too short to be an ICO" % path)
    reserved, kind, count = struct.unpack_from("<HHH", blob, 0)
    if reserved != 0:
        raise ValueError("%s: ICONDIR.reserved is %d, must be 0" % (path, reserved))
    if kind != 1:
        raise ValueError("%s: ICONDIR.type is %d, must be 1 (icon)" % (path, kind))
    if count == 0:
        raise ValueError("%s: contains no images" % path)
    if len(blob) < 6 + 16 * count:
        raise ValueError("%s: directory is truncated" % path)

    found = []
    for i in range(count):
        (w, h, colours, res, planes, bits, nbytes, off) = struct.unpack_from(
            "<BBBBHHII", blob, 6 + 16 * i
        )
        width = 256 if w == 0 else w
        height = 256 if h == 0 else h
        if width != height:
            raise ValueError("%s: entry %d is %dx%d, not square" % (path, i, width, height))
        if res != 0:
            raise ValueError("%s: entry %d reserved byte is %d" % (path, i, res))
        if off + nbytes > len(blob):
            raise ValueError("%s: entry %d payload runs past end of file" % (path, i))
        if nbytes == 0:
            raise ValueError("%s: entry %d has an empty payload" % (path, i))
        payload = blob[off:off + nbytes]
        if payload[:8] == b"\x89PNG\r\n\x1a\n":
            png_w, png_h = struct.unpack_from(">II", payload, 16)
            if (png_w, png_h) != (width, height):
                raise ValueError(
                    "%s: entry %d claims %dx%d but the PNG is %dx%d"
                    % (path, i, width, height, png_w, png_h)
                )
        else:
            bi_size, bi_w, bi_h, _bi_planes, bi_bits = struct.unpack_from("<IiiHH", payload, 0)
            if bi_size != 40:
                raise ValueError("%s: entry %d BITMAPINFOHEADER size is %d" % (path, i, bi_size))
            if bi_w != width or bi_h != height * 2:
                raise ValueError(
                    "%s: entry %d DIB is %dx%d, expected %dx%d (height doubled for the AND mask)"
                    % (path, i, bi_w, bi_h, width, height * 2)
                )
            if bi_bits != 32:
                raise ValueError("%s: entry %d is %d bpp, expected 32" % (path, i, bi_bits))
            mask_stride = ((width + 31) // 32) * 4
            expected = 40 + width * height * 4 + mask_stride * height
            if len(payload) != expected:
                raise ValueError(
                    "%s: entry %d payload is %d bytes, expected %d"
                    % (path, i, len(payload), expected)
                )
        del colours, planes, bits
        found.append(width)

    if len(set(found)) != len(found):
        raise ValueError("%s: duplicate sizes %r" % (path, found))
    return tuple(sorted(found))


def main(argv=None):
    """Command-line entry point. Returns a process exit code."""
    default = os.path.join(os.path.dirname(os.path.abspath(__file__)), "wmlst.ico")
    parser = argparse.ArgumentParser(
        prog="make_icon.py",
        description="Generate the WMLST multi-size Windows icon using only the stdlib.",
    )
    parser.add_argument("-o", "--output", default=default, help="output .ico path")
    parser.add_argument("--verify", action="store_true", help="only validate an existing .ico")
    parser.add_argument("--png", metavar="DIR", help="also write each size as a .png into DIR")
    ns = parser.parse_args(argv)

    if ns.verify:
        sizes = verify_ico(ns.output)
        sys.stdout.write("%s: valid ICO, sizes %s\n" % (ns.output, ", ".join(map(str, sizes))))
        return 0

    blob = build_ico()
    directory = os.path.dirname(os.path.abspath(ns.output))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(ns.output, "wb") as fh:
        fh.write(blob)

    if ns.png:
        if not os.path.isdir(ns.png):
            os.makedirs(ns.png)
        art = render_master(MASTER)
        for size in ICON_SIZES:
            rgba = unpremultiply(downsample(art, MASTER, size))
            with open(os.path.join(ns.png, "wmlst-%d.png" % size), "wb") as fh:
                fh.write(encode_png(rgba, size))

    sizes = verify_ico(ns.output)
    sys.stdout.write(
        "wrote %s (%d bytes), sizes %s\n"
        % (ns.output, len(blob), ", ".join(map(str, sizes)))
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
