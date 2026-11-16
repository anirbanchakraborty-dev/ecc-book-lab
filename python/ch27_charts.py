"""Write the five computed figures of chapter 27 of Error Correction Coding as SVG files.

    python3 python/ch27_charts.py --out figures

Every drawn value comes from the functions in ch27_bounds_on_codes.py, so a
figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet:
real text in IBM Plex, colors from the page's tokens with paper fallbacks,
and no element that the site's figure check would refuse.
"""

import argparse
import math
import os

from ch27_bounds_on_codes import (
    CODE_5_2_3, D_RANGE, DATA_WIDTHS, EXHIBITED_K, ROWS_7_4_3, dec_bits, elias_k, griesmer_k,
    hamming_k, largest_k, plotkin_best, sec_bits, show, singleton_k, span, varshamov_greedy,
    varshamov_k, weight)

W, H = 640, 400
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


class LinAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo, hi, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (v - self.lo) / (self.hi - self.lo)


def text(x: float, y: float, s: str, anchor: str = "start", fill: str = "", mono: bool = False,
         size: int = 0, bold: bool = False) -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    if bold:
        attrs.append('font-weight="600"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def ring(x: float, y: float, stroke: str = "currentColor", r: float = 5) -> str:
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="none" stroke="{stroke}"'
            f' stroke-width="1.5"/>')


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0,
         dash: str = "") -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def square(x: float, y: float, fill: str, s: float = 7) -> str:
    return rect(x - s / 2, y - s / 2, s, s, fill, "none", 0)


def diamond(x: float, y: float, fill: str, s: float = 5) -> str:
    pts = f"{fmt(x)},{fmt(y - s)} {fmt(x + s)},{fmt(y)} {fmt(x)},{fmt(y + s)} {fmt(x - s)},{fmt(y)}"
    return f'<polygon points="{pts}" fill="{fill}"/>'


def triangle(x: float, y: float, fill: str, s: float = 5) -> str:
    pts = f"{fmt(x)},{fmt(y - s)} {fmt(x + s)},{fmt(y + s * 0.8)} {fmt(x - s)},{fmt(y + s * 0.8)}"
    return f'<polygon points="{pts}" fill="{fill}"/>'


def nabla(x: float, y: float, fill: str, s: float = 5) -> str:
    pts = f"{fmt(x)},{fmt(y + s)} {fmt(x + s)},{fmt(y - s * 0.8)} {fmt(x - s)},{fmt(y - s * 0.8)}"
    return f'<polygon points="{pts}" fill="{fill}"/>'


def polyline(points: list, stroke: str, width: float = 1.5, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return ('<polyline points="' + " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
            + f'" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>')


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


# ------------------------------------------------------------------ figure 2

def singleton_deletion() -> str:
    """The 16 words of the (7,4,3) Hamming code with the last two positions deleted."""
    n, keep = 7, 5
    words = span(ROWS_7_4_3)
    prefixes = {w & (2 ** keep - 1) for w in words}
    if len(words) != 16 or len(prefixes) != 16 or len(prefixes) > 2 ** keep:
        raise SystemExit("deleting two positions of the (7,4,3) code merged two codewords")
    ordered = sorted(words, key=lambda w: show(w, n)[:keep])
    w_, h_ = 640, 352
    cell, row_h, top = 28, 26, 80
    blocks = [(96, ordered[:8]), (380, ordered[8:])]
    body = []
    for x0, block in blocks:
        cut = x0 + keep * cell
        body.append(rect(cut + 2, 40, 2 * cell - 4, 236, BG2, TEXT3, 1, 3, "3 3"))
        body.append(text(cut + cell, 30, "deleted", "middle", TEXT3, size=11))
        for j in range(n):
            body.append(text(x0 + j * cell + cell / 2, 56, str(j + 1), "middle", TEXT3, size=10))
        for i, w in enumerate(block):
            bits = show(w, n)
            y = top + i * row_h + 4
            for j, b in enumerate(bits):
                fill = TEXT3 if j >= keep else ""
                body.append(text(x0 + j * cell + cell / 2, y, b, "middle", fill, mono=True, size=13))
    body.append(text(w_ / 2, 308, "The first five positions still tell all 16 codewords apart,", "middle"))
    body.append(text(w_ / 2, 328, "so 16 words fit among the 2^5 = 32 words of length 5.", "middle"))
    return svg("The 16 codewords of the (7,4,3) Hamming code with their last two positions deleted, "
               "all still distinct", body, w_, h_)


# ------------------------------------------------------------------ figure 3

def plotkin_columns() -> str:
    """The (5,2,3) code, its column weights and both sides of the double count."""
    n, m, d = 5, len(CODE_5_2_3), 3
    cols = [sum(w >> j & 1 for w in CODE_5_2_3) for j in range(n)]
    contrib = [2 * c * (m - c) for c in cols]
    pair_sum = sum(weight(a ^ b) for a in CODE_5_2_3 for b in CODE_5_2_3)
    low, high = m * (m - 1) * d, n * m * m // 2
    if not (low <= pair_sum == sum(contrib) <= high):
        raise SystemExit("the Plotkin double count of the (5,2,3) code does not close")
    w_, h_ = 640, 320
    x0, pitch = 120, 44
    body = []
    body.append(text(96, 52, "position", "end", TEXT3, size=10))
    for j in range(n):
        body.append(text(x0 + j * pitch, 52, str(j + 1), "middle", TEXT3, size=10))
    for i, w in enumerate(CODE_5_2_3):
        y = 84 + i * 28
        body.append(text(96, y, f"c{i + 1}", "end", TEXT3, size=12))
        for j, b in enumerate(show(w, n)):
            body.append(text(x0 + j * pitch, y, b, "middle", mono=True, size=13))
    body.append(line(x0 - 20, 182, x0 + (n - 1) * pitch + 20, 182, "currentColor", 1))
    body.append(text(96, 208, "w", "end", size=12))
    body.append(text(96, 236, "2w(M - w)", "end", size=12))
    for j in range(n):
        body.append(text(x0 + j * pitch, 208, str(cols[j]), "middle", mono=True, size=13))
        body.append(text(x0 + j * pitch, 236, str(contrib[j]), "middle", GREEN, mono=True, size=13))
    body.append(text(96, 268, "column total", "end", size=12))
    body.append(text(x0 + 2 * pitch, 268, " + ".join(str(c) for c in contrib) + f" = {sum(contrib)}",
                     "middle", GREEN, mono=True, size=12))
    px, py, pw, ph = 372, 60, 244, 196
    body.append(rect(px, py, pw, ph, BG2, "none", 0, 3))
    body.append(text(px + 16, py + 26, "The distance sum, two ways", size=12, bold=True))
    body.append(text(px + 16, py + 54, "pairs at distance at least 3", size=12))
    body.append(text(px + 16, py + 74, f"M(M - 1)d = {low}", mono=True, size=12))
    body.append(text(px + 16, py + 102, "column by column, exactly", size=12))
    body.append(text(px + 16, py + 122, f"sum of 2w(M - w) = {pair_sum}", GREEN, mono=True, size=12))
    body.append(text(px + 16, py + 150, "no column beats w = M/2", size=12))
    body.append(text(px + 16, py + 170, f"n M^2 / 2 = {high}", mono=True, size=12))
    body.append(text(w_ / 2, 304, f"{low} is at most {pair_sum}, which is at most {high}, "
                                  f"for M = {m} words of length {n}", "middle", size=12))
    return svg("The four codewords of the (5,2,3) code, the ones in each column and the two sides "
               "of the Plotkin double count", body, w_, h_)


# ------------------------------------------------------------------ figure 5

def varshamov_columns() -> str:
    """The greedy choice of parity-check columns at d = 3 with three check bits."""
    r, d = 3, 3
    steps, chosen = varshamov_greedy(r, d)
    if len(chosen) != 7 or len(steps) != 8:
        raise SystemExit("the greedy construction at d = 3, r = 3 did not stop at length 7")
    w_, h_ = 640, 376
    x0, cw, ch, top, pitch = 96, 52, 24, 64, 28
    body = []
    for c in range(2 ** r):
        bits = "".join(str(c >> (r - 1 - i) & 1) for i in range(r))
        body.append(text(x0 + c * cw + cw / 2, 52, bits, "middle", mono=True, size=12))
    for s, forbidden in enumerate(steps):
        y = top + s * pitch
        body.append(text(88, y + 16, f"step {s + 1}", "end", size=12))
        new = chosen[s] if s < len(chosen) else None
        for c in range(2 ** r):
            x = x0 + c * cw + 4
            if c == new:
                body.append(rect(x, y, cw - 8, ch, GREEN, GREEN, 1.5, 3))
                body.append(text(x + (cw - 8) / 2, y + 16, "take", "middle", BG2, size=11))
            elif c in forbidden:
                body.append(rect(x, y, cw - 8, ch, BG2, TEXT3, 1, 3, "3 3"))
                body.append(text(x + (cw - 8) / 2, y + 16, "x", "middle", TEXT3, mono=True, size=11))
            else:
                body.append(rect(x, y, cw - 8, ch, "none", "currentColor", 1, 3))
        note = f"{len(forbidden)} forbidden" if new is not None else "none left"
        body.append(text(x0 + 8 * cw + 12, y + 16, note, fill=TEXT3 if new is None else "", size=12))
    ly = 312
    body.append(rect(96, ly - 12, 20, 14, GREEN, GREEN, 1.5, 3))
    body.append(text(124, ly, "taken at this step", size=12))
    body.append(rect(260, ly - 12, 20, 14, BG2, TEXT3, 1, 3, "3 3"))
    body.append(text(288, ly, "a sum of at most d - 2 = 1 earlier column", size=12))
    body.append(rect(96, ly + 14, 20, 14, "none", "currentColor", 1, 3))
    body.append(text(124, ly + 26, "still free", size=12))
    body.append(text(w_ / 2, ly + 52, "Seven columns fit, and they are the columns of the (7,4) Hamming code.",
                     "middle", size=12))
    return svg("Choosing three-bit parity-check columns one at a time for distance 3, "
               "until none is left after seven", body, w_, h_)


# ------------------------------------------------------------------ figure 6

def bounds_at_15() -> str:
    """Largest dimension k at length 15 from every bound, against the codes the chapter lists."""
    w_, h_ = 640, 420
    x0, x1, y0, y1 = 84, 604, 24, 344
    ax, ay = LinAxis(2, 15, x0, x1), LinAxis(0, 15, y1, y0 + 20)
    body = [f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}" fill="none" stroke="currentColor"'
            ' stroke-width="1"/>']
    for d in D_RANGE:
        body.append(line(ax(d), y1, ax(d), y1 + 5))
        body.append(text(ax(d), y1 + 19, str(d), "middle", mono=True, size=11))
    for k in range(0, 16):
        body.append(line(x0 - 5, ay(k), x0, ay(k)))
        if k % 3 == 0:
            body.append(text(x0 - 9, ay(k) + 4, str(k), "end", mono=True, size=11))
            if k:
                body.append(line(x0, ay(k), x1, ay(k), TEXT3, 0.5, "3 3"))
    body.append(text((x0 + x1) / 2, h_ - 30, "minimum distance d, at length n = 15", "middle"))
    body.append(f'<text transform="translate(24 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">'
                'largest dimension k</text>')
    series = [
        ("Singleton", [(d, singleton_k(15, d)) for d in D_RANGE], TEXT3, None, 0, "3 3"),
        ("Hamming", [(d, hamming_k(15, d)) for d in D_RANGE], CYAN, square, -6, ""),
        ("Griesmer", [(d, griesmer_k(15, d)) for d in D_RANGE], VIOLET, diamond, -2, ""),
        ("Elias", [(d, elias_k(15, d)) for d in D_RANGE], RED, triangle, 2, ""),
        ("Plotkin", [(d, largest_k(plotkin_best(15, d))) for d in D_RANGE if d >= 8], ORANGE, nabla, 6, ""),
        ("Varshamov", [(d, varshamov_k(15, d)) for d in D_RANGE], GREEN, dot, 0, ""),
    ]
    for name, pts, color, mark, dx, dash in series:
        xy = [(ax(d) + dx, ay(k)) for d, k in pts]
        body.append(polyline(xy, color, 1.5, dash))
        if mark is not None:
            for x, y in xy:
                body.append(mark(x, y, color, 3.5) if mark is dot else mark(x, y, color))
    for d in D_RANGE:
        k = EXHIBITED_K[d]
        upper = min(singleton_k(15, d), hamming_k(15, d), griesmer_k(15, d), elias_k(15, d))
        if not varshamov_k(15, d) <= k <= upper:
            raise SystemExit(f"the listed code at d = {d} is outside the bounds")
        body.append(ring(ax(d), ay(k), "currentColor", 6))
    lx, ly = 396, 46
    body.append(rect(lx - 16, ly - 14, 222, 140, BG2, "none", 0, 3))
    legend = [
        (lambda x, y: line(x - 10, y, x + 10, y, TEXT3, 1.5, "3 3"), "Singleton, k = n - d + 1"),
        (lambda x, y: square(x, y, CYAN), "Hamming, sphere packing"),
        (lambda x, y: diamond(x, y, VIOLET), "Griesmer, linear codes"),
        (lambda x, y: triangle(x, y, RED), "Elias"),
        (lambda x, y: nabla(x, y, ORANGE), "Plotkin, for d at least 8"),
        (lambda x, y: dot(x, y, GREEN, 3.5), "Varshamov, a code exists"),
        (lambda x, y: ring(x, y, "currentColor", 6), "a code listed in the text"),
    ]
    for i, (mark, label) in enumerate(legend):
        y = ly + i * 19
        body.append(mark(lx, y))
        body.append(text(lx + 18, y + 4, label, size=12))
    return svg("The largest dimension k allowed by each bound at length 15, for distances 2 to 15, "
               "with the codes that exist", body, w_, h_)


# ------------------------------------------------------------------ figure 7

def check_bits() -> str:
    """Check bits for single and double error correction against the data width."""
    w_, h_ = 640, 400
    x0, x1, y0, y1 = 84, 560, 24, 324
    ax, ay = LinAxis(2, 8, x0, x1), LinAxis(0, 20, y1, y0)
    body = [f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}" fill="none" stroke="currentColor"'
            ' stroke-width="1"/>']
    for k in DATA_WIDTHS:
        x = ax(math.log2(k))
        body.append(line(x, y1, x, y1 + 5))
        body.append(text(x, y1 + 19, str(k), "middle", mono=True, size=11))
    for r in range(0, 21, 4):
        body.append(line(x0 - 5, ay(r), x0, ay(r)))
        body.append(text(x0 - 9, ay(r) + 4, str(r), "end", mono=True, size=11))
        if r:
            body.append(line(x0, ay(r), x1, ay(r), TEXT3, 0.5, "3 3"))
    body.append(text((x0 + x1) / 2, h_ - 30, "data bits k, on a doubling scale", "middle"))
    body.append(f'<text transform="translate(24 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">'
                'check bits r</text>')
    x64 = ax(6)
    body.append(line(x64, y1, x64, ay(13), TEXT3, 1, "3 3"))
    sec = [sec_bits(k) for k in DATA_WIDTHS]
    secded = [r + 1 for r in sec]
    dec = [dec_bits(k) for k in DATA_WIDTHS]
    series = [("SEC", sec, CYAN, square, ""), ("SECDED", secded, GREEN, dot, ""),
              ("DEC", dec, ORANGE, diamond, "3 3")]
    for name, values, color, mark, dash in series:
        xy = [(ax(math.log2(k)), ay(r)) for k, r in zip(DATA_WIDTHS, values)]
        body.append(polyline(xy, color, 1.5, dash))
        for x, y in xy:
            body.append(mark(x, y, color, 3.5) if mark is dot else mark(x, y, color))
        body.append(text(xy[-1][0] + 12, xy[-1][1] + 4, name, fill=color, size=12))
    at64 = DATA_WIDTHS.index(64)
    lx, ly = 112, 44
    body.append(rect(lx - 18, ly - 16, 372, 76, BG2, "none", 0, 3))
    body.append(square(lx, ly, CYAN))
    body.append(text(lx + 14, ly + 4, "SEC, smallest r with 2^r at least k + r + 1", size=12))
    body.append(dot(lx, ly + 18, GREEN, 3.5))
    body.append(text(lx + 14, ly + 22, "SECDED, one check bit more than SEC", size=12))
    body.append(diamond(lx, ly + 36, ORANGE))
    body.append(text(lx + 14, ly + 40, "DEC, the sphere-packing floor for radius 2", size=12))
    body.append(text(lx - 4, ly + 56, f"At 64 data bits: {sec[at64]}, {secded[at64]} and at least {dec[at64]}",
                     size=12, bold=True))
    return svg("Check bits needed for single error correction, single error correction with double "
               "error detection, and double error correction, against the data width", body, w_, h_)


CHARTS = {
    "fig-27-2-singleton-deletion.svg": singleton_deletion,
    "fig-27-3-plotkin-columns.svg": plotkin_columns,
    "fig-27-5-varshamov-columns.svg": varshamov_columns,
    "fig-27-6-bounds-at-15.svg": bounds_at_15,
    "fig-27-7-check-bits.svg": check_bits,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
