"""Write the four computed figures of chapter 19 of Error Correction Coding as SVG files.

    python3 python/ch19_charts.py --out figures

Every drawn value comes from the functions in
ch19_hamming_codes_and_perfect_codes.py, so a figure is corrected by fixing
those functions and running this script again, never by editing its output.
The SVG follows the book's figure style sheet: real text in IBM Plex, colors
from the page's tokens with paper fallbacks, and no element that the site's
figure check would refuse.
"""

import argparse
import math
import os

from ch19_hamming_codes_and_perfect_codes import (
    C7, H15, RX15, check_positions, flip, index_of, syndrome, word)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
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


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0,
         dash: str = "") -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, stroke: str, label: str) -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.5, 3), text(x + 22, y, label, size=12)]


def hamming_15_11_columns() -> str:
    """H of the (15,11) code, column j holding j in binary, check columns marked."""
    n, m = 15, 4
    checks = set(check_positions(n))
    cell, x0, y0 = 30, 140, 50
    body = []
    body.append(text(x0 - 14, y0 - 12, "position", "end", TEXT3, size=11))
    for j in range(1, n + 1):
        cx = x0 + (j - 1) * cell
        body.append(text(cx + cell / 2, y0 - 12, str(j), "middle", GREEN if j in checks else "",
                         mono=True, size=11))
    for i in range(m):
        cy = y0 + i * cell
        body.append(text(x0 - 14, cy + cell / 2 + 4, f"{2 ** (m - 1 - i)}s bit", "end", size=12))
        for j in range(1, n + 1):
            cx = x0 + (j - 1) * cell
            v = H15[i][j - 1]
            if j in checks:
                body.append(rect(cx, cy, cell, cell, BG2, GREEN, 1.5, 3))
                body.append(text(cx + cell / 2, cy + cell / 2 + 5, str(v), "middle", GREEN, mono=True, size=13))
            else:
                body.append(rect(cx, cy, cell, cell, "none", TEXT3, 1, 3))
                body.append(text(cx + cell / 2, cy + cell / 2 + 5, str(v), "middle", mono=True, size=13))
    ly = y0 + m * cell + 22
    body.append(text(x0 - 14, ly, "bit", "end", TEXT3, size=11))
    data_count = 0
    for j in range(1, n + 1):
        cx = x0 + (j - 1) * cell
        if j in checks:
            body.append(text(cx + cell / 2, ly, f"c{j}", "middle", GREEN, mono=True, size=11))
        else:
            data_count += 1
            body.append(text(cx + cell / 2, ly, f"d{data_count}", "middle", TEXT3, mono=True, size=11))
    body += legend_box(x0, ly + 36, BG2, GREEN, "check position: one 1 in its column, a power of two")
    body += legend_box(x0, ly + 60, "none", TEXT3, "data position: two or more 1s in its column")
    return svg("The parity-check matrix of the (15,11) Hamming code, column j holding j in binary", body,
               640, int(ly + 80))


def syndrome_points_at_error() -> str:
    """The received word of the worked (15,11) example and the four parity checks on it."""
    n, m = 15, 4
    bad = index_of(syndrome(H15, RX15))
    cell, x0, y0, gap = 28, 150, 44, 10
    body = []
    band_x = x0 + (bad - 1) * cell
    rows_top = y0 + cell + 2 * gap
    rows_bottom = rows_top + m * (cell + 6) - 6
    body.append(rect(band_x, y0 - 4, cell, rows_bottom - y0 + 8, BG2, CYAN, 1.5, 4, "3 3"))
    body.append(text(x0 - 14, y0 - 12, "position", "end", TEXT3, size=11))
    for j in range(1, n + 1):
        body.append(text(x0 + (j - 1) * cell + cell / 2, y0 - 12, str(j), "middle",
                         CYAN if j == bad else "", mono=True, size=11))
    body.append(text(x0 - 14, y0 + cell / 2 + 4, "received r", "end", size=12))
    for j in range(1, n + 1):
        cx = x0 + (j - 1) * cell
        v = RX15[j - 1]
        if j == bad:
            body.append(rect(cx + 3, y0, cell - 6, cell, BG2, RED, 2, 3))
            body.append(text(cx + cell / 2, y0 + cell / 2 + 5, str(v), "middle", RED, mono=True, size=13))
        else:
            body.append(rect(cx + 3, y0, cell - 6, cell, "none", TEXT3, 1, 3))
            body.append(text(cx + cell / 2, y0 + cell / 2 + 5, str(v), "middle", mono=True, size=13))
    sx = x0 + n * cell + 24
    body.append(text(sx + 14, rows_top - 10, "parity", "middle", TEXT3, size=11))
    s = syndrome(H15, RX15)
    for i in range(m):
        cy = rows_top + i * (cell + 6)
        body.append(text(x0 - 14, cy + cell / 2 + 4, f"check on the {2 ** (m - 1 - i)}s bit", "end", size=12))
        for j in range(1, n + 1):
            if H15[i][j - 1]:
                cx = x0 + (j - 1) * cell
                body.append(rect(cx + 3, cy + 3, cell - 6, cell - 6, "none", TEXT3, 1, 3))
                body.append(text(cx + cell / 2, cy + cell / 2 + 5, str(RX15[j - 1]), "middle", mono=True, size=12))
        if s[i]:
            body.append(rect(sx, cy, cell, cell, BG2, GREEN, 1.5, 3))
            body.append(text(sx + cell / 2, cy + cell / 2 + 5, "1", "middle", GREEN, mono=True, size=13))
        else:
            body.append(rect(sx, cy, cell, cell, "none", TEXT3, 1, 3))
            body.append(text(sx + cell / 2, cy + cell / 2 + 5, "0", "middle", mono=True, size=13))
    ty = rows_bottom + 34
    body.append(text(x0, ty, f"Read down the parity column: s = {word(s)}, and binary {word(s)} is "
                     f"8 + 2 + 1 = {bad}.", size=12))
    body.append(text(x0, ty + 20, f"The checks that fail are exactly the checks that include position {bad}.",
                     TEXT3, size=12))
    return svg(f"The four parity checks on the received word, whose failures spell position {bad} in binary",
               body, 640, int(ty + 36))


def balls_tiling_length_7() -> str:
    """The 16 balls of radius one around the (7,4) codewords, 128 distinct words."""
    pw, ph, gap, x0, y0 = 146, 92, 6, 19, 12
    body = []
    seen = set()
    for idx, c in enumerate(C7):
        px = x0 + (idx % 4) * (pw + gap)
        py = y0 + (idx // 4) * (ph + gap)
        body.append(rect(px, py, pw, ph, "none", TEXT3, 1, 4))
        body.append(rect(px + 6, py + 5, pw - 12, 20, BG2, GREEN, 1.5, 3))
        body.append(text(px + pw / 2, py + 19, word(c), "middle", GREEN, mono=True, size=12, bold=True))
        seen.add(word(c))
        for j in range(1, 8):
            v = flip(c, j)
            seen.add(word(v))
            col, row = (0, j - 1) if j <= 4 else (1, j - 5)
            tx = px + 14 + col * 66
            ty = py + 40 + row * 14
            bits = word(v)
            spans = "".join(
                f'<tspan fill="{ORANGE}" font-weight="600">{b}</tspan>' if k == j - 1 else b
                for k, b in enumerate(bits))
            body.append(f'<text x="{fmt(tx)}" y="{fmt(ty)}" font-family="{MONO}" font-size="10">{spans}</text>')
    if len(seen) != 128:
        raise SystemExit(f"expected 128 distinct words, found {len(seen)}")
    ly = y0 + 4 * (ph + gap) + 14
    body += legend_box(x0, ly, BG2, GREEN, "codeword, the center of its ball")
    body.append(f'<text x="{x0 + 260}" y="{ly}" font-size="12">its seven neighbors, the changed bit in '
                f'<tspan fill="{ORANGE}" font-weight="600">orange</tspan></text>')
    return svg("The 16 balls of radius one around the (7,4) Hamming codewords, covering all 128 words of length 7",
               body, 640, int(ly + 14))


def rate_and_bound() -> str:
    """Rate of the best binary single-error-correcting code against the sphere-packing bound."""
    nmax = 255
    ax = LinAxis(1.5, math.log2(nmax + 1), X0, X1)
    ay = LinAxis(0.2, 1.0, Y1, Y0)
    body = list(frame_rate())
    for m in range(2, 9):
        n = 2 ** m - 1
        px = ax(m)
        body.append(line(px, Y1, px, Y1 + 4))
        body.append(text(px, Y1 + 18, str(n), "middle", mono=True, size=11))
    for v in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        py = ay(v)
        body.append(line(X0 - 4, py, X0, py))
        body.append(text(X0 - 8, py + 4, f"{v:.1f}", "end", mono=True, size=11))
    bound = [(ax(math.log2(n + 1)), ay(1 - math.log2(n + 1) / n)) for n in range(3, nmax + 1)]
    body.append('<polyline points="' + " ".join(f"{fmt(x)},{fmt(y)}" for x, y in bound)
                + '" fill="none" stroke="currentColor" stroke-width="1.5" stroke-dasharray="3 3"/>')
    stairs = []
    for n in range(3, nmax + 1):
        m = math.ceil(math.log2(n + 1))
        stairs.append((ax(math.log2(n + 1)), ay((n - m) / n)))
    body.append('<polyline points="' + " ".join(f"{fmt(x)},{fmt(y)}" for x, y in stairs)
                + f'" fill="none" stroke="{CYAN}" stroke-width="1.5"/>')
    for m in range(2, 9):
        n, k = 2 ** m - 1, 2 ** m - 1 - m
        px, py = ax(m), ay(k / n)
        if abs((1 - math.log2(n + 1) / n) - k / n) > 1e-12:
            raise SystemExit(f"Hamming code of length {n} does not meet the bound")
        body.append(dot(px, py, GREEN))
        label = f"({n},{k})"
        body.append(text(px - 6, py - 10, label, "end", GREEN, mono=True, size=10))
    lx, ly = 330, 196
    body.append(line(lx, ly, lx + 24, ly, "currentColor", 1.5, "3 3"))
    body.append(text(lx + 32, ly + 4, "sphere-packing bound, 1 - log2(n + 1)/n", size=12))
    body.append(line(lx, ly + 24, lx + 24, ly + 24, CYAN, 1.5))
    body.append(text(lx + 32, ly + 28, "best linear code with distance 3", size=12))
    body.append(dot(lx + 12, ly + 48, GREEN))
    body.append(text(lx + 32, ly + 52, "Hamming codes, on the bound", size=12))
    return svg("The rate of single-error-correcting codes against the sphere-packing bound, lengths 3 to 255", body)


def frame_rate() -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 22, "code length n, plotted at log2(n + 1)", "middle"),
        f'<text transform="translate(24 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">rate k/n</text>',
    ]


CHARTS = {
    "fig-19-1-hamming-15-11-columns.svg": hamming_15_11_columns,
    "fig-19-2-syndrome-points-at-error.svg": syndrome_points_at_error,
    "fig-19-3-balls-tiling-length-7.svg": balls_tiling_length_7,
    "fig-19-7-rate-and-bound.svg": rate_and_bound,
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
