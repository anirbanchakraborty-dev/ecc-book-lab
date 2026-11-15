"""Write the five computed figures of chapter 13 of Error Correction Coding as SVG files.

    python3 python/ch13_charts.py
    python3 python/ch13_charts.py --out figures

Every value drawn comes from the functions in ch13_integers_and_modular_arithmetic.py:
the squares that tile the 1147 by 899 rectangle are the divisions of Euclid's
algorithm, the remainder curves are the remainders of the same algorithm, the
arrows of the permutation are the products 3k mod 7, the grid holds each number
below 15 at its two remainders, and the two multiplication tables are computed
cell by cell. A figure is corrected by fixing those functions and running this
script again, never by editing its output. Inside the hub that holds both this
repository and the book's website, the default output folder is the site's
figure folder for the book. In a clone of this repository alone it is a folder
named figures in the working directory. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import hypot, log10

from ch13_integers_and_modular_arithmetic import euclid_steps, mult_table

W = 640
GREEN = "var(--green,#00614a)"
RED = "var(--red,#a1211f)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker>'
         '<marker id="ahg" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="var(--green,#00614a)"/></marker></defs>')


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def text(x: float, y: float, s: str, anchor: str = "start", fill: str = "", mono: bool = False,
         size: int = 0, weight_: str = "") -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    if weight_:
        attrs.append(f'font-weight="{weight_}"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", arrow: bool = False) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{"ahg" if stroke == GREEN else "ah"})"' if arrow else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{m}/>')


def rect(x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, dash: str = "", rx: float = 0) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}"'
            f' stroke="{stroke}" stroke-width="{width}"{d}/>')


def circle(x: float, y: float, r: float, fill: str = BG2, stroke: str = "currentColor", width: float = 1.5) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def dot(x: float, y: float, fill: str = "currentColor", r: float = 3.5) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, h: int, body: list, w: int = W) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


# Euclid's algorithm as a tiling ------------------------------------------------------

def tiling(a: int, b: int) -> list:
    """The squares Euclid's algorithm cuts from an a by b rectangle, as [x, y, side, step] in units.

    Each division q = r_(i-1) // r_i cuts q squares of side r_i from the strip
    that is left, alternately along its width and along its height, so the
    squares of the last division have side gcd(a, b) and fill what remains.
    """
    squares, x, y, w, h = [], 0, 0, a, b
    for step, (_, q, side, _) in enumerate(euclid_steps(a, b), start=1):
        for k in range(q):
            if w >= h:
                squares.append([x + k * side, y, side, step])
            else:
                squares.append([x, y + k * side, side, step])
        if w >= h:
            x, w = x + q * side, w - q * side
        else:
            y, h = y + q * side, h - q * side
    return squares


def euclid_tiling() -> str:
    a, b, scale = 1147, 899, 0.5
    x0, y0 = (W - a * scale) / 2, 52
    body = [text(W / 2, 30, "a = 1147 across, b = 899 down", "middle", mono=True, size=12)]
    last = max(s[3] for s in tiling(a, b))
    for x, y, side, step in tiling(a, b):
        px, py, ps = x0 + x * scale, y0 + y * scale, side * scale
        final = step == last
        body.append(rect(px, py, ps, ps, fill=BG2, stroke=GREEN if final else "currentColor",
                         width=1.5 if final else 1))
        size = 14 if ps > 100 else 12 if ps > 60 else 10 if ps > 25 else 9
        body.append(text(px + ps / 2, py + ps / 2 + size * 0.35, str(side), "middle",
                         fill=GREEN if final else "", mono=True, size=size, weight_="600" if final else ""))
    body.append(rect(x0, y0, a * scale, b * scale, width=1.5))
    h = round(y0 + b * scale + 40)
    body.append(text(x0, h - 14, "Squares of one size come from one division. Two green squares of side 31 end it.",
                     fill=TEXT3, size=12))
    return svg("Euclid's algorithm on 1147 and 899 drawn as squares cut from a rectangle", h, body)


# The remainders against the halving bound --------------------------------------------

def remainders(a: int, b: int) -> list:
    return [a, b] + [s[3] for s in euclid_steps(a, b) if s[3] > 0]


def euclid_remainders() -> str:
    h = 400
    x0, x1, y0, y1 = 84, 604, 40, 330
    n_max, lo, hi = 18, 1, 10000

    def px(i: float) -> float:
        return x0 + (x1 - x0) * i / n_max

    def py(v: float) -> float:
        return y1 - (y1 - y0) * (log10(v) - log10(lo)) / (log10(hi) - log10(lo))

    body = [ARROW]
    for k in range(5):
        v = 10 ** k
        body.append(line(x0, py(v), x1, py(v), TEXT3, 0.5) if k else line(x0, py(v), x1, py(v), width=1))
        body.append(text(x0 - 8, py(v) + 4, str(v), "end", mono=True, size=11))
    body.append(line(x0, y0, x0, y1, width=1))
    for i in range(0, n_max + 1, 2):
        body.append(line(px(i), y1, px(i), y1 + 5, width=1))
        body.append(text(px(i), y1 + 19, str(i), "middle", mono=True, size=11))
    body.append(text((x0 + x1) / 2, y1 + 40, "index i of the remainder r_i", "middle", size=12))
    body.append(text(x0 - 8, y0 - 18, "r_i", "end", mono=True, size=12))

    runs = [(4181, 2584, GREEN), (1147, 899, "currentColor")]
    for a, b, color in runs:
        rs = remainders(a, b)
        env = [(i, b / 2 ** ((i - 1) // 2)) for i in range(1, n_max + 1)]
        body.append(polyline([(px(i), py(v)) for i, v in env], TEXT3, 1, "3 3"))
        body.append(polyline([(px(i), py(v)) for i, v in enumerate(rs)], color, 1.5))
        body.extend(dot(px(i), py(v), color) for i, v in enumerate(rs))
    body.append(text(px(15) + 16, y1 - 8, "gcd 1 after 17 divisions", "end", fill=GREEN, size=12))
    body.append(text(px(6) - 8, py(31) + 20, "gcd 31 after 6 divisions", "end", size=12))
    lx, ly = 360, 62
    body += [line(lx, ly - 4, lx + 24, ly - 4, GREEN, 2), text(lx + 32, ly, "4181, 2584 (Fibonacci)", size=12),
             line(lx, ly + 16, lx + 24, ly + 16, "currentColor", 2), text(lx + 32, ly + 20, "1147, 899", size=12),
             line(lx, ly + 36, lx + 24, ly + 36, TEXT3, 1, "3 3"),
             text(lx + 32, ly + 40, "halving bound, b halved every two steps", size=12)]
    return svg("Remainders of Euclid's algorithm on a log scale against the halving bound", h, body)


# Multiplying by 3 permutes the nonzero residues modulo 7 -----------------------------

def fermat_permutation() -> str:
    h = 230
    ytop, ybot, r = 62, 172, 16

    def x(k: int) -> float:
        return 150 + (k - 1) * 80

    body = [ARROW, text(24, ytop + 4, "k", mono=True, size=13), text(24, ybot + 4, "3k mod 7", mono=True, size=13)]
    for k in range(1, 7):
        target = 3 * k % 7
        color = GREEN if target == 1 else "currentColor"
        dx, dy = x(target) - x(k), ybot - ytop
        d = hypot(dx, dy)
        sx, sy = x(k) + dx * r / d, ytop + dy * r / d
        ex, ey = x(target) - dx * (r + 3) / d, ybot - dy * (r + 3) / d
        body.append(line(sx, sy, ex, ey, color, 1.5 if target == 1 else 1, arrow=True))
    for k in range(1, 7):
        body.append(circle(x(k), ytop, r, stroke=GREEN if k == 5 else "currentColor"))
        body.append(text(x(k), ytop + 5, str(k), "middle", mono=True, size=13))
        body.append(circle(x(k), ybot, r, stroke=GREEN if k == 1 else "currentColor"))
        body.append(text(x(k), ybot + 5, str(k), "middle", mono=True, size=13))
    body.append(text(W / 2, 212, "Every residue from 1 to 6 is hit once. The green arrow says 3 times 5 is 1 mod 7.",
                     "middle", fill=TEXT3, size=12))
    return svg("Multiplication by 3 sends each nonzero residue modulo 7 to a different one", h, body)


# The Chinese remainder grid ----------------------------------------------------------

def crt_grid() -> str:
    m, n, target = 3, 5, 8
    cw, ch, x0, y0 = 70, 44, 200, 76
    h = y0 + m * ch + 52
    body = [text(x0 + n * cw / 2, 30, "x mod 5", "middle", mono=True, size=12),
            text(x0 - 40, y0 - 12, "x mod 3", "middle", mono=True, size=12)]
    for j in range(n):
        body.append(text(x0 + j * cw + cw / 2, y0 - 12, str(j), "middle", fill=TEXT3, mono=True, size=12))
    for i in range(m):
        body.append(text(x0 - 40, y0 + i * ch + ch / 2 + 4, str(i), "middle", fill=TEXT3, mono=True, size=12))
    for xv in range(m * n):
        i, j = xv % m, xv % n
        cx, cy = x0 + j * cw, y0 + i * ch
        hit = xv == target
        body.append(rect(cx, cy, cw, ch, fill=BG2 if hit else "none", stroke=GREEN if hit else "currentColor",
                         width=2 if hit else 1))
        body.append(text(cx + cw / 2, cy + ch / 2 + 5, str(xv), "middle", fill=GREEN if hit else "",
                         mono=True, size=14, weight_="600" if hit else ""))
    body.append(text(W / 2, h - 14, "Each of the 15 cells holds exactly one x from 0 to 14. Row 2, column 3 holds 8.",
                     "middle", fill=TEXT3, size=12))
    return svg("The numbers 0 to 14 placed by their remainders modulo 3 and modulo 5", h, body)


# The multiplication tables of Z_7 and Z_6 --------------------------------------------

def table(n: int, x0: float, y0: float, c: float) -> list:
    t = mult_table(n)
    body = [text(x0 + (n + 1) * c / 2, y0 - 14, f"multiplication mod {n}", "middle", size=13, weight_="600"),
            text(x0 + c / 2, y0 + c / 2 + 4, "x", "middle", fill=TEXT3, mono=True, size=12)]
    for k in range(n):
        body.append(text(x0 + (k + 1.5) * c, y0 + c / 2 + 4, str(k), "middle", fill=TEXT3, mono=True, size=12))
        body.append(text(x0 + c / 2, y0 + (k + 1.5) * c + 4, str(k), "middle", fill=TEXT3, mono=True, size=12))
    body.append(line(x0, y0 + c, x0 + (n + 1) * c, y0 + c, width=1))
    body.append(line(x0 + c, y0, x0 + c, y0 + (n + 1) * c, width=1))
    for a in range(n):
        for b in range(n):
            v = t[a][b]
            cx, cy = x0 + (b + 1) * c, y0 + (a + 1) * c
            one = v == 1
            zero_div = v == 0 and a and b
            if one:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill=BG2, stroke=GREEN, width=1.5, rx=3))
            elif zero_div:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3))
            fill = GREEN if one else RED if zero_div else ""
            body.append(text(cx + c / 2, cy + c / 2 + 4, str(v), "middle", fill=fill, mono=True, size=12,
                             weight_="600" if one or zero_div else ""))
    return body


def multiplication_tables() -> str:
    c = 30
    body = table(7, 40, 48, c) + table(6, 350, 48, c)
    y = 48 + 8 * c + 34
    body += [rect(60, y - 13, 18, 18, fill=BG2, stroke=GREEN, width=1.5, rx=3),
             text(88, y, "product 1: the row and the column are inverses", size=12),
             rect(60, y + 13, 18, 18, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3),
             text(88, y + 26, "product 0 from two nonzero factors: a zero divisor", size=12)]
    return svg("Multiplication tables modulo 7 and modulo 6, with inverses and zero divisors marked", y + 48, body)


CHARTS = {
    "fig-13-2-euclid-tiling.svg": euclid_tiling,
    "fig-13-3-euclid-remainders.svg": euclid_remainders,
    "fig-13-5-fermat-permutation.svg": fermat_permutation,
    "fig-13-6-crt-grid.svg": crt_grid,
    "fig-13-7-multiplication-tables.svg": multiplication_tables,
}


def default_out() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    site = os.path.join(here, "..", "..", "..", "Website-Building-Projects", "dranirbanchakraborty.com",
                        "public", "figures", "error-correction-coding")
    return os.path.normpath(site) if os.path.isdir(site) else "figures"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=default_out(), help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
