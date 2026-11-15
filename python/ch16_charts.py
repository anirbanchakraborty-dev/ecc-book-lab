"""Write the four computed figures of chapter 16 of Error Correction Coding as SVG files.

    python3 python/ch16_charts.py --out figures

Every drawn value comes from the functions in
ch16_matrices_and_elimination_over_finite_fields.py, so a figure is corrected by
fixing those functions and running this script again, never by editing its
output. The SVG follows the book's figure style sheet: real text in IBM Plex,
colors from the page's tokens with paper fallbacks, and no element that the
site's figure check would refuse.
"""

import argparse
import os

from ch16_matrices_and_elimination_over_finite_fields import (
    A, M, augment, coset_table, identity, invertible_fraction, rref)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" '
         'markerHeight="5" orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker></defs>')


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
         dash: str = "", opacity: str = "", arrow: bool = False) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    m = ' marker-end="url(#ah)"' if arrow else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}{m}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, stroke: str, label: str) -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.5, 3), text(x + 22, y, label, size=12)]


def op_text(step: tuple) -> str:
    if step[0] == "swap":
        return f"swap R{step[1] + 1} and R{step[2] + 1}"
    return f"add R{step[3] + 1} to R{step[1] + 1}"


def grid(x: float, y: float, rows: list, cell: float, pivots: list, split: int = 0,
         size: int = 13) -> list:
    """A matrix of bits, one square cell per entry, with the given (row, column) cells in green.

    When split is set, a gap of 10 units opens after that many columns and a rule
    is drawn in it, as in an augmented matrix.
    """
    out = []
    gap = 10 if split else 0
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            cx = x + j * cell + (gap if split and j >= split else 0)
            cy = y + i * cell
            if (i, j) in pivots:
                out.append(rect(cx, cy, cell, cell, BG2, GREEN, 1.5, 3))
                out.append(text(cx + cell / 2, cy + cell / 2 + 5, str(v), "middle", GREEN, mono=True, size=size))
            else:
                out.append(rect(cx, cy, cell, cell, "none", TEXT3, 1, 3))
                out.append(text(cx + cell / 2, cy + cell / 2 + 5, str(v), "middle", mono=True, size=size))
    if split:
        rx = x + split * cell + gap / 2
        out.append(line(rx, y - 4, rx, y + len(rows) * cell + 4, "currentColor", 1.5))
    return out


# The charts --------------------------------------------------------------------------

def elimination_of_a() -> str:
    """The 3-by-4 matrix A and the matrix after each pivot column is cleared."""
    _, pivots, _, stages = rref(A)
    panels = [("A", A, [])] + [(f"column {c + 1} cleared", R, done) for c, done, R in stages]
    cell, left, step = 26, 40, 152
    top = 52
    body = [ARROW]
    for k, (title, rows, done) in enumerate(panels):
        x = left + k * step
        body.append(text(x + 52, 22, title, "middle", size=12))
        piv = [(i, pivots[i]) for i in range(k)]
        last = k == len(panels) - 1
        for j in range(4):
            free = last and j not in pivots
            body.append(text(x + j * cell + cell / 2, top - 8, str(j + 1), "middle",
                             ORANGE if free else TEXT3, mono=True, size=10))
        body += grid(x, top, rows, cell, piv)
        if last:
            for j in (c for c in range(4) if c not in pivots):
                body.append(rect(x + j * cell - 3, top - 3, cell + 6, 3 * cell + 6, "none", ORANGE, 1.5, 3))
        for t, s in enumerate(done):
            body.append(text(x + 52, top + 3 * cell + 22 + 15 * t, op_text(s), "middle", size=11))
        if k < len(panels) - 1:
            y = top + 1.5 * cell
            body.append(line(x + 4 * cell + 8, y, x + step - 10, y, "currentColor", 1.5, arrow=True))
    ly = 204
    body += legend_box(left, ly, BG2, GREEN, "pivot, a leading 1 alone in its column")
    body += legend_box(left + 330, ly, "none", ORANGE, "free column, no pivot")
    return svg("Gauss-Jordan elimination of the 3-by-4 matrix A, one pivot column at a time", body, W, 220)


def solutions_as_cosets() -> str:
    """For each of the eight right sides b, the two solutions of A x = b."""
    table = coset_table(A)
    bw, bh, gx, gy = 128, 64, 16, 16
    left, top = 40, 16
    body = []
    for k, (b, sols) in enumerate(table):
        x = left + (k % 4) * (bw + gx)
        y = top + (k // 4) * (bh + gy)
        if b == "000":
            body.append(rect(x, y, bw, bh, BG2, GREEN, 1.5, 3))
            fill = GREEN
        else:
            body.append(rect(x, y, bw, bh, "none", "currentColor", 1, 3))
            fill = ""
        body.append(text(x + bw / 2, y + 22, f"b = {b}", "middle", TEXT3 if not fill else fill, mono=True, size=12))
        body.append(text(x + bw / 2 - 28, y + 48, sols[0], "middle", fill, mono=True, size=14))
        body.append(text(x + bw / 2 + 28, y + 48, sols[1], "middle", fill, mono=True, size=14))
    ly = top + 2 * bh + gy + 34
    body += legend_box(left, ly, BG2, GREEN, "b = 000, whose two solutions are the null space of A")
    body.append(rect(left, ly + 11, 14, 14, "none", "currentColor", 1, 3))
    body.append(text(left + 22, ly + 22, "any other b, whose solutions are the null space shifted by one solution", size=12))
    body.append(text(left, ly + 48, "In every box the two solutions differ by 1110, the one nonzero word of the null space.",
                     size=12))
    return svg("The two solutions of A x = b for each of the eight right sides b", body, W, ly + 62)


def inverse_by_elimination() -> str:
    """The augmented matrix [M | I] at the start and after each pivot column is cleared."""
    _, pivots, _, stages = rref(augment(M, identity(3)))
    panels = [("[M | I] at the start", augment(M, identity(3)), [])]
    panels += [(f"column {c + 1} cleared", R, done) for c, done, R in stages]
    cell = 24
    xs, tops = [40, 350], [52, 204]
    body = []
    for k, (title, rows, done) in enumerate(panels):
        x, top = xs[k % 2], tops[k // 2]
        body.append(text(x, top - 18, title, size=12))
        body += grid(x, top, rows, cell, [(i, pivots[i]) for i in range(k)], split=3)
        right = x + 3 * cell + 10
        if k == 0:
            body.append(text(x + 1.5 * cell, top + 3 * cell + 18, "M", "middle", TEXT3, size=12))
            body.append(text(right + 1.5 * cell, top + 3 * cell + 18, "I", "middle", TEXT3, size=12))
        if k == len(panels) - 1:
            body.append(rect(right - 2, top - 2, 3 * cell + 4, 3 * cell + 4, "none", CYAN, 1.5, 3))
            body.append(text(x + 1.5 * cell, top + 3 * cell + 18, "I", "middle", GREEN, size=12))
            body.append(text(right + 1.5 * cell, top + 3 * cell + 18, "inverse of M", "middle", CYAN, size=12))
        for t, s in enumerate(done):
            body.append(text(right + 3 * cell + 24, top + 20 + 16 * t, op_text(s), size=11))
    ly = 330
    body += legend_box(40, ly, BG2, GREEN, "pivot")
    body += legend_box(160, ly, "none", CYAN, "the inverse of M, read off the right half at the end")
    return svg("Inverting the 3-by-3 matrix M by row-reducing the augmented matrix [M | I]", body, W, 346)


def invertible_fraction_chart() -> str:
    """The fraction of n-by-n binary matrices that are invertible, for n from 1 to 10."""
    ax = LinAxis(0.5, 10.5, X0, X1)
    ay = LinAxis(0.2, 0.55, Y1, Y0)
    body = frame("matrix size n", "fraction invertible")
    for n in range(1, 11):
        px = ax(n)
        body.append(line(px, Y1, px, Y1 + 4))
        body.append(text(px, Y1 + 20, str(n), "middle", mono=True, size=11))
    for v in (0.2, 0.3, 0.4, 0.5):
        py = ay(v)
        body.append(line(X0 - 4, py, X0, py))
        body.append(text(X0 - 8, py + 4, f"{v:.2f}", "end", mono=True, size=11))
    limit = invertible_fraction(200)
    body.append(line(X0, ay(limit), X1, ay(limit), "currentColor", 1.5, "3 3"))
    body.append(line(X0, ay(0.25), X1, ay(0.25), ORANGE, 1.5, "3 3"))
    for n in range(1, 11):
        f = invertible_fraction(n)
        body.append(dot(ax(n), ay(f), GREEN))
        body.append(text(ax(n), ay(f) - 10, f"{f:.3f}", "middle", GREEN, mono=True, size=10))
    lx = 330
    body.append(dot(lx + 12, 50, GREEN))
    body.append(text(lx + 32, 54, "fraction invertible at size n", size=12))
    body.append(line(lx, 74, lx + 24, 74, "currentColor", 1.5, "3 3"))
    body.append(text(lx + 32, 78, f"the limit, {limit:.6f}", size=12))
    body.append(line(lx, 98, lx + 24, 98, ORANGE, 1.5, "3 3"))
    body.append(text(lx + 32, 102, "the bound 1/4", size=12))
    return svg("The fraction of n-by-n binary matrices that are invertible, falling toward 0.2888", body)


CHARTS = {
    "fig-16-2-elimination-of-a.svg": elimination_of_a,
    "fig-16-4-solutions-as-cosets.svg": solutions_as_cosets,
    "fig-16-5-inverse-by-elimination.svg": inverse_by_elimination,
    "fig-16-6-invertible-fraction.svg": invertible_fraction_chart,
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
