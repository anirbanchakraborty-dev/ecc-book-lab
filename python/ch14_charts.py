"""Write the two computed figures of chapter 14 of Error Correction Coding as SVG files.

    python3 python/ch14_charts.py --out figures

Every drawn value comes from the functions in ch14_groups_and_cosets.py, so a
figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import os

from ch14_groups_and_cosets import C5, order_of_two, primes_below, show, standard_array, syndrome, weight

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
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
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


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


def x_ticks(ax, values: list) -> list:
    out = []
    for v in values:
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, str(v), "middle", mono=True, size=11))
    return out


def y_ticks(ax, values: list) -> list:
    out = []
    for v in values:
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
        out.append(text(X0 - 8, py + 4, str(v), "end", mono=True, size=11))
    return out


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "") -> list:
    return [line(x, y - 4, x + 24, y - 4, color, 1.5, dash), text(x + 32, y, label, size=12)]


def legend_dot(x: float, y: float, color: str, label: str) -> list:
    return [dot(x + 12, y - 4, color), text(x + 32, y, label, size=12)]


# The charts --------------------------------------------------------------------------

def cosets_of_c5() -> str:
    """The eight cosets of C_5, one row each, with the leader first and the syndrome last."""
    rows = standard_array(C5, 5)
    cell_w, cell_h, gap = 92, 26, 6
    left, top = 67, 58
    xs = [left + i * (cell_w + gap) for i in range(4)]
    syn_x = xs[3] + cell_w + 40
    body = []
    heads = ["leader e"] + [f"e XOR {show(c, 5)}" for c in C5[1:]]
    for x, head in zip(xs, heads):
        body.append(text(x + cell_w / 2, top - 30, head, "middle", size=12))
    body.append(text(syn_x + 40, top - 30, "syndrome", "middle", size=12))
    for r, row in enumerate(rows):
        y = top - 8 + r * (cell_h + gap)
        least = min(weight(x) for x in row)
        for i, x in enumerate(row):
            if i == 0:
                body.append(rect(xs[i], y, cell_w, cell_h, BG2, GREEN, 1.5, 3))
                fill = GREEN
            else:
                body.append(rect(xs[i], y, cell_w, cell_h, "none", "currentColor", 1, 3))
                fill = ORANGE if weight(x) == least and least > 0 else ""
            body.append(text(xs[i] + cell_w / 2, y + 18, show(x, 5), "middle", fill, mono=True, size=13))
        body.append(rect(syn_x, y, 80, cell_h, "none", CYAN, 1, 3))
        body.append(text(syn_x + 40, y + 18, show(syndrome(row[0]), 3), "middle", CYAN, mono=True, size=13))
    ly =top - 8 + len(rows) * (cell_h + gap) + 22
    body.append(rect(left, ly - 12, 14, 14, BG2, GREEN, 1.5, 3))
    body.append(text(left + 22, ly, "coset leader, a word of least weight", size=12))
    tie = next(x for row in rows for x in row[1:] if weight(x) == min(weight(y) for y in row) and weight(x) > 0)
    body.append(text(left + 300, ly, show(tie, 5), fill=ORANGE, mono=True, size=12))
    body.append(text(left + 348, ly, "another word of the same least weight", size=12))
    height = ly + 18
    return svg("The eight cosets of the code C5 with their leaders and syndromes", body, W, height)


def order_of_two_chart() -> str:
    """The order of 2 modulo each odd prime below 100, against p - 1 and (p - 1) / 2."""
    primes = [p for p in primes_below(100) if p > 2]
    ax = LinAxis(0, 100, X0, X1)
    ay = LinAxis(0, 100, Y1, Y0)
    body = frame("prime p", "order of 2 modulo p")
    body += x_ticks(ax, [0, 20, 40, 60, 80, 100])
    body += y_ticks(ax=ay, values=[0, 20, 40, 60, 80, 100])
    body.append(line(ax(2), ay(1), ax(100), ay(99), "currentColor", 1.5, "3 3"))
    body.append(line(ax(2), ay(0.5), ax(100), ay(49.5), TEXT3, 1.5, "3 3"))
    for p in primes:
        o = order_of_two(p)
        color = GREEN if o == p - 1 else CYAN if 2 * o == p - 1 else ORANGE
        body.append(dot(ax(p), ay(o), color))
    lx, ly = X0 + 16, Y0 + 20
    body += legend_line(lx, ly, "currentColor", "p - 1, the order of the group", "3 3")
    body += legend_line(lx, ly + 20, TEXT3, "(p - 1) / 2", "3 3")
    body += legend_dot(lx, ly + 40, GREEN, "order p - 1: 2 generates every nonzero residue")
    body += legend_dot(lx, ly + 60, CYAN, "order (p - 1) / 2")
    body += legend_dot(lx, ly + 80, ORANGE, "order a smaller divisor of p - 1")
    return svg("The order of 2 modulo each odd prime below 100", body)


CHARTS = {
    "fig-14-6-cosets-of-c5.svg": cosets_of_c5,
    "fig-14-7-order-of-two.svg": order_of_two_chart,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="directory for the SVG files")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
