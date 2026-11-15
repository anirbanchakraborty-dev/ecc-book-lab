"""Write the computed figure of chapter 17 of Error Correction Coding as an SVG file.

    python3 python/ch17_charts.py --out figures

Every drawn value comes from the functions in ch17_linear_codes.py, so the
figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet:
real text in IBM Plex, colors from the page's tokens with paper fallbacks, and
no element that the site's figure check would refuse.
"""

import argparse
import os

from ch17_linear_codes import CODE16, HAMMING, code16_weights, hamming_weights, min_distance

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 70, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
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
         size: int = 0) -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}/>')


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(24 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, stroke: str, label: str) -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.5, 3), text(x + 22, y, label, size=12)]


# The chart ---------------------------------------------------------------------------

def weight_distribution() -> str:
    """Codewords of each weight in the (7,4) Hamming code and in the code of Matrices and Elimination."""
    ham, other = hamming_weights(), code16_weights()
    ax = LinAxis(-0.5, 7.5, X0, X1)
    ay = LinAxis(0, 8, Y1, Y0)
    body = frame("weight of the codeword", "number of codewords")
    for v in range(0, 9, 2):
        py = ay(v)
        body.append(line(X0 - 4, py, X0, py))
        body.append(text(X0 - 8, py + 4, str(v), "end", mono=True, size=11))
        if v:
            body.append(line(X0, py, X1, py, TEXT3, 0.5, "3 3"))
    bar = 22
    for w in range(8):
        cx = ax(w)
        body.append(line(cx, Y1, cx, Y1 + 4))
        body.append(text(cx, Y1 + 20, str(w), "middle", mono=True, size=11))
        for value, dx, fill, stroke, colour in ((ham[w], -bar - 1, GREEN, GREEN, GREEN),
                                                (other[w], 1, BG2, TEXT3, TEXT3)):
            if value:
                top = ay(value)
                body.append(rect(cx + dx, top, bar, Y1 - top, fill, stroke, 1.5, 2))
                body.append(text(cx + dx + bar / 2, top - 6, str(value), "middle", colour, mono=True, size=11))
    d_ham, d_other = min_distance(HAMMING), min_distance(CODE16)
    body += legend_box(84, 30, GREEN, GREEN, f"(7,4) Hamming code, d = {d_ham}")
    body += legend_box(310, 30, BG2, TEXT3, f"code of Matrices and Elimination, d = {d_other}")
    numbers = "zero one two three four five six seven eight".split()

    def listing(items: list) -> str:
        return ", ".join(items[:-1]) + " and " + items[-1]

    ham_parts = listing([f"{numbers[ham[w]]} of weight {w}" for w in range(8) if ham[w]])
    other_ws = listing([str(w) for w in range(8) if other[w]])
    other_counts = listing([str(other[w]) for w in range(8) if other[w]])
    ham_parts = ham_parts.replace("one of weight 0", "one codeword of weight 0", 1)
    return svg("Weight distributions of two binary linear codes of length 7 and dimension 4. The (7,4) Hamming "
               f"code has {ham_parts}, so its smallest nonzero weight, and its minimum distance, is {d_ham}. "
               "The code built from the systematic matrix of Matrices and Elimination over Finite Fields has "
               f"weights {other_ws} with counts {other_counts}, so its minimum distance is {d_other}.", body)


CHARTS = {
    "fig-17-6-weight-distribution.svg": weight_distribution,
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
