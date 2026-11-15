"""Write the two computed figures of chapter 22 of Error Correction Coding as SVG files.

    python3 python/ch22_charts.py --out figures

Every drawn value comes from the functions in
ch22_arithmetic_codes_and_checking_logic.py, so a figure is corrected by
fixing those functions and running this script again, never by editing its
output. The SVG follows the book's figure style sheet: real text in IBM Plex,
colors from the page's tokens with paper fallbacks, and no element that the
site's figure check would refuse.
"""

import argparse
import os

from ch22_arithmetic_codes_and_checking_logic import adder_experiment, correcting

W, H = 640, 400
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


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def smallest_correcting(k: int):
    """The smallest odd A whose code corrects every single arithmetic error in k data bits."""
    a = 3
    while True:
        m = (a * ((1 << k) - 1)).bit_length()
        if a % 2 == 1 and correcting(a, m):
            return a, m
        a += 1


def smallest_a_chart() -> str:
    x0, x1, y0, y1 = 84, 604, 28, 316
    xa = LinAxis(0, 32, x0, x1)
    ya = LinAxis(0, 90, y1, y0)
    body = []
    for v in range(0, 91, 15):
        y = ya(v)
        body.append(line(x0 - 5, y, x0, y, "currentColor", 1))
        body.append(text(x0 - 10, y + 4, str(v), "end", mono=True, size=11))
    for k in range(0, 33, 4):
        x = xa(k)
        body.append(line(x, y1, x, y1 + 5, "currentColor", 1))
        body.append(text(x, y1 + 20, str(k), "middle", mono=True, size=11))
    body.append(line(x0, y1, x1, y1, "currentColor", 1.5))
    body.append(line(x0, y0, x0, y1, "currentColor", 1.5))
    body.append(text((x0 + x1) / 2, y1 + 44, "data bits k", "middle", size=12))
    body.append(f'<text x="24" y="{fmt((y0 + y1) / 2)}" text-anchor="middle" font-size="12" '
                f'transform="rotate(-90 24 {fmt((y0 + y1) / 2)})">smallest correcting A</text>')
    marked = {4, 8, 12, 16, 32}
    points = [(k, *smallest_correcting(k)) for k in range(1, 33)]
    path = " ".join(f"{'M' if i == 0 else 'L'}{fmt(xa(k))} {fmt(ya(a))}" for i, (k, a, _) in enumerate(points))
    body.append(f'<path d="{path}" fill="none" stroke="{TEXT3}" stroke-width="1"/>')
    for k, a, m in points:
        if k in marked:
            body.append(dot(xa(k), ya(a), GREEN, 5))
        else:
            body.append(dot(xa(k), ya(a), TEXT3, 3))
    for k, a, m in points:
        if k not in marked:
            continue
        label = f"A = {a}, {m - k} check bits"
        if k == 32:
            body.append(text(xa(k) - 10, ya(a) - 12, label, "end", mono=True, size=11))
        else:
            body.append(text(xa(k) + 10, ya(a) + 18, label, mono=True, size=11))
    body.append(dot(xa(1) + 14, 46, GREEN, 5))
    body.append(text(xa(1) + 26, 50, "values the chapter states", size=12))
    body.append(dot(xa(1) + 14, 68, TEXT3, 3))
    body.append(text(xa(1) + 26, 72, "the remaining k from 1 to 32", size=12))
    return svg("The smallest A that corrects every single arithmetic error, for data words of 1 to 32 bits", body)


def checks_caught_chart() -> str:
    rows = adder_experiment()
    wrong = sum(r[1] for r in rows)
    caught = [
        ("mod-3 residue check", sum(r[2] for r in rows), GREEN),
        ("parity, shared carries", sum(r[3] for r in rows), ORANGE),
        ("parity, duplicated carry chain", sum(r[4] for r in rows), ORANGE),
        ("duplicated chain and compare", sum(r[5] for r in rows), GREEN),
        ("Berger check prediction", sum(r[6] for r in rows), GREEN),
    ]
    x0, x1 = 236, 560
    xa = LinAxis(0, wrong, x0, x1)
    top, step, bar = 40, 52, 28
    body = []
    axis_y = top + step * len(caught) - 10
    body.append(line(x0, axis_y, x1, axis_y, "currentColor", 1.5))
    for v in range(0, wrong + 1, 192):
        x = xa(v)
        body.append(line(x, axis_y, x, axis_y + 5, "currentColor", 1))
        body.append(text(x, axis_y + 20, str(v), "middle", mono=True, size=11))
    for i, (label, n, color) in enumerate(caught):
        y = top + i * step
        body.append(text(x0 - 12, y + bar / 2 + 4, label, "end", size=12))
        if n:
            body.append(rect(x0, y, xa(n) - x0, bar, color, color, 1, 2))
        body.append(text(max(xa(n), x0) + 8, y + bar / 2 + 4, f"{n} of {wrong}", mono=True, size=11))
    body.append(line(x0, top - 10, x0, axis_y, "currentColor", 1.5))
    body.append(text((x0 + x1) / 2, axis_y + 42, "wrong results caught", "middle", size=12))
    body.append(text((x0 + x1) / 2, axis_y + 62,
                     "4-bit adder, carries c1 to c3 stuck at 0 or 1, all 256 operand pairs", "middle",
                     fill=TEXT3, size=11))
    return svg(f"Of {wrong} wrong sums from a 4-bit adder with one stuck carry, how many each check catches", body,
               W, 360)


CHARTS = {
    "fig-22-2-smallest-correcting-a.svg": smallest_a_chart,
    "fig-22-5-checks-caught.svg": checks_caught_chart,
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
