"""Write the two computed figures of chapter 20 of Error Correction Coding as SVG files.

    python3 python/ch20_charts.py --out figures

Every drawn value comes from the code constructions in
ch20_secded_and_hsiao_codes.py, so a figure is corrected by fixing those
functions and running this script again, never by editing its output. The SVG
follows the book's figure style sheet: real text in IBM Plex, colors from the
page's tokens with paper fallbacks, and no element that the site's figure check
would refuse.
"""

import argparse
import os

from ch20_secded_and_hsiao_codes import (
    extended_hamming, hsiao_8_4, hsiao_13_8, hsiao_39_32, hsiao_72_64,
    hsiao_72_64_dictionary, row_loads, triple_outcomes)

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


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


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", width: float = 1.5,
         rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    s = f' stroke="{stroke}" stroke-width="{width}"' if stroke != "none" else ' stroke="none"'
    return f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}"{s}/>'


def svg(title: str, body: list, w: int, h: int) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, label: str) -> list:
    return [rect(x, y - 11, 14, 14, fill, rx=3), text(x + 22, y, label, size=12)]


def row_load_chart() -> str:
    """Ones in each row of three (72,64) SEC-DED check matrices."""
    h = 400
    x0, x1, y0, y1 = 64, 568, 60, 324
    top = 80
    series = [
        ("Extended Hamming (72,64)", ORANGE, row_loads(extended_hamming(72, 7), 8)),
        ("Hsiao, dictionary order", CYAN, row_loads(hsiao_72_64_dictionary(), 8)),
        ("Hsiao, balanced", GREEN, row_loads(hsiao_72_64(), 8)),
    ]

    def y(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = []
    # legend
    for (label, color, _), lx in zip(series, (64, 277, 483)):
        body += legend_box(lx, 30, color, label)
    # depth bands of a two-input XOR tree, drawn first so the bars cover them
    for v in (16, 32, 64):
        body.append(line(x0, y(v), x1, y(v), TEXT3, 1, "3 3"))
    for lo, hi, label in ((16, 32, "5 levels"), (32, 64, "6 levels"), (64, 80, "7 levels")):
        body.append(text(x1 + 8, (y(lo) + y(hi)) / 2 + 4, label, fill=TEXT3, size=11))
    # bars
    group = (x1 - x0) / 8
    bar = 15
    for row in range(8):
        gx = x0 + row * group + (group - 3 * bar) / 2
        for k, (_, color, loads) in enumerate(series):
            v = loads[row]
            body.append(rect(gx + k * bar, y(v), bar, y(0) - y(v), color))
        body.append(text(x0 + (row + 0.5) * group, y1 + 18, f"row {row}", "middle", mono=True, size=11))
    # axes
    body.append(f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}" fill="none" stroke="currentColor" stroke-width="1"/>')
    for v in (0, 16, 32, 48, 64, 80):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, str(v), "end", mono=True, size=11))
    body.append(text((x0 + x1) / 2, h - 18, "Row of the check matrix, one syndrome bit each", "middle"))
    body.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">Ones in the row</text>')
    return svg("Ones in each row of three (72,64) SEC-DED check matrices", body, W, h)


def triple_chart() -> str:
    """Share of three-bit errors miscorrected and flagged, for six SEC-DED codes."""
    codes = [
        ("Hsiao (8,4)", hsiao_8_4()),
        ("Hsiao (13,8)", hsiao_13_8()),
        ("Hsiao (39,32)", hsiao_39_32()),
        ("Extended Hamming (39,32)", extended_hamming(39, 6)),
        ("Hsiao (72,64)", hsiao_72_64()),
        ("Extended Hamming (72,64)", extended_hamming(72, 7)),
    ]
    bx0, bx1 = 196, 516
    row_h, y_first = 38, 76
    h = y_first + row_h * len(codes) + 64
    body = []
    body += legend_box(196, 30, RED, "miscorrected, a fourth bit flipped")
    body += legend_box(484, 30, GREEN, "flagged")
    body.append(text(bx1 + 10, y_first - 10, "miscorrected", fill=TEXT3, size=11))
    for i, (name, cols) in enumerate(codes):
        mis, det = triple_outcomes(cols)
        total = mis + det
        top = y_first + i * row_h
        split = bx0 + (bx1 - bx0) * mis / total
        body.append(text(bx0 - 10, top + 16, name, "end", size=12))
        body.append(rect(bx0, top + 2, split - bx0, 20, RED))
        if det:
            body.append(rect(split, top + 2, bx1 - split, 20, GREEN))
        body.append(text(bx1 + 10, top + 16, f"{100 * mis / total:.1f}%", mono=True, size=12))
    base = y_first + row_h * len(codes)
    body.append(line(bx0, base, bx1, base))
    for p in (0, 25, 50, 75, 100):
        px = bx0 + (bx1 - bx0) * p / 100
        body.append(line(px, base, px, base + 4))
        body.append(text(px, base + 18, str(p), "middle", mono=True, size=11))
    body.append(text((bx0 + bx1) / 2, base + 42, "Share of all three-bit errors, percent", "middle"))
    return svg("What a SEC-DED decoder does with a three-bit error", body, W, h)


CHARTS = {
    "fig-20-5-row-loads.svg": row_load_chart,
    "fig-20-6-triple-errors.svg": triple_chart,
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
