"""Write the two computed figures of chapter 29 of Error Correction Coding as SVG files.

    python3 python/ch29_charts.py --out figures

Every drawn value comes from ch29_case_study_ecc_in_processor_caches.py: the
OpenTitan (39,32) check matrix is rebuilt from the encoder masks transcribed
there, from hw/ip/prim/rtl/prim_secded_39_32_enc.sv at OpenTitan commit
4ced702ca207ad420fbae5226ac8685a17ebef7b, and the scrubbing curves come from
the same Poisson sum the chapter's numbers use. A figure is corrected by fixing
those functions and running this script again, never by editing its output. The
SVG follows the book's figure style sheet: real text in IBM Plex, colors from
the page's tokens with paper fallbacks, and no element that the site's figure
check would refuse.
"""

import argparse
import math
import os

from ch20_secded_and_hsiao_codes import hsiao_39_32
from ch29_case_study_ecc_in_processor_caches import (
    MASKS_39, at_least_two, ot_columns, row_loads, weight)

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
BG2 = "var(--bg-2,#f3eee3)"
TEXT3 = "var(--text-3,#5a5448)"
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


# ---------------------------------------------------------------- scrubbing

WORDS = 4096
YEAR = 8760
RATES = [  # upsets per word per hour, the FIT label, the color
    (1e-6, "1,000 FIT", GREEN),
    (1e-7, "100 FIT", CYAN),
    (1e-8, "10 FIT", ORANGE),
]


def words_lost(lam: float, t: float) -> float:
    """Expected words a year that collect two or more upsets between scrubs."""
    return WORDS * (YEAR / t) * at_least_two(lam * t)


def scrub_chart() -> str:
    """Words lost a year against the scrub interval, log-log, for three upset rates."""
    h = 420
    x0, x1, y0, y1 = 84, 524, 44, 344
    lx0, lx1 = 0.0, 3.0      # log10 of 1 h and 1,000 h
    ly0, ly1 = -9.0, -1.0    # log10 of the words-lost axis

    def px(t: float) -> float:
        return x0 + (x1 - x0) * (math.log10(t) - lx0) / (lx1 - lx0)

    def py(v: float) -> float:
        return y1 - (y1 - y0) * (math.log10(v) - ly0) / (ly1 - ly0)

    body = [rect(x0, y0, x1 - x0, y1 - y0, "none", "currentColor", 1)]
    # decade grid and labels
    for e in range(int(ly0), int(ly1) + 1):
        y = py(10.0 ** e)
        if ly0 < e < ly1:
            body.append(line(x0, y, x1, y, TEXT3, 0.5, "3 3"))
        body.append(line(x0 - 4, y, x0, y))
        body.append(f'<text x="{fmt(x0 - 8)}" y="{fmt(y + 4)}" text-anchor="end" font-family="{MONO}"'
                    f' font-size="11">10<tspan dy="-5" font-size="9">{e}</tspan></text>')
    for t, label in ((1, "1"), (10, "10"), (100, "100"), (1000, "1,000")):
        x = px(t)
        body.append(line(x, y1, x, y1 + 4))
        body.append(text(x, y1 + 18, label, "middle", mono=True, size=11))
    # the three intervals the chapter tabulates
    for t, label in ((24, "day"), (168, "week"), (720, "month")):
        x = px(t)
        body.append(line(x, y0, x, y1, TEXT3, 1, "3 3"))
        body.append(text(x, y0 - 8, label, "middle", fill=TEXT3, size=11))
    # one curve per rate, sampled every tenth of a decade
    samples = [10.0 ** (i / 10) for i in range(0, 31)]
    for lam, label, color in RATES:
        pts = " ".join(f"{fmt(px(t))},{fmt(py(words_lost(lam, t)))}" for t in samples)
        body.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>')
        end = py(words_lost(lam, 1000))
        body.append(text(x1 + 8, end + 4, label, fill=color, size=12, bold=True))
    # the four points the chapter's table gives, on the 1,000 FIT curve
    for t in (1, 24, 168, 720):
        body.append(f'<circle cx="{fmt(px(t))}" cy="{fmt(py(words_lost(1e-6, t)))}" r="3.5"'
                    f' fill="{GREEN}" stroke="currentColor" stroke-width="1"/>')
    body.append(text(x1 + 8, py(words_lost(1e-6, 1000)) + 20, "per word", fill=TEXT3, size=11))
    body.append(text((x0 + x1) / 2, y1 + 44, "Scrub interval T, hours, log scale", "middle"))
    body.append(f'<text x="22" y="{fmt((y0 + y1) / 2)}" text-anchor="middle"'
                f' transform="rotate(-90 22 {fmt((y0 + y1) / 2)})">'
                f"Words a year with two or more upsets</text>")
    body.append(text(x0, h - 14, "4,096 words of 72 bits, one year of 8,760 hours, Poisson upsets",
                     fill=TEXT3, size=11))
    return svg("Expected words a year that collect two or more upsets between scrubs, against the scrub "
               "interval T from 1 to 1,000 hours, on log scales, for upset rates of 1,000, 100 and 10 FIT "
               "per 72-bit word in a 4,096-word array. Each curve is a straight line of slope one, so "
               "doubling the interval doubles the loss, and each tenfold cut in the rate cuts the loss a "
               "hundredfold. Dots mark the hourly, daily, weekly and monthly intervals on the 1,000 FIT "
               "curve.", body, W, h)


# ---------------------------------------------------------------- the OpenTitan matrix

def matrix_chart() -> str:
    """The OpenTitan Hsiao (39,32) check matrix, one cell per entry."""
    k, r = 32, 7
    cols = ot_columns(MASKS_39, k)
    book = set(hsiao_39_32()[:k])
    loads = row_loads(cols, r)
    cell, gap = 12, 10
    gx = 82
    gy = 60
    h = 304

    def cx(j: int) -> float:
        return gx + j * cell + (gap if j >= k else 0)

    body = []
    right = cx(k + r - 1) + cell
    body.append(text(gx + k * cell / 2, gy - 28, "data bits 0 to 31", "middle", size=12))
    body.append(text(cx(k) + r * cell / 2, gy - 28, "check bits", "middle", size=12))
    for j in range(0, k + r):
        if j % 4 == 0 or j == k + r - 1:
            body.append(text(cx(j) + cell / 2, gy - 8, str(j), "middle", mono=True, size=9))
    for i in range(r):
        y = gy + i * (cell + 4)
        body.append(text(gx - 10, y + 10, f"row {i}", "end", size=11))
        for j, c in enumerate(cols):
            on = (c >> i) & 1
            if on:
                fill = GREEN if j < k else CYAN
                body.append(rect(cx(j) + 1, y + 1, cell - 2, cell - 2, fill, rx=2))
            else:
                body.append(rect(cx(j) + 1, y + 1, cell - 2, cell - 2, BG2, TEXT3, 0.5, rx=2))
        body.append(text(right + 14, y + 10, str(loads[i]), mono=True, size=11))
    body.append(text(right + 14, gy - 8, "ones", size=11))
    wy = gy + r * (cell + 4) + 14
    body.append(text(gx - 10, wy, "weight", "end", size=11))
    for j, c in enumerate(cols):
        body.append(text(cx(j) + cell / 2, wy, str(weight(c)), "middle", mono=True, size=9))
    # the data columns the book's matrix does not use
    my = wy + 10
    only = [j for j, c in enumerate(cols[:k]) if c not in book]
    for j in only:
        x = cx(j) + cell / 2
        body.append(f'<path d="M{fmt(x)} {fmt(my)}L{fmt(x - 5)} {fmt(my + 9)}H{fmt(x + 5)}z"'
                    f' fill="{ORANGE}"/>')
    body.append(text(gx, my + 28, "Triangles mark data columns "
                     f"{', '.join(str(j) for j in only[:-1])} and {only[-1]}, which the book's (39,32) matrix lacks.",
                     size=11))
    body.append(text(gx, my + 44, "Every column has odd weight, and no two columns are equal.", size=11))
    # legend
    ly = h - 38
    body.append(rect(gx, ly - 10, 10, 10, GREEN, rx=2))
    body.append(text(gx + 16, ly, "a 1 in a data column", size=11))
    body.append(rect(gx + 170, ly - 10, 10, 10, CYAN, rx=2))
    body.append(text(gx + 186, ly, "a 1 in a check column", size=11))
    body.append(rect(gx + 350, ly - 10, 10, 10, BG2, TEXT3, 0.5, rx=2))
    body.append(text(gx + 366, ly, "a 0", size=11))
    body.append(text(gx, h - 14, "Rebuilt from prim_secded_39_32_enc.sv, OpenTitan commit 4ced702ca207",
                     fill=TEXT3, size=11))
    # every count in the title comes from the masks, so it follows any change to them
    words = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]

    def listed(xs: list) -> str:
        xs = [str(x) for x in xs]
        return xs[0] if len(xs) == 1 else f"{', '.join(xs[:-1])} and {xs[-1]}"

    data_w = listed([words[w] for w in sorted({weight(c) for c in cols[:k]})])
    check_w = listed([words[w] for w in sorted({weight(c) for c in cols[k:]})])
    return svg(f"The check matrix of OpenTitan's Hsiao (39,32) code, rebuilt from the {words[r]} encoder masks, "
               f"with {words[r]} rows, {k} data columns of weight {data_w} and {words[r]} check columns of "
               f"weight {check_w}. The rows hold {listed(loads)} ones. "
               f"{words[len(only)].capitalize()} data columns, marked with triangles, are "
               "columns the book's (39,32) matrix does not use.", body, W, h)


CHARTS = {
    "fig-29-4-scrub-interval.svg": scrub_chart,
    "fig-29-7-opentitan-39-32.svg": matrix_chart,
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
