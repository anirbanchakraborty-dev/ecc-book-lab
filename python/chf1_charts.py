"""Write the two charts of front-matter chapter f1 of Error Correction Coding as SVG files.

    python3 python/chf1_charts.py --out figures

Every plotted value comes from chf1_ecc_before_you_start.py: the chapter counts
from its table of the book's Parts, and the binomial coefficients from math.comb.
A chart is corrected by fixing that table or this script and running it again,
never by editing its output. The SVG follows the book's figure style sheet:
real text in IBM Plex, colors from the page's tokens with paper fallbacks, and
no element that the site's figure check would refuse.
"""

import argparse
import os
from math import comb

from chf1_ecc_before_you_start import PARTS

GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
MONO = "IBM Plex Mono, monospace"

LEVEL_COLOR = {"front": TEXT3, "undergraduate": GREEN, "graduate": CYAN, "research": VIOLET, "back": TEXT3}


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def svg(w: int, h: int, title: str, body: list) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def y_axis(x0: int, x1: int, y_of, ticks: range, label: str, y0: int, y1: int) -> list:
    out = [f'<g stroke="currentColor" stroke-width="1" fill="none">'
           f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}"/></g>']
    for t in ticks:
        y = y_of(t)
        out.append(f'<line x1="{x0 - 4}" y1="{fmt(y)}" x2="{x0}" y2="{fmt(y)}" stroke="currentColor" stroke-width="1"/>')
        if t:
            out.append(f'<line x1="{x0}" y1="{fmt(y)}" x2="{x1}" y2="{fmt(y)}" stroke="{TEXT3}" stroke-width="1" stroke-opacity="0.35"/>')
        out.append(f'<text x="{x0 - 8}" y="{fmt(y + 4)}" text-anchor="end" font-family="{MONO}" font-size="11">{t}</text>')
    out.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">{label}</text>')
    return out


def book_map() -> str:
    """Chapters per Part, front matter to Reference, with each level shown by color and by a text label."""
    w, h = 720, 420
    x0, x1, y0, y1 = 70, 700, 30, 320
    top = 20
    pitch = (x1 - x0) / len(PARTS)
    bar = 24

    def y_of(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = y_axis(x0, x1, y_of, range(0, top + 1, 5), "chapters", y0, y1)
    short = {"Front matter": "F", "Reference": "R"}
    groups = []
    for i, (label, _, level, n) in enumerate(PARTS):
        cx = x0 + pitch * (i + 0.5)
        color = LEVEL_COLOR[level]
        body.append(f'<rect x="{fmt(cx - bar / 2)}" y="{fmt(y_of(n))}" width="{bar}" height="{fmt(y1 - y_of(n))}"'
                    f' rx="2" fill="{color}"/>')
        body.append(f'<text x="{fmt(cx)}" y="{fmt(y_of(n) - 6)}" text-anchor="middle" font-family="{MONO}"'
                    f' font-size="11">{n}</text>')
        tick = short.get(label, label.removeprefix("Part "))
        body.append(f'<text x="{fmt(cx)}" y="{y1 + 18}" text-anchor="middle" font-family="{MONO}"'
                    f' font-size="11">{tick}</text>')
        if groups and groups[-1][0] == level:
            groups[-1][2] = i
        else:
            groups.append([level, i, i])
    for level, first, last in groups:
        a = x0 + pitch * first + 4
        b = x0 + pitch * (last + 1) - 4
        color = LEVEL_COLOR[level]
        total = sum(n for _, _, lv, n in PARTS if lv == level)
        text = level if level in ("front", "back") else f"{level}, {total}"
        body.append(f'<polyline points="{fmt(a)},{y1 + 26} {fmt(a)},{y1 + 31} {fmt(b)},{y1 + 31} {fmt(b)},{y1 + 26}"'
                    f' fill="none" stroke="{color}" stroke-width="1.5"/>')
        body.append(f'<text x="{fmt((a + b) / 2)}" y="{y1 + 48}" text-anchor="middle" fill="{color}"'
                    f' font-weight="600">{text}</text>')
    body.append(f'<text x="{(x0 + x1) // 2}" y="{h - 14}" text-anchor="middle">'
                f'Part, with F for the front matter and R for the Reference back matter</text>')
    total = sum(n for *_, n in PARTS)
    return svg(w, h, f"Chapters in each Part of Error Correction Coding, {total} in all. Parts I to V are undergraduate "
               "with 71 chapters, Parts VI to XII graduate with 75, and Parts XIII to XV research with 36, beside 3 "
               "chapters of front matter and 3 of back matter.", body)


def binomial_row(n: int = 7) -> str:
    """The binomial coefficients n choose k for k = 0 to n, with their sum."""
    w, h = 640, 400
    x0, x1, y0, y1 = 84, 604, 24, 324
    top = 40
    row = [comb(n, k) for k in range(n + 1)]
    pitch = (x1 - x0) / len(row)
    bar = 32

    def y_of(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = y_axis(x0, x1, y_of, range(0, top + 1, 10), f"{n}-bit strings with exactly k ones", y0, y1)
    for k, v in enumerate(row):
        cx = x0 + pitch * (k + 0.5)
        body.append(f'<rect x="{fmt(cx - bar / 2)}" y="{fmt(y_of(v))}" width="{bar}" height="{fmt(y1 - y_of(v))}"'
                    f' rx="2" fill="{GREEN}"/>')
        body.append(f'<text x="{fmt(cx)}" y="{fmt(y_of(v) - 6)}" text-anchor="middle" font-family="{MONO}"'
                    f' font-size="11">{v}</text>')
        body.append(f'<line x1="{fmt(cx)}" y1="{y1}" x2="{fmt(cx)}" y2="{y1 + 4}" stroke="currentColor" stroke-width="1"/>')
        body.append(f'<text x="{fmt(cx)}" y="{y1 + 20}" text-anchor="middle" font-family="{MONO}" font-size="11">{k}</text>')
    half = len(row) // 2
    terms = [str(v) for v in row]
    lines = [" + ".join(terms[:half]), "+ " + " + ".join(terms[half:]), f"= {sum(row)}"]
    for j, line in enumerate(lines):
        body.append(f'<text x="{x0 + 10}" y="{y0 + 18 + 16 * j}" font-family="{MONO}" font-size="11">{line}</text>')
    body.append(f'<text x="{x0 + 10}" y="{y0 + 18 + 16 * len(lines)}" fill="{GREEN}" font-weight="600">'
                f'all {sum(row)} strings of {n} bits</text>')
    body.append(f'<text x="{(x0 + x1) // 2}" y="{h - 14}" text-anchor="middle">'
                f'k, the number of ones in a {n}-bit string</text>')
    return svg(w, h, f"The binomial coefficients {n} choose k for k from 0 to {n}: {', '.join(terms)}. Each counts the "
               f"{n}-bit strings with exactly k ones, and together they count all {sum(row)} strings of {n} bits.", body)


CHARTS = {
    "fig-f1-1-book-map.svg": book_map,
    "fig-f1-5-binomial-row.svg": binomial_row,
}


def default_out() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    site = os.path.join(here, "..", "..", "..", "Website-Building-Projects", "dranirbanchakraborty.com",
                        "public", "figures", "error-correction-coding")
    return os.path.normpath(site) if os.path.isdir(site) else "figures"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=default_out(), help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    for name, chart in CHARTS.items():
        path = os.path.join(out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(chart())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
