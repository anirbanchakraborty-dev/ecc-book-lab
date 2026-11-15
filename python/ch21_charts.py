"""Write the two computed figures of chapter 21 of Error Correction Coding as SVG files.

    python3 python/ch21_charts.py --out figures

Every drawn value comes from the functions in ch21_beyond_secded.py, so a
figure is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet:
real text in IBM Plex, colors from the page's tokens with paper fallbacks, and
no element that the site's figure check would refuse.
"""

import argparse
import os

from ch21_beyond_secded import (
    berger_bits, choice_rows, detects_up_to, f_bl1, f_bl2, f_mod, least_binary_symbol_bits,
    least_check_bits, least_daec_bits)

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
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


def check_bits_chart() -> str:
    """Check bits each fault model costs for 64 data bits, beside its counting minimum."""
    rows = choice_rows()
    # The counting minimum for each row, where the chapter derives one.
    minimum = [
        least_check_bits(64, 1, False),
        least_check_bits(64, 1, True),
        least_check_bits(64, 2, True),
        least_daec_bits(64),
        least_binary_symbol_bits(64, 4),
        berger_bits(64),
        None,
    ]
    bx0, bx1, top = 262, 566, 16
    row_h, y_first = 40, 70

    def x(v: float) -> float:
        return bx0 + (bx1 - bx0) * v / top

    body = []
    body.append(rect(262, 19, 14, 14, CYAN, rx=3))
    body.append(text(284, 30, "check bits of the code", size=12))
    body.append(line(452, 18, 452, 34, "currentColor", 2))
    body.append(text(462, 30, "counting minimum", size=12))
    for i, ((model, code, r), low) in enumerate(zip(rows, minimum)):
        y0 = y_first + i * row_h
        body.append(text(bx0 - 12, y0 + 9, model, "end", size=12))
        body.append(text(bx0 - 12, y0 + 24, code, "end", fill=TEXT3, mono=True, size=10))
        body.append(rect(bx0, y0 + 2, x(r) - bx0, 22, CYAN, rx=2))
        body.append(text(x(r) + 8, y0 + 18, str(r), mono=True, size=12))
        if low is not None:
            body.append(line(x(low), y0 - 3, x(low), y0 + 29, "currentColor", 2))
    base = y_first + row_h * len(rows)
    body.append(line(bx0, base, bx1, base))
    for v in range(0, top + 1, 4):
        body.append(line(x(v), base, x(v), base + 4))
        body.append(text(x(v), base + 18, str(v), "middle", mono=True, size=11))
    body.append(text((bx0 + bx1) / 2, base + 42, "Check bits beside 64 data bits", "middle"))
    h = base + 60
    return svg("Check bits each fault model costs for 64 data bits", body, W, h)


def bose_lin_chart() -> str:
    """Largest unidirectional error weight each zero-count code detects, against its check bits."""
    k = 64
    series = [
        ('zeros mod 2<tspan dy="-5" font-size="9">r</tspan>', ORANGE, {r: detects_up_to(f_mod(r), k) for r in range(2, 7)}),
        ("offset zero count", CYAN, {r: detects_up_to(f_bl1(r), k) for r in range(4, 8)}),
        ("Bose-Lin, second construction", GREEN, {r: detects_up_to(f_bl2(r), k) for r in range(5, 8)}),
    ]
    x0, x1, y0, y1 = 64, 604, 64, 316
    top = 64

    def y(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = []
    for (label, color, _), lx in zip(series, (64, 210, 418)):
        body += legend_box(lx, 30, color, label)
    for v in (16, 32, 48):
        body.append(line(x0, y(v), x1, y(v), TEXT3, 1, "3 3"))
    group = (x1 - x0) / 6
    bar = 22
    for g, r in enumerate(range(2, 8)):
        gx = x0 + g * group + (group - 3 * bar) / 2
        for s, (_, color, values) in enumerate(series):
            if r not in values:
                continue
            v = values[r]
            if v is None:
                continue
            bx = gx + s * bar
            body.append(rect(bx, y(v), bar - 2, y(0) - y(v), color))
            body.append(text(bx + (bar - 2) / 2, y(v) - 6, str(v), "middle", mono=True, size=11))
        body.append(text(x0 + (g + 0.5) * group, y1 + 18, f"r = {r}", "middle", mono=True, size=11))
    # Seven check bits counting zeros in full is the Berger code itself.
    assert detects_up_to(f_mod(7), k) is None
    bx = x0 + 5 * group + (group - 3 * bar) / 2
    body.append(rect(bx, y(64), bar - 2, y(0) - y(64), BG2, "currentColor", 1))
    body.append(text(bx - 6, y(64) + 14, "Berger: every weight", "end", size=11))
    body.append(f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}" fill="none" stroke="currentColor" stroke-width="1"/>')
    for v in (0, 16, 32, 48, 64):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, str(v), "end", mono=True, size=11))
    body.append(text((x0 + x1) / 2, 372, "Check bits r beside 64 data bits", "middle"))
    body.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">Largest weight always detected</text>')
    return svg("Unidirectional error weight detected against check bits, 64 data bits", body, W, 392)


CHARTS = {
    "fig-21-7-check-bits.svg": check_bits_chart,
    "fig-21-8-bose-lin.svg": bose_lin_chart,
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
