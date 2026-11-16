"""Write the four computed figures of chapter 31 of Error Correction Coding as SVG files.

    python3 python/ch31_charts.py --out figures

The chapter is "Lab: Measuring Miscorrection". Every measured value comes from
data/ch31/miscorrection_counts.json, which tb/miscorrection/run_miscorrection.py
wrote from the hardware campaign, and every exact value, interval and bound
from ch31_lab_measuring_miscorrection.py, so a figure is corrected by fixing
those and running this script again, never by editing its output. The SVG
follows the book's figure style sheet: real text in IBM Plex, colors from the
page's tokens with paper fallbacks, and no element that the site's figure check
would refuse. Each series is told apart by a shape or a label as well as by its
color.
"""

import argparse
import json
import os
from math import comb, log10, sqrt

from ch31_lab_measuring_miscorrection import (
    CHECKPOINTS, CODES, DATA, N, Z, Code, exact_split, half_chebyshev, half_hoeffding, layout,
    adjacent_outcomes, syndrome_counts, tracked, wilson, draws_normal, draws_hoeffding, draws_chebyshev)

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
DELTA = 0.05                              # the failure chance the two bounds are held to
STYLE = {"hsiao": ("circle", GREEN), "hsiao_dictionary": ("square", CYAN), "hamming": ("diamond", ORANGE)}
NAMES = {"hsiao": "balanced", "hsiao_dictionary": "dictionary", "hamming": "Hamming"}
LAYOUTS = ("natural", "stride5", "stride7", "stride11", "stride13", "stride29", "searched")


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


def polyline(points: list, stroke: str, width: float = 1.5, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def marker(shape: str, x: float, y: float, color: str, size: float = 4.5, hollow: bool = False) -> str:
    """A point drawn as a circle, a square, a diamond or a triangle, filled or hollow, so the shape names the series."""
    paint = f'fill="{BG2}" stroke="{color}" stroke-width="1.5"' if hollow else f'fill="{color}"'
    if shape == "circle":
        return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(size)}" {paint}/>'
    if shape == "square":
        return (f'<rect x="{fmt(x - size)}" y="{fmt(y - size)}" width="{fmt(2 * size)}"'
                f' height="{fmt(2 * size)}" {paint}/>')
    if shape == "diamond":
        s = size + 1
        pts = f"{fmt(x)},{fmt(y - s)} {fmt(x + s)},{fmt(y)} {fmt(x)},{fmt(y + s)} {fmt(x - s)},{fmt(y)}"
        return f'<polygon points="{pts}" {paint}/>'
    pts = f"{fmt(x)},{fmt(y - size - 1)} {fmt(x - size - 1)},{fmt(y + size)} {fmt(x + size + 1)},{fmt(y + size)}"
    return f'<polygon points="{pts}" {paint}/>'


def svg(title: str, body: list, w: int, h: int) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def axes(x0: float, y0: float, x1: float, y1: float) -> str:
    return (f'<polyline points="{fmt(x0)},{fmt(y0)} {fmt(x0)},{fmt(y1)} {fmt(x1)},{fmt(y1)}" fill="none"'
            ' stroke="currentColor" stroke-width="1"/>')


def y_label(x: float, y: float, s: str) -> str:
    return f'<text transform="translate({fmt(x)} {fmt(y)}) rotate(-90)" text-anchor="middle">{s}</text>'


def whisker_h(x_lo: float, x_hi: float, y: float, color: str) -> list:
    """A horizontal interval bar with end caps."""
    return [line(x_lo, y, x_hi, y, color, 1.5), line(x_lo, y - 5, x_lo, y + 5, color, 1.5),
            line(x_hi, y - 5, x_hi, y + 5, color, 1.5)]


def whisker_v(x: float, y_lo: float, y_hi: float, color: str) -> list:
    """A vertical interval bar with end caps."""
    return [line(x, y_lo, x, y_hi, color, 1.2), line(x - 3, y_lo, x + 3, y_lo, color, 1.2),
            line(x - 3, y_hi, x + 3, y_hi, color, 1.2)]


# ---------------------------------------------------------------- the numbers

def measured() -> dict:
    """The recorded hardware campaign, as the runner wrote it."""
    return json.loads(DATA.read_text())["codes"]


def recorded_draws() -> int:
    """The draws of each weight the recorded hardware campaign made."""
    return json.loads(DATA.read_text())["draws"]


def exact_shares() -> dict:
    """The exact share of the tracked outcome, for each code and weight 3 and 4."""
    shares = {}
    for name in CODES:
        code = Code(name)
        counts = syndrome_counts(code.cols, 4)
        shares[name] = {w: exact_split(code, counts, w)[tracked(w)] / comb(N, w) for w in (3, 4)}
    return shares


def pct(x: float, places: int = 3) -> str:
    return f"{100 * x:.{places}f}"


# ---------------------------------------------------------------- the charts

def measured_chart() -> str:
    """Measured rates with their 95% Wilson intervals beside the exact shares, both weights, all three matrices."""
    runs = measured()
    drawn = recorded_draws()
    exact = exact_shares()
    x0, x1 = 170, 610
    panels = ((3, "Triples, share miscorrected, percent", 24, 55, 80, 5, 0),
              (4, "Quadruples, share undetected, percent", 238, 0.7, 1.2, 0.1, 1))
    body = []
    for w, head, top, lo_v, hi_v, step, places in panels:
        def x(v: float) -> float:
            return x0 + (x1 - x0) * (100 * v - lo_v) / (hi_v - lo_v)

        body.append(text(20, top, head, bold=True))
        for i, name in enumerate(CODES):
            r = top + 40 + 46 * i
            shape, color = STYLE[name]
            m = runs[name]["campaign"][f"weight_{w}"][tracked(w)]
            lo, hi = wilson(m, drawn)
            ex = exact[name][w]
            body.append(text(x0 - 12, r + 4, NAMES[name], "end", fill=color))
            body.append(text(x0 + 4, r - 14, f"measured {pct(m / drawn)}, interval {pct(lo)} to {pct(hi)},"
                             f" exact {pct(ex)}", fill=TEXT3, mono=True, size=11))
            body.append(line(x(ex), r - 8, x(ex), r + 10, "currentColor", 1, "3 3"))
            body.append(marker("triangle", x(ex), r + 14, "currentColor", 3.5, hollow=True))
            body.extend(whisker_h(x(lo), x(hi), r, color))
            body.append(marker(shape, x(m / drawn), r, color, 4))
        base = top + 40 + 46 * 2 + 28
        body.append(line(x0, base, x1, base))
        k = 0
        while lo_v + k * step <= hi_v + 1e-9:
            v = lo_v + k * step
            xv = x0 + (x1 - x0) * (v - lo_v) / (hi_v - lo_v)
            body.append(line(xv, base, xv, base + 4))
            body.append(text(xv, base + 17, f"{v:.{places}f}", "middle", mono=True, size=11))
            k += 1
    ky = 448
    body.extend(whisker_h(40, 64, ky - 4, "currentColor"))
    body.append(marker("circle", 52, ky - 4, "currentColor", 4))
    body.append(text(74, ky, f"measured on {drawn:,} draws, with its 95% Wilson interval", size=11))
    body.append(line(392, ky - 14, 392, ky + 2, "currentColor", 1, "3 3"))
    body.append(marker("triangle", 392, ky + 4, "currentColor", 3.5, hollow=True))
    body.append(text(404, ky, "exact share, every pattern counted", size=11))
    return svg("Measured miscorrection and undetected rates beside the exact shares", body, W, 466)


def convergence_chart() -> str:
    """The running share at each checkpoint with its Wilson interval, triples and quadruples."""
    runs = measured()
    exact = exact_shares()
    x0, x1 = 70, 470
    lo_e, hi_e = 2, 5
    pad = 0.12                            # a margin, in decades, so no point sits on an axis

    def x(n: float) -> float:
        return x0 + (x1 - x0) * (log10(n) - lo_e + pad) / (hi_e - lo_e + 2 * pad)

    offsets = {"hsiao": -6, "hsiao_dictionary": 0, "hamming": 6}
    panels = ((3, "Triples, running share miscorrected, percent", 24, 40, 90, 10, 100),
              (4, "Quadruples, running share undetected, percent", 290, 0.4, 2.8, 0.4, 1000))
    body = []
    for w, head, top, lo_v, hi_v, step, first in panels:
        y0, y1 = top + 22, top + 192

        def y(v: float) -> float:
            return y1 - (y1 - y0) * (100 * v - lo_v) / (hi_v - lo_v)

        body.append(text(20, top, head, bold=True))
        for name in CODES:
            shape, color = STYLE[name]
            body.append(line(x0, y(exact[name][w]), x1, y(exact[name][w]), color, 1, "3 3"))
        for name in CODES:
            shape, color = STYLE[name]
            for n, m in zip(CHECKPOINTS, runs[name]["running"][f"weight_{w}"]):
                if n < first:
                    continue
                lo, hi = wilson(m, n)
                xn = x(n) + offsets[name]
                body.extend(whisker_v(xn, y(min(hi, hi_v / 100)), y(lo), color))
                body.append(marker(shape, xn, y(m / n), color, 3.5, hollow=not lo <= exact[name][w] <= hi))
        body.append(axes(x0, y0, x1, y1))
        k = 0
        while lo_v + k * step <= hi_v + 1e-9:
            v = lo_v + k * step
            yv = y1 - (y1 - y0) * (v - lo_v) / (hi_v - lo_v)
            body.append(line(x0 - 4, yv, x0, yv))
            body.append(text(x0 - 8, yv + 4, f"{v:.{1 if step < 1 else 0}f}", "end", mono=True, size=11))
            k += 1
        for e in range(lo_e, hi_e + 1):
            body.append(line(x(10**e), y1, x(10**e), y1 + 4))
            body.append(text(x(10**e), y1 + 17, f"{10**e:,}", "middle", mono=True, size=11))
        if w == 3:
            yh, yd, ym = y(exact["hsiao"][3]), y(exact["hsiao_dictionary"][3]), y(exact["hamming"][3])
            body.append(text(x1 + 12, ym + 4, f"Hamming, exact {pct(exact['hamming'][3])}", fill=ORANGE, size=11))
            body.append(text(x1 + 12, yd - 4, f"dictionary, exact {pct(exact['hsiao_dictionary'][3])}",
                             fill=CYAN, size=11))
            body.append(text(x1 + 12, yh + 12, f"balanced, exact {pct(exact['hsiao'][3])}", fill=GREEN, size=11))
        else:
            yh, ym = y(exact["hsiao"][4]), y(exact["hamming"][4])
            body.append(text(x1 + 12, ym + 4, f"Hamming, exact {pct(exact['hamming'][4])}", fill=ORANGE, size=11))
            body.append(text(x1 + 12, yh + 4, f"balanced, exact {pct(exact['hsiao'][4])}", fill=GREEN, size=11))
            body.append(text(x1 + 12, yh + 18, f"dictionary, exact {pct(exact['hsiao_dictionary'][4])}",
                             fill=CYAN, size=11))
    body.append(text((x0 + x1) / 2, 528, "Draws so far, on a logarithmic scale", "middle"))
    ky = 552
    body.append(marker("circle", 76, ky - 4, "currentColor", 3.5))
    body.append(text(86, ky, "interval contains the exact share", size=11))
    body.append(marker("circle", 290, ky - 4, "currentColor", 3.5, hollow=True))
    body.append(text(300, ky, "interval misses it", size=11))
    body.append(line(420, ky - 4, 444, ky - 4, "currentColor", 1, "3 3"))
    body.append(text(452, ky, "exact share", size=11))
    return svg("Running shares of the campaign with their intervals, against the exact shares", body, W, 566)


def adjacent_chart() -> str:
    """Miscorrected windows of three adjacent cells, of 70, in seven layouts and under a random order.

    The three layouts the bench drove, natural, stride 11 and searched, are read
    from the recorded campaign. The other four, and the expected count under a
    random order, come from the lab script, which checks that its own counts for
    the first three equal the recorded ones. No search runs for extended
    Hamming, whose floor is 46 windows in any layout, so its searched bar repeats
    stride 11 and carries a footnote. Each bar carries its series' marker shape
    as well as its color, so color is never the only signal.
    """
    runs = measured()
    exact = exact_shares()
    windows = N - 2
    groups = []
    for kind in LAYOUTS:
        row = []
        for name in CODES:
            bench = runs[name]["adjacent"]
            if kind in bench:
                row.append(bench[kind]["miscorrected"])
            else:
                row.append(adjacent_outcomes(Code(name), layout(kind, name))["miscorrected"])
        groups.append((kind, row, False))
    groups.append(("random", [windows * exact[name][3] for name in CODES], True))
    x0, x1, y0, y1 = 110, 620, 56, 286
    top = 70
    slot = (x1 - x0) / len(groups)
    bar, gap = 15, 3

    def y(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = []
    for v in range(0, top + 1, 10):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, str(v), "end", mono=True, size=11))
    for g, (kind, row, hollow) in enumerate(groups):
        gx = x0 + g * slot + (slot - 3 * bar - 2 * gap) / 2
        for i, name in enumerate(CODES):
            shape, color = STYLE[name]
            bx = gx + i * (bar + gap)
            height = max(y(0) - y(row[i]), 0)
            if height > 0:
                if hollow:
                    body.append(rect(bx + 0.75, y(row[i]) + 0.75, bar - 1.5, height - 0.75, BG2, color, 1.5))
                else:
                    body.append(rect(bx, y(row[i]), bar, height, color))
                if height >= 20:
                    body.append(marker(shape, bx + bar / 2, y1 - 9, color if hollow else BG2, 3.5))
        label = kind.replace("stride", "stride ") if kind.startswith("stride") else kind
        cx = x0 + (g + 0.5) * slot
        body.append(text(cx, y1 + 17, label, "middle", size=11))
        if hollow:
            body.append(text(cx, y1 + 31, "expected", "middle", size=11, fill=TEXT3))
        for i, name in enumerate(CODES):
            value = f"{row[i]:.2f}" if hollow else str(row[i])
            if kind == "searched" and name == "hamming":
                value += "*"
            body.append(text(cx, y1 + 52 + 15 * i, value, "middle", mono=True, size=10))
    for i, name in enumerate(CODES):
        body.append(text(x0 - 8, y1 + 52 + 15 * i, NAMES[name], "end", fill=STYLE[name][1], size=11))
    body.append(line(x0, y(46), x1, y(46), "currentColor", 1, "3 3"))
    body.append(text(x0 + 3 * slot + 6, y(46) - 6, "Hamming floor in any layout, 46 windows", size=11))
    body.append(axes(x0, y0 - 6, x1, y1))
    body.append(y_label(62, (y0 + y1) / 2, "Windows miscorrected, of 70"))
    keys = (("balanced Hsiao", GREEN), ("dictionary Hsiao", CYAN), ("extended Hamming", ORANGE))
    kx = x0
    for (label, color), name in zip(keys, CODES):
        body.append(rect(kx, 14, 12, 16, color))
        body.append(marker(STYLE[name][0], kx + 6, 22, BG2, 3.5))
        body.append(text(kx + 18, 26, label, size=11))
        kx += 128
    body.append(rect(kx + 0.75, 16.75, 10.5, 10.5, BG2, "currentColor", 1.5))
    body.append(text(kx + 18, 26, "hollow, a random order", size=11))
    body.append(text(40, 392, "* No search runs for extended Hamming, so its searched bar repeats stride 11.",
                     size=11, fill=TEXT3))
    return svg("Miscorrected windows of three adjacent cells in seven layouts", body, W, 402)


def draws_chart() -> str:
    """Half-width of a 95% interval against the number of draws, from the normal rule and from two bounds."""
    x0, x1, y0, y1 = 80, 470, 40, 340
    lo_x, hi_x = 2, 6
    lo_y, hi_y = -2, 1.5

    def x(n: float) -> float:
        return x0 + (x1 - x0) * (log10(n) - lo_x) / (hi_x - lo_x)

    def y(points: float) -> float:
        return y1 - (y1 - y0) * (log10(points) - lo_y) / (hi_y - lo_y)

    shares = exact_shares()["hsiao"]
    p3, p4 = shares[3], shares[4]           # the balanced matrix's exact shares
    drawn = recorded_draws()
    rules = (("Chebyshev, any rate", RED, lambda n: 100 * half_chebyshev(n, DELTA), -2, "diamond"),
             ("Hoeffding, any rate", VIOLET, lambda n: 100 * half_hoeffding(n, DELTA), 0, "triangle"),
             (f"normal, rate {100 * p3:.2f}%", GREEN, lambda n: 100 * Z * sqrt(p3 * (1 - p3) / n), 2, "circle"),
             (f"normal, rate {100 * p4:.3f}%", CYAN, lambda n: 100 * Z * sqrt(p4 * (1 - p4) / n), 0, "square"))
    body = []
    body.append(line(x0, y(1), x1, y(1), TEXT3, 1, "3 3"))
    body.append(text(x1 - 4, y(1) - 6, "one point", "end", fill=TEXT3, size=11))
    body.append(line(x(drawn), y0, x(drawn), y1, TEXT3, 1, "3 3"))
    body.append(text(x(drawn), y0 - 8, f"{drawn:,} draws", "middle", fill=TEXT3, size=11))
    for label, color, half, nudge, shape in rules:
        pts = []
        steps = 80
        for i in range(steps + 1):
            n = 10 ** (lo_x + (hi_x - lo_x) * i / steps)
            pts.append((x(n), y(half(n))))
        body.append(polyline(pts, color, 1.5))
        body.append(marker(shape, x(drawn), y(half(drawn)), color, 3.5))
        body.append(text(x1 + 12, y(half(10 ** hi_x)) + 4 + nudge, label, fill=color, size=11))
    crossings = ((draws_chebyshev(0.01, DELTA), RED, "diamond"), (draws_hoeffding(0.01, DELTA), VIOLET, "triangle"),
                 (draws_normal(p3, 0.01), GREEN, "circle"))
    for n, color, shape in crossings:
        body.append(marker(shape, x(n), y(1), color, 3.5, hollow=True))
    body.append(axes(x0, y0, x1, y1))
    for e in range(lo_x, hi_x + 1):
        body.append(line(x(10**e), y1, x(10**e), y1 + 4))
        body.append(text(x(10**e), y1 + 17, f"{10**e:,}", "middle", mono=True, size=11))
    for v in (0.01, 0.1, 1, 10):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, f"{v:g}", "end", mono=True, size=11))
    body.append(text((x0 + x1) / 2, y1 + 40, "Draws, on a logarithmic scale", "middle"))
    body.append(y_label(22, (y0 + y1) / 2, "Half-width, percentage points"))
    return svg("Half-width of a 95 percent interval against the number of draws", body, W, 400)


CHARTS = {
    "fig-31-3-measured-vs-exact.svg": measured_chart,
    "fig-31-4-convergence.svg": convergence_chart,
    "fig-31-6-adjacent-layouts.svg": adjacent_chart,
    "fig-31-7-draws-needed.svg": draws_chart,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="folder to write the SVG files to (default figures)")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
