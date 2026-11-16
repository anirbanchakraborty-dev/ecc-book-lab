"""Write the four computed figures of chapter 28 of Error Correction Coding as SVG files.

    python3 python/ch28_charts.py --out figures

Every drawn value comes from ch28_verifying_ecc_hardware.py, the decoder model
and the seeded campaigns included, so a figure is corrected by fixing those
functions and running this script again, never by editing its output. The SVG
follows the book's figure style sheet: real text in IBM Plex, colors from the
page's tokens with paper fallbacks, and no element that the site's figure check
would refuse. Each series is told apart by a shape or a label as well as by its
color.
"""

import argparse
import os
import random
from math import comb, log10

from ch28_verifying_ecc_hardware import (
    K, N, RATE, SEED, WORDS, YEAR, draws, envelope, envelope_size, harmonic, wrong_cases)

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
CASES = WORDS * len(envelope())          # 23,552
CAMPAIGNS_PER_POINT = 10_000             # campaigns behind each seeded point, as in the chapter's run
T_GRID = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)


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


def marker(shape: str, x: float, y: float, color: str, size: float = 4.5) -> str:
    """A filled point drawn as a circle, a square, a diamond or a triangle, so the shape names the series."""
    if shape == "circle":
        return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(size)}" fill="{color}"/>'
    if shape == "square":
        return rect(x - size, y - size, 2 * size, 2 * size, color)
    if shape == "diamond":
        s = size + 1
        pts = f"{fmt(x)},{fmt(y - s)} {fmt(x + s)},{fmt(y)} {fmt(x)},{fmt(y + s)} {fmt(x - s)},{fmt(y)}"
        return f'<polygon points="{pts}" fill="{color}"/>'
    pts = f"{fmt(x)},{fmt(y - size - 1)} {fmt(x - size - 1)},{fmt(y + size)} {fmt(x + size + 1)},{fmt(y + size)}"
    return f'<polygon points="{pts}" fill="{color}"/>'


def power(e: int) -> str:
    """Ten to the power e, with the exponent raised."""
    return f'10<tspan dy="-5" font-size="8">{e}</tspan>'


def svg(title: str, body: list, w: int, h: int) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def axes(x0: float, y0: float, x1: float, y1: float) -> str:
    return (f'<polyline points="{fmt(x0)},{fmt(y0)} {fmt(x0)},{fmt(y1)} {fmt(x1)},{fmt(y1)}" fill="none"'
            ' stroke="currentColor" stroke-width="1"/>')


def y_label(x: float, y: float, s: str) -> str:
    return f'<text transform="translate({fmt(x)} {fmt(y)}) rotate(-90)" text-anchor="middle">{s}</text>'


# ---------------------------------------------------------------- the charts

def weight_chart() -> str:
    """How many of the 8,192 error patterns of the (13,8) code have each weight."""
    h = 400
    x0, x1, y0, y1 = 70, 610, 92, 330
    top = 1800
    counts = [comb(N, w) for w in range(N + 1)]

    def y(v: float) -> float:
        return y1 - (y1 - y0) * v / top

    body = []
    slot = (x1 - x0) / (N + 1)
    bar = slot * 0.62
    for v in (500, 1000, 1500):
        body.append(line(x0, y(v), x1, y(v), TEXT3, 1, "3 3"))
    for w, c in enumerate(counts):
        bx = x0 + w * slot + (slot - bar) / 2
        color = GREEN if w <= 2 else (ORANGE if w == 3 else TEXT3)
        height = max(y(0) - y(c), 1.5)
        body.append(rect(bx, y(0) - height, bar, height, color))
        body.append(text(bx + bar / 2, y(c) - 6, f"{c:,}", "middle", mono=True, size=10))
        body.append(text(bx + bar / 2, y1 + 17, str(w), "middle", mono=True, size=11))
    # brackets naming the two regions, so the color is not the only signal
    inside = sum(counts[:3])
    ex = x0 + 3 * slot
    body.append(polyline([(x0 + 4, 64), (x0 + 4, 58), (ex - 4, 58), (ex - 4, 64)], GREEN, 1.5))
    body.append(text((x0 + ex) / 2, 50, f"inside the envelope, {inside} patterns", "middle", fill=GREEN, size=12))
    body.append(polyline([(ex + 4, 64), (ex + 4, 58), (x1 - 4, 58), (x1 - 4, 64)], TEXT3, 1.5))
    body.append(text((ex + x1) / 2, 50, f"outside the envelope, {2**N - inside:,} patterns", "middle",
                     fill=TEXT3, size=12))
    body.append(text(x0 + 3.5 * slot, y(counts[3]) - 22, "triples", "middle", fill=ORANGE, size=11))
    body.append(axes(x0, y0 - 10, x1, y1))
    for v in (0, 500, 1000, 1500):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, f"{v:,}", "end", mono=True, size=11))
    body.append(text((x0 + x1) / 2, h - 30, "Weight of the error pattern, the number of bits flipped", "middle"))
    body.append(y_label(20, (y0 + y1) / 2, "Patterns of that weight"))
    return svg("The 8,192 error patterns of the Hsiao (13,8) code, by weight", body, W, h)


def case_count_chart() -> str:
    """Cases a complete check needs, for four Hsiao codes, on a logarithmic axis."""
    codes = [(13, 8), (22, 16), (39, 32), (72, 64)]
    series = [
        ("envelope on one codeword", "circle", CYAN, lambda n, k: envelope_size(n)),
        ("envelope, every data word", "square", GREEN, lambda n, k: 2**k * envelope_size(n)),
        ("every pattern, every data word", "triangle", RED, lambda n, k: 2**k * 2**n),
    ]
    h = 450
    x0, x1, y0, y1 = 76, 470, 80, 390
    top = 42

    def y(v: float) -> float:
        return y1 - (y1 - y0) * log10(v) / top

    def x(i: int) -> float:
        return x0 + 40 + i * (x1 - x0 - 80) / (len(codes) - 1)

    body = []
    # what a checker running 10^9 cases a second gets through, an assumed rate
    marks = [("one second", RATE), ("one day", RATE * 86400), ("one year", RATE * YEAR),
             ("a million years", RATE * YEAR * 10**6)]
    for label, v in marks:
        body.append(line(x0, y(v), x1, y(v), TEXT3, 1, "3 3"))
        body.append(text(x1 + 6, y(v) + 4, label, fill=TEXT3, size=11))
    body.append(text(x1 + 6, y(marks[-1][1]) - 30, "time taken at an", fill=TEXT3, size=11))
    body.append(text(x1 + 6, y(marks[-1][1]) - 16, f"assumed {power(9)} a second", fill=TEXT3, size=11))
    for label, shape, color, f in series:
        pts = [(x(i), y(f(n, k))) for i, (n, k) in enumerate(codes)]
        body.append(polyline(pts, color, 1.5))
        for px, py in pts:
            body.append(marker(shape, px, py, color))
    for i, (n, k) in enumerate(codes):
        body.append(text(x(i), y1 + 18, f"({n},{k})", "middle", mono=True, size=11))
    body.append(axes(x0, y0, x1, y1))
    for e in range(0, top + 1, 6):
        body.append(line(x0 - 4, y(10**e), x0, y(10**e)))
        body.append(text(x0 - 8, y(10**e) + 4, power(e), "end", mono=True, size=11))
    for (label, shape, color, _), (lx, ly) in zip(series, ((84, 22), (300, 22), (84, 44))):
        body.append(marker(shape, lx + 6, ly - 4, color))
        body.append(text(lx + 18, ly, label, size=12))
    body.append(text((x0 + x1) / 2, h - 24, "Hsiao SEC-DED code (n,k)", "middle"))
    body.append(y_label(20, (y0 + y1) / 2, "Cases to check"))
    return svg("How many cases a complete check of a SEC-DED decoder needs", body, W, h)


def coupon_chart() -> str:
    """Distinct envelope patterns a seeded campaign has drawn, against the expected count."""
    n = len(envelope())
    t_max = 700
    rng = random.Random(SEED)
    ks = draws(rng, t_max)
    seen = set()
    path = []
    all_at = None
    for t, k in enumerate(ks, start=1):
        seen.add(k % n)
        path.append((t, len(seen)))
        if all_at is None and len(seen) == n:
            all_at = t
    h = 400
    x0, x1, y0, y1 = 70, 600, 60, 330

    def x(t: float) -> float:
        return x0 + (x1 - x0) * t / t_max

    def y(v: float) -> float:
        return y1 - (y1 - y0) * v / 100

    body = []
    body.append(line(x0, y(n), x1, y(n), TEXT3, 1, "3 3"))
    body.append(text(x1 - 4, y(n) - 6, f"all {n} patterns", "end", fill=TEXT3, size=11))
    expected = [(x(t), y(n * (1 - (1 - 1 / n) ** t))) for t in range(0, t_max + 1, 5)]
    body.append(polyline(expected, CYAN, 1.5, "3 3"))
    seeded = [(x(0), y(0))]
    for t, s in path:
        seeded.append((x(t), y(s)))
    body.append(polyline(seeded, GREEN, 1.5))
    mean = n * harmonic(n)
    body.append(line(x(mean), y(0), x(mean), y(n), CYAN, 1))
    body.append(text(x(mean) + 6, y(30), f"{mean:.1f}", fill=CYAN, mono=True, size=11))
    body.append(text(x(mean) + 6, y(30) + 14, "draws expected", fill=CYAN, size=11))
    body.append(text(x(mean) + 6, y(30) + 28, "for all 92", fill=CYAN, size=11))
    body.append(marker("circle", x(all_at), y(n), GREEN))
    body.append(text(x(all_at) - 8, y(n) - 8, f"seed 28: all {n} by draw {all_at}", "end", fill=GREEN, size=11))
    body.append(axes(x0, y0, x1, y1))
    for t in range(0, t_max + 1, 100):
        body.append(line(x(t), y1, x(t), y1 + 4))
        body.append(text(x(t), y1 + 18, str(t), "middle", mono=True, size=11))
    for v in (0, 20, 40, 60, 80, 100):
        body.append(line(x0 - 4, y(v), x0, y(v)))
        body.append(text(x0 - 8, y(v) + 4, str(v), "end", mono=True, size=11))
    body.append(line(x0 + 10, 30, x0 + 40, 30, GREEN, 1.5))
    body.append(text(x0 + 48, 34, "seeded campaign, seed 28", size=12))
    body.append(line(x0 + 250, 30, x0 + 280, 30, CYAN, 1.5, "3 3"))
    body.append(text(x0 + 288, 34, f'expected, 92 (1 - (1 - 1/92)<tspan dy="-5" font-size="8">T</tspan><tspan dy="5">)</tspan>', size=12))
    body.append(text((x0 + x1) / 2, h - 30, "Draws T, each a uniform case from the 23,552", "middle"))
    body.append(y_label(20, (y0 + y1) / 2, "Distinct patterns drawn"))
    return svg("How many of the 92 envelope patterns a random campaign has drawn", body, W, h)


def seeded_miss_fractions() -> dict:
    """For each mutant, the fraction of seeded campaigns of each size that drew no wrong case.

    Each size starts again from random.Random(SEED), so the point at T = 100 is the
    chapter's own run of 10,000 campaigns of 100 draws.
    """
    wrong = {m: {k for k, bad in enumerate(wrong_cases(m)) if bad} for m in (1, 2, 3, 4)}
    missed = {m: [] for m in wrong}
    for t in T_GRID:
        rng = random.Random(SEED)
        counts = {m: 0 for m in wrong}
        for _ in range(CAMPAIGNS_PER_POINT):
            batch = set(draws(rng, t))
            for m, bad in wrong.items():
                if batch.isdisjoint(bad):
                    counts[m] += 1
        for m in wrong:
            missed[m].append(counts[m] / CAMPAIGNS_PER_POINT)
    return missed


def miss_chart() -> str:
    """The chance a campaign of T draws misses a bug with m wrong cases, closed form and seeded."""
    sizes = {m: sum(wrong_cases(m)) for m in (1, 2, 3, 4)}
    seeded = seeded_miss_fractions()
    series = ((1, "diamond", RED), (2, "square", ORANGE), (3, "circle", GREEN), (4, "triangle", VIOLET))
    h = 420
    x0, x1, y0, y1 = 70, 520, 60, 340
    lo, hi = 0, 4

    def x(t: float) -> float:
        return x0 + (x1 - x0) * (log10(t) - lo) / (hi - lo)

    def y(p: float) -> float:
        return y1 - (y1 - y0) * p

    body = []
    for p in (0.25, 0.5, 0.75):
        body.append(line(x0, y(p), x1, y(p), TEXT3, 1, "3 3"))
    for m, color, dash in [(sizes[k], c, "") for k, _, c in series] + [(1, TEXT3, "3 3")]:
        pts = []
        steps = 160
        for i in range(steps + 1):
            t = 10 ** (lo + (hi - lo) * i / steps)
            pts.append((x(t), y((1 - m / CASES) ** t)))
        body.append(polyline(pts, color, 1.5, dash))
    for mutant, shape, color in series:
        for t, p in zip(T_GRID, seeded[mutant]):
            if p > 0:
                body.append(marker(shape, x(t), y(p), color, 4))
    # a key for the four mutants, in the empty band between the 0.5 and 0.75 guides
    for i, (mutant, shape, color) in enumerate(series):
        ky = y(0.75) + 18 + 16 * i
        body.append(line(x(250), ky - 4, x(250) + 22, ky - 4, color, 1.5))
        body.append(marker(shape, x(250) + 11, ky - 4, color, 4))
        body.append(text(x(250) + 30, ky, f"mutant {mutant}, m = {sizes[mutant]:,}", fill=color, size=11))
    body.append(text(x1 + 6, y((1 - 1 / CASES) ** (10**hi)) + 4, "one wrong", fill=TEXT3, size=11))
    body.append(text(x1 + 6, y((1 - 1 / CASES) ** (10**hi)) + 18, "case, m = 1", fill=TEXT3, size=11))
    body.append(axes(x0, y0, x1, y1))
    for e in range(lo, hi + 1):
        body.append(line(x(10**e), y1, x(10**e), y1 + 4))
        body.append(text(x(10**e), y1 + 18, f"{10**e:,}", "middle", mono=True, size=11))
    for p in (0, 0.25, 0.5, 0.75, 1):
        body.append(line(x0 - 4, y(p), x0, y(p)))
        body.append(text(x0 - 8, y(p) + 4, f"{p:.2f}", "end", mono=True, size=11))
    body.append(line(x0 + 10, 30, x0 + 40, 30, "currentColor", 1.5))
    body.append(text(x0 + 48, 34, f'(1 - m/23,552)<tspan dy="-5" font-size="8">T</tspan>', mono=True, size=12))
    body.append(marker("circle", x0 + 220, 30, "currentColor", 4))
    body.append(text(x0 + 232, 34, f"share of {CAMPAIGNS_PER_POINT:,} seeded campaigns that missed", size=12))
    body.append(text((x0 + x1) / 2, h - 30, "Campaign size T, draws, on a logarithmic axis", "middle"))
    body.append(y_label(20, (y0 + y1) / 2, "Chance the campaign misses the bug"))
    return svg("The chance a random campaign misses a planted bug", body, W, h)


CHARTS = {
    "fig-28-2-weight-classes.svg": weight_chart,
    "fig-28-3-case-counts.svg": case_count_chart,
    "fig-28-5-coupon-collector.svg": coupon_chart,
    "fig-28-6-miss-probability.svg": miss_chart,
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
