"""Write the four computed figures of chapter 11 of Error Correction Coding as SVG files.

    python3 python/ch11_charts.py
    python3 python/ch11_charts.py --out figures

Every plotted value comes from the functions in
ch11_project_channel_simulator.py, with the same seeds and trial counts the
chapter names, so a figure is corrected by fixing those functions and running
this script again, never by editing its output. Inside the hub that holds both
this repository and the book's website, the default output folder is the site's
figure folder for the book. In a clone of this repository alone it is a folder
named figures in the working directory. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import exp, log10, pi, sqrt

from ch11_project_channel_simulator import (BOX_MULLER_COUNT, BOX_MULLER_SEED, COVERAGE_EXACT, COVERAGE_LONG,
                                            SWEEP_DB, SWEEP_TRIALS, awgn_sweep, box_muller_samples,
                                            coverage_counts, hard_repetition_awgn_exact, run_experiments,
                                            uncoded_bpsk_exact, wilson)

W = 640
X0, X1 = 84, 604
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
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


class LogAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = log10(lo), log10(hi), px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (log10(v) - self.lo) / (self.hi - self.lo)


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
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def ring(x: float, y: float, stroke: str, r: float = 4) -> str:
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{BG2}" stroke="{stroke}"'
            ' stroke-width="1.5"/>')


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


class Plot:
    """A plot area of the standard width, with its own top, bottom and canvas height."""

    def __init__(self, y0: int, y1: int, h: int):
        self.y0, self.y1, self.h = y0, y1, h

    def frame(self, xlabel: str, ylabel: str) -> list:
        return [
            f'<polyline points="{X0},{self.y0} {X0},{self.y1} {X1},{self.y1}" fill="none" stroke="currentColor"'
            ' stroke-width="1"/>',
            text((X0 + X1) / 2, self.h - 14, xlabel, "middle"),
            f'<text transform="translate(20 {(self.y0 + self.y1) // 2}) rotate(-90)" text-anchor="middle">'
            f'{ylabel}</text>',
        ]

    def x_ticks(self, ax, values: list, labels: list) -> list:
        out = []
        for v, lab in zip(values, labels):
            px = ax(v)
            out.append(line(px, self.y1, px, self.y1 + 4))
            out.append(text(px, self.y1 + 20, lab, "middle", mono=True, size=11))
        return out

    def y_ticks(self, ax, values: list, labels: list) -> list:
        out = []
        for v, lab in zip(values, labels):
            py = ax(v)
            out.append(line(X0 - 4, py, X0, py))
            out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
            out.append(text(X0 - 8, py + 4, lab, "end", mono=True, size=11))
        return out

    def svg(self, title: str, body: list) -> str:
        return svg(title, self.h, body)


def svg(title: str, h: int, body: list) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {h}" width="{W}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "", width: float = 2) -> list:
    return [line(x, y - 4, x + 24, y - 4, color, width, dash), text(x + 32, y, label, size=12)]


# The charts --------------------------------------------------------------------------

BIN = 0.25


def box_muller_histogram() -> str:
    """Histogram of the chapter's 100,000 Box-Muller samples against the standard normal density."""
    zs = box_muller_samples(BOX_MULLER_COUNT, BOX_MULLER_SEED)
    edges = [-4 + BIN * i for i in range(33)]
    counts = [0] * 32
    for z in zs:
        if -4 <= z < 4:
            counts[int((z + 4) // BIN)] += 1
    outside = BOX_MULLER_COUNT - sum(counts)
    plot = Plot(30, 300, 350)
    xa = LinAxis(-4, 4, X0, X1)
    ya = LinAxis(0, 0.45, plot.y1, plot.y0)
    body = plot.y_ticks(ya, [0, 0.1, 0.2, 0.3, 0.4], ["0", "0.1", "0.2", "0.3", "0.4"])
    for i, c in enumerate(counts):
        height = c / (BOX_MULLER_COUNT * BIN)
        x, y = xa(edges[i]), ya(height)
        body.append(f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(xa(edges[i + 1]) - x)}" '
                    f'height="{fmt(plot.y1 - y)}" fill="{BG2}" stroke="{TEXT3}" stroke-width="1"/>')
    xs = [-4 + 8 * i / 200 for i in range(201)]
    body.append(polyline([(xa(x), ya(exp(-x * x / 2) / sqrt(2 * pi))) for x in xs], GREEN, 2))
    body += plot.x_ticks(xa, list(range(-4, 5)), [str(v) for v in range(-4, 5)])
    body += plot.frame("sample value z", "fraction of samples per unit of z")
    body.append(f'<rect x="430" y="40" width="20" height="12" fill="{BG2}" stroke="{TEXT3}" stroke-width="1"/>')
    body.append(text(458, 50, "samples, bins of 0.25", size=12))
    body += legend_line(428, 76, GREEN, "standard normal density")
    peak = max(range(32), key=lambda i: counts[i])
    title = (f"Histogram of {BOX_MULLER_COUNT:,} Box-Muller samples made from random.Random({BOX_MULLER_SEED}), in "
             f"32 bins of width 0.25 from -4 to 4, each bar scaled to count divided by 25,000 so its area is the "
             f"fraction of samples in the bin. The green curve is the standard normal density, which peaks at "
             f"{1 / sqrt(2 * pi):.4f}. The tallest bar, {edges[peak]:g} to {edges[peak + 1]:g}, holds "
             f"{counts[peak]:,} samples, and {outside} samples fall outside the plotted range.")
    return plot.svg(title, body)


def interval_coverage() -> str:
    """The first 50 Wilson intervals of the coverage experiment, misses in red."""
    runs = 50
    counts = coverage_counts(COVERAGE_LONG, runs)
    plot = Plot(60, 320, 370)
    xa = LinAxis(0.5, runs + 0.5, X0, X1)
    ya = LinAxis(0.005, 0.065, plot.y1, plot.y0)
    body = plot.y_ticks(ya, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
                        ["0.01", "0.02", "0.03", "0.04", "0.05", "0.06"])
    body.append(line(X0, ya(COVERAGE_EXACT), X1, ya(COVERAGE_EXACT), GREEN, 2))
    misses = []
    for i, k in enumerate(counts):
        lo, hi = wilson(k, COVERAGE_LONG)
        hit = lo <= COVERAGE_EXACT <= hi
        color = "currentColor" if hit else RED
        x = xa(i + 1)
        body.append(line(x, ya(lo), x, ya(hi), color, 2 if not hit else 1.5))
        body.append(dot(x, ya(k / COVERAGE_LONG), color, 2.5))
        if not hit:
            misses.append((i + 1, k))
    body += plot.x_ticks(xa, [1, 10, 20, 30, 40, 50], ["1", "10", "20", "30", "40", "50"])
    body += plot.frame("run number, which is also its seed, with 1,000 words in each run", "word error rate")
    body += legend_line(96, 34, GREEN, "exact value 0.028")
    body += legend_line(290, 34, "currentColor", "interval covers", width=1.5)
    body += legend_line(460, 34, RED, "interval misses")
    miss_text = " and ".join(f"run {r} with {k} failures" for r, k in misses)
    title = (f"The 95 percent Wilson intervals of the first {runs} runs of the coverage experiment, each run "
             f"sending 1,000 words of the three-fold repetition code over a binary symmetric channel with p = 0.1 "
             f"from random.Random(seed) with seeds 1 to {runs}. A dot marks each run's estimate and a vertical bar "
             f"its interval, and the green line marks the exact word error rate 0.028. {len(misses)} of the {runs} "
             f"intervals miss it, {miss_text}, drawn in red.")
    return plot.svg(title, body)


def awgn_check() -> str:
    """The Gaussian sweep: simulated points with Wilson bars against the closed forms."""
    sweep = awgn_sweep()
    plot = Plot(30, 300, 350)
    xa = LinAxis(0, 7, X0 + 16, X1 - 16)
    ya = LogAxis(1e-4, 1, plot.y1, plot.y0)
    body = plot.y_ticks(ya, [1e-4, 1e-3, 1e-2, 1e-1, 1], ["1e-4", "1e-3", "0.01", "0.1", "1"])
    dbs = [7 * i / 140 for i in range(141)]
    body.append(polyline([(xa(d), ya(uncoded_bpsk_exact(d))) for d in dbs], GREEN, 2))
    body.append(polyline([(xa(d), ya(hard_repetition_awgn_exact(3, d))) for d in dbs], ORANGE, 2))
    for db, ku, kh, ks in sweep:
        for k, dx, kind in ((ku, -6, "uncoded"), (ks, 6, "soft"), (kh, 0, "hard")):
            lo, hi = wilson(k, SWEEP_TRIALS)
            x = xa(db) + dx
            color = {"uncoded": "currentColor", "soft": GREEN, "hard": ORANGE}[kind]
            body.append(line(x, ya(lo), x, ya(hi), color, 1.5))
            body.append(dot(x, ya(k / SWEEP_TRIALS), color, 3) if kind == "uncoded"
                        else ring(x, ya(k / SWEEP_TRIALS), color, 3.5))
    body += plot.x_ticks(xa, SWEEP_DB, [str(d) for d in SWEEP_DB])
    body += plot.frame("energy per data bit over noise density, Eb/N0 in dB", "bit error rate")
    body += legend_line(100, 200, GREEN, "closed form, uncoded and soft R3")
    body += legend_line(100, 220, ORANGE, "closed form, hard R3")
    body += [dot(112, 236, "currentColor", 3), text(132, 240, "simulated uncoded", size=12)]
    body += [ring(112, 256, GREEN, 3.5), text(132, 260, "simulated soft R3", size=12)]
    body += [ring(112, 276, ORANGE, 3.5), text(132, 280, "simulated hard R3", size=12)]
    k7 = sweep[-1]
    title = (f"Bit error rate of BPSK over the Gaussian channel from 0 to 7 dB, uncoded and with the three-fold "
             f"repetition code decoded hard and soft, on a logarithmic axis. Each point is {SWEEP_TRIALS:,} words "
             f"simulated with seed 100 plus the value in dB, drawn with its 95 percent Wilson interval. The green "
             f"curve is the closed form Q(sqrt(2Eb/N0)), which uncoded and soft decoding share, and the orange "
             f"curve is the closed form for hard decoding, which lies above it at every point. At 7 dB the "
             f"simulation counted {k7[1]} uncoded, {k7[3]} soft and {k7[2]} hard failures.")
    return plot.svg(title, body)


SHORT = {
    "R3 on BSC(0.1)": "R3 on BSC(0.1)",
    "R3 on BEC(0.2)": "R3 on BEC(0.2)",
    "parity (8,7) on BSC(0.01), undetected": "parity, undetected, BSC(0.01)",
    "parity (8,7) on BEC(0.1), unfilled": "parity, unfilled, BEC(0.1)",
    "uncoded BPSK at 4 dB": "uncoded BPSK, 4 dB",
    "hard R3 at 4 dB": "hard R3, 4 dB",
    "soft R3 at 4 dB": "soft R3, 4 dB",
}


def check_summary() -> str:
    """A forest plot: each experiment's estimate and Wilson interval divided by its closed form."""
    rows = run_experiments()
    top, step = 62, 40
    h = top + step * len(rows) + 46
    xa = LinAxis(0.8, 1.2, 280, 604)
    axis_y = top + step * (len(rows) - 1) + 26
    body = [text(16, 30, "experiment", size=12), text(442, 30, "simulated divided by exact", "middle", size=12)]
    for v in (0.8, 0.9, 1.0, 1.1, 1.2):
        x = xa(v)
        body.append(line(x, 44, x, axis_y, TEXT3, 1, opacity="0.35"))
        body.append(line(x, axis_y, x, axis_y + 4))
        body.append(text(x, axis_y + 20, fmt(v) if v != 1.0 else "1", "middle", mono=True, size=11))
    body.append(line(280, axis_y, 604, axis_y))
    body.append(line(xa(1.0), 44, xa(1.0), axis_y, GREEN, 2))
    misses = []
    for i, (name, seed, trials, k, lo, hi, exact) in enumerate(rows):
        y = top + step * i
        hit = lo <= exact <= hi
        color = "currentColor" if hit else RED
        body.append(text(16, y + 1, SHORT[name], size=12))
        body.append(text(16, y + 15, f"seed {seed}, N = {trials:,}", fill=TEXT3, mono=True, size=10))
        body.append(line(xa(lo / exact), y, xa(hi / exact), y, color, 2))
        body.append(line(xa(lo / exact), y - 5, xa(lo / exact), y + 5, color, 1.5))
        body.append(line(xa(hi / exact), y - 5, xa(hi / exact), y + 5, color, 1.5))
        body.append(dot(xa(k / trials / exact), y, color, 3.5))
        if not hit:
            body.append(text(600, y + 4, "misses", "end", RED, size=12))
            misses.append((name, trials, k, exact))
    body.append(text(442, h - 12, "1 means the simulation agrees with the closed form", "middle", size=12))
    name, trials, k, exact = misses[0]
    mean = trials * exact
    z = (k - mean) / sqrt(trials * exact * (1 - exact))
    title = (f"The chapter's {len(rows)} seeded simulations, each checked against its closed form. Each row shows "
             f"the estimate K/N divided by the exact value, with its 95 percent Wilson interval divided the same "
             f"way, so agreement means the bar crosses the green line at 1. {len(rows) - len(misses)} intervals cross it. The one "
             f"that misses is {name}, {k} failures in {trials:,} words where {mean:.0f} were expected, "
             f"{abs(z):.2f} standard deviations low.")
    return svg(title, h, body)


CHARTS = {
    "fig-11-4-box-muller-histogram.svg": box_muller_histogram,
    "fig-11-6-interval-coverage.svg": interval_coverage,
    "fig-11-7-awgn-check.svg": awgn_check,
    "fig-11-8-check-summary.svg": check_summary,
}


def default_out() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    site = os.path.join(here, "..", "..", "..", "Website-Building-Projects", "dranirbanchakraborty.com", "public",
                        "figures", "error-correction-coding")
    return os.path.normpath(site) if os.path.isdir(site) else "figures"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=default_out(), help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
