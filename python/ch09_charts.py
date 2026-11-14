"""Write the four computed figures of chapter 9 of Error Correction Coding as SVG files.

    python3 python/ch09_charts.py --out figures

Every plotted value comes from the functions in
ch09_error_control_strategies_and_performance.py, so a figure is corrected by
fixing those functions and running this script again, never by editing its
output. The Monte Carlo figures replay the same seeded runs that the companion
script checks. The SVG follows the book's figure style sheet: real text in IBM
Plex, colors from the page's tokens with paper fallbacks, and no element that the
site's figure check would refuse.
"""

import argparse
import os
from math import log10

from ch09_error_control_strategies_and_performance import (MAIN_P, MAIN_TRIALS, eta_gbn, eta_sr, eta_sw,
                                                           frame_error, hard_rep_ber, hybrid_failure, main_run,
                                                           mc_points, required_repetition, required_uncoded,
                                                           uncoded_ber, wilson, word_error)

W = 640
X0, X1 = 84, 604
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


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def polygon(points: list, fill: str) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polygon points="{pts}" fill="{fill}" stroke="none"/>'


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
        head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {self.h}" width="{W}" height="{self.h}">'
        return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "", width: float = 2) -> list:
    return [line(x, y - 4, x + 24, y - 4, color, width, dash), text(x + 32, y, label, size=12)]


def log_grid(lo: float, hi: float, count: int) -> list:
    a, b = log10(lo), log10(hi)
    return [10 ** (a + (b - a) * i / count) for i in range(count + 1)]


# The charts --------------------------------------------------------------------------

LINK_N, LINK_A, LINK_RATE, HYBRID_RATE, HYBRID_T = 1000, 10, 0.968, 0.938, 3


def protocol_curves(p: float) -> dict:
    P = frame_error(LINK_N, p)
    return {
        "sw": eta_sw(P, LINK_A, LINK_RATE),
        "gbn": eta_gbn(P, 1 + LINK_A, LINK_RATE),
        "sr": eta_sr(P, LINK_RATE),
        "hybrid": (1 - hybrid_failure(LINK_N, HYBRID_T, p)) * HYBRID_RATE,
    }


def arq_throughput() -> str:
    """Throughput efficiency of the three ARQ protocols and type-I hybrid ARQ against the bit error probability."""
    plot = Plot(64, 364, 440)
    x = LogAxis(1e-6, 1e-2, X0, X1)
    y = LinAxis(0, 1, plot.y1, plot.y0)
    body = plot.frame("bit error probability p, log scale", "throughput efficiency")
    body += plot.y_ticks(y, [0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    body += plot.x_ticks(x, [1e-6, 1e-5, 1e-4, 1e-3, 1e-2], ["0.000001", "0.00001", "0.0001", "0.001", "0.01"])
    body.append(line(x(1e-4), plot.y0, x(1e-4), plot.y1, TEXT3, 1, "3 3"))
    grid = log_grid(1e-6, 1e-2, 200)
    curves = [protocol_curves(p_) for p_ in grid]
    styles = [("sw", ORANGE, "", "stop-and-wait"), ("gbn", CYAN, "3 3", "go-back-N"),
              ("sr", GREEN, "", "selective repeat"), ("hybrid", "currentColor", "6 3", "type-I hybrid ARQ")]
    at = protocol_curves(1e-4)
    for key, color, dash, _ in styles:
        body.append(polyline([(x(p_), y(c[key])) for p_, c in zip(grid, curves)], color, 2, dash))
        if key != "hybrid":
            body.append(dot(x(1e-4), y(at[key]), color, 4))
    for i, (key, color, dash, label) in enumerate(styles):
        lx = X0 + 260 * (i % 2)
        ly = 24 + 22 * (i // 2)
        body += legend_line(lx, ly, color, f"{label}, {at[key]:.3f} at p = 0.0001", dash)
    return plot.svg("Throughput efficiency against the bit error probability p, from 0.000001 to 0.01 on a log "
                    "scale, for 1000-bit frames with 32 check bits on a link where the round trip lasts ten frame "
                    f"times. At p = 0.0001, marked by a dashed vertical line, stop-and-wait reaches {at['sw']:.3f}, "
                    f"go-back-N {at['gbn']:.3f} and selective repeat {at['sr']:.3f}. Stop-and-wait never rises "
                    f"above {protocol_curves(1e-6)['sw']:.3f}. A type-I hybrid scheme that corrects up to three "
                    f"errors per frame holds near {HYBRID_RATE} until p passes 0.001, then falls, while the three "
                    "plain protocols fall toward zero above about p = 0.001.", body)


def coding_gain_chart() -> str:
    """Bit error rate against Eb/N0 for uncoded BPSK and hard-decision R3 and R5."""
    plot = Plot(24, 324, 400)
    x = LinAxis(0, 14, X0, X1)
    ymin = 1e-9
    y = LogAxis(ymin, 1, plot.y1, plot.y0)
    body = plot.frame("energy per information bit Eb/N0, dB", "bit error rate, log scale")
    body += plot.y_ticks(y, [1e-9, 1e-6, 1e-3, 1], ["1e-9", "1e-6", "0.001", "1"])
    body += plot.x_ticks(x, list(range(0, 15, 2)), [str(v) for v in range(0, 15, 2)])
    body.append(line(X0, y(1e-6), X1, y(1e-6), TEXT3, 1, "3 3"))
    grid = [i * 0.05 for i in range(281)]
    unc = required_uncoded(1e-6)[1]
    r3 = required_repetition(3, 1e-6)[3]
    r5 = required_repetition(5, 1e-6)[3]
    styles = [(uncoded_ber, GREEN, "", unc, "uncoded BPSK and soft-decision R3"),
              (lambda e: hard_rep_ber(3, e), ORANGE, "", r3, "hard-decision R3"),
              (lambda e: hard_rep_ber(5, e), CYAN, "3 3", r5, "hard-decision R5")]
    for f, color, dash, at, _ in styles:
        pts = [(x(e), y(f(e))) for e in grid if f(e) >= ymin]
        body.append(polyline(pts, color, 2, dash))
        body.append(dot(x(at), y(1e-6), color, 4))
    for i, (_, color, dash, at, label) in enumerate(styles):
        body += legend_line(100, 250 + 22 * i, color, f"{label}, {at:.2f} dB", dash)
    return plot.svg("Bit error rate on a log scale against the energy per information bit Eb/N0 from 0 to 14 dB, "
                    "for uncoded BPSK, which soft-decision R3 matches exactly, and for hard-decision R3 and R5. A "
                    f"dashed line marks a bit error rate of one in a million, which uncoded BPSK reaches at {unc:.2f} "
                    f"dB, hard-decision R3 at {r3:.2f} dB and hard-decision R5 at {r5:.2f} dB, so both repetition "
                    "curves lie to the right of the uncoded one and their coding gain is negative.", body)


def running_estimate() -> str:
    """The running estimate of the R3 word error rate with its Wilson interval, from 100 to 10,000 trials."""
    outcomes = main_run()
    exact = word_error(3, 1, MAIN_P)
    plot = Plot(24, 324, 400)
    x = LogAxis(100, MAIN_TRIALS, X0, X1)
    y = LinAxis(0, 0.06, plot.y1, plot.y0)
    body = plot.frame("number of trials N, log scale", "estimated word error rate")
    body += plot.y_ticks(y, [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
                         ["0", "0.01", "0.02", "0.03", "0.04", "0.05", "0.06"])
    body += plot.x_ticks(x, [100, 300, 1000, 3000, 10000], ["100", "300", "1000", "3000", "10000"])
    running, k = [], 0
    for n, o in enumerate(outcomes, start=1):
        k += o
        running.append(k)
    ns = sorted({int(round(v)) for v in log_grid(100, MAIN_TRIALS, 400)})
    upper = [(x(n), y(wilson(running[n - 1], n)[1])) for n in ns]
    lower = [(x(n), y(wilson(running[n - 1], n)[0])) for n in ns]
    body.append(polygon(upper + lower[::-1], BG2))
    body.append(polyline(upper, CYAN, 1.5, "3 3"))
    body.append(polyline(lower, CYAN, 1.5, "3 3"))
    body.append(line(X0, y(exact), X1, y(exact), "currentColor", 1, "6 3"))
    body.append(polyline([(x(n), y(running[n - 1] / n)) for n in ns], GREEN, 2))
    final = running[-1] / MAIN_TRIALS
    body.append(dot(x(MAIN_TRIALS), y(final), GREEN, 4))
    entries = [(GREEN, "", 2, "running estimate K/N"), (CYAN, "3 3", 1.5, "95% Wilson interval"),
               ("currentColor", "6 3", 1, f"exact value {exact:.3f}")]
    for i, (color, dash, width, label) in enumerate(entries):
        body += legend_line(330, 48 + 22 * i, color, label, dash, width)
    lo, hi = wilson(running[-1], MAIN_TRIALS)
    return plot.svg("Running estimate of the word error rate of R3 on a binary symmetric channel with p = 0.1, "
                    "after each number of trials from 100 to 10000 on a log scale, with its 95 percent Wilson "
                    f"interval shaded around it and the exact value {exact:.3f} as a dashed line. After 100 trials "
                    f"the estimate is {running[99] / 100:.2f} and the interval runs from "
                    f"{wilson(running[99], 100)[0]:.4f} to {wilson(running[99], 100)[1]:.4f}. After 10000 trials "
                    f"the estimate is {final:.4f} and the interval has narrowed to {lo:.4f} to {hi:.4f}, which "
                    "contains the exact value.", body)


def mc_intervals() -> str:
    """Seeded estimates of the R3 word error rate at five p, with Wilson intervals, over the exact curve."""
    rows = mc_points()
    plot = Plot(24, 324, 400)
    x = LogAxis(0.01, 0.5, X0, X1)
    y = LogAxis(1e-4, 1, plot.y1, plot.y0)
    body = plot.frame("bit error probability p, log scale", "word error rate of R3, log scale")
    body += plot.y_ticks(y, [1e-4, 1e-3, 1e-2, 1e-1, 1], ["0.0001", "0.001", "0.01", "0.1", "1"])
    body += plot.x_ticks(x, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5], ["0.01", "0.02", "0.05", "0.1", "0.2", "0.5"])
    grid = log_grid(0.01, 0.5, 200)
    body.append(polyline([(x(p), y(word_error(3, 1, p))) for p in grid], "currentColor", 1.5, "6 3"))
    for p, _, k, n, lo, hi in rows:
        px = x(p)
        body.append(line(px, y(lo), px, y(hi), GREEN, 1.5))
        body.append(line(px - 4, y(lo), px + 4, y(lo), GREEN, 1.5))
        body.append(line(px - 4, y(hi), px + 4, y(hi), GREEN, 1.5))
        body.append(dot(px, y(k / n), GREEN, 4))
        body.append(text(px + 8, y(k / n) + 14, f"K = {k}", mono=True, size=11))
    body += legend_line(100, 48, "currentColor", "exact word error rate", "6 3", 1.5)
    body.append(dot(112, 66, GREEN, 4))
    body.append(text(132, 70, "estimate from 10000 trials", size=12))
    body.append(line(112, 82, 112, 98, GREEN, 1.5))
    body.append(text(132, 94, "95% Wilson interval", size=12))
    summary = ", ".join(f"at p = {p} the estimate is {k}/{n} against the exact {ex:.4g}"
                        for p, ex, k, n, _, _ in rows)
    return plot.svg("Monte Carlo estimates of the word error rate of R3 at p = 0.02, 0.05, 0.1, 0.2 and 0.3, each "
                    "from 10000 seeded trials and drawn with its 95 percent Wilson interval, over the exact rate "
                    f"as a dashed curve, on log scales. In turn, {summary}. Every interval contains the exact "
                    "value, and the intervals are widest, relative to the rate, where the rate is smallest.", body)


CHARTS = {
    "fig-9-3-arq-throughput.svg": arq_throughput,
    "fig-9-5-coding-gain.svg": coding_gain_chart,
    "fig-9-6-running-estimate.svg": running_estimate,
    "fig-9-7-mc-intervals.svg": mc_intervals,
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
