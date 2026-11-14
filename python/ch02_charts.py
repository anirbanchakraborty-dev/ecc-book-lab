"""Write the six charts of chapter 2 of Error Correction Coding as SVG files.

    python3 python/ch02_charts.py --out figures

Every plotted value comes from the functions in ch02_probability_for_coding.py,
so a chart is corrected by fixing those functions and running this script
again, never by editing its output. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from fractions import Fraction
from math import exp, log, log10, sqrt

from ch02_probability_for_coding import (binom_pmf, binom_tail, chebyshev, chernoff, markov, phi, Q,
                                         q_half, q_lower, q_upper, running_average)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
MINUS = "−"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


class LogAxis:
    def __init__(self, lo_exp: int, hi_exp: int, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo_exp, hi_exp, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (log10(v) - self.lo) / (self.hi - self.lo)


class LinAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo, hi, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (v - self.lo) / (self.hi - self.lo)


def text(x: float, y: float, s: str, anchor: str = "start", fill: str = "", mono: bool = False,
         size: int = 0, weight: str = "") -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    if weight:
        attrs.append(f'font-weight="{weight}"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def pow10(x: float, y: float, e: int, anchor: str) -> str:
    """A tick label 10^e written as real text, with the exponent raised by a tspan."""
    if e == 0:
        return text(x, y, "1", anchor, mono=True, size=11)
    exp_ = f"{MINUS}{-e}" if e < 0 else str(e)
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" text-anchor="{anchor}" font-family="{MONO}" font-size="11">'
            f'10<tspan dy="-5" font-size="8">{exp_}</tspan></text>')


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def lin_x_ticks(ax: LinAxis, values: list, labels: list, grid: bool = False) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        if grid:
            out.append(line(px, Y0, px, Y1, TEXT3, 1, opacity="0.35"))
        out.append(text(px, Y1 + 20, lab, "middle", mono=True, size=11))
    return out


def lin_y_ticks(ax: LinAxis, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
        out.append(text(X0 - 8, py + 4, lab, "end", mono=True, size=11))
    return out


def log_y_ticks(ax: LogAxis, labelled: range, every: range) -> list:
    out = []
    for e in every:
        y = ax(10.0**e)
        out.append(line(X0 - 4, y, X0, y))
        if e in labelled:
            out.append(line(X0, y, X1, y, TEXT3, 1, opacity="0.35"))
            out.append(pow10(X0 - 8, y + 4, e, "end"))
    return out


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def binomial_pmf_chart(n: int = 20, p: Fraction = Fraction(1, 10)) -> str:
    """Bars of P(X = k) for the number of flips among n bits, with the mode and the tail marked."""
    x = LinAxis(-0.5, 10.5, X0, X1)
    y = LinAxis(0, 0.30, Y1, Y0)
    body = frame(f"number of flipped bits k among n = {n} bits, each flipped with p = {float(p)}",
                 "probability P(X = k)")
    body += lin_y_ticks(y, [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
                        ["0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"])
    body += lin_x_ticks(x, list(range(11)), [str(k) for k in range(11)])
    for k in range(11):
        v = float(binom_pmf(n, k, p))
        top, half = y(v), 14
        if k == 2:
            fill, stroke, width = BG2, GREEN, 2
        elif k >= 5:
            fill, stroke, width = BG2, ORANGE, 1.5
        else:
            fill, stroke, width = BG2, "currentColor", 1
        if v * 300 / 0.30 >= 1:
            body.append(f'<rect x="{fmt(x(k) - half)}" y="{fmt(top)}" width="{2 * half}" height="{fmt(Y1 - top)}"'
                        f' fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')
        if k <= 6:
            body.append(text(x(k), top - 6, f"{v:.3f}", "middle", mono=True, size=10))
    body.append(text(x(2) + 24, y(0.262), "mode and mean: k = np = 2", fill=GREEN, weight="600"))
    tail = float(binom_tail(n, 5, p))
    body.append(text(x(5) - 14, y(0.135), "orange bars: at least 5 flips,", fill=ORANGE))
    body.append(text(x(5) - 14, y(0.135) + 16, f"total probability {tail:.4f}", fill=ORANGE))
    return svg(f"Bar chart of the binomial probabilities P(X = k) for n = {n} and p = {float(p)}. The bars rise to "
               f"a peak of 0.285 at k = 2, the mean np, and fall away quickly, so that five or more flips have "
               f"total probability {tail:.4f}.", body)


def dominating_functions(a: float = 2.0, mu: float = 1.0, s: float = 1.0) -> str:
    """The step 1{x >= a} with the three functions that lie above it, one per tail bound."""
    x = LinAxis(0, 4, X0, X1)
    y = LinAxis(-0.1, 2.5, Y1, Y0)
    body = frame("value x of a nonnegative random variable  (a = 2, mean 1)", "height of each function")
    body += lin_y_ticks(y, [0, 0.5, 1, 1.5, 2, 2.5], ["0", "0.5", "1", "1.5", "2", "2.5"])
    body += lin_x_ticks(x, [0, 1, 2, 3, 4], ["0", "1", "2", "3", "4"])
    xs = [4 * i / 400 for i in range(401)]
    ymax = 2.5

    def clipped(f):
        return [(x(v), y(f(v))) for v in xs if f(v) <= ymax]

    body.append(polyline(clipped(lambda v: v / a), CYAN, 1.5))
    body.append(polyline(clipped(lambda v: ((v - mu) / (a - mu)) ** 2), ORANGE, 1.5))
    body.append(polyline(clipped(lambda v: exp(s * (v - a))), VIOLET, 1.5))
    body.append(polyline([(x(0), y(0)), (x(a), y(0))], GREEN, 2))
    body.append(polyline([(x(a), y(1)), (x(4), y(1))], GREEN, 2))
    body.append(line(x(a), y(0), x(a), y(1), GREEN, 1, "3 3"))
    body.append(f'<circle cx="{fmt(x(a))}" cy="{fmt(y(1))}" r="3" fill="{GREEN}"/>')
    body.append(text(x(2.1), y(1) + 20, "step: 1 when x is at least a", fill=GREEN, weight="600"))
    x_cheb = mu + (a - mu) * sqrt(ymax)
    x_cher = a + log(ymax) / s
    body.append(text(x(x_cheb) - 8, Y0 + 14, "Chebyshev", "end", fill=ORANGE, weight="600"))
    body.append(text(x(x_cher) + 8, Y0 + 14, "Chernoff", fill=VIOLET, weight="600"))
    body.append(text(x(3.5), y(3.5 / a) - 12, "Markov", "end", fill=CYAN, weight="600"))
    return svg("The step function that is 0 below a = 2 and 1 from a on, drawn with the three functions that lie "
               "above it everywhere, namely the line x over a used by Markov, the parabola used by Chebyshev, and the "
               "exponential used by Chernoff.", body)


def running_average_chart(p: float = 0.1, n: int = 10000, seeds: tuple = (1, 2, 3)) -> str:
    """The fraction of flipped bits after k bits, for three simulated runs, inside a two-sd band."""
    x = LogAxis(1, 4, X0, X1)
    y = LinAxis(0, 0.4, Y1, Y0)
    body = frame("number of bits sent k (log scale)", "fraction of the k bits that flipped")
    body += lin_y_ticks(y, [0, 0.1, 0.2, 0.3, 0.4], ["0", "0.1", "0.2", "0.3", "0.4"])
    for e in range(1, 5):
        px = x(10.0**e)
        body.append(line(px, Y1, px, Y1 + 4))
        body.append(pow10(px, Y1 + 20, e, "middle"))
    ks = sorted({round(10 ** (1 + 3 * i / 300)) for i in range(301)})
    upper = [(x(k), y(min(0.4, p + 2 * sqrt(p * (1 - p) / k)))) for k in ks]
    lower = [(x(k), y(max(0.0, p - 2 * sqrt(p * (1 - p) / k)))) for k in ks]
    body.append(polyline(upper, TEXT3, 1.5, "3 3"))
    body.append(polyline(lower, TEXT3, 1.5, "3 3"))
    for seed, color in zip(seeds, (CYAN, ORANGE, VIOLET)):
        run = running_average(seed, n, p)
        body.append(polyline([(x(k), y(min(0.4, run[k - 1]))) for k in ks], color, 1.2))
    body.append(line(X0, y(p), X1, y(p), GREEN, 1.5))
    body.append(text(X0 + 6, y(p) + 18, "p = 0.1", fill=GREEN, weight="600"))
    k_lab = 1500
    body.append(text(x(k_lab), y(p + 2 * sqrt(p * (1 - p) / k_lab)) - 14, "p plus two standard deviations",
                     "middle", fill=TEXT3))
    body.append(text(x(k_lab), y(p - 2 * sqrt(p * (1 - p) / k_lab)) + 26, "p minus two standard deviations",
                     "middle", fill=TEXT3))
    body.append(text(X1 - 4, Y0 + 14, "three simulated runs", "end", fill=TEXT3))
    return svg(f"The fraction of flipped bits among the first k bits, for three simulated runs of a channel with "
               f"p = {p}, for k from 10 to {n}. The runs wander early and settle toward 0.1, inside a band that "
               f"narrows as one over the square root of k.", body)


def gaussian_error_region(sigma: float = 0.5) -> str:
    """Densities of the received voltage for +1 and -1 sent, with the error area for +1 shaded."""
    x = LinAxis(-3, 3, X0, X1)
    y = LinAxis(0, 1.2, Y1, Y0)
    body = frame("received voltage r", "probability density")
    body += lin_y_ticks(y, [0, 0.4, 0.8, 1.2], ["0", "0.4", "0.8", "1.2"])
    body += lin_x_ticks(x, [-3, -2, -1, 0, 1, 2, 3], [f"{MINUS}3", f"{MINUS}2", f"{MINUS}1", "0", "1", "2", "3"])

    def dens(v: float, centre: float) -> float:
        return phi((v - centre) / sigma) / sigma

    xs = [-3 + 6 * i / 600 for i in range(601)]
    shade = [(x(v), y(dens(v, 1))) for v in xs if v <= 0]
    pts = " ".join(f"{fmt(a)},{fmt(b)}" for a, b in shade)
    body.append(f'<polygon points="{fmt(x(-3))},{Y1} {pts} {fmt(x(0))},{Y1}" fill="{RED}" fill-opacity="0.3"'
                f' stroke="none"/>')
    body.append(polyline([(x(v), y(dens(v, -1))) for v in xs], TEXT3, 1.5, "3 3"))
    body.append(polyline([(x(v), y(dens(v, 1))) for v in xs], GREEN, 2))
    body.append(line(x(0), Y0 + 26, x(0), Y1, "currentColor", 1, "3 3"))
    peak = dens(1, 1)
    body.append(text(x(1), y(peak) - 8, "sent +1", "middle", fill=GREEN, weight="600"))
    body.append(text(x(-1), y(peak) - 8, f"sent {MINUS}1", "middle", fill=TEXT3, weight="600"))
    body.append(text(x(0), Y0 + 16, "threshold 0", "middle"))
    body.append(text(x(0) - 8, Y0 + 44, "decide 0", "end"))
    body.append(text(x(0) + 8, Y0 + 44, "decide 1", "start"))
    body.append(text(X0 + 6, Y0 + 16, "red area: +1 sent,", fill=RED))
    body.append(text(X0 + 6, Y0 + 32, "received below 0,", fill=RED))
    body.append(text(X0 + 6, Y0 + 48, f"probability Q(2) = {Q(1 / sigma):.4f}", fill=RED))
    return svg(f"Two bell-shaped densities of the received voltage, centered on +1 and on minus 1, each with noise "
               f"standard deviation {sigma}. The part of the +1 curve below the threshold at 0 is shaded, and its "
               f"area, Q(2) = {Q(1 / sigma):.4f}, is the probability of a decision error.", body)


def q_bounds_chart() -> str:
    """Q(x) on a log scale with its two-sided bounds and the simple bound exp(-x^2/2)/2."""
    x = LinAxis(0, 6, X0, X1)
    y = LogAxis(-10, 0, Y1, Y0)
    body = frame("x", "tail probability")
    body += log_y_ticks(y, range(-10, 1, 2), range(-10, 1))
    body += lin_x_ticks(x, list(range(7)), [str(k) for k in range(7)], grid=True)
    xs = [6 * i / 600 for i in range(1, 601)]
    body.append(polyline([(x(v), y(q_half(v))) for v in xs], TEXT3, 1.5, "3 3"))
    body.append(polyline([(x(v), y(q_lower(v))) for v in xs if v >= 0.02], ORANGE, 1.5))
    body.append(polyline([(x(v), y(q_upper(v))) for v in xs if q_upper(v) <= 1], CYAN, 1.5))
    body.append(polyline([(x(v), y(Q(v))) for v in xs], GREEN, 2))
    legend = [(GREEN, 2, "", "Q(x), the exact tail"),
              (CYAN, 1.5, "", "upper bound phi(x) / x"),
              (ORANGE, 1.5, "", "lower bound x phi(x) / (1 + x squared)"),
              (TEXT3, 1.5, "3 3", "bound one half of exp(" + MINUS + "x squared / 2)")]
    body.append(f'<rect x="96" y="194" width="300" height="104" rx="3" fill="{BG2}" stroke="none"/>')
    for i, (color, width, dash, label) in enumerate(legend):
        ly = 214 + 22 * i
        body.append(line(104, ly - 4, 132, ly - 4, color, width, dash))
        body.append(text(140, ly, label, fill=color if color != TEXT3 else TEXT3))
    return svg("Q(x) for x from 0 to 6 on a logarithmic scale, falling from one half to about one in a billion, "
               "squeezed between an upper and a lower bound that meet as x grows, with the looser bound one half "
               "of exp of minus x squared over 2 above all three.", body)


def tail_bounds_chart(n: int = 1000, p: Fraction = Fraction(1, 100)) -> str:
    """Markov, Chebyshev and Chernoff bounds on P(X >= k) against the exact binomial tail."""
    mean, var = n * float(p), n * float(p) * (1 - float(p))
    x = LinAxis(10, 40, X0, X1)
    y = LogAxis(-14, 0, Y1, Y0)
    body = frame(f"threshold k: at least k errors among n = {n} bits at p = {float(p)}", "probability")
    body += log_y_ticks(y, range(-14, 1, 2), range(-14, 1))
    body += lin_x_ticks(x, [10, 15, 20, 25, 30, 35, 40], ["10", "15", "20", "25", "30", "35", "40"], grid=True)
    ks = list(range(11, 41))
    curves = {
        "Markov": ([min(1.0, markov(mean, k)) for k in ks], TEXT3, 1.5, "3 3"),
        "Chebyshev": ([min(1.0, chebyshev(var, k - mean)) for k in ks], ORANGE, 1.5, ""),
        "Chernoff": ([chernoff(n, k / n, float(p)) for k in ks], CYAN, 1.5, ""),
        "exact": ([float(binom_tail(n, k, p)) for k in ks], GREEN, 2, ""),
    }
    for name, (vals, color, width, dash) in curves.items():
        body.append(polyline([(x(k), y(v)) for k, v in zip(ks, vals)], color, width, dash))
    # Each label sits on the empty side of its curve: above Markov and Chebyshev, above and to the
    # right of the falling Chernoff curve, and below and to the left of the falling exact tail.
    place = {"Markov": (38, 0, -8, "middle", "Markov"), "Chebyshev": (38, 0, -8, "middle", "Chebyshev"),
             "Chernoff": (30, 6, -6, "start", "Chernoff"), "exact": (30, -6, 16, "end", "exact tail")}
    for name, (k, dx, dy, anchor, label) in place.items():
        vals, color = curves[name][0], curves[name][1]
        body.append(text(x(k) + dx, y(vals[ks.index(k)]) + dy, label, anchor, fill=color, weight="600"))
    body.append(line(x(20), Y0, x(20), Y1, "currentColor", 1, "3 3"))
    body.append(text(x(20) - 6, Y1 - 8, "k = 20", "end", mono=True, size=11))
    return svg(f"Bounds on the probability of at least k errors among {n} bits with p = {float(p)}, for k from 11 "
               "to 40 on a logarithmic scale. Markov and Chebyshev stay near one in ten or above, while the "
               "Chernoff bound falls steeply and runs parallel to the exact tail a little above it.", body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    charts = (("fig-2-2-binomial-pmf.svg", binomial_pmf_chart()),
              ("fig-2-4-dominating-functions.svg", dominating_functions()),
              ("fig-2-5-running-average.svg", running_average_chart()),
              ("fig-2-6-gaussian-error-region.svg", gaussian_error_region()),
              ("fig-2-7-q-bounds.svg", q_bounds_chart()),
              ("fig-2-9-tail-bounds.svg", tail_bounds_chart()))
    for name, content in charts:
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
