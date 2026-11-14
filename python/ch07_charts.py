"""Write the six charts of chapter 7 of Error Correction Coding as SVG files.

    python3 python/ch07_charts.py --out figures

Every plotted value comes from the functions in ch07_capacity_and_shannons_theorem.py,
so a chart is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import os
from math import log, log10

from ch07_capacity_and_shannons_theorem import (binary_entropy, bpsk_hard_capacity, bpsk_soft_capacity, bsc_capacity,
                                                db, ebn0_boundary, gaussian_capacity, normal_approximation,
                                                q_inverse, sphere_decoder_bound, typical_probability_float)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
MINUS = "−"


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


def square(x: float, y: float, fill: str, half: float = 4) -> str:
    return (f'<rect x="{fmt(x - half)}" y="{fmt(y - half)}" width="{fmt(2 * half)}" height="{fmt(2 * half)}" '
            f'fill="{fill}" stroke="{fill}" stroke-width="1"/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def x_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        out.append(text(px, Y1 + 20, lab, "middle", mono=True, size=11))
    return out


def y_ticks(ax, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        py = ax(v)
        out.append(line(X0 - 4, py, X0, py))
        out.append(line(X0, py, X1, py, TEXT3, 1, opacity="0.35"))
        out.append(text(X0 - 8, py + 4, lab, "end", mono=True, size=11))
    return out


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "") -> list:
    return [line(x, y - 4, x + 24, y - 4, color, 2, dash), text(x + 32, y, label, size=12)]


def legend_mark(x: float, y: float, color: str, label: str, shape: str = "dot") -> list:
    mark = dot(x + 12, y - 4, color, 4) if shape == "dot" else square(x + 12, y - 4, color)
    return [mark, text(x + 32, y, label, size=12)]


def neg(label: str) -> str:
    """Write a leading hyphen as the typographic minus sign."""
    return MINUS + label[1:] if label.startswith("-") else label


# The charts --------------------------------------------------------------------------

def binary_capacity() -> str:
    """Capacity of BSC(p) and BEC(epsilon) against the channel parameter."""
    x = LinAxis(0, 1, X0, X1)
    y = LinAxis(0, 1, Y1, Y0)
    body = frame("crossover probability p, or erasure probability epsilon", "capacity in bits per channel use")
    body += y_ticks(y, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body += x_ticks(x, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body.append(line(X0, y(0.5), X1, y(0.5), TEXT3, 1.5, "3 3"))
    grid = [i / 400 for i in range(401)]
    body.append(polyline([(x(p), y(bsc_capacity(p))) for p in grid], GREEN, 2))
    body.append(polyline([(x(0), y(1)), (x(1), y(0))], CYAN, 2))
    body.append(dot(x(0.5), y(0.5), CYAN, 4))
    body.append(dot(x(0.110028), y(0.5), GREEN, 4))
    body.append(dot(x(0.889972), y(0.5), GREEN, 4))
    c = bsc_capacity(0.1)
    body.append(square(x(0.1), y(c), GREEN, 3.5))
    body.append(text(x(0.1) + 10, y(c) - 14, f"BSC(0.1): {c:.4f}", mono=True, size=11))
    body.append(text(420, y(0.5) - 6, "0.5 bit per use", size=12, fill=TEXT3))
    body += legend_line(290, 50, GREEN, f"binary symmetric, 1 {MINUS} h(p)")
    body += legend_line(290, 72, CYAN, f"binary erasure, 1 {MINUS} epsilon")
    return svg("Capacity in bits per channel use against the channel parameter from 0 to 1. The binary symmetric "
               "channel gives 1 minus h(p), a U-shaped curve that is 1 at both ends and 0 at p = 0.5. The binary "
               f"erasure channel gives the straight line 1 minus epsilon. A dashed line at 0.5 crosses the "
               f"symmetric curve at p = 0.110 and p = 0.890 and the erasure line at epsilon = 0.5. A square marks "
               f"the binary symmetric channel with p = 0.1 at {c:.4f}.", body)


def gaussian_vs_bpsk() -> str:
    """Capacity per real dimension of the Gaussian channel with Gaussian, soft BPSK and hard BPSK inputs."""
    x = LinAxis(-10, 20, X0, X1)
    y = LinAxis(0, 3.5, Y1, Y0)
    body = frame("signal-to-noise ratio P/N in dB", "bits per real dimension")
    body += y_ticks(y, [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5], ["0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5"])
    body += x_ticks(x, [-10, -5, 0, 5, 10, 15, 20], [neg(str(v)) for v in (-10, -5, 0, 5, 10, 15, 20)])
    grid = [-10 + i * 0.5 for i in range(61)]
    snr = [10 ** (d / 10) for d in grid]
    body.append(polyline([(x(d), y(gaussian_capacity(s))) for d, s in zip(grid, snr)], GREEN, 2))
    body.append(polyline([(x(d), y(bpsk_soft_capacity(s))) for d, s in zip(grid, snr)], CYAN, 2))
    body.append(polyline([(x(d), y(bpsk_hard_capacity(s))) for d, s in zip(grid, snr)], ORANGE, 2, "3 3"))
    body += legend_line(110, 50, GREEN, "Gaussian input, half log2(1 + P/N)")
    body += legend_line(110, 72, CYAN, "BPSK input, soft decisions")
    body += legend_line(110, 94, ORANGE, "BPSK input, hard decisions", "3 3")
    g20 = gaussian_capacity(100)
    return svg("Bits per real dimension against signal-to-noise ratio from minus 10 to 20 dB. The Gaussian-input "
               f"capacity half log2(1 + P/N) rises without limit, reaching 0.5 at 0 dB and {g20:.2f} at 20 dB. BPSK "
               "with soft decisions follows it closely below about 0 dB and then levels off at 1 bit. BPSK with "
               "hard decisions lies below both and also levels off at 1 bit.", body)


def bandwidth_power() -> str:
    """Spectral efficiency against the smallest Eb/N0, with the Shannon limit and uncoded BPSK."""
    x = LinAxis(-4, 16, X0, X1)
    y = LogAxis(0.01, 10, Y1, Y0)
    body = frame("Eb/N0 in dB", "spectral efficiency eta in bit/s/Hz (log scale)")
    body += y_ticks(y, [0.01, 0.1, 1, 10], ["0.01", "0.1", "1", "10"])
    body += x_ticks(x, [-4, 0, 4, 8, 12, 16], [neg(str(v)) for v in (-4, 0, 4, 8, 12, 16)])
    limit = db(log(2))
    body.append(line(x(limit), Y0, x(limit), Y1, TEXT3, 1.5, "3 3"))
    etas = [10 ** (-2 + i * 0.01) for i in range(301)]
    pts = [(x(db(ebn0_boundary(eta))), y(eta)) for eta in etas if db(ebn0_boundary(eta)) <= 16]
    body.append(polyline(pts, GREEN, 2))
    bpsk = db(q_inverse(1e-5) ** 2 / 2)
    bound2 = db(ebn0_boundary(2))
    body.append(line(x(bound2), y(2), x(bpsk), y(2), TEXT3, 1.5, "3 3"))
    body.append(square(x(bpsk), y(2), ORANGE, 4))
    body.append(text((x(bound2) + x(bpsk)) / 2, y(2) - 8, f"{bpsk - bound2:.2f} dB", "middle", mono=True, size=11))
    body.append(text(x(bpsk) - 8, y(2) + 22, "uncoded BPSK", "end", size=12))
    body.append(text(x(bpsk) - 8, y(2) + 36, "bit error 1e-5", "end", size=12))
    body.append(text(x(limit) + 6, 300, f"Shannon limit {MINUS}1.59 dB", size=12))
    body.append(text(x(limit) + 14, 56, "impossible", size=12, fill=TEXT3))
    body.append(text(470, 250, "achievable", size=12, fill=TEXT3))
    return svg("Spectral efficiency eta in bits per second per hertz on a log scale from 0.01 to 10 against Eb/N0 "
               "in dB. The boundary curve (2 to the eta minus 1) over eta separates the impossible region on the "
               f"left from the achievable region on the right. It passes 1.76 dB at eta = 2 and 0 dB at eta = 1, "
               f"and approaches the dashed vertical Shannon limit at minus 1.59 dB as eta falls. A square marks "
               f"uncoded BPSK at bit error 1e-5, {bpsk:.2f} dB at eta = 2, {bpsk - bound2:.2f} dB to the right "
               "of the boundary.", body)


def typical_set() -> str:
    """Probability of the typical set against block length for two windows."""
    x = LinAxis(0, 1000, X0, X1)
    y = LinAxis(0, 1, Y1, Y0)
    body = frame("block length n", "probability of the typical set")
    body += y_ticks(y, [0, 0.2, 0.4, 0.6, 0.8, 1], ["0", "0.2", "0.4", "0.6", "0.8", "1"])
    body += x_ticks(x, [0, 200, 400, 600, 800, 1000], ["0", "200", "400", "600", "800", "1000"])
    ns = list(range(10, 1001, 10))
    wide = [typical_probability_float(n, 0.1, 0.05) for n in ns]
    narrow = [typical_probability_float(n, 0.1, 0.02) for n in ns]
    body.append(polyline([(x(n), y(v)) for n, v in zip(ns, wide)], GREEN, 2))
    body.append(polyline([(x(n), y(v)) for n, v in zip(ns, narrow)], CYAN, 2))
    body += legend_line(330, 250, GREEN, "ones fraction within 0.05 of p")
    body += legend_line(330, 272, CYAN, "ones fraction within 0.02 of p")
    return svg("Probability that a block from a binary source with p = 0.1 is typical, against block length n from "
               "10 to 1000. With the window of 0.05 the probability starts at 0.39 at n = 10, passes "
               f"{typical_probability_float(100, 0.1, 0.05):.2f} at n = 100 and is above 0.999 by n = 400. With "
               f"the narrower window of 0.02 it rises more slowly, to {narrow[-1]:.2f} at n = 1000.", body)


def random_coding() -> str:
    """Union bound on the average error of random codes over BSC(0.1), radius optimized, against n."""
    x = LinAxis(0, 4000, X0, X1)
    y = LogAxis(1e-16, 1, Y1, Y0)
    body = frame("block length n", "average block error, upper bound (log scale)")
    body += y_ticks(y, [1e-16, 1e-12, 1e-8, 1e-4, 1], ["1e-16", "1e-12", "1e-8", "1e-4", "1"])
    body += x_ticks(x, [0, 1000, 2000, 3000, 4000], ["0", "1000", "2000", "3000", "4000"])
    ns = list(range(100, 4001, 50))
    styles = ((0.3, GREEN, ""), (0.4, CYAN, ""), (0.45, ORANGE, "3 3"))
    for rate, color, dash in styles:
        pts = []
        for n in ns:
            v = sphere_decoder_bound(n, rate, 0.1)
            if v < 1e-16:
                break
            pts.append((x(n), y(v)))
        body.append(polyline(pts, color, 2, dash))
    walked = sphere_decoder_bound(1000, 0.4, 0.1)
    body.append(dot(x(1000), y(walked), CYAN, 4))
    body.append(text(x(1000) + 10, y(walked) - 8, f"n = 1000: {walked:.1e}".replace("e-0", "e-"), mono=True, size=11))
    body += legend_line(380, 46, GREEN, "rate 0.30")
    body += legend_line(380, 66, CYAN, "rate 0.40")
    body += legend_line(380, 86, ORANGE, "rate 0.45", "3 3")
    body.append(text(290, 290, f"capacity of BSC(0.1): {bsc_capacity(0.1):.4f}", size=12, fill=TEXT3))
    return svg("Upper bound on the average block error of a random code over the binary symmetric channel with "
               "p = 0.1, against block length n from 100 to 4000, on a log scale from 1e-16 to 1. Each curve uses "
               "the best decoding radius. At rate 0.30 the bound falls below 1e-16 near n = 1400, at rate 0.40 it "
               f"is {walked:.1e} at n = 1000, and at rate 0.45, closer to the capacity 0.5310, it falls more slowly, "
               f"to {sphere_decoder_bound(4000, 0.45, 0.1):.1e} at n = 4000.".replace("e-0", "e-"), body)


def finite_length() -> str:
    """Normal approximation to the best rate at block length n for BSC(0.1)."""
    x = LogAxis(100, 10000, X0, X1)
    y = LinAxis(0, 0.6, Y1, Y0)
    body = frame("block length n (log scale)", "best rate in bits per channel use")
    body += y_ticks(y, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6"])
    body += x_ticks(x, [100, 1000, 10000], ["100", "1000", "10000"])
    c = bsc_capacity(0.1)
    body.append(line(X0, y(c), X1, y(c), TEXT3, 1.5, "3 3"))
    ns = [round(10 ** (2 + i / 50)) for i in range(101)]
    body.append(polyline([(x(n), y(normal_approximation(n, 1e-3, 0.1))) for n in ns], GREEN, 2))
    body.append(polyline([(x(n), y(normal_approximation(n, 1e-6, 0.1))) for n in ns], CYAN, 2))
    body += legend_line(330, 250, GREEN, "block error 1e-3")
    body += legend_line(330, 272, CYAN, "block error 1e-6")
    body += legend_line(330, 294, TEXT3, f"capacity {c:.4f}", "3 3")
    return svg("Normal approximation to the best rate of a code for the binary symmetric channel with p = 0.1, "
               "against block length n from 100 to 10000 on a log scale. For block error 1e-3 the rate is "
               f"{normal_approximation(100, 1e-3, 0.1):.4f} at n = 100, {normal_approximation(1000, 1e-3, 0.1):.4f} "
               f"at n = 1000 and {normal_approximation(10000, 1e-3, 0.1):.4f} at n = 10000. For block error 1e-6 "
               f"it is {normal_approximation(100, 1e-6, 0.1):.4f}, {normal_approximation(1000, 1e-6, 0.1):.4f} and "
               f"{normal_approximation(10000, 1e-6, 0.1):.4f}. Both approach the dashed capacity line at {c:.4f}.",
               body)


CHARTS = {
    "fig-7-3-binary-capacity.svg": binary_capacity,
    "fig-7-5-gaussian-capacity.svg": gaussian_vs_bpsk,
    "fig-7-6-bandwidth-power.svg": bandwidth_power,
    "fig-7-7-typical-set.svg": typical_set,
    "fig-7-8-random-coding.svg": random_coding,
    "fig-7-9-finite-length.svg": finite_length,
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
