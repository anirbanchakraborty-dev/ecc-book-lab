"""Write the four computed figures of chapter 10 of Error Correction Coding as SVG files.

    python3 python/ch10_charts.py --out figures

Every plotted value comes from the functions in
ch10_case_study_codes_in_the_wild.py, so a figure is corrected by fixing those
functions and running this script again, never by editing its output. The
Voyager figure plots the data rates that the companion script checks against
the table it was taken from, and the compact disc figure plots the frame-by-frame
delay model that the script checks against its closed form. The SVG follows the
book's figure style sheet: real text in IBM Plex, colors from the page's tokens
with paper fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import log10

from ch10_case_study_codes_in_the_wild import (VOYAGER_TABLE, cd_max_hits, hamming74_block, repetition3_block,
                                               repetition5_word_error, rm15_word_error, uncoded_block)

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
        head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {self.h}" width="{W}" height="{self.h}">'
        return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "", width: float = 2) -> list:
    return [line(x, y - 4, x + 24, y - 4, color, width, dash), text(x + 32, y, label, size=12)]


def log_grid(lo: float, hi: float, count: int) -> list:
    a, b = log10(lo), log10(hi)
    return [10 ** (a + (b - a) * i / count) for i in range(count + 1)]


def power_label(e: int) -> str:
    return "1" if e == 0 else f"10^{e}" if e > 0 else f"1e-{-e}"


# The charts --------------------------------------------------------------------------

def hamming_block_error() -> str:
    plot = Plot(30, 300, 350)
    xa = LogAxis(1e-4, 1e-1, X0, X1)
    ya = LogAxis(1e-8, 1, plot.y1, plot.y0)
    ps = log_grid(1e-4, 1e-1, 120)
    body = plot.y_ticks(ya, [10.0 ** e for e in range(-8, 1, 2)], [power_label(e) for e in range(-8, 1, 2)])
    body += plot.x_ticks(xa, [1e-4, 1e-3, 1e-2, 1e-1], ["0.0001", "0.001", "0.01", "0.1"])
    body += plot.frame("bit error probability p of the binary symmetric channel",
                       "probability the 4 data bits arrive wrong")
    curves = [
        (lambda p: uncoded_block(p, 4), TEXT3, "3 3", "4 bits sent bare"),
        (hamming74_block, GREEN, "", "Hamming (7,4), 7 bits"),
        (repetition3_block, CYAN, "", "each bit sent three times, 12 bits"),
    ]
    for f, color, dash, _ in curves:
        body.append(polyline([(xa(p), ya(f(p))) for p in ps], color, 2, dash))
    for f, color, _, _ in curves:
        body.append(dot(xa(0.01), ya(f(0.01)), color))
    for i, (_, color, dash, label) in enumerate(curves):
        body += legend_line(330, 210 + 22 * i, color, label, dash)
    title = (f"Probability that four data bits are delivered wrong over a binary symmetric channel, for p from "
             f"0.0001 to 0.1 on logarithmic axes. Sent bare the block fails with probability "
             f"{uncoded_block(0.01, 4):.4f} at p = 0.01. The Hamming (7,4) code, which corrects any single error "
             f"in 7 bits, fails with {hamming74_block(0.01):.5f}, and three-fold repetition, which spends 12 bits, "
             f"fails with {repetition3_block(0.01):.5f}. Both coded curves fall with slope 2, since each needs two "
             f"errors in the wrong place.")
    return plot.svg(title, body)


def mariner_word_error() -> str:
    plot = Plot(30, 300, 350)
    xa = LogAxis(0.01, 0.2, X0, X1)
    ya = LogAxis(1e-10, 1, plot.y1, plot.y0)
    ps = log_grid(0.01, 0.2, 120)
    body = plot.y_ticks(ya, [10.0 ** e for e in range(-10, 1, 2)], [power_label(e) for e in range(-10, 1, 2)])
    body += plot.x_ticks(xa, [0.01, 0.02, 0.05, 0.1, 0.2], ["0.01", "0.02", "0.05", "0.1", "0.2"])
    body += plot.frame("bit error probability p of the binary symmetric channel",
                       "probability a 6-bit word is wrong")
    curves = [
        (lambda p: uncoded_block(p, 6), TEXT3, "3 3", "6 bits sent bare"),
        (repetition5_word_error, CYAN, "", "each bit sent five times, 30 bits"),
        (rm15_word_error, GREEN, "", "Reed–Muller (32,6), up to 7 errors fixed"),
    ]
    for f, color, dash, _ in curves:
        body.append(polyline([(xa(p), ya(f(p))) for p in ps], color, 2, dash))
    for f, color, _, _ in curves:
        body.append(dot(xa(0.05), ya(f(0.05)), color))
    for i, (_, color, dash, label) in enumerate(curves):
        body += legend_line(290, 222 + 22 * i, color, label, dash)
    title = (f"Probability that a 6-bit word is decoded wrong over a binary symmetric channel, for p from 0.01 "
             f"to 0.2 on logarithmic axes. At p = 0.05 the bare 6 bits fail with probability "
             f"{uncoded_block(0.05, 6):.4f}, five-fold repetition in 30 bits fails with "
             f"{repetition5_word_error(0.05):.5f}, and the 32-bit Reed–Muller code of the Mariner probes, "
             f"decoded to correct up to 7 errors, fails with {rm15_word_error(0.05):.6f}. Its curve is the "
             f"steepest because 8 errors are needed to defeat it.")
    return plot.svg(title, body)


def voyager_data_rates() -> str:
    plot = Plot(30, 300, 350)
    names = [row[0] for row in VOYAGER_TABLE]
    xa = LinAxis(-0.5, len(names) - 0.5, X0, X1)
    ya = LogAxis(1000, 200_000, plot.y1, plot.y0)
    ticks = [1000, 3000, 10_000, 30_000, 100_000]
    body = plot.y_ticks(ya, ticks, ["1,000", "3,000", "10,000", "30,000", "100,000"])
    body += plot.x_ticks(xa, list(range(len(names))), names)
    body += plot.frame("planet at encounter", "downlink data rate, bits per second")
    body.append(polyline([(xa(i), ya(row[2])) for i, row in enumerate(VOYAGER_TABLE)], TEXT3, 1.5, "3 3"))
    body.append(polyline([(xa(i), ya(row[3])) for i, row in enumerate(VOYAGER_TABLE)], GREEN, 2))
    for i, row in enumerate(VOYAGER_TABLE):
        body.append(ring(xa(i), ya(row[2]), TEXT3))
        body.append(dot(xa(i), ya(row[3]), GREEN))
        if i:
            body.append(text(xa(i) + 10, ya(row[3]) - 8, f"{row[3]:,}", fill=GREEN, mono=True, size=11))
            body.append(text(xa(i) - 10, ya(row[2]) + 18, f"{row[2]:,}", "end", TEXT3, mono=True, size=11))
    body.append(text(xa(0) + 10, ya(VOYAGER_TABLE[0][3]) - 8, f"{VOYAGER_TABLE[0][3]:,}", mono=True, size=11))
    body += legend_line(110, 262, TEXT3, "expected from the drop in received power", "3 3", 1.5)
    body += legend_line(110, 284, GREEN, "achieved by Voyager", "", 2)
    pairs = ", ".join(f"{row[0]} had {row[2]:,} expected and {row[3]:,} achieved" for row in VOYAGER_TABLE[1:])
    title = (f"Voyager's downlink data rates at its four planetary encounters on a logarithmic axis. Jupiter set "
             f"the reference of {VOYAGER_TABLE[0][3]:,} bits per second. The inverse-square fall of received "
             f"power with distance predicts the dashed curve, and the solid curve is the maximum rate achieved. "
             f"Saturn used the same code and antennas as Jupiter, and its gain came from the rates available and the "
             f"margins chosen. At Uranus and Neptune the Reed–Solomon code and an enlarged, arrayed ground system both "
             f"contributed. In bits per second, {pairs}.")
    return plot.svg(title, body)


def cd_burst_erasures() -> str:
    plot = Plot(30, 300, 350)
    xa = LinAxis(0, 40, X0, X1)
    ya = LinAxis(0, 10, plot.y1, plot.y0)
    body = plot.y_ticks(ya, list(range(0, 11, 2)), [str(v) for v in range(0, 11, 2)])
    body += plot.x_ticks(xa, list(range(0, 41, 8)), [str(v) for v in range(0, 41, 8)])
    body += plot.frame("length of the burst, in consecutive frames lost from the disc",
                       "bytes lost from one C2 codeword")
    pts = []
    for b in range(0, 41):
        h = cd_max_hits(b)
        if pts:
            pts.append((xa(b), pts[-1][1]))
        pts.append((xa(b), ya(h)))
    body.append(polyline(pts, GREEN, 2))
    body.append(line(X0, ya(4), X1, ya(4), ORANGE, 1.5, "3 3"))
    body.append(text(X1 - 4, ya(4) - 8, "C2 fills at most 4 erasures", "end", ORANGE, size=12))
    body.append(line(xa(16), ya(4), xa(16), plot.y1, TEXT3, 1, "3 3"))
    body.append(dot(xa(16), ya(4), GREEN))
    body.append(text(xa(16) + 8, ya(4) + 22, "16 frames, 2.18 ms", size=12))
    title = (f"Number of bytes of one C2 codeword that a burst of lost frames removes, against the burst length "
             f"from 0 to 40 frames. Because the delay lines put 4 frames between consecutive bytes of a C2 "
             f"codeword, a burst of B frames removes at most ceil(B/4) bytes, a staircase. C2 fills up to 4 erasures, so "
             f"the longest burst it survives is 16 frames, {cd_max_hits(16)} bytes, about 2.18 milliseconds of "
             f"playing time, while a 17-frame burst can remove {cd_max_hits(17)}.")
    return plot.svg(title, body)


CHARTS = {
    "fig-10-2-hamming-block-error.svg": hamming_block_error,
    "fig-10-3-mariner-word-error.svg": mariner_word_error,
    "fig-10-5-voyager-data-rates.svg": voyager_data_rates,
    "fig-10-7-cd-burst-erasures.svg": cd_burst_erasures,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="figures", help="folder to write the SVG files into")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name, make in CHARTS.items():
        path = os.path.join(args.out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(make())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
