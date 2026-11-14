"""Write the eight charts of chapter 4 of Error Correction Coding as SVG files.

    python3 python/ch04_charts.py --out figures

Every plotted value comes from the functions in ch04_information_and_entropy.py,
so a chart is corrected by fixing those functions and running this script
again, never by editing its output. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import log2

from ch04_information_and_entropy import (DISTRIBUTIONS, P_FLIP, RUNNING_Q, bec_mi, bsc_joint, bsc_mi,
                                          bsc_output_one, entropy, h, surprise, z_mi)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
MINUS = "−"
ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5" '
         'orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker>'
         '<marker id="ar" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5" '
         f'orient="auto"><path d="M0 0L10 5L0 10z" fill="{RED}"/></marker></defs>')


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


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


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "", arrow: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    m = f' marker-end="url(#{arrow})"' if arrow else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}{m}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def rect(x: float, y: float, w: float, h_: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, rx: float = 3) -> str:
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h_)}" rx="{fmt(rx)}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def lin_x_ticks(ax: LinAxis, values: list, labels: list) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
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


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def unit_ticks() -> tuple:
    return [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"]


# The charts --------------------------------------------------------------------------

def surprise_curve() -> str:
    """s(p) = log2(1/p), from p = 1/128 (seven bits) to p = 1 (no surprise), with six marked outcomes."""
    x = LinAxis(0, 1, X0, X1)
    y = LinAxis(0, 7, Y1, Y0)
    body = frame("probability p of the outcome", "surprise in bits")
    body += lin_y_ticks(y, list(range(8)), [str(k) for k in range(8)])
    body += lin_x_ticks(x, *unit_ticks())
    pts = [2.0 ** (-7 + 7 * i / 400) for i in range(401)]
    body.append(polyline([(x(p), y(surprise(p))) for p in pts], "currentColor", 2))
    marks = ((1 / 52, "card, 1/52", "start"), (0.1, "rare coin side, 0.1", "start"),
             (1 / 6, "die face, 1/6", "start"), (1 / 4, "1/4", "start"), (1 / 2, "fair coin, 1/2", "start"),
             (0.9, "0.9", "middle"))
    for p, name, anchor in marks:
        px, py = x(p), y(surprise(p))
        body.append(dot(px, py, GREEN))
        label = f"{name}: {surprise(p):.3f}".replace(": 1.000", ": 1").replace(": 2.000", ": 2")
        if anchor == "middle":
            body.append(text(px, py - 16, label, "middle", mono=True, size=11))
        else:
            body.append(text(px + 10, py - 8, label, mono=True, size=11))
    return svg("A curve of surprise in bits against probability. It falls from seven bits at probability 1 in 128 "
               "to zero at probability 1. Marked points: a card, probability 1 in 52, 5.700 bits; the rare side of a "
               "biased coin, 0.1, 3.322 bits; a die face, 2.585 bits; one of four, 2 bits; a fair coin, 1 bit; the "
               "common side of a biased coin, 0.9, 0.152 bits.", body)


def distribution_entropies() -> str:
    """Four distributions as bars on one probability scale, each with its entropy."""
    w, h_ = 640, 410
    names = {"fair coin": ["heads", "tails"], "biased coin": ["heads", "tails"],
             "fair die": [str(k) for k in range(1, 7)], "loaded die": [str(k) for k in range(1, 7)]}
    labels = {1 / 6: "1/6", 0.5: "0.5", 0.9: "0.9", 0.1: "0.1"}
    body = []
    for idx, (name, probs) in enumerate(DISTRIBUTIONS):
        left = 70 + 300 * (idx % 2)
        top = 24 + 198 * (idx // 2)
        base = top + 150
        width = 220
        hval = entropy(probs)
        hs = "1 bit" if abs(hval - 1) < 1e-12 else f"{hval:.3f} bits"
        body.append(text(left, top + 14, name, weight="600"))
        body.append(text(left + width, top + 14, f"H = {hs}", "end", fill=GREEN, mono=True, size=11))
        body.append(line(left, base, left + width, base))
        slot = width / len(probs)
        for k, p in enumerate(probs):
            bw = slot * 0.6
            bx = left + slot * k + (slot - bw) / 2
            bh = 120 * p
            body.append(rect(bx, base - bh, bw, bh, BG2, "currentColor", 1, 2))
            lab = min(labels.items(), key=lambda kv: abs(kv[0] - p))[1]
            body.append(text(bx + bw / 2, base - bh - 4, lab, "middle", mono=True, size=10))
            body.append(text(bx + bw / 2, base + 16, names[name][k], "middle", size=11, fill=TEXT3))
    return svg("Four bar charts on one probability scale. A fair coin, two bars of 0.5, has entropy 1 bit. A biased "
               "coin, bars of 0.9 and 0.1, has 0.469 bits. A fair die, six bars of 1/6, has 2.585 bits. A loaded "
               "die, one bar of 0.5 and five of 0.1, has 2.161 bits.", body, w, h_)


def binary_entropy() -> str:
    """h(p) from 0 to 1, with its peak at 1/2 and the symmetric pair at 0.1 and 0.9."""
    x = LinAxis(0, 1, X0, X1)
    y = LinAxis(0, 1.1, Y1, Y0)
    body = frame("probability p", "h(p) in bits")
    body += lin_y_ticks(y, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body += lin_x_ticks(x, *unit_ticks())
    body.append(line(x(0.5), y(1), x(0.5), Y1, TEXT3, 1, "3 3"))
    body.append(polyline([(x(i / 500), y(h(i / 500))) for i in range(501)], "currentColor", 2))
    body.append(dot(x(0.5), y(1), GREEN))
    body.append(text(x(0.5), y(1) - 12, "h(1/2) = 1", "middle", fill=GREEN, mono=True, size=11, weight="600"))
    body.append(dot(x(0.1), y(h(0.1))))
    body.append(text(x(0.1) + 10, y(h(0.1)) + 18, f"h(0.1) = {h(0.1):.3f}", mono=True, size=11))
    body.append(dot(x(0.9), y(h(0.9))))
    body.append(text(x(0.9) - 10, y(h(0.9)) + 18, f"h(0.9) = {h(0.9):.3f}", "end", mono=True, size=11))
    body.append(dot(x(0.25), y(h(0.25))))
    body.append(text(x(0.25) + 10, y(h(0.25)) + 18, f"h(0.25) = {h(0.25):.3f}", mono=True, size=11))
    return svg("The binary entropy function h of p, an arch from 0 at p equal to 0, up to 1 bit at p equal to one "
               f"half, and back to 0 at p equal to 1. It is symmetric: h(0.1) and h(0.9) are both {h(0.1):.3f}. "
               f"h(0.25) is {h(0.25):.3f}.", body)


def jensen_chord() -> str:
    """log2 x with the chord from (1, 0) to (4, 2): the chord's midpoint lies under the curve."""
    x = LinAxis(0, 5, X0, X1)
    y = LinAxis(-1.5, 2.5, Y1, Y0)
    body = frame("x", "log2 x")
    body += lin_y_ticks(y, [-1, 0, 1, 2], [f"{MINUS}1", "0", "1", "2"])
    body += lin_x_ticks(x, [0, 1, 2, 3, 4, 5], ["0", "1", "2", "3", "4", "5"])
    xs = [0.4 + 4.6 * i / 400 for i in range(401)]
    body.append(polyline([(x(v), y(log2(v))) for v in xs], "currentColor", 2))
    body.append(text(x(5) - 44, y(log2(5)) - 10, "log2 x", "end"))
    body.append(line(x(1), y(0), x(4), y(2), GREEN, 2))
    body.append(line(x(2.5), Y1, x(2.5), y(log2(2.5)), TEXT3, 1, "3 3"))
    body.append(dot(x(1), y(0)))
    body.append(dot(x(4), y(2)))
    body.append(text(x(1) + 8, y(0) + 18, "(1, 0)", mono=True, size=11))
    body.append(text(x(4) + 10, y(2) + 14, "(4, 2)", mono=True, size=11))
    body.append(dot(x(2.5), y(log2(2.5))))
    body.append(text(x(2.5) - 10, y(log2(2.5)) - 8, f"log2 E[X] = {log2(2.5):.3f}", "end", mono=True, size=11))
    body.append(dot(x(2.5), y(1), GREEN))
    body.append(text(x(2.5) + 10, y(1) + 16, "E[log2 X] = 1", fill=GREEN, mono=True, size=11, weight="600"))
    return svg("The curve log2 x with a green chord from the point (1, 0) to the point (4, 2). Above x equal to 2.5, "
               "the chord's midpoint sits at height 1, the average of log2 1 and log2 4, while the curve sits higher, "
               f"at log2 2.5, about {log2(2.5):.3f}.", body)


def joint_table() -> str:
    """The running example: the input distribution, BSC(0.1), and the joint table with its margins."""
    w, h_ = 700, 300
    q, p = RUNNING_Q, P_FLIP
    J = bsc_joint(q, p)
    body = [ARROW]
    body.append(text(150, 44, "sent X", "middle", fill=TEXT3, size=11))
    body.append(text(370, 44, "received Y", "middle", fill=TEXT3, size=11))
    for cx in (150, 370):
        for cy, bit in ((90, "0"), (230, "1")):
            body.append(f'<circle cx="{cx}" cy="{cy}" r="16" fill="{BG2}" stroke="currentColor" stroke-width="1.5"/>')
            body.append(text(cx, cy + 5, bit, "middle", mono=True, size=13))
    body.append(text(124, 94, f"P(X = 0) = {1 - q:g}", "end", mono=True, size=11))
    body.append(text(124, 234, f"P(X = 1) = {q:g}", "end", mono=True, size=11))
    body.append(line(168, 90, 350, 90, arrow="ah", width=1.5))
    body.append(line(168, 230, 350, 230, arrow="ah", width=1.5))
    body.append(line(165, 98, 351, 222, RED, 1.5, arrow="ar"))
    body.append(line(165, 222, 351, 98, RED, 1.5, arrow="ar"))
    body.append(text(260, 81, f"{1 - p:g}", "middle", mono=True, size=12))
    body.append(text(260, 250, f"{1 - p:g}", "middle", mono=True, size=12))
    body.append(text(222, 164, f"{p:g}", "middle", fill=RED, mono=True, size=12, weight="600"))
    body.append(text(294, 164, f"{p:g}", "middle", fill=RED, mono=True, size=12, weight="600"))
    # The table: a row-label column, then Y = 0, Y = 1 and the X margin.
    cols = (500, 562, 624)
    cw, rows, rh = 62, (90, 135, 180), 45
    body.append(text(593, 44, "joint probabilities", "middle", fill=TEXT3, size=11))
    for cx, lab in zip(cols, ("Y = 0", "Y = 1", "P(X)")):
        body.append(text(cx + cw / 2, 80, lab, "middle", mono=True, size=11))
    for ry, lab in zip(rows, ("X = 0", "X = 1", "P(Y)")):
        body.append(text(466, ry + rh / 2 + 4, lab, "middle", mono=True, size=11))
    py0, py1 = J[(0, 0)] + J[(1, 0)], J[(0, 1)] + J[(1, 1)]
    cells = ((J[(0, 0)], J[(0, 1)], 1 - q), (J[(1, 0)], J[(1, 1)], q), (py0, py1, 1.0))
    for r, (ry, row) in enumerate(zip(rows, cells)):
        for c, (cx, v) in enumerate(zip(cols, row)):
            joint = r < 2 and c < 2
            body.append(rect(cx, ry, cw, rh, BG2 if joint else "none", "currentColor" if joint else TEXT3, 1, 2))
            body.append(text(cx + cw / 2, ry + rh / 2 + 4, f"{v:g}", "middle", mono=True, size=12,
                             fill="" if joint else TEXT3))
    hxy = entropy(J.values())
    body.append(text(593, 254, f"H(X, Y) = {hxy:.4f} bits", "middle", fill=GREEN, mono=True, size=12, weight="600"))
    body.append(text(593, 276, f"H(X) = {h(q):.4f}   H(Y) = {h(py1):.4f}", "middle", fill=TEXT3, mono=True,
                     size=11))
    return svg(f"A binary symmetric channel with flip probability {p:g}, fed 0 with probability {1 - q:g} and 1 with "
               f"probability {q:g}, beside the joint table it produces. The four joint probabilities are "
               f"{J[(0, 0)]:g} and {J[(0, 1)]:g} in the row X = 0 and {J[(1, 0)]:g} and {J[(1, 1)]:g} in the row "
               f"X = 1. The column sums are {py0:g} and {py1:g}. The joint entropy is {hxy:.4f} bits.", body, w, h_)


def information_diagram() -> str:
    """Two overlapping circles for H(X) and H(Y), with the running example's values in each region."""
    w, h_ = 640, 350
    q, p = RUNNING_Q, P_FLIP
    J = bsc_joint(q, p)
    hx, hy, hxy = h(q), h(bsc_output_one(q, p)), entropy(J.values())
    hx_y, hy_x, mi = hxy - hy, hxy - hx, hx + hy - hxy
    body = []
    body.append(f'<circle cx="250" cy="170" r="120" fill="{BG2}" fill-opacity="0.6" stroke="{CYAN}" stroke-width="2"/>')
    body.append(f'<circle cx="390" cy="170" r="125" fill="{BG2}" fill-opacity="0.6" stroke="{ORANGE}" stroke-width="2"/>')
    body.append(text(200, 36, f"H(X) = {hx:.4f}", "middle", fill=CYAN, mono=True, size=12, weight="600"))
    body.append(text(440, 36, f"H(Y) = {hy:.4f}", "middle", fill=ORANGE, mono=True, size=12, weight="600"))
    body.append(text(185, 166, "H(X | Y)", "middle"))
    body.append(text(185, 184, f"{hx_y:.4f}", "middle", mono=True))
    body.append(text(317, 166, "I(X; Y)", "middle", fill=GREEN, weight="600"))
    body.append(text(317, 184, f"{mi:.4f}", "middle", fill=GREEN, mono=True, weight="600"))
    body.append(text(452, 166, "H(Y | X)", "middle"))
    body.append(text(452, 184, f"{hy_x:.4f}", "middle", mono=True))
    body.append(line(130, 312, 515, 312))
    body.append(line(130, 304, 130, 320))
    body.append(line(515, 304, 515, 320))
    body.append(text(322.5, 336, f"whole area: H(X, Y) = {hxy:.4f}", "middle", mono=True, size=12))
    return svg("An information diagram: two overlapping circles, H(X) on the left and H(Y) on the right. The left "
               f"crescent is H(X given Y) = {hx_y:.4f} bits, the overlap is the mutual information I(X; Y) = "
               f"{mi:.4f} bits, the right crescent is H(Y given X) = {hy_x:.4f} bits, and the whole shaded area is "
               f"the joint entropy H(X, Y) = {hxy:.4f} bits.", body, w, h_)


def bsc_mutual_information() -> str:
    """I(X; Y) against the flip probability p, for the uniform input and for q = 1/4."""
    x = LinAxis(0, 1, X0, X1)
    y = LinAxis(0, 1.1, Y1, Y0)
    body = frame("flip probability p", "I(X; Y) in bits")
    body += lin_y_ticks(y, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body += lin_x_ticks(x, *unit_ticks())
    ps = [i / 500 for i in range(501)]
    body.append(polyline([(x(p), y(bsc_mi(0.5, p))) for p in ps], GREEN, 2))
    body.append(polyline([(x(p), y(bsc_mi(RUNNING_Q, p))) for p in ps], "currentColor", 1.5, "3 3"))
    p = P_FLIP
    body.append(dot(x(p), y(bsc_mi(0.5, p)), GREEN))
    body.append(text(x(p) + 10, y(bsc_mi(0.5, p)) - 7, f"{bsc_mi(0.5, p):.3f}", fill=GREEN, mono=True, size=11,
                     weight="600"))
    body.append(dot(x(p), y(bsc_mi(RUNNING_Q, p))))
    body.append(text(x(p) - 8, y(bsc_mi(RUNNING_Q, p)) + 22, f"{bsc_mi(RUNNING_Q, p):.3f}", "end", mono=True,
                     size=11))
    body.append(dot(x(0.5), Y1, "currentColor"))
    body.append(line(x(0.5), 286, x(0.5), Y1 - 6, TEXT3, 1, "3 3"))
    body.append(text(x(0.5), 280, "p = 1/2: nothing gets through", "middle", size=11, fill=TEXT3))
    lx = 250
    body.append(line(lx, 50, lx + 20, 50, GREEN, 2))
    body.append(text(lx + 26, 54, "uniform input, q = 1/2", size=11))
    body.append(line(lx, 70, lx + 20, 70, "currentColor", 1.5, "3 3"))
    body.append(text(lx + 26, 74, "q = 1/4", size=11))
    return svg("Mutual information of the binary symmetric channel against its flip probability p. The green curve, "
               "for a uniform input, is 1 minus h(p): 1 bit at p equal to 0 and at p equal to 1, and 0 at p equal to "
               f"one half. The dashed curve, for an input that is 1 a quarter of the time, lies below it everywhere "
               f"except at p equal to one half, where both are 0. At p = {P_FLIP:g} the two give {bsc_mi(0.5, P_FLIP):.3f} and "
               f"{bsc_mi(RUNNING_Q, P_FLIP):.3f} bits.", body)


def input_distribution() -> str:
    """I(X; Y) against the input probability q for three channels, each with noise parameter 0.1."""
    x = LinAxis(0, 1, X0, X1)
    y = LinAxis(0, 1, Y1, Y0)
    p = P_FLIP
    body = frame("input probability q = P(X = 1)", "I(X; Y) in bits")
    body += lin_y_ticks(y, [0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"])
    body += lin_x_ticks(x, *unit_ticks())
    qs = [i / 500 for i in range(501)]
    body.append(polyline([(x(q), y(bec_mi(q, p))) for q in qs], TEXT3, 1.5, "3 3"))
    body.append(polyline([(x(q), y(z_mi(q, p))) for q in qs], GREEN, 2))
    body.append(polyline([(x(q), y(bsc_mi(q, p))) for q in qs], "currentColor", 2))
    best = max(range(1, 1000), key=lambda k: z_mi(k / 1000, p)) / 1000
    body.append(dot(x(0.5), y(bec_mi(0.5, p)), TEXT3))
    body.append(text(x(0.5), y(bec_mi(0.5, p)) - 10, f"{bec_mi(0.5, p):g}", "middle", fill=TEXT3, mono=True,
                     size=11))
    body.append(dot(x(best), y(z_mi(best, p)), GREEN))
    body.append(text(x(best), y(z_mi(best, p)) + 33, f"{z_mi(best, p):.4f} at q = {best:g}", "middle", fill=GREEN,
                     mono=True, size=11, weight="600"))
    body.append(dot(x(0.5), y(bsc_mi(0.5, p))))
    body.append(text(x(0.5), y(bsc_mi(0.5, p)) + 22, f"{bsc_mi(0.5, p):.3f}", "middle", mono=True, size=11))
    lx = 262
    body.append(line(lx, 250, lx + 20, 250, GREEN, 2))
    body.append(text(lx + 26, 254, f"Z-channel Z({p:g})", size=11))
    body.append(line(lx, 270, lx + 20, 270, "currentColor", 2))
    body.append(text(lx + 26, 274, f"binary symmetric BSC({p:g})", size=11))
    body.append(line(lx, 290, lx + 20, 290, TEXT3, 1.5, "3 3"))
    body.append(text(lx + 26, 294, f"binary erasure BEC({p:g})", size=11))
    return svg("Mutual information against the input probability q for three channels with noise parameter 0.1. "
               f"The erasure channel peaks at {bec_mi(0.5, p):g} bits and the binary symmetric channel at "
               f"{bsc_mi(0.5, p):.3f} bits, both at q equal to one half. The Z-channel, whose noise touches only "
               f"the 1s, peaks at {z_mi(best, p):.4f} bits at q = {best:g}, to the left of one half.", body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    charts = (("fig-4-1-surprise-curve.svg", surprise_curve()),
              ("fig-4-2-distribution-entropies.svg", distribution_entropies()),
              ("fig-4-3-binary-entropy.svg", binary_entropy()),
              ("fig-4-4-jensen-chord.svg", jensen_chord()),
              ("fig-4-5-joint-table.svg", joint_table()),
              ("fig-4-6-information-diagram.svg", information_diagram()),
              ("fig-4-7-bsc-mutual-information.svg", bsc_mutual_information()),
              ("fig-4-8-input-distribution.svg", input_distribution()))
    for name, content in charts:
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
