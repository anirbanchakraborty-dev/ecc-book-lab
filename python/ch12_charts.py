"""Write the four computed figures of chapter 12 of Error Correction Coding as SVG files.

    python3 python/ch12_charts.py
    python3 python/ch12_charts.py --out figures

The two tree diagrams are laid out from the same recursive split that
hdl/parity/xor_tree.sv uses, and the waveform and the coverage chart come from the
functions in ch12_lab_first_hardware_checker.py, with the seed, crossover
probabilities and word counts the cocotb tests use, so a figure is corrected by
fixing those functions and running this script again, never by editing its
output. Inside the hub that holds both this repository and the book's website,
the default output folder is the site's figure folder for the book. In a clone of
this repository alone it is a folder named figures in the working directory. The
SVG follows the book's figure style sheet: real text in IBM Plex, colors from the
page's tokens with paper fallbacks, and no element that the site's figure check
would refuse.
"""

import argparse
import os
from math import log10

from ch12_lab_first_hardware_checker import (P_VALUES, STEP_NS, WAVE_P, WAVE_WORDS, bsc_stream, codeword,
                                             coverage, flag, run_bsc, weight, wilson)

W = 640
X0, X1 = 84, 604
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker>'
         '<marker id="ahg" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="var(--green,#00614a)"/></marker></defs>')


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
         size: int = 0, weight_: str = "") -> str:
    attrs = [f'x="{fmt(x)}"', f'y="{fmt(y)}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if fill:
        attrs.append(f'fill="{fill}"')
    if mono:
        attrs.append(f'font-family="{MONO}"')
    if size:
        attrs.append(f'font-size="{size}"')
    if weight_:
        attrs.append(f'font-weight="{weight_}"')
    return f'<text {" ".join(attrs)}>{s}</text>'


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def polyline(points: list, stroke: str, width: float, dash: str = "", arrow: bool = False) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{"ahg" if stroke == GREEN else "ah"})"' if arrow else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}{m}/>'


def svg(title: str, h: int, body: list, w: int = W) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_line(x: float, y: float, color: str, label: str, dash: str = "", width: float = 2) -> list:
    return [line(x, y - 4, x + 24, y - 4, color, width, dash), text(x + 32, y, label, size=12)]


# The XOR trees -----------------------------------------------------------------------

def split_tree(lo: int, hi: int):
    """The tree xor_tree.sv builds over bits lo to hi - 1: a leaf index, or a (low, high) pair."""
    if hi - lo == 1:
        return lo
    mid = lo + (hi - lo) // 2
    return (split_tree(lo, mid), split_tree(mid, hi))


def height(node) -> int:
    return 0 if isinstance(node, int) else 1 + max(height(node[0]), height(node[1]))


def deepest_path(node) -> set:
    """The gates and the leaf on one longest input-to-output path, ties broken toward the high half."""
    path = set()
    while not isinstance(node, int):
        path.add(id(node))
        low, high = node
        node = high if height(high) >= height(low) else low
    path.add(("leaf", node))
    return path


def xor_gate(cx: float, cy: float, stroke: str, width: float) -> list:
    """The ANSI/IEEE Std 91 distinctive-shape XOR gate, 40 wide and 32 tall, centred on (cx, cy).

    The body is the OR shape, a concave back and a pointed front, and the second
    curve behind it is what marks it as exclusive OR. Its output is at cx + 20.
    """
    left, h = cx - 20, 16
    body = (f'<path d="M{fmt(left)} {fmt(cy - h)}Q{fmt(left + 10)} {fmt(cy)} {fmt(left)} {fmt(cy + h)}'
            f'Q{fmt(left + 29)} {fmt(cy + h)} {fmt(cx + 20)} {fmt(cy)}'
            f'Q{fmt(left + 29)} {fmt(cy - h)} {fmt(left)} {fmt(cy - h)}z"'
            f' fill="{BG2}" stroke="{stroke}" stroke-width="{width}"/>')
    back = (f'<path d="M{fmt(left - 5)} {fmt(cy - h)}Q{fmt(left + 5)} {fmt(cy)} {fmt(left - 5)} {fmt(cy + h)}"'
            f' fill="none" stroke="{stroke}" stroke-width="{width}"/>')
    return [body, back]


def tree_figure(width_in: int, prefix: str, out_name: str, title: str) -> str:
    gap, top, leaf_x, dx, gate_w = 34, 40, 64, 130, 40
    root = split_tree(0, width_in)
    levels = height(root)
    on_path = deepest_path(root)
    body = [ARROW]
    wires, boxes, labels = [], [], []
    count = [0]

    def place(node):
        """Return the output point of node, drawing it and everything below it."""
        if isinstance(node, int):
            y = top + node * gap
            labels.append(text(leaf_x - 8, y + 4, f"{prefix}{node}", "end", mono=True, size=12))
            return leaf_x, y, ("leaf", node) in on_path
        (lx, ly, lp), (hx, hy, hp) = place(node[0]), place(node[1])
        cx = leaf_x + 36 + height(node) * dx
        cy = (ly + hy) / 2
        hot = id(node) in on_path
        for sx, sy, child_hot, ty in ((lx, ly, lp, cy - 8), (hx, hy, hp, cy + 8)):
            color = GREEN if hot and child_hot else "currentColor"
            bend = cx - gate_w / 2 - 18
            # The wire stops on the gate's back curve, which sits about 1 px left of the body here.
            wires.append(polyline([(sx, sy), (bend, sy), (bend, ty), (cx - gate_w / 2 - 1, ty)], color,
                                  2 if color == GREEN else 1.2))
        boxes.extend(xor_gate(cx, cy, GREEN if hot else "currentColor", 2 if hot else 1.5))
        count[0] += 1
        return cx + gate_w / 2, cy, hot

    ox, oy, _ = place(root)
    wires.append(polyline([(ox, oy), (ox + 46, oy)], GREEN, 2, arrow=True))
    labels.append(text(ox + 54, oy + 4, out_name, mono=True, size=12))
    bottom = top + (width_in - 1) * gap + 34
    for level in range(1, levels + 1):
        labels.append(text(leaf_x + 36 + level * dx, bottom, f"level {level}", "middle", TEXT3, size=11))
    labels += legend_line(leaf_x - 40, bottom + 22, GREEN, f"a deepest path, {levels} gates")
    labels.append(text(W - 16, bottom + 22, f"{count[0]} gates in all", "end", size=12))
    body += wires + boxes + labels
    return svg(title, bottom + 36, body)


def generator_tree() -> str:
    return tree_figure(
        7, "d", "parity",
        "The parity tree inside parity_gen for K = 7 data bits, d0 to d6, built by the split rule of xor_tree. "
        "The seven bits split into a low group of three, d0 to d2, and a high group of four, d3 to d6. The low "
        "three split again into d0 alone and the pair d1, d2, and the high four into the pairs d3, d4 and d5, "
        "d6. Each split adds one XOR gate, six in all, and the longest path from an input to the parity output, "
        "drawn in green, passes through three gates.")


def checker_tree() -> str:
    return tree_figure(
        8, "r", "error",
        "The checking tree inside parity_check for a received word of n = 8 bits, r0 to r7. The split rule "
        "halves the eight bits into r0 to r3 and r4 to r7, then into four pairs. Four gates combine the pairs, "
        "two combine those results, and one more makes the error flag, seven gates in all, and every input "
        "reaches the output through exactly three gates, so every path is a deepest path. The green one is the "
        "path reached by always taking the high half.")


# The waveform -------------------------------------------------------------------------

def waveform() -> str:
    w, tx0, seg = 760, 96, 54
    k = 7
    stream = bsc_stream(WAVE_P, k)
    words = []
    for _ in range(WAVE_WORDS):
        data, mask = next(stream)
        c = codeword(data, k)
        words.append((data, c, mask, c ^ mask, flag(c ^ mask)))
    rows = ("data", "codeword", "error_mask", "flips", "received", "error", "verdict")
    yc = [58 + 32 * r for r in range(len(rows))]
    h = yc[-1] + 26
    body = [text(12, 24, "time, ns", fill=TEXT3, size=11)]
    for i in range(WAVE_WORDS + 1):
        x = tx0 + i * seg
        body.append(line(x, 30, x, 34))
        body.append(text(x, 24, str(i * STEP_NS), "middle", mono=True, size=11))
        body.append(line(x, 36, x, h - 8, TEXT3, 1, dash="2 3", opacity="0.5"))
    for r, name in enumerate(rows):
        derived = name in ("flips", "verdict")
        body.append(text(12, yc[r] + 4, name, fill=TEXT3 if derived else "", mono=not derived, size=12))
    hi_lo = []
    for i, (data, c, mask, received, f) in enumerate(words):
        xa, xb = tx0 + i * seg, tx0 + (i + 1) * seg
        mid = (xa + xb) / 2
        outcome = GREEN if f else ORANGE
        for r, value, width_bits in ((0, data, 7), (1, c, 8), (2, mask, 8), (4, received, 8)):
            y = yc[r]
            pts = [(xa, y), (xa + 4, y - 10), (xb - 4, y - 10), (xb, y), (xb - 4, y + 10), (xa + 4, y + 10)]
            body.append(f'<polygon points="{" ".join(f"{fmt(px)},{fmt(py)}" for px, py in pts)}"'
                        f' fill="{BG2}" stroke="currentColor" stroke-width="1"/>')
            body.append(text(mid, y + 3, f"{value:0{width_bits}b}", "middle", mono=True, size=9))
        body.append(text(mid, yc[3] + 4, str(weight(mask)), "middle", outcome, mono=True, size=11))
        body.append(text(mid, yc[6] + 4, "caught" if f else "slipped", "middle", outcome, size=10))
        level = yc[5] - 10 if f else yc[5] + 10
        hi_lo += [(xa, level), (xb, level)]
    body.append(polyline(hi_lo, "currentColor", 1.5))
    return svg(
        "The first twelve words of the p = 0.2 fault-injection run with seed 12, one word every 10 ns, as the "
        "Verilator simulation records them. Each row of boxes is a bus, its value written in binary with the "
        "highest bit on the left: the seven data bits, the eight-bit codeword whose top bit is the parity bit, "
        "the error mask, and the received word, which is the codeword with the masked bits flipped. The flips "
        "row counts the 1s in the mask. The error line rises in the five words whose mask flips an odd number "
        "of bits, which the checker catches, and stays low in the seven words with an even number of flips, "
        "which slip past it although every one of them is damaged.",
        h, body, w)


# Measured coverage --------------------------------------------------------------------

def measured_coverage() -> str:
    h = 446
    ax = LogAxis(0.005, 0.6, X0, X1)
    top0, top1 = 36, 196
    bot0, bot1 = 262, 392
    ay = LinAxis(0.45, 1.0, top1, top0)
    ad = LinAxis(-0.010, 0.008, bot1, bot0)
    body = []
    for y0, y1 in ((top0, top1), (bot0, bot1)):
        body.append(f'<polyline points="{X0},{y0} {X0},{y1} {X1},{y1}" fill="none" stroke="currentColor"'
                    ' stroke-width="1"/>')
    for v in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        py = ay(v)
        body += [line(X0 - 4, py, X0, py), line(X0, py, X1, py, TEXT3, 1, opacity="0.35"),
                 text(X0 - 8, py + 4, f"{v:.1f}", "end", mono=True, size=11)]
    for v in (-0.008, -0.004, 0.0, 0.004, 0.008):
        py = ad(v)
        body += [line(X0 - 4, py, X0, py), line(X0, py, X1, py, TEXT3, 1, opacity="0.35"),
                 text(X0 - 8, py + 4, f"{v:.3f}" if v else "0", "end", mono=True, size=11)]
    for v, lab in ((0.01, "0.01"), (0.02, "0.02"), (0.05, "0.05"), (0.1, "0.1"), (0.2, "0.2"), (0.5, "0.5")):
        px = ax(v)
        body += [line(px, top1, px, top1 + 4), line(px, bot1, px, bot1 + 4),
                 text(px, bot1 + 20, lab, "middle", mono=True, size=11)]
    body.append(text((X0 + X1) / 2, h - 14, "crossover probability p, on a log scale", "middle"))
    body.append(f'<text transform="translate(20 {(top0 + top1) // 2}) rotate(-90)" text-anchor="middle">'
                'detection coverage</text>')
    body.append(f'<text transform="translate(20 {(bot0 + bot1) // 2}) rotate(-90)" text-anchor="middle">'
                'difference</text>')
    curve = []
    steps = 200
    for i in range(steps + 1):
        p = 10 ** (log10(0.005) + (log10(0.6) - log10(0.005)) * i / steps)
        curve.append((ax(p), ay(coverage(p, 8))))
    body.append(polyline(curve, CYAN, 2))
    body.append(line(X0, ad(0), X1, ad(0), CYAN, 2, dash="6 4"))
    for p in P_VALUES:
        damaged, _, slipped = run_bsc(p)
        caught = damaged - slipped
        lo, hi = wilson(caught, damaged)
        c = coverage(p, 8)
        miss = not lo <= c <= hi
        color = RED if miss else "currentColor"
        px = ax(p)
        body.append(dot(px, ay(caught / damaged), "currentColor", 3.5))
        body.append(line(px, ad(lo - c), px, ad(hi - c), color, 2))
        body.append(dot(px, ad(caught / damaged - c), color, 3.5))
    body += legend_line(96, 152, CYAN, "closed form, n = 8")
    body += [dot(108, 168, "currentColor", 3.5), text(128, 172, "measured, 100,000 words each", size=12)]
    body += [line(96, 230, 96, 216, "currentColor", 2), dot(96, 223, "currentColor", 3.5),
             text(108, 228, "measured minus closed form, with its Wilson 95% interval", size=12)]
    body += [line(96, 250, 96, 236, RED, 2), dot(96, 243, RED, 3.5),
             text(108, 248, "an interval that misses the closed form", size=12)]
    return svg(
        "Measured detection coverage of the (8, 7) single parity check code against the closed form, from the "
        "Verilator simulation with seed 12 and 100,000 words at each of eight crossover probabilities. The upper "
        "panel shows the fraction of damaged words the checker caught beside the closed-form curve, which falls "
        "from near 1 at small p toward 128/255 at p = 0.5. The lower panel magnifies the gap. Each point is the "
        "measured coverage minus the closed form, with its Wilson 95% interval. Six intervals contain zero, and "
        "the two at p = 0.1 and p = 0.3, drawn in red, miss it.",
        h, body)


CHARTS = {
    "fig-12-3-generator-tree.svg": generator_tree,
    "fig-12-4-checker-tree.svg": checker_tree,
    "fig-12-6-waveform.svg": waveform,
    "fig-12-7-measured-coverage.svg": measured_coverage,
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
