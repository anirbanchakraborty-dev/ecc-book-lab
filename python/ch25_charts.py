"""Write the five computed figures of chapter 25 of Error Correction Coding as SVG files.

    python3 python/ch25_charts.py --out figures

Every drawn value comes from the functions in ch25_golay_codes.py, so a figure
is corrected by fixing those functions and running this script again, never by
editing its output. The SVG follows the book's figure style sheet: real text in
IBM Plex, colors from the page's tokens with paper fallbacks, and no element
that the site's figure check would refuse.
"""

import argparse
import os

from ch25_golay_codes import (
    B, BROW, CODE23, CODE24, GROWS, column_weights, decoder_cases, gate_counts, pop)

W, H = 640, 400
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"

ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker></defs>')


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


class LinAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo, hi, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (v - self.lo) / (self.hi - self.lo)


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


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0,
         dash: str = "") -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def polyline(points, stroke: str = "currentColor", width: float = 1.5, arrow: bool = False) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    m = ' marker-end="url(#ah)"' if arrow else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{m}/>'


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, stroke: str, label: str, dash: str = "") -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.5, 3, dash), text(x + 22, y, label, size=12)]


def weight_counts(code) -> dict:
    counts = {}
    for c in code:
        w = pop(c)
        counts[w] = counts.get(w, 0) + 1
    return dict(sorted(counts.items()))


# ---------------------------------------------------------------- the generator matrix

def generator_grid() -> str:
    """G = [I | B] as a 12 by 24 grid of cells, the 1s of B in green."""
    cell, x0, y0, gap = 20, 96, 64, 8

    def cx(j: int) -> float:
        return x0 + j * cell + (gap if j >= 12 else 0)

    body = []
    body.append(text(cx(0) + 6 * cell, 30, "message part, the identity", "middle", size=12, bold=True))
    body.append(text(cx(12) + 6 * cell, 30, "check part, the matrix B", "middle", GREEN, size=12, bold=True))
    body.append(text(x0 - 8, 56, "column", "end", TEXT3, size=11))
    for j in range(24):
        body.append(text(cx(j) + cell / 2, 56, str(j), "middle", TEXT3, mono=True, size=9))
    body.append(text(612, 56, "weight", "middle", TEXT3, size=11))
    for i in range(12):
        y = y0 + i * cell
        body.append(text(x0 - 8, y + 14, f"row {i}", "end", size=12))
        g = GROWS[i]
        for j in range(24):
            bit = (g >> (23 - j)) & 1
            if bit and j < 12:
                body.append(rect(cx(j) + 1, y + 1, cell - 2, cell - 2, TEXT3, TEXT3, 1, 2))
            elif bit:
                body.append(rect(cx(j) + 1, y + 1, cell - 2, cell - 2, GREEN, GREEN, 1, 2))
            else:
                body.append(rect(cx(j) + 1, y + 1, cell - 2, cell - 2, "none", TEXT3, 1, 2))
        body.append(text(612, y + 14, str(pop(g)), "middle", mono=True, size=12))
    body.append(rect(cx(12) - 2, y0 - 2, 11 * cell + 2, 11 * cell + 2, "none", CYAN, 1.5, 3, "3 3"))
    ly = y0 + 12 * cell + 34
    body += legend_box(x0, ly, TEXT3, TEXT3, "a 1 of the identity")
    body += legend_box(x0 + 170, ly, GREEN, GREEN, "a 1 of B")
    body += legend_box(x0 + 290, ly, "none", TEXT3, "a 0")
    body += legend_box(x0, ly + 24, "none", CYAN, "the 11 by 11 block drawn from the squares modulo 11", "3 3")
    body.append(text(x0 + 22, ly + 46, "column 23 and row 11 border it, with a 0 where they meet", size=12))
    return svg("The generator matrix of the extended binary Golay code as a 12 by 24 grid, the identity on "
               "the left and the matrix B on the right, every row of weight 8 or 12", body, 640, ly + 62)


# ---------------------------------------------------------------- the row overlaps

def row_overlaps() -> str:
    """Row 0 of B against rows 1 to 5 and row 11, shared 1s outlined in green."""
    cell, x0, gap = 28, 130, 8

    def cx(j: int) -> float:
        return x0 + j * cell + (gap if j == 11 else 0)

    def draw_row(y: float, row: int, ref: int) -> list:
        out = []
        for j in range(12):
            bit = B[row][j]
            shared = bit and B[ref][j] and row != ref
            if shared:
                out.append(rect(cx(j) + 2, y, cell - 4, cell, BG2, GREEN, 2, 3))
                out.append(text(cx(j) + cell / 2, y + 19, "1", "middle", GREEN, mono=True, size=13, bold=True))
            else:
                out.append(rect(cx(j) + 2, y, cell - 4, cell, "none", TEXT3, 1, 3))
                out.append(text(cx(j) + cell / 2, y + 19, str(bit), "middle", "" if bit else TEXT3,
                                mono=True, size=13))
        return out

    body = []
    body.append(text(x0 - 10, 30, "column", "end", TEXT3, size=11))
    for j in range(12):
        body.append(text(cx(j) + cell / 2, 30, str(j), "middle", TEXT3, mono=True, size=11))
    body.append(text(500, 30, "shared 1s", "start", TEXT3, size=11))
    body.append(text(590, 30, "shift", "start", TEXT3, size=11))
    body.append(text(x0 - 10, 63, "row 0", "end", size=12, bold=True))
    body += draw_row(44, 0, 0)
    body.append(text(500, 63, "the reference row", "start", TEXT3, size=12))
    body.append(line(x0 - 50, 88, 630, 88, TEXT3, 1, "3 3"))
    rows = [1, 2, 3, 4, 5, 11]
    for k, row in enumerate(rows):
        y = 100 + 36 * k
        body.append(text(x0 - 10, y + 19, f"row {row}", "end", size=12))
        body += draw_row(y, row, 0)
        inside = sum(1 for j in range(11) if B[row][j] and B[0][j])
        last = B[row][11] & B[0][11]
        body.append(text(500, y + 19, f"{inside} + {last} = {inside + last}", "start", GREEN, mono=True, size=12))
        body.append(text(590, y + 19, str(row) if row < 11 else "border", "start", TEXT3, mono=row < 11,
                         size=12))
    note_y = 100 + 36 * len(rows) + 24
    body.append(text(x0 - 50, note_y, "Each count adds the shared 1s in columns 0 to 10 to those in column 11.",
                     size=12))
    body.append(text(x0 - 50, note_y + 18, "Every count is even, so B times B is 0 off the diagonal, modulo 2.",
                     size=12))
    body += legend_box(x0 - 50, note_y + 46, BG2, GREEN, "a 1 the row shares with row 0")
    return svg("Rows 1 to 5 and row 11 of the matrix B compared with row 0: rows 1 to 5 share three 1s in the "
               "first eleven columns and one in the last, row 11 shares six", body, 640, note_y + 62)


# ---------------------------------------------------------------- the weight distributions

def weight_distributions() -> str:
    """Bars of A_w for the (24,12,8) and (23,12,7) codes on one count axis."""
    w_, h_ = 760, 380
    y0, y1 = 50, 320
    ya = LinAxis(0, 3000, y1, y0)
    body = []
    for v in range(0, 3001, 500):
        body.append(text(62, ya(v) + 4, str(v), "end", TEXT3, mono=True, size=11))
    body.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle" '
                f'font-size="12">codewords of weight w</text>')
    panels = [(70, 320, weight_counts(CODE24), GREEN, "extended code, length 24", 30),
              (380, 740, weight_counts(CODE23), CYAN, "perfect code, length 23", 28)]
    for px0, px1, counts, color, title, bar in panels:
        for v in range(500, 3001, 500):
            body.append(line(px0, ya(v), px1, ya(v), TEXT3, 0.5, "3 3"))
        body.append(line(px0, y1, px1, y1, "currentColor", 1))
        body.append(line(px0, y0, px0, y1, "currentColor", 1))
        body.append(text((px0 + px1) / 2, 30, title, "middle", color, size=12, bold=True))
        slot = (px1 - px0) / len(counts)
        for k, (w, a) in enumerate(counts.items()):
            mid = px0 + slot * (k + 0.5)
            top = ya(a)
            body.append(rect(mid - bar / 2, top, bar, y1 - top, color, color, 1, 0))
            body.append(text(mid, top - 6, str(a), "middle", color, mono=True, size=11))
            body.append(text(mid, y1 + 16, str(w), "middle", mono=True, size=11))
        body.append(text((px0 + px1) / 2, y1 + 38, "weight w", "middle", TEXT3, size=12))
    return svg("Weight distributions of the extended (24,12,8) and perfect (23,12,7) binary Golay codes: "
               "1, 759, 2576, 759, 1 and 1, 253, 506, 1288, 1288, 506, 253, 1", body, w_, h_)


# ---------------------------------------------------------------- the decoder outcomes

def decoder_outcomes() -> str:
    """One 100 percent bar per error weight, split by the decoder step that resolved it.

    Each segment carries its step letter as text, so the colour only repeats what the
    label says: green for a corrected pattern, red for a flagged one."""
    table = decoder_cases()
    x0, x1, bh = 110, 530, 38
    xa = LinAxis(0, 1, x0, x1)
    letters = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    body = []
    body.append(text(x0, 46, "share of the error patterns of each weight resolved by each step", "start",
                     TEXT3, size=12))
    body.append(text(540, 46, "patterns", "start", TEXT3, size=12))
    for k, w in enumerate(range(1, 5)):
        y = 70 + 58 * k
        counts = table[w]
        total = sum(counts)
        body.append(text(x0 - 10, y + bh / 2 + 4, f"weight {w}", "end", size=12))
        left = 0
        for c, n in enumerate(counts):
            if n == 0:
                continue
            color = RED if c == 4 else GREEN
            a, b = xa(left / total), xa((left + n) / total)
            body.append(rect(a + 1, y, b - a - 2, bh, BG2, color, 1.5, 3))
            body.append(text((a + b) / 2, y + 15, letters[c], "middle", color, mono=True, size=11, bold=True))
            body.append(text((a + b) / 2, y + 30, str(n), "middle", "", mono=True, size=11))
            left += n
        verdict = "flagged" if counts[4] == total else "corrected"
        body.append(text(540, y + bh / 2 + 4, f"{total} {verdict}", "start", mono=False, size=12))
    ty = 70 + 58 * 3 + bh + 4
    for p in range(0, 101, 25):
        x = xa(p / 100)
        body.append(line(x, ty, x, ty + 6, TEXT3, 1))
        body.append(text(x, ty + 20, f"{p}%", "middle", TEXT3, mono=True, size=11))
    ly = ty + 50
    body += legend_box(x0, ly, BG2, GREEN, "steps (a) to (d), the error corrected")
    body += legend_box(x0 + 270, ly, BG2, RED, "step (e), the word flagged")
    return svg("How the five steps of the Golay decoder resolve every error pattern of weight 1 to 4, each bar "
               "segment labelled with its step and its count: all 2324 patterns of weight 1 to 3 are corrected by "
               "steps (a) to (d) and all 10626 of weight 4 are flagged by step (e)",
               body, 640, ly + 18)


# ---------------------------------------------------------------- the datapath

def xor_gate(cx: float, cy: float) -> list:
    """An ANSI distinctive-shape XOR gate 40 wide and 32 tall, output at cx + 20."""
    left, h = cx - 20, 16
    body = (f'<path d="M{fmt(left)} {fmt(cy - h)}Q{fmt(left + 10)} {fmt(cy)} {fmt(left)} {fmt(cy + h)}'
            f'Q{fmt(left + 29)} {fmt(cy + h)} {fmt(cx + 20)} {fmt(cy)}Q{fmt(left + 29)} {fmt(cy - h)} '
            f'{fmt(left)} {fmt(cy - h)}z" fill="{BG2}" stroke="currentColor" stroke-width="1.5"/>')
    back = (f'<path d="M{fmt(left - 5)} {fmt(cy - h)}Q{fmt(left + 5)} {fmt(cy)} {fmt(left - 5)} {fmt(cy + h)}" '
            f'fill="none" stroke="currentColor" stroke-width="1.5"/>')
    return [body, back]


def block(x: float, y: float, w: float, h: float, lines: list, stroke: str = "currentColor") -> list:
    out = [rect(x, y, w, h, BG2, stroke, 1.5, 3)]
    step = h / (len(lines) + 1)
    for k, (s, mono) in enumerate(lines):
        out.append(text(x + w / 2, y + step * (k + 1) + 4, s, "middle", "" if k == 0 else TEXT3, mono=mono,
                        size=12 if k == 0 else 11, bold=k == 0))
    return out


def decoder_datapath() -> str:
    """The XOR tree of syndrome bit 0 beside a block diagram of the whole decoder."""
    stats = gate_counts()
    cols = column_weights()
    bank = 12 * stats["chain3"] + stats["chain4"]
    inputs = ["r0"] + [f"r{12 + k}" for k in range(12) if B[k][0]]
    body = [ARROW]
    body.append(text(20, 40, f"syndrome bit 0: {len(inputs)} inputs, {len(inputs) - 1} XOR gates, 3 levels",
                     "start", size=11, bold=True))
    ys = [75, 91, 127, 143, 179, 195, 231, 247]
    for name, y in zip(inputs, ys):
        body.append(text(52, y + 4, name, "end", mono=True, size=11))
        body.append(line(58, y, 109, y, "currentColor", 1.5))
    level1 = [83, 135, 187, 239]
    for cy in level1:
        body += xor_gate(130, cy)
    level2 = [109, 213]
    for k, cy in enumerate(level2):
        for src, ty in ((level1[2 * k], cy - 8), (level1[2 * k + 1], cy + 8)):
            body.append(polyline([(150, src), (170, src), (170, ty), (189, ty)]))
        body += xor_gate(210, cy)
    for src, ty in ((level2[0], 153), (level2[1], 169)):
        body.append(polyline([(230, src), (250, src), (250, ty), (269, ty)]))
    body += xor_gate(290, 161)
    body.append(polyline([(310, 161), (346, 161)], arrow=True))
    body.append(text(354, 165, "s0", "start", mono=True, size=12, bold=True))
    notes = [f"Bits s1 to s10 have the same shape. Bit s11 has {cols[11] + 1}",
             f"inputs and needs {cols[11]} gates in {stats['syndrome_depth']} levels, so the 12 trees",
             f"hold {stats['syndrome_xor']} two-input XOR gates in all."]
    for k, s in enumerate(notes):
        body.append(text(20, 296 + 17 * k, s, "start", TEXT3, size=11))
    body.append(line(392, 30, 392, 455, TEXT3, 1, "3 3"))
    # the block diagram on the right
    c1, c2, bw = 490, 680, 140
    body.append(text(c1, 44, "r, 24 bits", "middle", mono=True, size=11))
    body.append(polyline([(c1, 50), (c1, 66)], arrow=True))
    body += block(c1 - bw / 2, 68, bw, 50, [("syndrome network", False),
                                           (f"{stats['syndrome_xor']} XOR2, {stats['syndrome_depth']} levels",
                                            True)])
    body.append(polyline([(c1 + bw / 2, 93), (c2 - bw / 2 - 2, 93)], arrow=True))
    body.append(text((c1 + c2) / 2, 86, "s", "middle", mono=True, size=11))
    body += block(c2 - bw / 2, 68, bw, 50, [("s B network", False),
                                           (f"{stats['sb_xor']} XOR2, {stats['sb_depth']} levels", True)])
    body.append(polyline([(c1, 118), (c1, 166)], arrow=True))
    body.append(text(c1 + 8, 146, "s", "start", mono=True, size=11))
    body.append(polyline([(c2, 118), (c2, 166)], arrow=True))
    body.append(text(c2 + 8, 146, "s B", "start", mono=True, size=11))
    for cx, name, word in ((c1, "test bank 1", "s"), (c2, "test bank 2", "s B")):
        body += block(cx - bw / 2, 168, bw, 66, [(name, False), (f"13 weight tests on {word}", False),
                                                (f"{bank} AND/OR, {stats['shared_inverters']} NOT", True)])
        body.append(polyline([(cx - 35, 234), (cx - 35, 286)], arrow=True))
        body.append(text(cx - 43, 264, "13 flags", "end", TEXT3, size=11))
        body.append(polyline([(cx + 35, 234), (cx + 35, 286)], arrow=True))
        body.append(text(cx + 27, 264, "13 words", "end", TEXT3, size=11))
    body += block(c1 - bw / 2 + 10, 288, c2 - c1 + bw - 20, 58,
                  [("priority select", False), ("the first test that fires picks e from", False),
                   ("s, s + bi, s B or s B + bi", True)], GREEN)
    mid = (c1 + c2) / 2
    body.append(polyline([(mid, 346), (mid, 372)], arrow=True))
    body.append(text(mid + 8, 364, "e", "start", mono=True, size=11))
    body += block(mid - 100, 374, 200, 44, [("correction", False),
                                            (f"r + e, {stats['correction_xor']} XOR2", True)])
    body.append(text(mid - 146, 400, "r", "end", mono=True, size=11))
    body.append(polyline([(mid - 140, 396), (mid - 102, 396)], arrow=True))
    body.append(polyline([(mid, 418), (mid, 440)], arrow=True))
    body.append(text(mid, 458, "corrected word, or the flag", "middle", size=12))
    return svg("The decoder of the extended binary Golay code in hardware: the XOR tree for syndrome bit 0, and "
               "the syndrome network, the s B network, two banks of weight tests whose flags and candidate words "
               "feed a priority select, and a row of XOR gates that adds the chosen error pattern e to r",
               body, 760, 470)


CHARTS = {
    "fig-25-2-generator-grid.svg": generator_grid,
    "fig-25-3-row-overlaps.svg": row_overlaps,
    "fig-25-4-weight-distributions.svg": weight_distributions,
    "fig-25-6-decoder-outcomes.svg": decoder_outcomes,
    "fig-25-7-decoder-datapath.svg": decoder_datapath,
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
