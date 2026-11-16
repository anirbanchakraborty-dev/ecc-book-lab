"""Write the three computed figures of the project chapter on a SEC-DED encoder and decoder as SVG files.

    python3 python/ch30_charts.py --out figures

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction
Coding. The check matrices are drawn from the constructions in
ch20_secded_and_hsiao_codes.py, and every gate count, depth and percentage is
read from data/ch30/secded72_gates.json, the file that synth/synth_secded72.py
writes from the Yosys netlist. A figure is corrected by fixing those sources and
running this script again, never by editing its output. The SVG follows the
book's figure style sheet: real text in IBM Plex, colors from the page's tokens
with paper fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import json
import os
from pathlib import Path

from ch20_secded_and_hsiao_codes import extended_hamming, hsiao_72_64

LAB = Path(__file__).resolve().parent.parent
GATES = LAB / "data" / "ch30" / "secded72_gates.json"

W = 640
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
RED = "var(--red,#a1211f)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


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


def svg(title: str, body: list, w: int, h: int) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, label: str, stroke: str = "none") -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.2, rx=3), text(x + 22, y, label, size=12)]


def measured() -> dict:
    return {c["code"]: c for c in json.loads(GATES.read_text(encoding="utf-8"))["codes"]}


def row_loads(cols: list, r: int) -> list:
    return [sum((c >> i) & 1 for c in cols) for i in range(r)]


def hand_counts(code: dict) -> dict:
    """The chapter's hand counts, from the rows and check equations the data file records."""
    n, k, r = code["n"], code["k"], code["r"]
    syndrome = sum(w - 1 for w in code["row_weights"])
    encoder = sum(f - 1 for f in code["check_fan_in"])
    base = syndrome + n * (r - 1) + k + (n - 1) + 1
    if code["code"] == "hamming":
        flag0, flag1 = base + (r - 1), base + (r - 2)
    else:
        flag0, flag1 = base + (r - 1), base + 2 * (r - 1)
    return {"syndrome": syndrome, "encoder": encoder, "flag0": flag0, "flag1": flag1}


def matrix_chart() -> str:
    """The two (72,64) check matrices, one square per 1, in the bit order of chapter 20."""
    codes = [
        ("Hsiao (72,64), balanced", hsiao_72_64(), list(range(64, 72))),
        ("Extended Hamming (72,64)", extended_hamming(72, 7), [0, 1, 2, 4, 8, 16, 32, 64]),
    ]
    h = 352
    x0, pitch, cell = 52, 7, 6
    row_h = 11
    body = [text(x0 + 72 * pitch + 12, 22, "ones", fill=TEXT3, size=11)]
    for panel, (label, cols, checks) in enumerate(codes):
        top = 34 + panel * 150
        loads = row_loads(cols, 8)
        body.append(text(x0, top - 10, f"{label}, {sum(loads)} ones, largest row {max(loads)}", size=13,
                         bold=True))
        for i in range(8):
            y = top + i * row_h
            body.append(rect(x0, y, 72 * pitch - 1, row_h - 2, BG2))
            body.append(text(x0 - 8, y + 8, f"s{i}", "end", mono=True, size=10))
            for j, c in enumerate(cols):
                if (c >> i) & 1:
                    fill = GREEN if j in checks else "currentColor"
                    body.append(rect(x0 + j * pitch, y, cell, row_h - 2, fill))
            body.append(text(x0 + 72 * pitch + 12, y + 8, str(loads[i]), mono=True, size=10))
        base = top + 8 * row_h + 2
        for j in checks:
            cx = x0 + j * pitch + cell / 2
            body.append(f'<path d="M{fmt(cx)} {fmt(base)}l3.5 6h-7z" fill="{GREEN}"/>')
        for j in (0, 8, 16, 24, 32, 40, 48, 56, 64, 71):
            cx = x0 + j * pitch + cell / 2
            body.append(text(cx, base + 22, str(j), "middle", mono=True, size=10, fill=TEXT3))
    body.append(text(x0, h - 30, "Column = bit position 0 to 71 of the stored word, row = syndrome bit.",
                     size=12))
    body.append(f'<path d="M{x0 + 3.5} {h - 18}l3.5 6h-7z" fill="{GREEN}"/>')
    body.append(text(x0 + 14, h - 12, "marks a check-bit position, whose 1s are drawn in green.", size=12))
    return svg("The two (72,64) check matrices, one square for each 1", body, W, h)


def gate_chart() -> str:
    """Gates counted by hand against gates in the Yosys netlist, for both (72,64) codes."""
    m = measured()
    blocks = [
        ("Syndrome stage", "syndrome", lambda c: c["syndrome_whole"]["two_input"]),
        ("Encoder", "encoder", lambda c: c["encoder_whole"]["two_input"]),
        ("Decoder, no-match flag", "flag0", lambda c: c["decoder_whole_flag0"]["two_input"]),
        ("Decoder, parity flag", "flag1", lambda c: c["decoder_whole_flag1"]["two_input"]),
    ]
    codes = [("Hsiao", m["hsiao"]), ("extended Hamming", m["hamming"])]
    x0, x1 = 200, 560
    top = 1000

    def x(v: float) -> float:
        return x0 + (x1 - x0) * v / top

    body = []
    body += legend_box(x0, 26, BG2, "counted by hand, two-input gates", "currentColor")
    body += legend_box(x0 + 236, 26, GREEN, "Yosys and ABC, whole block")
    y = 62
    for title, key, got in blocks:
        body.append(text(16, y + 2, title, size=12, bold=True))
        y += 10
        for name, code in codes:
            hand = hand_counts(code)[key]
            meas = got(code)
            body.append(text(x0 - 8, y + 15, name, "end", size=12))
            body.append(rect(x0, y + 2, x(hand) - x0, 10, BG2, "currentColor", 1.2))
            body.append(text(x(hand) + 6, y + 11, str(hand), mono=True, size=11))
            body.append(rect(x0, y + 15, x(meas) - x0, 10, GREEN))
            body.append(text(x(meas) + 6, y + 24, str(meas), mono=True, size=11))
            y += 34
        y += 8
    body.append(f'<polyline points="{x0},44 {x0},398 {x1},398" fill="none" stroke="currentColor" stroke-width="1"/>')
    for v in (0, 200, 400, 600, 800, 1000):
        body.append(line(x(v), 398, x(v), 403))
        body.append(text(x(v), 416, str(v), "middle", mono=True, size=11))
    body.append(text((x0 + x1) / 2, 436, "Two-input gates, with inverters not counted", "middle"))
    return svg("Gates counted by hand and gates in the synthesized netlist, for both (72,64) codes", body,
               W, 452)


def saving_chart() -> str:
    """What Hsiao saves or costs against extended Hamming, as a percent of extended Hamming's figure."""
    m = measured()
    h, e = m["hsiao"], m["hamming"]

    def trees(code, key):
        return sum(t["gates"] for t in code[key])

    def depth(code):
        return max(t["depth"] for t in code["syndrome_trees"])

    rows = [
        ("separate trees", "syndrome XOR gates", trees(h, "syndrome_trees"), trees(e, "syndrome_trees")),
        ("separate trees", "encoder XOR gates", trees(h, "encoder_trees"), trees(e, "encoder_trees")),
        ("separate trees", "syndrome levels", depth(h), depth(e)),
        ("one whole block", "syndrome gates", h["syndrome_whole"]["two_input"], e["syndrome_whole"]["two_input"]),
        ("one whole block", "encoder gates", h["encoder_whole"]["two_input"], e["encoder_whole"]["two_input"]),
        ("one whole block", "decoder gates, no-match flag", h["decoder_whole_flag0"]["two_input"],
         e["decoder_whole_flag0"]["two_input"]),
        ("one whole block", "decoder gates, parity flag", h["decoder_whole_flag1"]["two_input"],
         e["decoder_whole_flag1"]["two_input"]),
        ("one whole block", "levels to the corrected data", h["decoder_whole_flag0"]["data_depth"],
         e["decoder_whole_flag0"]["data_depth"]),
    ]
    mid, span = 380, 5.5          # the zero line, and pixels per percent
    body = []
    body += legend_box(232, 26, GREEN, "Hsiao needs fewer")
    body += legend_box(392, 26, ORANGE, "Hsiao needs more")
    y = 62
    group = ""
    for kind, label, hv, ev in rows:
        if kind != group:
            body.append(text(16, y + 6, f"Measured as {kind}", size=12, bold=True))
            group = kind
            y += 16
        saving = 100 * (ev - hv) / ev
        word = f"{abs(saving):.1f}% {'fewer' if saving >= 0 else 'more'}"
        body.append(text(16, y + 12, label, size=12))
        if saving >= 0:
            body.append(rect(mid, y + 2, saving * span, 14, GREEN))
            body.append(text(mid + saving * span + 6, y + 13, word, mono=True, size=11))
        else:
            body.append(rect(mid + saving * span, y + 2, -saving * span, 14, ORANGE))
            body.append(text(mid + saving * span - 6, y + 13, word, "end", mono=True, size=11))
        y += 30
        if label == "syndrome levels":
            y += 6
    body.append(line(mid, 44, mid, 368, "currentColor", 1.5))
    body.append(line(mid - 30 * span, 368, mid + 30 * span, 368))
    for v in (-30, -20, -10, 0, 10, 20, 30):
        body.append(line(mid + v * span, 368, mid + v * span, 373))
        body.append(text(mid + v * span, 386, f"{abs(v)}", "middle", mono=True, size=11))
    body.append(text(mid - 15 * span, 406, "Hsiao more", "middle", size=12, fill=TEXT3))
    body.append(text(mid + 15 * span, 406, "Hsiao fewer", "middle", size=12, fill=TEXT3))
    body.append(text(W / 2, 428, "Percent of the extended Hamming figure", "middle"))
    return svg("What the Hsiao code saves or costs against extended Hamming, measured", body, W, 444)


CHARTS = {
    "fig-30-2-check-matrices.svg": matrix_chart,
    "fig-30-6-gate-counts.svg": gate_chart,
    "fig-30-7-hsiao-saving.svg": saving_chart,
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
