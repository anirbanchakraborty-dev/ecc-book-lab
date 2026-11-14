"""Write the four charts of front-matter chapter f2 of Error Correction Coding as SVG files.

    python3 python/chf2_charts.py --out figures

Every plotted value is computed here or taken from chf2_ecc_lab_bench.py: the
bits of 100 and of the two numbers made from it, the error left when 1/10 is cut
to k binary places, the spacing of doubles from math.ulp, and the mod-2 product
of a message with the Hamming generator matrix G. A chart is corrected by fixing
this script and running it again, never by editing its output. The SVG follows
the book's figure style sheet: real text in IBM Plex, colors from the page's
tokens with paper fallbacks, and no element that the site's figure check would
refuse.
"""

import argparse
import os
from fractions import Fraction
from math import floor, log2, ulp

from chf2_ecc_lab_bench import G

GREEN = "var(--green,#00614a)"
VIOLET = "var(--violet,#6b3fa0)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def svg(w: int, h: int, title: str, body: list) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def bits_of_100() -> str:
    """Three bytes: 100, 100 shifted right by 2, and 100 with its top bit flipped by an XOR."""
    w, h = 660, 300
    cell, left, top = 40, 196, 62
    rows = [("100", 100, "decimal 100"),
            ("100 >> 2", 100 >> 2, f"decimal {100 >> 2}"),
            ("100 ^ 0b10000000", 100 ^ 0b10000000, f"decimal {100 ^ 0b10000000}")]
    body = [f'<g font-family="{MONO}" font-size="11" fill="{TEXT3}" text-anchor="middle">']
    for i in range(8):
        cx = left + cell * i + cell / 2
        body.append(f'<text x="{fmt(cx)}" y="28">bit {7 - i}</text>')
        body.append(f'<text x="{fmt(cx)}" y="44">{2 ** (7 - i)}</text>')
    body.append("</g>")
    for r, (label, value, note) in enumerate(rows):
        y = top + r * 76
        bits = format(value, "08b")
        body.append(f'<text x="{left - 14}" y="{y + 25}" text-anchor="end" font-family="{MONO}" font-size="13">{label}</text>')
        for i, b in enumerate(bits):
            x = left + cell * i
            changed = r == 2 and bits[i] != format(100, "08b")[i]
            fill = f'fill="{GREEN}" fill-opacity="0.18"' if b == "1" else f'fill="{BG2}"'
            stroke = VIOLET if changed else "currentColor"
            width = 3 if changed else 1.5
            body.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" {fill} stroke="{stroke}" stroke-width="{width}"/>')
            color = GREEN if b == "1" else "currentColor"
            body.append(f'<text x="{x + cell // 2}" y="{y + 26}" text-anchor="middle" font-family="{MONO}" font-size="16"'
                        f' font-weight="600" fill="{color}">{b}</text>')
        body.append(f'<text x="{left + 8 * cell + 16}" y="{y + 25}" font-size="13">{note}</text>')
    body.append(f'<text x="{w // 2}" y="{h - 14}" text-anchor="middle" font-size="12" fill="{TEXT3}">'
                f'Tinted cells hold a 1. The violet outline marks the one bit the XOR changed.</text>')
    return svg(w, h, "The byte 01100100 that holds 100, the byte 00011001 left after shifting it right by two places, "
               "which holds 25, and the byte 11100100 left after an XOR with 10000000 flips its top bit, which holds 228. "
               "Each column is labeled with its bit position and its weight.", body)


def tenth_errors(kmax: int = 24) -> list:
    """The exact error left when 1/10 is cut to k binary places, for k = 1 to kmax."""
    tenth = Fraction(1, 10)
    return [tenth - Fraction(floor(tenth * 2 ** k), 2 ** k) for k in range(1, kmax + 1)]


def tenth_in_binary() -> str:
    """log2 of the truncation error of 1/10 against the bits kept, beside the bound 2^-k."""
    w, h = 700, 420
    x0, x1, y0, y1 = 78, 670, 24, 350
    kmax = 24
    lo, hi = -28, 0
    errors = tenth_errors(kmax)

    def x_of(k: float) -> float:
        return x0 + (x1 - x0) * (k - 0.5) / kmax

    def y_of(v: float) -> float:
        return y0 + (y1 - y0) * (hi - v) / (hi - lo)

    body = [f'<g stroke="currentColor" stroke-width="1" fill="none"><polyline points="{x0},{y0} {x0},{y1} {x1},{y1}"/></g>']
    for t in range(lo, hi + 1, 4):
        y = y_of(t)
        body.append(f'<line x1="{x0 - 4}" y1="{fmt(y)}" x2="{x0}" y2="{fmt(y)}" stroke="currentColor" stroke-width="1"/>')
        if t != lo:
            body.append(f'<line x1="{x0}" y1="{fmt(y)}" x2="{x1}" y2="{fmt(y)}" stroke="{TEXT3}" stroke-width="1" stroke-opacity="0.35"/>')
        body.append(f'<text x="{x0 - 8}" y="{fmt(y + 4)}" text-anchor="end" font-family="{MONO}" font-size="11">{t}</text>')
    for k in range(4, kmax + 1, 4):
        x = x_of(k)
        body.append(f'<line x1="{fmt(x)}" y1="{y1}" x2="{fmt(x)}" y2="{y1 + 4}" stroke="currentColor" stroke-width="1"/>')
        body.append(f'<text x="{fmt(x)}" y="{y1 + 18}" text-anchor="middle" font-family="{MONO}" font-size="11">{k}</text>')
    body.append(f'<text x="{(x0 + x1) // 2}" y="{y1 + 40}" text-anchor="middle">binary places kept, k</text>')
    body.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">log2 of the error</text>')
    bound = " ".join(f"{fmt(x_of(k))},{fmt(y_of(-k))}" for k in range(1, kmax + 1))
    body.append(f'<polyline points="{bound}" fill="none" stroke="{TEXT3}" stroke-width="1.5" stroke-dasharray="5 4"/>')
    for k, e in enumerate(errors, start=1):
        body.append(f'<circle cx="{fmt(x_of(k))}" cy="{fmt(y_of(log2(e)))}" r="4" fill="{GREEN}"/>')
    lx, ly = x_of(13), y_of(-3)
    body.append(f'<line x1="{fmt(lx)}" y1="{fmt(ly)}" x2="{fmt(lx + 28)}" y2="{fmt(ly)}" stroke="{TEXT3}" stroke-width="1.5" stroke-dasharray="5 4"/>')
    body.append(f'<text x="{fmt(lx + 36)}" y="{fmt(ly + 4)}" font-size="12">the bound 2^-k</text>')
    body.append(f'<circle cx="{fmt(lx + 14)}" cy="{fmt(ly + 22)}" r="4" fill="{GREEN}"/>')
    body.append(f'<text x="{fmt(lx + 36)}" y="{fmt(ly + 26)}" font-size="12">1/10 minus its first k bits</text>')
    return svg(w, h, "The error left when 1/10 is cut to its first k binary places, on a base-2 logarithmic scale, for k from "
               "1 to 24. Every point lies below the dashed bound 2 to the minus k and no lower than a fifth of it, and none "
               "reaches zero, because 1/10 has no finite binary expansion.", body)


def ulp_staircase() -> str:
    """log2 of the gap between consecutive doubles against log2 of x, a staircase that climbs one step per binade."""
    w, h = 700, 420
    x0, x1, y0, y1 = 78, 670, 24, 350
    elo, ehi = -6, 12
    ylo, yhi = -60, -38

    def x_of(v: float) -> float:
        return x0 + (x1 - x0) * (v - elo) / (ehi - elo)

    def y_of(v: float) -> float:
        return y0 + (y1 - y0) * (yhi - v) / (yhi - ylo)

    body = [f'<g stroke="currentColor" stroke-width="1" fill="none"><polyline points="{x0},{y0} {x0},{y1} {x1},{y1}"/></g>']
    for t in range(ylo, yhi + 1, 4):
        y = y_of(t)
        body.append(f'<line x1="{x0 - 4}" y1="{fmt(y)}" x2="{x0}" y2="{fmt(y)}" stroke="currentColor" stroke-width="1"/>')
        if t != ylo:
            body.append(f'<line x1="{x0}" y1="{fmt(y)}" x2="{x1}" y2="{fmt(y)}" stroke="{TEXT3}" stroke-width="1" stroke-opacity="0.35"/>')
        body.append(f'<text x="{x0 - 8}" y="{fmt(y + 4)}" text-anchor="end" font-family="{MONO}" font-size="11">{t}</text>')
    for e in range(elo, ehi + 1, 2):
        x = x_of(e)
        body.append(f'<line x1="{fmt(x)}" y1="{y1}" x2="{fmt(x)}" y2="{y1 + 4}" stroke="currentColor" stroke-width="1"/>')
        body.append(f'<text x="{fmt(x)}" y="{y1 + 18}" text-anchor="middle" font-family="{MONO}" font-size="11">{e}</text>')
    body.append(f'<text x="{(x0 + x1) // 2}" y="{y1 + 40}" text-anchor="middle">log2 of x</text>')
    body.append(f'<text transform="translate(20 {(y0 + y1) // 2}) rotate(-90)" text-anchor="middle">log2 of the spacing at x</text>')
    pts = []
    for e in range(elo, ehi):
        g = log2(ulp(2.0 ** e))
        pts += [f"{fmt(x_of(e))},{fmt(y_of(g))}", f"{fmt(x_of(e + 1))},{fmt(y_of(g))}"]
    body.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{GREEN}" stroke-width="2.5"/>')
    for x, name in [(0.1, "0.1"), (1.0, "1.0"), (1000.0, "1000.0")]:
        px, py = x_of(log2(x)), y_of(log2(ulp(x)))
        body.append(f'<circle cx="{fmt(px)}" cy="{fmt(py)}" r="5" fill="{VIOLET}"/>')
        anchor, dx, dy = ("end", -10, -12) if x > 100 else ("start", 12, 22)
        body.append(f'<text x="{fmt(px + dx)}" y="{fmt(py + dy)}" text-anchor="{anchor}" font-family="{MONO}" font-size="12">'
                    f'ulp({name}) = 2^{int(log2(ulp(x)))}</text>')
    return svg(w, h, "The gap between consecutive double-precision numbers, plotted as log2 of the gap against log2 of x for x "
               "from 2 to the minus 6 up to 2 to the 12. The gap is constant across each interval from one power of two "
               "to the next and doubles at every power of two. The marked points are 0.1 with a gap of 2 to the minus 56, "
               "1.0 with 2 to the minus 52, and 1000.0 with 2 to the minus 43.", body)


def mod2_product() -> str:
    """The message 1011 picks rows of G, the picked rows are summed column by column, and the sums are reduced mod 2."""
    m = [1, 0, 1, 1]
    w, h = 640, 360
    cell, left, top = 40, 220, 40
    body = [f'<g font-family="{MONO}" font-size="11" fill="{TEXT3}" text-anchor="middle">']
    for j in range(7):
        body.append(f'<text x="{left + cell * j + cell // 2}" y="{top - 10}">c{j}</text>')
    body.append("</g>")
    for i, row in enumerate(G):
        y = top + cell * i
        picked = m[i] == 1
        body.append(f'<text x="{left - 150}" y="{y + 25}" text-anchor="end" font-family="{MONO}" font-size="13">m{i} = {m[i]}</text>')
        body.append(f'<text x="{left - 14}" y="{y + 25}" text-anchor="end" font-size="12"'
                    f'{"" if picked else f" fill={chr(34)}{TEXT3}{chr(34)}"}>{"row " + str(i) + " kept" if picked else "row " + str(i) + " dropped"}</text>')
        for j, b in enumerate(row):
            x = left + cell * j
            fill = f'fill="{GREEN}" fill-opacity="0.18"' if (picked and b) else f'fill="{BG2}"'
            opacity = "" if picked else ' stroke-opacity="0.4"'
            body.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" {fill} stroke="currentColor" stroke-width="1.5"{opacity}/>')
            color = GREEN if (picked and b) else ("currentColor" if picked else TEXT3)
            body.append(f'<text x="{x + cell // 2}" y="{y + 26}" text-anchor="middle" font-family="{MONO}" font-size="16" fill="{color}">{b}</text>')
    sums = [sum(G[i][j] for i in range(4) if m[i]) for j in range(7)]
    code = [s % 2 for s in sums]
    ys, yc = top + cell * 4 + 34, top + cell * 4 + 96
    body.append(f'<line x1="{left}" y1="{ys - 14}" x2="{left + 7 * cell}" y2="{ys - 14}" stroke="currentColor" stroke-width="2"/>')
    body.append(f'<text x="{left - 14}" y="{ys + 25}" text-anchor="end" font-size="12">column sums</text>')
    body.append(f'<text x="{left - 14}" y="{yc + 25}" text-anchor="end" font-size="12">sums mod 2</text>')
    for j in range(7):
        x = left + cell * j
        body.append(f'<rect x="{x}" y="{ys}" width="{cell}" height="{cell}" fill="{BG2}" stroke="currentColor" stroke-width="1.5"/>')
        body.append(f'<text x="{x + cell // 2}" y="{ys + 26}" text-anchor="middle" font-family="{MONO}" font-size="16">{sums[j]}</text>')
        fill = f'fill="{GREEN}" fill-opacity="0.18"' if code[j] else f'fill="{BG2}"'
        body.append(f'<rect x="{x}" y="{yc}" width="{cell}" height="{cell}" {fill} stroke="{VIOLET}" stroke-width="2"/>')
        color = GREEN if code[j] else "currentColor"
        body.append(f'<text x="{x + cell // 2}" y="{yc + 26}" text-anchor="middle" font-family="{MONO}" font-size="16"'
                    f' font-weight="600" fill="{color}">{code[j]}</text>')
    body.append(f'<text x="{left + 7 * cell + 16}" y="{yc + 25}" font-family="{MONO}" font-size="13">c = {"".join(map(str, code))}</text>')
    return svg(w, h, "The product of the message 1011 with the four-row generator matrix G of a Hamming code. The message "
               "keeps rows 0, 2 and 3 and drops row 1, the kept rows are added column by column to give 1 0 1 1 2 3 2, "
               f"and reducing each sum mod 2 gives the codeword {''.join(map(str, code))}.", body)


CHARTS = {
    "fig-f2-5-bits-of-100.svg": bits_of_100,
    "fig-f2-6-tenth-in-binary.svg": tenth_in_binary,
    "fig-f2-7-ulp-staircase.svg": ulp_staircase,
    "fig-f2-8-mod2-product.svg": mod2_product,
}


def default_out() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    site = os.path.join(here, "..", "..", "..", "Website-Building-Projects", "dranirbanchakraborty.com",
                        "public", "figures", "error-correction-coding")
    return os.path.normpath(site) if os.path.isdir(site) else "figures"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=default_out(), help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    for name, chart in CHARTS.items():
        path = os.path.join(out, name)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(chart())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
