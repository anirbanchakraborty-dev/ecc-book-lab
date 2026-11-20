"""Write the two computed figures of chapter 34 of Error Correction Coding as SVG files.

    python3 python/ch34_charts.py
    python3 python/ch34_charts.py --out figures

Every value drawn comes from the functions in ch34_ideals_and_quotient_rings.py:
the multiplication tables, units, zero divisors and nilpotent elements of the
four quotients of GF(2)[x] by a polynomial of degree 2, with the matchings that
identify two of them with rings from Rings and Fields, and the number of units
in GF(2)[x]/(p(x)) for every p of degree 1 to 6. A figure is corrected by fixing
those functions and running this script again, never by editing its output.
Inside the hub that holds both this repository and the book's website, the
default output folder is the site's figure folder for the book. In a clone of
this repository alone it is a folder named figures in the working directory.
The SVG follows the book's figure style sheet: real text in IBM Plex, colors
from the page's tokens with paper fallbacks, and no element that the site's
figure check would refuse.
"""

import argparse
import os

from ch34_ideals_and_quotient_rings import (FOUR_MODULI, bitwise_ring, f4_ring, is_irreducible2, matchings,
                                            nilpotent_list, ring_mod2, unit_counts, unit_list,
                                            zero_divisor_list)

W = 640
GREEN = "var(--green,#00614a)"
RED = "var(--red,#a1211f)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


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
         dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}/>')


def rect(x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}" fill="{fill}"'
            f' stroke="{stroke}" stroke-width="{width}"{d}/>')


def circle(x: float, y: float, r: float, fill: str = BG2, stroke: str = "currentColor", width: float = 1.5) -> str:
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}" stroke="{stroke}"'
            f' stroke-width="{width}"/>')


def polygon(points: list, fill: str = "none", stroke: str = "currentColor", width: float = 1.5) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def svg(title: str, h: int, body: list, w: int = W) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def rich(pieces: list) -> str:
    """Text pieces, each either plain or raised as an exponent, as tspans."""
    out, raised, lowered = [], False, False
    for s, up in pieces:
        if up == "sub":
            out.append(f'<tspan dy="3" font-size="9">{s}</tspan>')
            lowered = True
            continue
        if lowered:
            out.append(f'<tspan dy="-3" dx="3">{s.lstrip()}</tspan>')
            lowered = False
            continue
        if up and not raised:
            out.append(f'<tspan dy="-5" font-size="9">{s}</tspan>')
            raised = True
        elif not up and raised:
            out.append(f'<tspan dy="5" dx="3">{s}</tspan>')
            raised = False
        else:
            out.append(s)
    return "".join(out)


def triangle(cx: float, cy: float, fill: str = ORANGE) -> str:
    return polygon([(cx - 4.5, cy + 3.5), (cx + 4.5, cy + 3.5), (cx, cy - 4)], fill=fill, stroke=fill, width=1)


def diamond(cx: float, cy: float, fill: str = RED, half: float = 4.5) -> str:
    return polygon([(cx, cy - half), (cx + half, cy), (cx, cy + half), (cx - half, cy)], fill=fill, stroke=fill,
                   width=1)


# The four quotients with four elements ------------------------------------------

MODULUS_PIECES = {
    0b111: [("x", False), ("2", True), (" + x + 1", False)],
    0b110: [("x", False), ("2", True), (" + x", False)],
    0b101: [("x", False), ("2", True), (" + 1", False)],
    0b100: [("x", False), ("2", True)],
}


def four_quotients() -> str:
    names = ["00", "01", "10", "11"]
    f4 = matchings(ring_mod2(0b111), f4_ring())
    bitwise = matchings(ring_mod2(0b110), bitwise_ring())
    swap = matchings(ring_mod2(0b101), ring_mod2(0b100))
    if len(f4) != 2 or len(bitwise) != 2 or len(swap) != 1:
        raise SystemExit("the matchings differ from the ones the chapter states")
    notes = {
        0b111: [[("a field: every nonzero word is a unit, and", False)],
                [("with a = 10 and b = 11 it is the F", False), ("4", "sub"), (" of Rings and Fields", False)]],
        0b110: [[("x", False), ("2", True), (" = x and (x + 1)", False), ("2", True), (" = x + 1", False)],
                [("the 2-bit words, once 01 is sent to 11", False)]],
        0b101: [[("(x + 1)", False), ("2", True), (" = 0, so 11 is nilpotent", False)],
                [("substituting x + 1 for x gives the ring at right", False)]],
        0b100: [[("x", False), ("2", True), (" = 0, so 10 is nilpotent", False)],
                [("the same ring as the one at left, relabeled", False)]],
    }
    body = [text(W / 2, 24, "The four rings GF(2)[x]/(p(x)) with p of degree 2", "middle", size=13, weight_="600"),
            text(W / 2, 44, "Each word b1 b0 stands for b1 x + b0. Addition is XOR in all four, so only products"
                 " are shown.", "middle", fill=TEXT3, size=12)]
    cw, rh, hw = 46, 22, 42
    for idx, p in enumerate(FOUR_MODULI):
        add, mul, one = ring_mod2(p)
        zd = set(zero_divisor_list(mul))
        nil = set(nilpotent_list(mul))
        units = unit_list(mul, one)
        x0 = 28 + (idx % 2) * 312
        y0 = 70 + (idx // 2) * 214
        body.append(text(x0, y0 + 12, rich([("modulo ", False)] + MODULUS_PIECES[p]), size=12, weight_="600"))
        ty = y0 + 24
        body.append(rect(x0 + hw, ty, 4 * cw, rh, fill=BG2, stroke="none"))
        body.append(rect(x0, ty + rh, hw, 4 * rh, fill=BG2, stroke="none"))
        body.append(text(x0 + hw / 2, ty + 15, "times", "middle", fill=TEXT3, size=11))
        for j in range(4):
            cx = x0 + hw + j * cw + cw / 2
            body.append(text(cx - (5 if j in zd else 0), ty + 15, names[j], "middle", mono=True, size=12,
                             weight_="600"))
            if j in nil:
                body.append(diamond(cx + 15, ty + 11, half=4))
            elif j in zd:
                body.append(triangle(cx + 15, ty + 11))
        for i in range(4):
            ry = ty + rh * (i + 1)
            body.append(text(x0 + hw / 2 - (5 if i in zd else 0), ry + 15, names[i], "middle", mono=True, size=12,
                             weight_="600"))
            if i in nil:
                body.append(diamond(x0 + hw / 2 + 13, ry + 11, half=4))
            elif i in zd:
                body.append(triangle(x0 + hw / 2 + 13, ry + 11))
            for j in range(4):
                cx = x0 + hw + j * cw + cw / 2
                v = mul[i][j]
                hit = v == 0 and i != 0 and j != 0
                if hit:
                    body.append(rect(cx - 17, ry + 3, 34, rh - 6, fill="none", stroke=RED, width=1.2, dash="3 2"))
                body.append(text(cx, ry + 15, names[v], "middle", mono=True, size=12, fill=RED if hit else ""))
        body.append(rect(x0, ty, hw + 4 * cw, 5 * rh, fill="none", stroke=TEXT3, width=1))
        body.append(line(x0 + hw, ty, x0 + hw, ty + 5 * rh, stroke=TEXT3))
        body.append(line(x0, ty + rh, x0 + hw + 4 * cw, ty + rh, stroke=TEXT3))
        unit_words = " ".join(names[u] for u in units)
        body.append(text(x0, ty + 5 * rh + 18, f"units: {unit_words}", size=11))
        for k, pieces in enumerate(notes[p]):
            body.append(text(x0, ty + 5 * rh + 34 + 15 * k, rich(pieces), fill=TEXT3, size=11))
    ly = 70 + 2 * 214 + 6
    body += [triangle(46, ly - 4), text(56, ly, "zero divisor", size=12),
             diamond(160, ly - 4, half=4), text(170, ly, "nilpotent zero divisor (a power is 00)", size=12),
             rect(424, ly - 13, 30, 16, fill="none", stroke=RED, width=1.2, dash="3 2"),
             text(460, ly, "00 from nonzero words", size=12)]
    return svg("The multiplication tables of the four quotients of GF(2)[x] by x^2 + x + 1, x^2 + x, x^2 + 1 and x^2,"
               " with zero divisors and nilpotent elements marked", ly + 20, body)


# Units in every quotient by a polynomial of degree 1 to 6 ---------------------------

def unit_count_chart() -> str:
    top = 6
    counts = unit_counts(top)
    px0, px1, py0, scale = 84, 624, 450, 6.0
    colw = (px1 - px0) / top

    def py(v: float) -> float:
        return py0 - v * scale

    body = [line(px0, py0, px1, py0), line(px0, py0, px0, py(64) - 10)]
    for v in (0, 16, 32, 48, 64):
        if v:
            body.append(line(px0 + 1, py(v), px1, py(v), stroke=BG2, width=1))
        body.append(line(px0 - 4, py(v), px0, py(v)))
        body.append(text(px0 - 8, py(v) + 4, str(v), "end", fill=TEXT3, mono=True, size=11))
    body.append(text(px0 - 52, py(64) - 22, "units in GF(2)[x]/(p(x))", size=12))
    for m in range(1, top + 1):
        cx = px0 + (m - 0.5) * colw
        body.append(line(cx, py0, cx, py0 + 4))
        body.append(text(cx, py0 + 18, str(m), "middle", fill=TEXT3, mono=True, size=11))
        full = (1 << m) - 1
        body.append(line(cx - colw / 2 + 6, py(full), cx + colw / 2 - 6, py(full), stroke=GREEN, width=1, dash="4 3"))
        rows: dict[int, list[bool]] = {}
        for p, u in zip(range(1 << m, 1 << (m + 1)), counts[m - 1]):
            rows.setdefault(u, []).append(is_irreducible2(p))
        irreducible = 0
        for u, flags in sorted(rows.items()):
            if any(flags) and not all(flags):
                raise SystemExit(f"degree {m}: irreducible and reducible moduli share the count {u}")
            k = len(flags)
            step = 6.4
            for i, irr in enumerate(flags):
                x = cx + (i - (k - 1) / 2) * step
                if irr:
                    body.append(diamond(x, py(u), fill=GREEN, half=3.4))
                    irreducible += 1
                else:
                    body.append(circle(x, py(u), 2.4, fill=BG2, stroke=CYAN, width=1.2))
        body.append(text(cx, py(full) - 12, f"{irreducible} of {1 << m}", "middle", fill=GREEN, size=11,
                         weight_="600"))
    body.append(text((px0 + px1) / 2, py0 + 38, "m, the degree of the modulus p(x)", "middle", size=12))
    ly = 26
    body += [diamond(110, ly - 4, fill=GREEN, half=3.4), text(120, ly, "irreducible modulus", size=12),
             circle(262, ly - 4, 2.4, fill=BG2, stroke=CYAN, width=1.2), text(272, ly, "reducible modulus", size=12),
             line(404, ly - 4, 428, ly - 4, stroke=GREEN, width=1, dash="4 3"),
             text(436, ly, rich([("2", False), ("m", True), (" - 1, every nonzero element", False)]), size=12)]
    return svg("The number of units in GF(2)[x]/(p(x)) for every modulus p of degree 1 to 6, with the irreducible"
               " moduli reaching 2^m - 1", py0 + 56, body)


CHARTS = {
    "fig-34-4-four-element-quotients.svg": four_quotients,
    "fig-34-7-unit-counts.svg": unit_count_chart,
}


def default_out() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    site = os.path.join(here, "..", "..", "..", "Website-Building-Projects", "dranirbanchakraborty.com",
                        "public", "figures", "error-correction-coding")
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
