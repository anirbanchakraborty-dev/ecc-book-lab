"""Write the five computed figures of chapter 32 of Error Correction Coding as SVG files.

    python3 python/ch32_charts.py
    python3 python/ch32_charts.py --out figures

Every value drawn comes from the functions in ch32_rings_and_fields.py: the
multiplication table of the integers modulo 12 and the marks on its units and
zero divisors, the arrows of the two multiplication maps, the counts of units
and zero divisors of each ring of residues from 2 to 40, the addition and
multiplication tables of the three rings with four elements, and the sorting
of the sizes from 2 to 64 into primes, higher prime powers and the rest. A
figure is corrected by fixing those functions and running this script again,
never by editing its output. Inside the hub that holds both this repository
and the book's website, the default output folder is the site's figure folder
for the book. In a clone of this repository alone it is a folder named figures
in the working directory. The SVG follows the book's figure style sheet: real
text in IBM Plex, colors from the page's tokens with paper fallbacks, and no
element that the site's figure check would refuse.
"""

import argparse
import os
from math import hypot

from ch32_rings_and_fields import (F4_NAMES, is_prime, is_prime_power, ring_f4, ring_words, ring_zn,
                                   times_map, unit_and_zero_divisor_counts, units, zero_divisors)

W = 640
GREEN = "var(--green,#00614a)"
RED = "var(--red,#a1211f)"
CYAN = "var(--cyan,#005f7a)"
TEXT3 = "var(--text-3,#5a5448)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
ARROW = ('<defs><marker id="ah" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker>'
         '<marker id="ahg" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="var(--green,#00614a)"/></marker>'
         '<marker id="ahr" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" markerHeight="5"'
         ' orient="auto"><path d="M0 0L10 5L0 10z" fill="var(--red,#a1211f)"/></marker></defs>')


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
         dash: str = "", arrow: bool = False) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    marker = "ahg" if stroke == GREEN else "ahr" if stroke == RED else "ah"
    m = f' marker-end="url(#{marker})"' if arrow else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{m}/>')


def rect(x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, dash: str = "", rx: float = 0) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}"'
            f' stroke="{stroke}" stroke-width="{width}"{d}/>')


def circle(x: float, y: float, r: float, fill: str = BG2, stroke: str = "currentColor", width: float = 1.5,
           dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}/>')


def polygon(points: list, fill: str = "none", stroke: str = "currentColor", width: float = 1.5) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def svg(title: str, h: int, body: list, w: int = W) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def unit_mark(x: float, y: float, r: float) -> str:
    """A unit's label sits in a green circle."""
    return circle(x, y, r, fill="none", stroke=GREEN, width=1.5)


def zero_divisor_mark(x: float, y: float, half: float) -> str:
    """A zero divisor's label sits in a red dashed square."""
    return rect(x - half, y - half, 2 * half, 2 * half, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3)


# The multiplication table of Z_12 ----------------------------------------------------

def z12_multiplication() -> str:
    n, c = 12, 30
    ring = ring_zn(n)
    t, us, zs = ring["mul"], set(units(ring)), set(zero_divisors(ring))
    x0, y0 = (W - (n + 1) * c) / 2, 56
    body = [text(W / 2, 30, "multiplication mod 12", "middle", size=13, weight_="600"),
            text(x0 + c / 2, y0 + c / 2 + 4, "x", "middle", fill=TEXT3, mono=True, size=12)]
    for k in range(n):
        for hx, hy in ((x0 + (k + 1.5) * c, y0 + c / 2), (x0 + c / 2, y0 + (k + 1.5) * c)):
            if k in us:
                body.append(unit_mark(hx, hy, 11))
            elif k in zs:
                body.append(zero_divisor_mark(hx, hy, 10))
            body.append(text(hx, hy + 4, str(k), "middle", fill=TEXT3, mono=True, size=12))
    body.append(line(x0, y0 + c, x0 + (n + 1) * c, y0 + c))
    body.append(line(x0 + c, y0, x0 + c, y0 + (n + 1) * c))
    for a in range(n):
        for b in range(n):
            v = t[a][b]
            cx, cy = x0 + (b + 1) * c, y0 + (a + 1) * c
            one = v == 1
            zero_product = v == 0 and a and b
            if one:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill=BG2, stroke=GREEN, width=1.5, rx=3))
            elif zero_product:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3))
            fill = GREEN if one else RED if zero_product else ""
            body.append(text(cx + c / 2, cy + c / 2 + 4, str(v), "middle", fill=fill, mono=True, size=12,
                             weight_="600" if one or zero_product else ""))
    y = y0 + (n + 1) * c + 34
    body += [unit_mark(69, y - 4, 9),
             text(88, y, "unit: its row holds a 1", size=12),
             zero_divisor_mark(69, y + 22, 8),
             text(88, y + 26, "zero divisor: its row holds a 0 off the zero column", size=12),
             rect(350, y - 13, 18, 18, fill=BG2, stroke=GREEN, width=1.5, rx=3),
             text(376, y, "product 1: an inverse pair", size=12),
             rect(350, y + 13, 18, 18, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3),
             text(376, y + 26, "product 0 of two nonzero factors", size=12),
             text(W / 2, y + 60, "Units 1, 5, 7, 11. Zero divisors 2, 3, 4, 6, 8, 9, 10. The element 0 is neither.",
                  "middle", fill=TEXT3, size=12)]
    return svg("Multiplication table modulo 12 with its units and zero divisors marked", y + 80, body)


# Multiplication by a fixed element as a map ------------------------------------------

def multiplication_maps() -> str:
    h, ytop, ybot, r = 262, 76, 176, 14
    panels = ((3, 7, 24, 76, 39, "times 3 mod 7"), (2, 6, 346, 392, 43, "times 2 mod 6"))
    body = [ARROW]
    for a, n, lab, x0, step, title in panels:
        image = times_map(a, n)
        hit = set(image)

        def x(k: int) -> float:
            return x0 + k * step

        mid = (x(0) + x(n - 1)) / 2
        body.append(text(mid, 30, title, "middle", size=13, weight_="600"))
        body.append(text(lab, ytop + 4, "x", mono=True, size=12, fill=TEXT3))
        body.append(text(lab, ybot + 4, f"{a}x", mono=True, size=12, fill=TEXT3))
        for k in range(n):
            target = image[k]
            collide = image.count(target) > 1 and k != image.index(target)
            color = GREEN if target == 1 else RED if collide else "currentColor"
            dx, dy = x(target) - x(k), ybot - ytop
            d = hypot(dx, dy)
            sx, sy = x(k) + dx * r / d, ytop + dy * r / d
            ex, ey = x(target) - dx * (r + 3) / d, ybot - dy * (r + 3) / d
            body.append(line(sx, sy, ex, ey, color, 1.5 if color != "currentColor" else 1,
                             dash="3 3" if collide else "", arrow=True))
        for k in range(n):
            body.append(circle(x(k), ytop, r))
            body.append(text(x(k), ytop + 5, str(k), "middle", mono=True, size=13))
            missed = k not in hit
            body.append(circle(x(k), ybot, r, fill="none" if missed else BG2, dash="3 3" if missed else ""))
            body.append(text(x(k), ybot + 5, str(k), "middle", mono=True, size=13, fill=TEXT3 if missed else ""))
        verdict = ("every output hit once, so 3 is a unit" if len(hit) == n
                   else "0, 2, 4 hit twice, 1, 3, 5 never")
        body.append(text(mid, 216, verdict, "middle", size=12))
    body.append(text(W / 2, 246, "Green arrow: 3 times 5 is 1 mod 7. Dashed red arrows: a second input landing on an"
                     " output already hit.", "middle", fill=TEXT3, size=11))
    return svg("Multiplication by 3 modulo 7 is one to one, multiplication by 2 modulo 6 is not", h, body)


# How many units and zero divisors each Z_n has ---------------------------------------

def unit_counts() -> str:
    lo, hi, top = 2, 40, 40
    px0, px1, py0, py1 = 70, 610, 372, 52

    def px(n: float) -> float:
        return px0 + (n - lo) * (px1 - px0) / (hi - lo)

    def py(v: float) -> float:
        return py0 - v * (py0 - py1) / top

    body = [line(px0, py0, px1 + 6, py0), line(px0, py0, px0, py1 - 6)]
    for v in range(0, top + 1, 10):
        body.append(line(px0 - 4, py(v), px0, py(v)))
        body.append(text(px0 - 8, py(v) + 4, str(v), "end", fill=TEXT3, mono=True, size=11))
    for n in (2, 5, 10, 15, 20, 25, 30, 35, 40):
        body.append(line(px(n), py0, px(n), py0 + 4))
        body.append(text(px(n), py0 + 18, str(n), "middle", fill=TEXT3, mono=True, size=11))
    body.append(text((px0 + px1) / 2, py0 + 38, "n, the modulus", "middle", size=12))
    body.append(text(px0 - 36, py1 - 16, "count", size=12))
    body.append(line(px(lo), py(lo - 1), px(hi), py(hi - 1), stroke=TEXT3, width=1, dash="3 3"))
    body.append(text(px(33), py(32) - 4, "n - 1, every nonzero element", "end", fill=TEXT3, size=11))
    for n in range(lo, hi + 1):
        u, z = unit_and_zero_divisor_counts(n)
        x = px(n)
        if is_prime(n):
            y = py(u)
            body.append(polygon([(x, y - 5.5), (x + 5.5, y), (x, y + 5.5), (x - 5.5, y)], fill=BG2, stroke=GREEN))
        else:
            body.append(circle(x, py(u), 4, fill=BG2, stroke=GREEN))
        if z:
            body.append(rect(x - 3.5, py(z) - 3.5, 7, 7, fill=RED, stroke=RED, width=1))
    ly = 22
    body += [circle(82, ly - 4, 4, fill=BG2, stroke=GREEN), text(92, ly, "units, n composite", size=12),
             polygon([(240, ly - 9.5), (245.5, ly - 4), (240, ly + 1.5), (234.5, ly - 4)], fill=BG2, stroke=GREEN),
             text(252, ly, "units, n prime", size=12),
             rect(376.5, ly - 7.5, 7, 7, fill=RED, stroke=RED, width=1), text(390, ly, "zero divisors", size=12),
             text(W / 2, py0 + 62, "Units and zero divisors always add up to n - 1. Only a prime modulus has no zero"
                  " divisor.", "middle", fill=TEXT3, size=12)]
    return svg("The number of units and of zero divisors in the integers modulo n, for n from 2 to 40",
               py0 + 80, body)


# The three rings with four elements ---------------------------------------------------

def four_tables(names: list, t: list, x0: float, y0: float, c: float, mark_ones: bool, unity: int = 1) -> list:
    n = len(names)
    body = [text(x0 + c / 2, y0 + c / 2 + 4, "x" if mark_ones else "+", "middle", fill=TEXT3, mono=True, size=12)]
    for k in range(n):
        body.append(text(x0 + (k + 1.5) * c, y0 + c / 2 + 4, names[k], "middle", fill=TEXT3, mono=True, size=12))
        body.append(text(x0 + c / 2, y0 + (k + 1.5) * c + 4, names[k], "middle", fill=TEXT3, mono=True, size=12))
    body.append(line(x0, y0 + c, x0 + (n + 1) * c, y0 + c))
    body.append(line(x0 + c, y0, x0 + c, y0 + (n + 1) * c))
    for a in range(n):
        for b in range(n):
            v = t[a][b]
            cx, cy = x0 + (b + 1) * c, y0 + (a + 1) * c
            one = mark_ones and v == unity
            zero_product = mark_ones and v == 0 and a and b
            zero_sum = not mark_ones and v == 0
            if one:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill=BG2, stroke=GREEN, width=1.5, rx=3))
            elif zero_product:
                body.append(rect(cx + 3, cy + 3, c - 6, c - 6, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3))
            elif zero_sum:
                body.append(circle(cx + c / 2, cy + c / 2, c / 2 - 4, fill="none", stroke=CYAN, width=1.5))
            fill = GREEN if one else RED if zero_product else CYAN if zero_sum else ""
            body.append(text(cx + c / 2, cy + c / 2 + 4, names[v], "middle", fill=fill, mono=True, size=12,
                             weight_="600" if fill else ""))
    return body


def four_element_rings() -> str:
    c = 34
    f4 = ring_f4()
    rings = (("integers mod 4", ["0", "1", "2", "3"], ring_zn(4), "2 times 2 is 0"),
             ("2-bit words, XOR and AND", ["00", "01", "10", "11"], ring_words(2), "01 AND 10 is 00"),
             ("the four-element field", list(F4_NAMES), f4, "no zero product"))
    xs = (28, 232, 436)
    ym, ya = 74, 274
    body = []
    for (title, names, ring, verdict), x0 in zip(rings, xs):
        mid = x0 + 2.5 * c
        body.append(text(mid, 30, title, "middle", size=13, weight_="600"))
        body.append(text(mid, 50, verdict, "middle", fill=TEXT3, size=12))
        body += four_tables(names, ring["mul"], x0, ym, c, True, ring["one"])
        body += four_tables(names, ring["add"], x0, ya, c, False)
    y = ya + 5 * c + 32
    body += [rect(40, y - 13, 18, 18, fill=BG2, stroke=GREEN, width=1.5, rx=3),
             text(66, y, "product equal to the ring's 1", size=12),
             rect(330, y - 13, 18, 18, fill="none", stroke=RED, width=1.5, dash="3 3", rx=3),
             text(356, y, "product 0 of nonzero factors", size=12),
             circle(49, y + 22, 9, fill="none", stroke=CYAN, width=1.5),
             text(66, y + 26, "sum 0: the row and the column are negatives of each other", size=12),
             text(W / 2, y + 58, "In the two right-hand rings every element is its own negative, so 1 + 1 = 0.",
                  "middle", fill=TEXT3, size=12)]
    return svg("Multiplication and addition tables of the three rings with four elements", y + 76, body)


# Which sizes a finite field can have --------------------------------------------------

def prime_and_exponent(n: int) -> tuple[int, int]:
    p = next(d for d in range(2, n + 1) if n % d == 0)
    m, r = 0, n
    while r % p == 0:
        r //= p
        m += 1
    return p, m


def field_sizes() -> str:
    cols, cw, ch = 9, 60, 46
    x0, y0 = (W - cols * cw) / 2, 30
    body = []
    for i, n in enumerate(range(2, 65)):
        row, col = divmod(i, cols)
        cx, cy = x0 + col * cw + cw / 2, y0 + row * ch + ch / 2
        if is_prime(n):
            body.append(rect(cx - 24, cy - 18, 48, 36, fill=BG2, stroke=GREEN, width=1.5, rx=3))
            body.append(text(cx, cy + 5, str(n), "middle", fill=GREEN, mono=True, size=14, weight_="600"))
        elif is_prime_power(n):
            p, m = prime_and_exponent(n)
            body.append(circle(cx, cy, 19, fill=BG2, stroke=CYAN))
            body.append(text(cx, cy, str(n), "middle", fill=CYAN, mono=True, size=13, weight_="600"))
            body.append(f'<text x="{fmt(cx)}" y="{fmt(cy + 13)}" text-anchor="middle" fill="{CYAN}"'
                        f' font-family="{MONO}" font-size="9">{p}<tspan dy="-3" font-size="8">{m}</tspan></text>')
        else:
            body.append(text(cx, cy + 5, str(n), "middle", fill=TEXT3, mono=True, size=13))
            body.append(line(cx - 14, cy + 9, cx + 14, cy - 9, stroke=TEXT3, width=1))
    y = y0 + 7 * ch + 34
    body += [rect(40, y - 14, 24, 18, fill=BG2, stroke=GREEN, width=1.5, rx=3),
             text(72, y, "prime: a field exists", size=12),
             circle(234, y - 5, 10, fill=BG2, stroke=CYAN),
             text(252, y, "higher prime power: a field exists", size=12),
             text(470, y, "6", "middle", fill=TEXT3, mono=True, size=13),
             line(462, y + 2, 478, y - 12, stroke=TEXT3, width=1),
             text(488, y, "struck out: no field", size=12),
             text(W / 2, y + 30, "From 2 to 64 there are 18 primes, 9 higher prime powers and 36 sizes no field can have.",
                  "middle", fill=TEXT3, size=12)]
    return svg("The sizes from 2 to 64, sorted into primes, higher prime powers and sizes no field can have",
               y + 48, body)


CHARTS = {
    "fig-32-1-z12-multiplication.svg": z12_multiplication,
    "fig-32-2-multiplication-maps.svg": multiplication_maps,
    "fig-32-3-unit-counts.svg": unit_counts,
    "fig-32-4-four-element-rings.svg": four_element_rings,
    "fig-32-8-field-sizes.svg": field_sizes,
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
