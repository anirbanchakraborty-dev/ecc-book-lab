"""Write the three computed figures of chapter 33 of Error Correction Coding as SVG files.

    python3 python/ch33_charts.py
    python3 python/ch33_charts.py --out figures

Every value drawn comes from the functions in ch33_polynomials_over_a_field.py:
the rows of the extended Euclidean table for x^5 + x^4 + 1 and x^4 + x^2 + 1
over GF(2), the sieve that strikes out the multiples of x, x + 1 and
x^2 + x + 1 among the polynomials of degree 1 to 5, and the number of
irreducible polynomials over GF(2) of each degree from 1 to 10. A figure is
corrected by fixing those functions and running this script again, never by
editing its output. Inside the hub that holds both this repository and the
book's website, the default output folder is the site's figure folder for the
book. In a clone of this repository alone it is a folder named figures in the
working directory. The SVG follows the book's figure style sheet: real text in
IBM Plex, colors from the page's tokens with paper fallbacks, and no element
that the site's figure check would refuse.
"""

import argparse
import os
from math import log2

from ch33_polynomials_over_a_field import (bits, clmul, deg2, divides2, irreducible_counts, sieve2, weight,
                                           xgcd2)

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


def polyline(points: list, stroke: str = "currentColor", width: float = 1.5, dash: str = "") -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def rect(x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = "currentColor",
         width: float = 1, rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}"'
            f' stroke="{stroke}" stroke-width="{width}"/>')


def circle(x: float, y: float, r: float, fill: str = BG2, stroke: str = "currentColor", width: float = 1.5) -> str:
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}" stroke="{stroke}"'
            f' stroke-width="{width}"/>')


def polygon(points: list, fill: str = "none", stroke: str = "currentColor", width: float = 1.5) -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>'


def svg(title: str, h: int, body: list, w: int = W) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


# Polynomials written with raised exponents ------------------------------------------

def poly_segments(f: int) -> list:
    """The polynomial f over GF(2) as (text, raised) pieces, highest power first."""
    if f == 0:
        return [("0", False)]
    out = []
    for k in range(deg2(f), -1, -1):
        if f >> k & 1:
            if out:
                out.append((" + ", False))
            if k == 0:
                out.append(("1", False))
            elif k == 1:
                out.append(("x", False))
            else:
                out += [("x", False), (str(k), True)]
    return out


def factor_segments(f: int) -> list:
    """f as a factor in a product: in parentheses when it has more than one term."""
    inner = poly_segments(f)
    return [("(", False), *inner, (")", False)] if weight(f) > 1 else inner


def rich(segments: list) -> str:
    """The inside of a text element, with raised pieces in a smaller size."""
    merged = []
    for s, up in segments:
        if merged and merged[-1][1] == up:
            merged[-1] = (merged[-1][0] + s, up)
        else:
            merged.append((s, up))
    out, raised = [], False
    for s, up in merged:
        if up:
            out.append(f'<tspan dy="-5" font-size="9">{s}</tspan>')
            raised = True
        elif raised:
            gap = ' dx="3"' if s.startswith(" ") else ""
            out.append(f'<tspan dy="5"{gap}>{s.lstrip()}</tspan>')
            raised = False
        else:
            out.append(s)
    return "".join(out)


def rich_width(segments: list, size: int = 12) -> float:
    return sum(len(s) * 0.6 * (9 if up else size) for s, up in segments)


# The extended Euclidean table ---------------------------------------------------------

def euclid_table() -> str:
    a, b = 0b110001, 0b10101
    rows = xgcd2(a, b)
    g = len(rows) - 2
    r, _, s, t = rows[g]
    if clmul(s, a) ^ clmul(t, b) != r:
        raise SystemExit("the gcd row does not satisfy s f + t g = r")
    cols = [30, 74, 264, 354, 482, 610]
    top0, rh = 92, 40
    body = [text(W / 2, 26, "The extended Euclidean algorithm over GF(2)", "middle", size=13, weight_="600"),
            text(W / 2, 48, rich([("f = ", False), *poly_segments(a), (f" ({bits(a)}) and g = ", False),
                                  *poly_segments(b), (f" ({bits(b)})", False)]), "middle", fill=TEXT3, size=12)]
    for k, head in enumerate(["i", "r", "q", "s", "t"]):
        x = (cols[0] + cols[1]) / 2 if k == 0 else cols[k] + 12
        body.append(text(x, top0 - 12, head, "middle" if k == 0 else "start", size=12, weight_="600"))
    body.append(line(cols[0], top0, cols[-1], top0))
    for i, (ri, qi, si, ti) in enumerate(rows):
        top = top0 + i * rh
        fill = TEXT3 if ri == 0 else ""
        if i == g:
            body.append(rect(cols[0] + 2, top + 3, cols[-1] - cols[0] - 4, rh - 4, stroke=GREEN, width=1.5, rx=3))
        body.append(text((cols[0] + cols[1]) / 2, top + 20, str(i), "middle", fill=fill, mono=True, size=12))
        for k, v in ((1, ri), (2, qi), (3, si), (4, ti)):
            if v < 0:
                continue
            body.append(text(cols[k] + 12, top + 20, rich(poly_segments(v)), fill=fill, size=12))
            body.append(text(cols[k] + 12, top + 33, bits(v), fill=TEXT3, mono=True, size=10))
    bottom = top0 + len(rows) * rh
    body.append(line(cols[0], bottom, cols[-1], bottom))
    bezout = [*factor_segments(s), *factor_segments(a), (" + ", False), *factor_segments(t),
              *factor_segments(b), (" = ", False), *poly_segments(r)]
    body.append(text(W / 2, bottom + 30, rich(bezout), "middle", size=13))
    body.append(text(W / 2, bottom + 54, f"Row {g} holds the last nonzero remainder, the gcd, with the s and t"
                     " that build it from f and g.", "middle", fill=TEXT3, size=12))
    return svg("The extended Euclidean table for x^5 + x^4 + 1 and x^4 + x^2 + 1 over GF(2), with the gcd row"
               " boxed", bottom + 74, body)


# The sieve for irreducible polynomials -----------------------------------------------

def sieve_class(f: int) -> str:
    """Which of x, x + 1 and x^2 + x + 1 strikes f out first, or kept."""
    if deg2(f) == 1:
        return "kept"
    if f & 1 == 0:
        return "x"
    if weight(f) % 2 == 0:
        return "x+1"
    if f != 0b111 and divides2(0b111, f):
        return "x^2+x+1"
    return "kept"


def triangle(cx: float, cy: float) -> str:
    return polygon([(cx - 4.5, cy + 3.5), (cx + 4.5, cy + 3.5), (cx, cy - 4)], fill=ORANGE, stroke=ORANGE, width=1)


def diamond(cx: float, cy: float) -> str:
    return polygon([(cx, cy - 4.5), (cx + 4.5, cy), (cx, cy + 4.5), (cx - 4.5, cy)], fill=RED, stroke=RED, width=1)


def irreducible_sieve() -> str:
    top_degree = 5
    kept_by_sieve = {f for f in range(2, 1 << (top_degree + 1)) if sieve_class(f) == "kept"}
    if kept_by_sieve != set(sieve2(top_degree)):
        raise SystemExit("the three-factor sieve and sieve2 disagree")
    counts = irreducible_counts(top_degree)
    x0, cw, rh, top = 108, 64, 38, 64
    body = [text(W / 2, 24, "A sieve for the irreducible polynomials over GF(2) of degree 1 to 5", "middle",
                 size=13, weight_="600"),
            text(W / 2, 44, "Past degree 1 a constant term 0 means a factor x, so only constant term 1 is shown.",
                 "middle", fill=TEXT3, size=12)]
    for d in range(1, top_degree + 1):
        polys = [f for f in range(1 << d, 1 << (d + 1)) if d == 1 or f & 1]
        kept = sum(1 for f in polys if sieve_class(f) == "kept")
        if kept != counts[d - 1]:
            raise SystemExit(f"degree {d}: the sieve keeps {kept}, the count is {counts[d - 1]}")
        body.append(text(20, top + 17, f"degree {d}", size=12))
        body.append(text(20, top + 31, f"{kept} kept", fill=GREEN, size=11))
        for start in range(0, len(polys), 8):
            for j, f in enumerate(polys[start:start + 8]):
                cx = x0 + j * cw + cw / 2
                kind = sieve_class(f)
                if kind == "kept":
                    bw = len(bits(f)) * 7.2 + 10
                    body.append(rect(cx - bw / 2, top + 4, bw, 18, fill=BG2, stroke=GREEN, width=1.5, rx=3))
                    body.append(text(cx, top + 17, bits(f), "middle", mono=True, size=12, weight_="600"))
                else:
                    body.append(text(cx, top + 17, bits(f), "middle", fill=TEXT3, mono=True, size=12))
                    body.append(triangle(cx, top + 29) if kind == "x+1" else diamond(cx, top + 29))
            top += rh
    ly = top + 26
    body += [triangle(46, ly - 4), text(56, ly, "multiple of x + 1 (even weight)", size=12),
             diamond(309, ly - 4),
             text(319, ly, rich([("multiple of x", False), ("2", True), (" + x + 1", False)]), size=12),
             rect(488, ly - 12, 18, 14, fill=BG2, stroke=GREEN, width=1.5, rx=3),
             text(514, ly, "irreducible", size=12)]
    listed = ", ".join(str(c) for c in counts[:-1]) + f" and {counts[-1]}"
    body.append(text(W / 2, ly + 28, f"The survivors number {listed} for degrees 1 to {top_degree}.", "middle",
                     fill=TEXT3, size=12))
    return svg("A sieve that strikes out the multiples of x, x + 1 and x^2 + x + 1 among the polynomials over GF(2)"
               " of degree 1 to 5", ly + 46, body)


# How many polynomials of each degree are irreducible --------------------------------

def irreducible_count_chart() -> str:
    top_m = 10
    counts = irreducible_counts(top_m)
    lo, hi = -1, 10
    px0, px1, py0, py1 = 80, 600, 380, 70

    def px(m: float) -> float:
        return px0 + (m - 0.5) * (px1 - px0) / top_m

    def py(v: float) -> float:
        return py0 - (log2(v) - lo) * (py0 - py1) / (hi - lo)

    body = [line(px0, py0, px1 + 10, py0), line(px0, py0, px0, py1 - 10)]
    for v in (1, 4, 16, 64, 256, 1024):
        body.append(line(px0 + 1, py(v), px1 + 10, py(v), stroke=BG2, width=1))
        body.append(line(px0 - 4, py(v), px0, py(v)))
        body.append(text(px0 - 8, py(v) + 4, str(v), "end", fill=TEXT3, mono=True, size=11))
    for m in range(1, top_m + 1):
        body.append(line(px(m), py0, px(m), py0 + 4))
        body.append(text(px(m), py0 + 18, str(m), "middle", fill=TEXT3, mono=True, size=11))
    body.append(text((px0 + px1) / 2, py0 + 38, "m, the degree", "middle", size=12))
    body.append(text(px0 - 44, py1 - 18, "how many (doubling scale)", size=12))
    body.append(polyline([(px(m), py(2 ** m)) for m in range(1, top_m + 1)], stroke=TEXT3, width=1.5))
    body.append(polyline([(px(m), py(2 ** m / m)) for m in range(1, top_m + 1)], stroke=CYAN, width=1.5,
                         dash="3 3"))
    for m, c in enumerate(counts, start=1):
        if m * c > 2 ** m:
            raise SystemExit(f"degree {m}: m I(m) exceeds 2^m")
        body.append(circle(px(m), py(c), 4, fill=BG2, stroke=GREEN))
        body.append(text(px(m), py(c) + 17, str(c), "middle", fill=GREEN, mono=True, size=11, weight_="600"))
    ratio = top_m * counts[-1] / 2 ** top_m
    body += [text(px0 + 20, 100, f"at degree {top_m}, {counts[-1]} of the {2 ** top_m} are irreducible,",
                  fill=TEXT3, size=12),
             text(px0 + 20, 118, f"and {top_m} × {counts[-1]} / {2 ** top_m} = {ratio:.3f}", fill=TEXT3,
                  size=12)]
    ly = 26
    body += [line(96, ly - 4, 120, ly - 4, stroke=TEXT3, width=1.5),
             text(128, ly, rich([("all 2", False), ("m", True), (" of degree m", False)]), size=12),
             circle(270, ly - 4, 4, fill=BG2, stroke=GREEN),
             text(282, ly, "irreducible ones, I(m)", size=12),
             line(444, ly - 4, 468, ly - 4, stroke=CYAN, width=1.5, dash="3 3"),
             text(476, ly, rich([("2", False), ("m", True), (" / m", False)]), size=12)]
    if abs(ratio - 0.967) > 0.0005:
        raise SystemExit(f"the ratio at degree {top_m} is {ratio:.4f}, not 0.967")
    return svg("The number of irreducible polynomials over GF(2) of each degree from 1 to 10 on a doubling scale,"
               " against 2^m and 2^m / m", py0 + 56, body)


CHARTS = {
    "fig-33-4-extended-euclid-table.svg": euclid_table,
    "fig-33-5-irreducible-sieve.svg": irreducible_sieve,
    "fig-33-6-irreducible-counts.svg": irreducible_count_chart,
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
