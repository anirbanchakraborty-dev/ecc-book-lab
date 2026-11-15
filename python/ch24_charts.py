"""Write the two computed figures of chapter 24 of Error Correction Coding as SVG files.

    python3 python/ch24_charts.py --out figures

Every drawn value comes from the codes built in
ch24_modifying_and_combining_codes.py, so a figure is corrected by fixing that
module and running this script again, never by editing its output. The SVG
follows the book's figure style sheet: real text in IBM Plex, colors from the
page's tokens with paper fallbacks, and no element that the site's figure check
would refuse.
"""

import argparse
import os

from ch24_modifying_and_combining_codes import (
    CHART_CODES, DIRECT8, UUV8, dimension, dmin, weight)

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
ORANGE = "var(--orange,#8a4b00)"
BG2 = "var(--bg-2,#f3eee3)"
MONO = "IBM Plex Mono, monospace"
DISTANCE_COLOR = {2: ORANGE, 3: CYAN, 4: GREEN}


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


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, width: float = 1.5, rx: float = 0) -> str:
    r = f' rx="{fmt(rx)}"' if rx else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}"{r} fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')


def svg(title: str, body: list, w: int = W, h_: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h_}" width="{w}" height="{h_}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def legend_box(x: float, y: float, fill: str, stroke: str, label: str) -> list:
    return [rect(x, y - 11, 14, 14, fill, stroke, 1.5, 3), text(x + 22, y, label, size=12)]


def counts_by_weight(C: list, n: int) -> list:
    counts = [0] * (n + 1)
    for c in C:
        counts[weight(c)] += 1
    return counts


def direct_sum_versus_u_u_plus_v() -> str:
    """Codewords at each weight in the direct sum and the (u | u + v) code of (4,3,2) and (4,1,4)."""
    n = len(DIRECT8[0])
    series = [(counts_by_weight(DIRECT8, n), ORANGE, -19, f"direct sum, d = {dmin(DIRECT8)}"),
              (counts_by_weight(UUV8, n), GREEN, 1, f"(u | u + v), d = {dmin(UUV8)}")]
    slot = (X1 - X0) / (n + 1)
    ay = LinAxis(0, 16, Y1, Y0)
    body = [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, "codeword weight", "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">number of codewords</text>',
    ]
    for v in (0, 4, 8, 12, 16):
        py = ay(v)
        body.append(line(X0 - 4, py, X0, py))
        body.append(text(X0 - 8, py + 4, str(v), "end", mono=True, size=11))
    for w in range(n + 1):
        cx = X0 + slot * (w + 0.5)
        body.append(line(cx, Y1, cx, Y1 + 4))
        body.append(text(cx, Y1 + 20, str(w), "middle", mono=True, size=11))
        for counts, color, dx, _ in series:
            if counts[w]:
                top = ay(counts[w])
                body.append(rect(cx + dx, top, 18, Y1 - top, color, color, 1, 2))
                body.append(text(cx + dx + 9, top - 6, str(counts[w]), "middle", color, mono=True, size=10))
    for i, (_, color, _, label) in enumerate(series):
        body += legend_box(440, 46 + 24 * i, color, color, label)
    return svg("Weight distributions of the direct sum and of the (u | u + v) code built from the "
               "(4,3,2) and (4,1,4) codes", body)


def rate_and_distance() -> str:
    """Rate and minimum distance of each code built in the chapter, sorted by rate."""
    rows = []
    for name, C in CHART_CODES:
        n, k, d = len(C[0]), dimension(C), dmin(C)
        rows.append((k / n, n, k, d, name))
    rows.sort(key=lambda r: -r[0])
    x0, x1, y0, pitch = 200, 580, 20, 32
    y1 = y0 + pitch * len(rows)
    h_ = y1 + 52
    ax = LinAxis(0, 1, x0, x1)
    body = [f'<polyline points="{x0},{y0} {x0},{y1} {x1},{y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
            text((x0 + x1) / 2, h_ - 12, "rate k/n", "middle")]
    for v in (0, 0.25, 0.5, 0.75, 1):
        px = ax(v)
        body.append(line(px, y1, px, y1 + 4))
        body.append(text(px, y1 + 18, f"{v:g}", "middle", mono=True, size=11))
    for i, (r, n, k, d, name) in enumerate(rows):
        cy = y0 + pitch * i + pitch / 2
        color = DISTANCE_COLOR[d]
        end = ax(r)
        body.append(f'<text x="{x0 - 10}" y="{fmt(cy + 4)}" text-anchor="end" font-size="12">'
                    f'<tspan font-family="{MONO}">({n},{k},{d})</tspan> {name}</text>')
        body.append(rect(x0, cy - 10, end - x0, 20, BG2, color, 1.5, 3))
        body.append(text(end + 8, cy + 4, f"{r:.3f}", mono=True, size=11))
        body.append(text(end + 52, cy + 4, f"d = {d}", fill=color, size=12, bold=True))
    return svg("Rate and minimum distance of nine codes built from the small codes of this chapter", body, W, h_)


CHARTS = {
    "fig-24-5-direct-sum-versus-u-u-plus-v.svg": direct_sum_versus_u_u_plus_v,
    "fig-24-7-rate-and-distance.svg": rate_and_distance,
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
