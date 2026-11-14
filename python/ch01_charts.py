"""Write the two charts of chapter 1 of Error Correction Coding as SVG files.

    python3 python/ch01_charts.py --out figures

Every plotted value comes from the functions in ch01_why_codes_exist.py, so a
chart is corrected by fixing those functions and running this script again,
never by editing its output. The SVG follows the book's figure style sheet:
real text in IBM Plex, colors from the page's tokens with paper fallbacks, and
no element that the site's figure check would refuse.
"""

import argparse
import os
from math import log10

from ch01_why_codes_exist import binary_entropy, rep_error

W, H = 640, 400
X0, X1, Y0, Y1 = 84, 604, 24, 324          # plot area in viewBox units
GREEN = "var(--green,#00614a)"
CYAN = "var(--cyan,#005f7a)"
TEXT3 = "var(--text-3,#5a5448)"
MONO = "IBM Plex Mono, monospace"
MINUS = "−"


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


class LogAxis:
    def __init__(self, lo_exp: int, hi_exp: int, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo_exp, hi_exp, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (log10(v) - self.lo) / (self.hi - self.lo)


class LinAxis:
    def __init__(self, lo: float, hi: float, px_lo: float, px_hi: float):
        self.lo, self.hi, self.a, self.b = lo, hi, px_lo, px_hi

    def __call__(self, v: float) -> float:
        return self.a + (self.b - self.a) * (v - self.lo) / (self.hi - self.lo)


def pow10(x: float, y: float, e: int, anchor: str) -> str:
    """A tick label 10^e written as real text, with the exponent raised by a tspan."""
    if e == 0:
        return f'<text x="{fmt(x)}" y="{fmt(y)}" text-anchor="{anchor}" font-family="{MONO}" font-size="11">1</text>'
    exp = f"{MINUS}{-e}" if e < 0 else str(e)
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" text-anchor="{anchor}" font-family="{MONO}" font-size="11">'
            f'10<tspan dy="-5" font-size="8">{exp}</tspan></text>')


def frame(xticks: str, yticks: str, xlabel: str, ylabel: str) -> list:
    return [
        f'<g stroke="currentColor" stroke-width="1" fill="none">'
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}"/></g>',
        xticks,
        yticks,
        f'<text x="{(X0 + X1) // 2}" y="{H - 14}" text-anchor="middle">{xlabel}</text>',
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def log_y_ticks(ax: LogAxis, labelled: range, every: range) -> str:
    out = []
    for e in every:
        y = ax(10.0**e)
        out.append(f'<line x1="{X0 - 4}" y1="{fmt(y)}" x2="{X0}" y2="{fmt(y)}" stroke="currentColor" stroke-width="1"/>')
        if e in labelled:
            out.append(f'<line x1="{X0}" y1="{fmt(y)}" x2="{X1}" y2="{fmt(y)}" stroke="{TEXT3}" stroke-width="1" stroke-opacity="0.35"/>')
            out.append(pow10(X0 - 8, y + 4, e, "end"))
    return "\n".join(out)


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def repetition_error_curves() -> str:
    """Decoded error probability of R1, R3 and R5 against the crossover probability p."""
    x = LogAxis(-3, 0, X0, X1)
    y = LogAxis(-9, 0, Y1, Y0)
    xt = []
    for e in range(-3, 1):
        px = x(10.0**e)
        xt.append(f'<line x1="{fmt(px)}" y1="{Y1}" x2="{fmt(px)}" y2="{Y1 + 4}" stroke="currentColor" stroke-width="1"/>')
        if e < 0:
            xt.append(f'<line x1="{fmt(px)}" y1="{Y0}" x2="{fmt(px)}" y2="{Y1}" stroke="{TEXT3}" stroke-width="1" stroke-opacity="0.35"/>')
        xt.append(pow10(px, Y1 + 20, e, "middle"))
    body = frame("\n".join(xt), log_y_ticks(y, range(-9, 1, 3), range(-9, 1)),
                 "crossover probability p of the channel",
                 "probability the decoded bit is wrong")
    ps = [10 ** (-3 + (log10(0.5) + 3) * i / 100) for i in range(101)]
    curves = [("no code", 1, TEXT3, 1.5, "3 3"), ("R3", 3, GREEN, 2, ""), ("R5", 5, CYAN, 1.5, "")]
    for name, n, color, width, dash in curves:
        body.append(polyline([(x(p), y(rep_error(p, n))) for p in ps], color, width, dash))
        lp = 1.6e-3
        body.append(f'<text x="{fmt(x(lp))}" y="{fmt(y(rep_error(lp, n)) + 17)}" fill="{color if n > 1 else TEXT3}"'
                    f' font-weight="600">{name}</text>')
    half = x(0.5)
    body.append(f'<line x1="{fmt(half)}" y1="{Y0}" x2="{fmt(half)}" y2="{Y1}" stroke="currentColor" stroke-width="1" stroke-dasharray="3 3"/>')
    body.append(f'<text x="{fmt(half - 6)}" y="{Y1 - 8}" text-anchor="end" font-family="{MONO}" font-size="11">p = 1/2</text>')
    return svg("Decoded bit error probability of no code, R3 and R5 against the crossover probability p, on log scales. "
               "R3 falls as 3p squared and R5 as 10p cubed, and all three meet at p = 1/2.", body)


def rate_versus_reliability(p: float = 0.1) -> str:
    """Rate against decoded error probability for R1 to R25 at p = 0.1, with the capacity line."""
    x = LinAxis(0, 1, X0, X1)
    y = LogAxis(-8, 0, Y1, Y0)
    xt = []
    for i in range(6):
        v = i / 5
        px = x(v)
        xt.append(f'<line x1="{fmt(px)}" y1="{Y1}" x2="{fmt(px)}" y2="{Y1 + 4}" stroke="currentColor" stroke-width="1"/>')
        xt.append(f'<text x="{fmt(px)}" y="{Y1 + 20}" text-anchor="middle" font-family="{MONO}" font-size="11">{fmt(v) if v else "0"}</text>')
    body = frame("\n".join(xt), log_y_ticks(y, range(-8, 1, 2), range(-8, 1)),
                 "rate: message bits per channel bit",
                 "probability the decoded bit is wrong")
    cap = 1 - binary_entropy(p)
    cx = x(cap)
    body.append(f'<line x1="{fmt(cx)}" y1="{Y0}" x2="{fmt(cx)}" y2="{Y1}" stroke="{GREEN}" stroke-width="2" stroke-dasharray="3 3"/>')
    body.append(f'<text x="{fmt(cx + 6)}" y="{Y0 + 16}" fill="{GREEN}" font-weight="600">capacity C = {cap:.3f}</text>')
    body.append(f'<text x="{fmt(cx + 6)}" y="{Y1 - 46}" fill="{GREEN}">Shannon: every rate left of this</text>')
    body.append(f'<text x="{fmt(cx + 6)}" y="{Y1 - 30}" fill="{GREEN}">line reaches any error target,</text>')
    body.append(f'<text x="{fmt(cx + 6)}" y="{Y1 - 14}" fill="{GREEN}">given long enough codes</text>')
    pts = [(n, x(1 / n), y(rep_error(p, n))) for n in range(1, 26, 2)]
    body.append(polyline([(px, py) for _, px, py in pts], CYAN, 1.5))
    body.append(f'<g fill="{CYAN}">' + "".join(f'<circle cx="{fmt(px)}" cy="{fmt(py)}" r="3"/>' for _, px, py in pts) + "</g>")
    where = {1: (-6, 18, "end"), 3: (6, 14, "start"), 5: (6, 14, "start"), 25: (0, 18, "middle")}
    for n, px, py in pts:
        if n in where:
            dx, dy, anchor = where[n]
            body.append(f'<text x="{fmt(px + dx)}" y="{fmt(py + dy)}" text-anchor="{anchor}" fill="{CYAN}" font-weight="600">R{n}</text>')
    return svg(f"Rate against decoded error probability for the repetition codes R1 to R25 on a channel with p = {p}. "
               f"The error falls only as the rate falls toward zero, while the capacity line at {cap:.3f} marks the rates "
               "Shannon's theorem says can be made as reliable as wanted.", body)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    for name, text in (("fig-1-6-repetition-error-curves.svg", repetition_error_curves()),
                       ("fig-1-7-rate-versus-reliability.svg", rate_versus_reliability())):
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
