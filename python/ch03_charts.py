"""Write the five charts of chapter 3 of Error Correction Coding as SVG files.

    python3 python/ch03_charts.py --out figures

Every plotted value comes from the functions in ch03_channels_signals_and_likelihoods.py,
so a chart is corrected by fixing those functions and running this script
again, never by editing its output. The SVG follows the book's figure style
sheet: real text in IBM Plex, colors from the page's tokens with paper
fallbacks, and no element that the site's figure check would refuse.
"""

import argparse
import os
from math import log, log10

from ch03_channels_signals_and_likelihoods import (BURST_BITS, BSC_SEED, GE_SEED, Q, bpsk_ber, bsc_simulate,
                                                   gauss_pdf, ge_simulate, hard_r3, runs, solve_db)

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
         'orient="auto"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker></defs>')


def fmt(v: float) -> str:
    return f"{v:.1f}".rstrip("0").rstrip(".")


def signed(v: float, places: int) -> str:
    s = f"{abs(v):.{places}f}"
    return f"{MINUS}{s}" if v < 0 else s


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


def pow10(x: float, y: float, e: int, anchor: str) -> str:
    """A tick label 10^e written as real text, with the exponent raised by a tspan."""
    if e == 0:
        return text(x, y, "1", anchor, mono=True, size=11)
    exp_ = f"{MINUS}{-e}" if e < 0 else str(e)
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" text-anchor="{anchor}" font-family="{MONO}" font-size="11">'
            f'10<tspan dy="-5" font-size="9">{exp_}</tspan></text>')


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = "currentColor", width: float = 1,
         dash: str = "", opacity: str = "", arrow: bool = False) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' stroke-opacity="{opacity}"' if opacity else ""
    m = ' marker-end="url(#ah)"' if arrow else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{stroke}"'
            f' stroke-width="{width}"{d}{o}{m}/>')


def dot(x: float, y: float, fill: str = "currentColor", r: float = 4) -> str:
    return f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{fill}"/>'


def frame(xlabel: str, ylabel: str) -> list:
    return [
        f'<polyline points="{X0},{Y0} {X0},{Y1} {X1},{Y1}" fill="none" stroke="currentColor" stroke-width="1"/>',
        text((X0 + X1) / 2, H - 14, xlabel, "middle"),
        f'<text transform="translate(20 {(Y0 + Y1) // 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>',
    ]


def lin_x_ticks(ax: LinAxis, values: list, labels: list, grid: bool = False) -> list:
    out = []
    for v, lab in zip(values, labels):
        px = ax(v)
        out.append(line(px, Y1, px, Y1 + 4))
        if grid:
            out.append(line(px, Y0, px, Y1, TEXT3, 1, opacity="0.35"))
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


def log_y_ticks(ax: LogAxis, labelled: range, every: range) -> list:
    out = []
    for e in every:
        y = ax(10.0**e)
        out.append(line(X0 - 4, y, X0, y))
        if e in labelled:
            out.append(line(X0, y, X1, y, TEXT3, 1, opacity="0.35"))
            out.append(pow10(X0 - 8, y + 4, e, "end"))
    return out


def polyline(points: list, stroke: str, width: float, dash: str = "") -> str:
    pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{pts}" fill="none" stroke="{stroke}" stroke-width="{width}"{d}/>'


def svg(title: str, body: list, w: int = W, h: int = H) -> str:
    head = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
    return "\n".join([head, f"<title>{title}</title>", *body, "</svg>", ""])


def bpsk_likelihoods(a: float = 1.0, sigma: float = 0.5, prior0: float = 0.9) -> str:
    """The two prior-weighted densities, with the ML threshold at 0 and the MAP threshold where they cross."""
    x = LinAxis(-3, 3, X0, X1)
    y = LinAxis(0, 0.8, Y1, Y0)
    tau = sigma**2 / (2 * a) * log((1 - prior0) / prior0)
    body = frame("received voltage y", "prior times density")
    body += lin_y_ticks(y, [0, 0.2, 0.4, 0.6, 0.8], ["0", "0.2", "0.4", "0.6", "0.8"])
    body += lin_x_ticks(x, [-3, -2, -1, 0, 1, 2, 3], [f"{MINUS}3", f"{MINUS}2", f"{MINUS}1", "0", "1", "2", "3"])
    xs = [-3 + 6 * i / 600 for i in range(601)]
    body.append(polyline([(x(v), y(prior0 * gauss_pdf(v, a, sigma))) for v in xs], "currentColor", 2))
    body.append(polyline([(x(v), y((1 - prior0) * gauss_pdf(v, -a, sigma))) for v in xs], TEXT3, 1.5))
    body.append(line(x(0), Y0 + 24, x(0), Y1, "currentColor", 1, "3 3"))
    body.append(line(x(tau), Y0 + 24, x(tau), Y1, GREEN, 2))
    body.append(text(x(0) + 6, Y0 + 14, "ML, y = 0"))
    body.append(text(x(tau) - 6, Y0 + 14, f"MAP, y = {signed(tau, 3)}", "end", fill=GREEN, weight="600"))
    body.append(text(x(-2), Y0 + 76, "decide 1", "middle"))
    body.append(text(x(2.2), Y0 + 76, "decide 0", "middle"))
    # Legend in the empty top-right corner, clear of the tall curve, which falls below it there.
    lx = 460
    body.append(line(lx, Y0 + 10, lx + 20, Y0 + 10, "currentColor", 2))
    body.append(text(lx + 26, Y0 + 14, f"{prior0} f(y | sent 0)", mono=True, size=11))
    body.append(line(lx, Y0 + 28, lx + 20, Y0 + 28, TEXT3, 1.5))
    body.append(text(lx + 26, Y0 + 32, f"{1 - prior0:.1f} f(y | sent 1)", mono=True, size=11, fill=TEXT3))
    return svg(f"Two prior-weighted bell curves over the received voltage. The tall one is centered at plus {a:g} "
               f"and weighted by {prior0}, the short one at minus {a:g} and weighted by {1 - prior0:.1f}. A dashed "
               f"line at 0 marks the ML threshold and a green line at minus {abs(tau):.3f} marks the MAP threshold, "
               "where the curves cross.", body)


def signal_space() -> str:
    """The plane of two samples with s, minus s, the received point and the decision boundary."""
    cx, cy, k = 320, 200, 120              # origin and pixels per unit, 1.25 units each side

    def px(u: float, v: float) -> tuple:
        return cx + k * u, cy - k * v

    s, y = (0.6, 0.8), (0.9, -0.1)
    d_plus = sum((a - b) ** 2 for a, b in zip(y, s))
    d_minus = sum((a + b) ** 2 for a, b in zip(y, s))
    body = [ARROW]
    body.append(line(*px(-1.25, 0), *px(1.3, 0), TEXT3, 1, arrow=True))
    body.append(line(*px(0, -1.25), *px(0, 1.3), TEXT3, 1, arrow=True))
    body.append(text(*px(1.32, 0.04), "first sample", fill=TEXT3))
    body.append(text(px(0, 1.3)[0] + 8, px(0, 1.3)[1] + 10, "second sample", fill=TEXT3))
    # The boundary runs along (0.8, -0.6), perpendicular to s, clipped to the square.
    body.append(line(*px(-1.25, 0.9375), *px(1.25, -0.9375), "currentColor", 1.5, "3 3"))
    body.append(text(*px(0.98, 0.42), "decide 0", "middle"))
    body.append(text(*px(-0.95, -0.38), "decide 1", "middle"))
    body.append(line(*px(0, 0), *px(*s), TEXT3, 1, "3 3"))
    body.append(line(*px(*y), *px(-s[0], -s[1]), ORANGE, 1.5, "3 3"))
    body.append(line(*px(*y), *px(*s), GREEN, 2))
    body.append(dot(*px(*s)))
    body.append(dot(*px(-s[0], -s[1])))
    body.append(dot(*px(*y), CYAN))
    body.append(text(px(*s)[0] + 8, px(*s)[1] - 6, "s (0.6, 0.8)", mono=True, size=11))
    body.append(text(px(-s[0], -s[1])[0] - 8, px(-s[0], -s[1])[1] + 14, f"{MINUS}s ({MINUS}0.6, {MINUS}0.8)",
                     "end", mono=True, size=11))
    body.append(text(px(*y)[0] + 8, px(*y)[1] + 16, f"y (0.9, {MINUS}0.1)", mono=True, size=11, fill=CYAN))
    mid = px((y[0] + s[0]) / 2, (y[1] + s[1]) / 2)
    body.append(text(mid[0] - 8, mid[1] + 2, f"{d_plus:.2f}", "end", mono=True, size=11, fill=GREEN, weight="600"))
    mid = px((y[0] - s[0]) / 2, (y[1] - s[1]) / 2)
    body.append(text(mid[0] + 8, mid[1] + 22, f"{d_minus:.2f}", "middle", mono=True, size=11, fill=ORANGE))
    return svg("The plane of two received samples. The signal points s at 0.6, 0.8 and minus s lie on opposite "
               "sides of a dashed boundary line through the origin, perpendicular to s. The received point y at "
               f"0.9, minus 0.1 is joined to s by a green segment of squared length {d_plus:.2f} and to minus s by "
               f"an orange dashed segment of squared length {d_minus:.2f}, so it is decided as 0.", body)


def bpsk_ber_chart() -> str:
    """Uncoded BPSK, which soft R3 equals, against hard-decision R3, with the 1e-5 crossings marked."""
    x = LinAxis(0, 12, X0, X1)
    y = LogAxis(-9, 0, Y1, Y0)
    body = frame("Eb/N0 (dB)", "bit error probability")
    body += log_y_ticks(y, range(-9, 1, 1), range(-9, 1))
    body += lin_x_ticks(x, list(range(0, 13, 2)), [str(v) for v in range(0, 13, 2)], grid=True)
    xs = [12 * i / 480 for i in range(481)]
    body.append(polyline([(x(v), y(hard_r3(v))) for v in xs], ORANGE, 2))
    body.append(polyline([(x(v), y(bpsk_ber(v))) for v in xs], GREEN, 2))
    target = 1e-5
    soft_db, hard_db = solve_db(bpsk_ber, target), solve_db(hard_r3, target)
    body.append(line(X0, y(target), X1, y(target), "currentColor", 1, "3 3"))
    body.append(f'<text x="{X0 + 6}" y="{fmt(y(target) - 6)}">target 10<tspan dy="-5" font-size="9">'
                f'{MINUS}5</tspan></text>')
    for db_, anchor, dx, color in ((soft_db, "end", -4, GREEN), (hard_db, "start", 4, ORANGE)):
        body.append(line(x(db_), y(target), x(db_), Y1, color, 1, "3 3"))
        body.append(dot(x(db_), y(target), color, 3.5))
        body.append(text(x(db_) + dx, Y1 - 8, f"{db_:.2f} dB", anchor, fill=color, mono=True, size=11))
    # The soft curve's label sits below and left of it, the hard curve's above and right: both are empty.
    v = 8
    body.append(text(x(v) - 8, y(bpsk_ber(v)) + 14, "uncoded = soft R3", "end", fill=GREEN, weight="600"))
    v = 9
    body.append(text(x(v) + 8, y(hard_r3(v)) - 8, "hard R3", "start", fill=ORANGE, weight="600"))
    return svg("Bit error probability on a logarithmic scale against Eb/N0 from 0 to 12 dB. The green curve is "
               "uncoded BPSK, which soft-decision R3 equals, and the orange curve is hard-decision R3, which lies "
               f"above it. A dashed line at 10 to the minus 5 meets the green curve at {soft_db:.2f} dB and the "
               f"orange curve at {hard_db:.2f} dB.", body)


def llr_versus_voltage(a: float = 1.0, sigma: float = 0.5) -> str:
    """The soft LLR 2ay/sigma^2 against the hard decision's two values, with the points where they agree."""
    c = 2 * a / sigma**2
    p = Q(a / sigma)
    hard = log((1 - p) / p)
    cross = hard / c
    x = LinAxis(-1.5, 1.5, X0, X1)
    y = LinAxis(-12, 12, Y1, Y0)
    body = frame("received voltage y", "LLR")
    body += lin_y_ticks(y, [-12, -8, -4, 0, 4, 8, 12], [f"{MINUS}12", f"{MINUS}8", f"{MINUS}4", "0", "4", "8", "12"])
    body += lin_x_ticks(x, [-1.5, -1, -0.5, 0, 0.5, 1, 1.5],
                        [f"{MINUS}1.5", f"{MINUS}1", f"{MINUS}0.5", "0", "0.5", "1", "1.5"], grid=True)
    body.append(line(x(-1.5), y(-hard), x(0), y(-hard), ORANGE, 2))
    body.append(line(x(0), y(hard), x(1.5), y(hard), ORANGE, 2))
    body.append(line(x(0), y(-hard), x(0), y(hard), ORANGE, 1, "3 3"))
    body.append(line(x(-1.5), y(-1.5 * c), x(1.5), y(1.5 * c), GREEN, 2))
    for sgn in (1, -1):
        body.append(line(x(sgn * cross), y(sgn * hard), x(sgn * cross), y(0), "currentColor", 1, "3 3"))
        body.append(dot(x(sgn * cross), y(sgn * hard)))
    body.append(text(x(cross) - 8, y(hard) - 8, f"{cross:.2f}", "end", mono=True, size=11))
    body.append(text(x(-cross) + 8, y(-hard) + 18, f"{MINUS}{cross:.2f}", "start", mono=True, size=11))
    body.append(text(x(1.2) - 10, y(1.2 * c), f"soft, L = {c:g}y", "end", fill=GREEN, weight="600"))
    body.append(text(x(1.45), y(hard) - 8, f"hard, +{hard:.2f}", "end", fill=ORANGE, weight="600"))
    body.append(text(x(-1.45), y(-hard) + 20, f"hard, {MINUS}{hard:.2f}", "start", fill=ORANGE, weight="600"))
    return svg(f"The LLR against the received voltage from minus 1.5 to 1.5. The soft LLR is a green straight line "
               f"through the origin with slope {c:g}. The hard LLR is an orange step, minus {hard:.2f} for negative "
               f"voltages and plus {hard:.2f} for positive ones. The two meet at plus and minus {cross:.2f}.", body)


def burst_pattern() -> str:
    """Error positions from the seeded Gilbert-Elliott and BSC simulators, bad spells shaded."""
    w, h = 640, 250
    left, right = 40, 600
    states, ge_errors = ge_simulate()
    bsc_errors = bsc_simulate()
    n = BURST_BITS

    def bx(i: float) -> float:
        return left + (right - left) * i / n

    body = []
    strips = ((40, f"Gilbert-Elliott channel, seed {GE_SEED}: {sum(ge_errors)} errors, bad spells shaded", ge_errors,
               True),
              (130, f"memoryless BSC(0.0109), seed {BSC_SEED}: {sum(bsc_errors)} errors", bsc_errors, False))
    for top, label, errors, shade in strips:
        body.append(text(left, top - 10, label))
        if shade:
            for start, length in runs(states):
                body.append(f'<rect x="{fmt(bx(start))}" y="{top}" width="{fmt(max(bx(start + length) - bx(start), 1))}"'
                            f' height="40" fill="{ORANGE}" fill-opacity="0.3" stroke="none"/>')
        body.append(f'<rect x="{left}" y="{top}" width="{right - left}" height="40" rx="3" fill="none" '
                    f'stroke="currentColor" stroke-width="1"/>')
        for i, e in enumerate(errors):
            if e:
                body.append(line(bx(i + 0.5), top + 6, bx(i + 0.5), top + 34, RED, 1.5))
    base = 178
    for i in range(0, n + 1, 500):
        body.append(line(bx(i), base, bx(i), base + 4))
        body.append(text(bx(i), base + 18, f"{i:,}", "middle", mono=True, size=11))
    body.append(text((left + right) / 2, base + 46, "bit position", "middle"))
    ge_spells = len(runs(states))
    return svg(f"Two strips of {n} bits with a red tick at every bit error. In the top strip, from the Gilbert-Elliott "
               f"channel, {sum(ge_errors)} errors crowd into {ge_spells} shaded spells in the bad state. In the bottom "
               f"strip, from a memoryless channel with the same average error rate, {sum(bsc_errors)} errors are "
               "scattered.", body, w, h)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="figures", help="directory to write the SVG files into")
    out = ap.parse_args().out
    os.makedirs(out, exist_ok=True)
    charts = (("fig-3-2-bpsk-likelihoods.svg", bpsk_likelihoods()),
              ("fig-3-3-signal-space.svg", signal_space()),
              ("fig-3-5-bpsk-ber.svg", bpsk_ber_chart()),
              ("fig-3-6-llr-versus-voltage.svg", llr_versus_voltage()),
              ("fig-3-8-burst-pattern.svg", burst_pattern()))
    for name, content in charts:
        with open(os.path.join(out, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
        print(f"wrote {os.path.join(out, name)}")


if __name__ == "__main__":
    main()
