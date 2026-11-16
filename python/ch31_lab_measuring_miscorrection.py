"""Reproduce every number in chapter 31 of Error Correction Coding, "Lab: Measuring Miscorrection".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch31_lab_measuring_miscorrection.py
    python3 python/ch31_lab_measuring_miscorrection.py --seed 7
    python3 python/ch31_lab_measuring_miscorrection.py --draws 5000
    python3 python/ch31_lab_measuring_miscorrection.py --weight 5
    python3 python/ch31_lab_measuring_miscorrection.py --stride 5

Run with no options, each check prints the value computed here beside the
value the chapter states, and the program stops with exit status 1 at the
first pair that disagrees. Whole numbers are compared exactly, and every
fractional value is printed with a fixed number of places, so its comparison
is exact too.

The three (72,64) check matrices are the ones "SECDED and Hsiao Codes" builds,
taken from ch20_secded_and_hsiao_codes.py and compared with the columns of
the model that generated the bench, tb/secded72/secded72_model.py, which the
hardware run records in data/ch31/miscorrection_counts.json. Every three-bit and four-bit error pattern
of each code is then classified, 59,640 and 1,028,790 of them, and the counts
of weight 5 and 6 come from a recursion over syndrome counts. The seeded
campaign draws the same error patterns as tb/miscorrection/test_miscorrection.py,
from random.Random(31) through its random() method alone, so its counts are the
counts the hardware must produce, and the program checks them against the data
file the hardware run wrote.

The options change the experiment and print what changes, without the checks
against the chapter: --seed and --draws rerun the software campaign, --weight 5
prints the exact split of five-bit errors beside a campaign of that weight, and
--stride S prints how many adjacent triples miscorrect when cell k of the word
holds codeword position S * k mod 72.

A column of a check matrix with r rows is held as an integer whose bit i is
the entry in row i. An error pattern is a list of codeword positions, and a
codeword position is the index of a column, so position 0 is the first column.
"""

import argparse
import json
import random
import sys
from itertools import combinations
from math import ceil, comb, exp, log, sqrt
from pathlib import Path

from ch20_secded_and_hsiao_codes import (
    extended_hamming, hsiao_72_64, hsiao_72_64_dictionary, require, same, section)

LAB = Path(__file__).resolve().parent.parent
DATA = LAB / "data" / "ch31" / "miscorrection_counts.json"

N = 72                      # codeword length
K = 64                      # data bits
SEED = 31                   # the campaign's seed
DRAWS = 100_000             # draws of each weight in the recorded campaign
Z = 1.96                    # the two-sided 95 percent point of the normal distribution
CHECKPOINTS = (100, 200, 500, 1000, 2000, 5000, 10_000, 20_000, 50_000, 100_000)

# The three codes, in the order the chapter lists them, by the names the
# hardware model in tb/secded72/secded72_model.py uses.
CODES = ("hsiao", "hsiao_dictionary", "hamming")
LABELS = {"hsiao": "Hsiao balanced", "hsiao_dictionary": "Hsiao dictionary", "hamming": "Extended Hamming"}

CORRECT, FLAGGED, MISCORRECTED, UNDETECTED = "correct", "flagged", "miscorrected", "undetected"
OUTCOMES = (CORRECT, FLAGGED, MISCORRECTED, UNDETECTED)


# ---------------------------------------------------------------- the codes

def columns(name: str) -> list[int]:
    """The 72 columns of one code's check matrix, from chapter 20's constructions."""
    if name == "hsiao":
        return hsiao_72_64()
    if name == "hsiao_dictionary":
        return hsiao_72_64_dictionary()
    return extended_hamming(N, 7)


def check_positions(name: str) -> list[int]:
    """The codeword positions that hold check bits; the rest hold the 64 data bits."""
    if name == "hamming":
        return [0, 1, 2, 4, 8, 16, 32, 64]
    return list(range(64, 72))


class Code:
    def __init__(self, name: str):
        self.name = name
        self.cols = columns(name)
        self.position_of = {c: j for j, c in enumerate(self.cols)}
        self.checks = check_positions(name)
        self.data_mask = sum(1 << p for p in range(N) if p not in self.checks)


def outcome(code: Code, positions) -> str:
    """What the decoder does with the error pattern that flips these positions.

    The syndrome is the XOR of the flipped columns. A zero syndrome leaves the
    word alone, a syndrome equal to column j flips bit j, and any other syndrome
    raises the uncorrectable flag. The outcome is then read from the data bits
    that leave the decoder: no flag and no data bit wrong is correct, no flag and
    some data bit wrong is a miscorrection if a bit was flipped and an undetected
    error if none was.
    """
    s = 0
    e = 0
    for p in positions:
        s ^= code.cols[p]
        e ^= 1 << p
    if s == 0:
        return UNDETECTED if e & code.data_mask else CORRECT
    j = code.position_of.get(s)
    if j is None:
        return FLAGGED
    return MISCORRECTED if (e ^ (1 << j)) & code.data_mask else CORRECT


# ---------------------------------------------------------------- exact counts

def enumerate_weight(code: Code, w: int) -> dict:
    """Classify every error pattern of weight w by listing them all."""
    counts = dict.fromkeys(OUTCOMES, 0)
    for positions in combinations(range(N), w):
        counts[outcome(code, positions)] += 1
    return counts


def syndrome_counts(cols, top: int) -> list[list[int]]:
    """counts[w][s], the number of weight-w error patterns whose syndrome is s, for w up to top.

    Each column in turn either joins a pattern or stays out of it, so the
    patterns of weight w that use the first t + 1 columns are those of weight w
    on the first t columns, plus those of weight w - 1 with column t added,
    whose syndromes are moved by an XOR with column t.
    """
    counts = [[0] * 256 for _ in range(top + 1)]
    counts[0][0] = 1
    for c in cols:
        for w in range(top, 0, -1):
            row, below = counts[w], counts[w - 1]
            for s in range(256):
                if below[s]:
                    row[s ^ c] += below[s]
    return counts


def matched(code: Code, counts_w: list[int]) -> int:
    """Patterns of one weight whose syndrome equals some column, from the syndrome counts."""
    return sum(counts_w[c] for c in code.cols)


def check_only(code: Code) -> list[int]:
    """Q_w, the number of codewords of weight w whose ones all sit in check positions.

    Only the 256 subsets of the eight check positions need trying. A pattern
    whose bits after decoding form such a codeword leaves the data correct.
    """
    q = [0] * (N + 1)
    for r in range(len(code.checks) + 1):
        for subset in combinations(code.checks, r):
            s = 0
            for p in subset:
                s ^= code.cols[p]
            if s == 0:
                q[r] += 1
    return q


def exact_split(code: Code, counts: list[list[int]], w: int) -> dict:
    """The four outcome counts of weight w, from the syndrome counts and the check-only codewords.

    A pattern with zero syndrome is itself a codeword, and it is undetected unless
    the codeword lies in the check positions alone. A pattern whose syndrome is
    column j has e + u_j as a codeword, of weight w + 1 when j is outside e and
    w - 1 when j is inside it, and the data come out correct only when that
    codeword lies in the check positions alone.
    """
    q = check_only(code)
    a_w = counts[w][0]
    m_w = matched(code, counts[w])
    q_up = (w + 1) * q[w + 1] if w + 1 <= N else 0
    q_down = (N - w + 1) * q[w - 1] if w >= 1 else 0
    return {CORRECT: q[w] + q_up + q_down, FLAGGED: comb(N, w) - a_w - m_w,
            MISCORRECTED: m_w - q_up - q_down, UNDETECTED: a_w - q[w]}


def dual_weights(code: Code) -> list[int]:
    """B_j, the number of the 256 codewords of the dual code with weight j.

    The dual code is spanned by the eight rows of the check matrix, and row i,
    read across the 72 positions, has a 1 wherever a column has bit i set.
    """
    rows = [sum(1 << j for j, c in enumerate(code.cols) if (c >> i) & 1) for i in range(8)]
    b = [0] * (N + 1)
    for subset in range(256):
        word = 0
        for i in range(8):
            if (subset >> i) & 1:
                word ^= rows[i]
        b[word.bit_count()] += 1
    return b


def krawtchouk(i: int, x: int, n: int = N) -> int:
    """The coefficient of y^i in (1 + y)^(n - x) (1 - y)^x."""
    return sum((-1) ** t * comb(x, t) * comb(n - x, i - t) for t in range(i + 1))


def macwilliams(b: list[int], i: int) -> int:
    """A_i of the code from the dual's weights, by the MacWilliams identity."""
    total = sum(b[x] * krawtchouk(i, x) for x in range(N + 1))
    assert total % 256 == 0
    return total // 256


# ---------------------------------------------------------------- bursts

def layout(kind: str, code_name: str = "") -> list[int]:
    """order[k] is the codeword position stored in cell k of the word."""
    if kind == "natural":
        return list(range(N))
    if kind.startswith("stride"):
        s = int(kind[6:])
        return [s * k % N for k in range(N)]
    if kind == "searched":
        # No order flags every window of extended Hamming, whose floor is 46
        # windows in any layout, so no search runs for it and stride 11, which
        # meets that floor, stands in.
        if code_name == "hamming":
            return layout("stride11")
        order, _ = search_layout(columns(code_name))
        return order
    raise ValueError(kind)


def adjacent_outcomes(code: Code, order: list[int]) -> dict:
    """Outcomes of the 70 triples of three neighboring cells, cells k, k + 1 and k + 2."""
    counts = dict.fromkeys(OUTCOMES, 0)
    for k in range(N - 2):
        counts[outcome(code, order[k:k + 3])] += 1
    return counts


def search_layout(cols) -> tuple:
    """A cell order in which no three neighboring cells miscorrect, found by depth-first search.

    Cells are filled from the left. A position may go into the next cell only if
    the XOR of its column with the columns in the two cells before it is not a
    column. When no unused position fits, the search backs up one cell and tries
    the next candidate there. Candidates are tried lowest position first, so the
    result is the first such order in that sense. Returns the order and the
    number of steps the search took.
    """
    present = set(cols)
    order, used, choice = [], [False] * N, [0] * N
    depth = steps = 0
    while True:
        steps += 1
        if depth == N:
            return order, steps
        j = choice[depth]
        while j < N and (used[j] or (depth >= 2 and cols[order[-2]] ^ cols[order[-1]] ^ cols[j] in present)):
            j += 1
        if j < N:
            choice[depth] = j + 1
            order.append(j)
            used[j] = True
            depth += 1
            if depth < N:
                choice[depth] = 0
        else:
            if depth == 0:
                return None, steps
            depth -= 1
            used[order.pop()] = False


# ---------------------------------------------------------------- the campaign

def draw(rng: random.Random, weight: int) -> tuple:
    """One case of the campaign: a 64-bit data word and `weight` distinct positions.

    Only rng.random() is used, the one method whose sequence Python's
    documentation promises across versions. Each half of the data word is a
    random() draw scaled to 32 bits, and each position is a draw scaled to 72
    and truncated, drawn again whenever it repeats one already chosen, so every
    set of `weight` positions is equally likely, to a few parts in 10**14.
    """
    data = int(rng.random() * 2 ** 32) | int(rng.random() * 2 ** 32) << 32
    positions = []
    while len(positions) < weight:
        p = int(rng.random() * N)
        if p not in positions:
            positions.append(p)
    return data, positions


def campaign_cases(seed: int, draws: int, weights) -> list:
    """Every case of a campaign, weight by weight, from one generator: (weight, data, positions)."""
    rng = random.Random(seed)
    return [(w, *draw(rng, w)) for w in weights for _ in range(draws)]


def adjacent_data(seed: int, count: int) -> list[int]:
    """The data words the hardware writes under the adjacent triples, from a generator of their own."""
    rng = random.Random(seed + 1)
    return [draw(rng, 0)[0] for _ in range(count)]


def tally(cases, code: Code) -> dict:
    """counts[w][outcome] over the campaign, plus the running count of the outcome each weight tracks."""
    counts, running = {}, {}
    for w, _, positions in cases:
        c = counts.setdefault(w, dict.fromkeys(OUTCOMES, 0))
        c[outcome(code, positions)] += 1
        n = sum(c.values())
        if n in CHECKPOINTS:
            running.setdefault(w, []).append(c[tracked(w)])
    return {"counts": counts, "running": running}


def tracked(w: int) -> str:
    """The outcome a weight's rate is about: miscorrection for odd weights, undetected for even."""
    return MISCORRECTED if w % 2 else UNDETECTED


# ---------------------------------------------------------------- reading a rate

def wilson(m: int, n: int, z: float = Z) -> tuple:
    """The Wilson score interval for a proportion, m successes in n trials."""
    ph = m / n
    center = ph + z * z / (2 * n)
    half = z * sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    scale = 1 + z * z / n
    return (center - half) / scale, (center + half) / scale


def ceil_whole(x: float) -> int:
    """The smallest whole number at least x, forgiving the rounding of a value that is whole on paper."""
    return ceil(x - 1e-9)


def draws_normal(p: float, h: float, z: float = Z) -> int:
    """Draws for a normal-approximation half-width of h at rate p."""
    return ceil_whole(z * z * p * (1 - p) / (h * h))


def draws_chebyshev(h: float, delta: float) -> int:
    """Draws that make P(|estimate - p| >= h) at most delta for every p, by Chebyshev's inequality."""
    return ceil_whole(1 / (4 * delta * h * h))


def draws_hoeffding(h: float, delta: float) -> int:
    """Draws that make P(|estimate - p| >= h) at most delta for every p, by Hoeffding's inequality."""
    return ceil_whole(log(2 / delta) / (2 * h * h))


def half_chebyshev(n: int, delta: float) -> float:
    return sqrt(1 / (4 * delta * n))


def half_hoeffding(n: int, delta: float) -> float:
    return sqrt(log(2 / delta) / (2 * n))


def pct(x: float, places: int = 3) -> str:
    return f"{100 * x:.{places}f}"


# ---------------------------------------------------------------- designer's number

WORDS = 4096                # words in the illustrative array
RATE = 1e-6                 # upsets per word per hour, for illustration
YEAR = 8760                 # hours in a year
SCRUBS = (("1 day", 24), ("30 days", 720), ("1 year", 8760))


DISTINCT = (N - 1) * (N - 2) / N ** 2     # chance three strikes land on three different cells


def silent_per_year(share: float, t_hours: float) -> float:
    """Silent corruptions a year when each word is decoded once per interval of t_hours.

    The number of upsets in one word over one interval is Poisson with mean
    mu = RATE * t_hours, so three of them arrive with probability
    mu^3 / 6 * exp(-mu). The three strikes land on three different cells with
    probability DISTINCT, and a share of those triples leaves the decoder as a
    wrong word with no flag.
    """
    mu = RATE * t_hours
    p3 = mu ** 3 / 6 * exp(-mu)
    return WORDS * (YEAR / t_hours) * p3 * DISTINCT * share


# ---------------------------------------------------------------- what the chapter states

TRIPLES, QUADRUPLES = comb(N, 3), comb(N, 4)

STATED_EXACT = {      # (triples miscorrected, triples flagged, quadruples undetected, quadruples flagged)
    "hsiao": (33568, 26072, 8392, 1020398),
    "hsiao_dictionary": (34032, 25608, 8508, 1020282),
    "hamming": (45304, 14336, 11326, 1017464),
}
STATED_SHARES = {     # (triples miscorrected, one place; three places; quadruples undetected)
    "hsiao": ("56.3", "56.284", "0.816"),
    "hsiao_dictionary": ("57.1", "57.062", "0.827"),
    "hamming": ("76.0", "75.962", "1.101"),
}
STATED_WEIGHT_6 = {"hsiao": 1216136, "hsiao_dictionary": 1215217, "hamming": 1446144}
STATED_WEIGHT_5 = {"hsiao": (7867472, "56.230"), "hsiao_dictionary": (7869846, "56.247"),
                   "hamming": (9447032, "67.520")}
STATED_CAMPAIGN = {   # weight: (count of the tracked outcome, flagged, percent, Wilson low, Wilson high)
    "hsiao": {3: (56222, 43778, "56.222", "55.914", "56.529"), 4: (864, 99136, "0.864", "0.808", "0.923")},
    "hsiao_dictionary": {3: (56984, 43016, "56.984", "56.677", "57.291"),
                         4: (829, 99171, "0.829", "0.775", "0.887")},
    "hamming": {3: (75975, 24025, "75.975", "75.709", "76.239"), 4: (1122, 98878, "1.122", "1.059", "1.189")},
}
STATED_ADJACENT = {   # miscorrected windows of 70 in each layout
    "natural": (35, 28, 70), "stride5": (49, 48, 65), "stride7": (43, 42, 53),
    "stride11": (39, 42, 46), "stride13": (39, 33, 46), "stride29": (33, 29, 46), "searched": (0, 0, 46),
}
STATED_EXPECTED_WINDOWS = {"hsiao": "39.40", "hsiao_dictionary": "39.94", "hamming": "53.17"}
STATED_SEARCH_STEPS = {"hsiao": 5039, "hsiao_dictionary": 8291}
STATED_SILENT = {     # silent corruptions a year at each scrub interval
    "hsiao": ("1.86e-09", "1.67e-06", "2.45e-04"),
    "hamming": ("2.51e-09", "2.26e-06", "3.31e-04"),
}
HARDWARE_LAYOUTS = ("natural", "stride11", "searched")


def label(name: str) -> str:
    return LABELS[name]


# ---------------------------------------------------------------- the checks

def check_matrices(data: dict) -> None:
    section("The three check matrices, against the model that generated the bench")
    for name in CODES:
        require(f"{label(name)}: columns equal the model's that generated the bench", columns(name) == data["codes"][name]["columns"])
    for name in ("hsiao", "hsiao_dictionary"):
        require(f"{label(name)}: every column has odd weight", all(c.bit_count() % 2 for c in columns(name)))
    require("Extended Hamming: every column has a 1 in the parity row",
            all((c >> 7) & 1 for c in columns("hamming")))
    for name in CODES:
        same(f"{label(name)}: codewords in the check positions alone", sum(check_only(Code(name))[1:]), 0)


def check_enumeration() -> dict:
    section(f"Every triple and every quadruple, enumerated ({TRIPLES} and {QUADRUPLES} of them)")
    exact = {}
    for name in CODES:
        code = Code(name)
        t, f = enumerate_weight(code, 3), enumerate_weight(code, 4)
        exact[name] = {3: t, 4: f}
        mis, flag3, und, flag4 = STATED_EXACT[name]
        one, three, quad = STATED_SHARES[name]
        same(f"{label(name)}: triples miscorrected", t[MISCORRECTED], mis)
        same(f"{label(name)}: triples flagged", t[FLAGGED], flag3)
        same(f"{label(name)}: triples correct or undetected", t[CORRECT] + t[UNDETECTED], 0)
        same(f"{label(name)}: percent of triples miscorrected", f"{100 * t[MISCORRECTED] / TRIPLES:.1f}", one)
        same(f"{label(name)}: the same percent to three places", pct(t[MISCORRECTED] / TRIPLES), three)
        same(f"{label(name)}: quadruples undetected", f[UNDETECTED], und)
        same(f"{label(name)}: quadruples flagged", f[FLAGGED], flag4)
        same(f"{label(name)}: quadruples correct or miscorrected", f[CORRECT] + f[MISCORRECTED], 0)
        same(f"{label(name)}: percent of quadruples undetected", pct(f[UNDETECTED] / QUADRUPLES), quad)
    return exact


def check_weight_distribution(exact: dict) -> None:
    section("Weight-4 codewords, the syndrome-count recursion and the MacWilliams identity")
    for name in CODES:
        code = Code(name)
        counts = syndrome_counts(code.cols, 6)
        a4 = counts[4][0]
        same(f"{label(name)}: A_4 from the recursion", a4, STATED_EXACT[name][2])
        same(f"{label(name)}: A_4 from the MacWilliams identity", macwilliams(dual_weights(code), 4), a4)
        require(f"{label(name)}: triples miscorrected equal 4 A_4", exact[name][3][MISCORRECTED] == 4 * a4)
        require(f"{label(name)}: recursion split equals enumeration, weights 3, 4",
                all(exact_split(code, counts, w) == exact[name][w] for w in (3, 4)))
        same(f"{label(name)}: A_6 from the recursion", counts[6][0], STATED_WEIGHT_6[name])
        same(f"{label(name)}: A_6 from the MacWilliams identity", macwilliams(dual_weights(code), 6), counts[6][0])
        five = exact_split(code, counts, 5)
        mis5, pct5 = STATED_WEIGHT_5[name]
        same(f"{label(name)}: weight-5 patterns miscorrected", five[MISCORRECTED], mis5)
        require(f"{label(name)}: weight 5 miscorrected = 68 A_4 + 6 A_6",
                five[MISCORRECTED] == 68 * a4 + 6 * counts[6][0])
        same(f"{label(name)}: weight 5 correct or undetected", five[CORRECT] + five[UNDETECTED], 0)
        same(f"{label(name)}: percent of weight 5 miscorrected", pct(five[MISCORRECTED] / comb(N, 5)), pct5)


def campaign_lines(seed: int, draws: int, weights, exact: dict | None, show: bool = True) -> dict:
    """Run the software campaign and print its counts, its intervals and, given exact shares, the verdicts."""
    cases = campaign_cases(seed, draws, weights)
    results = {}
    for name in CODES:
        code = Code(name)
        results[name] = tally(cases, code)
        for w in weights if show else ():
            c = results[name]["counts"][w]
            k = c[tracked(w)]
            lo, hi = wilson(k, draws)
            line = f"  {label(name):<17} weight {w}  {tracked(w):<12} {k:>7} of {draws}  {pct(k / draws)} percent" \
                   f"  95% interval {pct(lo)} to {pct(hi)}"
            if exact is not None:
                share = exact[name][w][tracked(w)] / comb(N, w)
                line += f"  exact {pct(share)}  {'inside' if lo <= share <= hi else 'OUTSIDE'}"
            print(line)
    return {"cases": cases, "results": results}


def check_campaign(exact: dict) -> dict:
    section(f"The seeded campaign: seed {SEED}, {DRAWS} draws of weight 3 and of weight 4")
    run = campaign_lines(SEED, DRAWS, (3, 4), exact, show=False)
    for name in CODES:
        for w in (3, 4):
            c = run["results"][name]["counts"][w]
            k = c[tracked(w)]
            count, flagged, percent, low, high = STATED_CAMPAIGN[name][w]
            lo, hi = wilson(k, DRAWS)
            share = exact[name][w][tracked(w)] / comb(N, w)
            same(f"{label(name)}, weight {w}: {tracked(w)}", k, count)
            same(f"{label(name)}, weight {w}: flagged", c[FLAGGED], flagged)
            same(f"{label(name)}, weight {w}: other outcomes", DRAWS - k - c[FLAGGED], 0)
            same(f"{label(name)}, weight {w}: percent {tracked(w)}", pct(k / DRAWS), percent)
            same(f"{label(name)}, weight {w}: Wilson interval", f"{pct(lo)} to {pct(hi)}", f"{low} to {high}")
            require(f"{label(name)}, weight {w}: interval holds the exact share", lo <= share <= hi)
    same("chance that six 95% intervals all hold their shares", f"{0.95 ** 6:.3f}", "0.735")
    return run


def paired_counts(cases) -> tuple:
    """Weight-3 draws miscorrected by the balanced matrix alone, and by the dictionary matrix alone."""
    h, d = Code("hsiao"), Code("hsiao_dictionary")
    only_h = only_d = 0
    for w, _, positions in cases:
        if w != 3:
            continue
        a, b = outcome(h, positions) == MISCORRECTED, outcome(d, positions) == MISCORRECTED
        only_h += a and not b
        only_d += b and not a
    return only_h, only_d


def check_hardware(data: dict, run: dict) -> None:
    section(f"The hardware campaign recorded in data/ch31/{DATA.name}")
    same("seed of the hardware campaign", data["seed"], SEED)
    same("draws of each weight", data["draws"], DRAWS)
    same("weights drawn", " ".join(str(w) for w in data["weights"]), "3 4")
    for name in CODES:
        record = data["codes"][name]
        ours = run["results"][name]
        require(f"{label(name)}: hardware counts equal this program's",
                all(record["campaign"][f"weight_{w}"] == ours["counts"][w] for w in (3, 4)))
        require(f"{label(name)}: hardware running counts equal this program's",
                all(record["running"][f"weight_{w}"] == ours["running"][w] for w in (3, 4)))
        same(f"{label(name)}: cases where the hardware differs from the model", record["model_mismatches"], 0)
        same(f"{label(name)}: cases where it differs from the prediction", record["prediction_mismatches"], 0)
        code = Code(name)
        require(f"{label(name)}: hardware adjacent windows equal this program's",
                all(record["adjacent"][kind] == adjacent_outcomes(code, layout(kind, name))
                    for kind in HARDWARE_LAYOUTS))
    only_h, only_d = paired_counts(run["cases"])
    same("weight-3 draws miscorrected by the balanced matrix alone", data["paired"]["hsiao_only"], only_h)
    same("weight-3 draws miscorrected by the dictionary matrix alone", data["paired"]["dictionary_only"], only_d)


def check_rates(run: dict) -> None:
    section("Reading a rate")
    p = STATED_EXACT["hsiao"][0] / TRIPLES
    same("standard error at the balanced share, 100000 draws, points", pct(sqrt(p * (1 - p) / DRAWS)), "0.157")
    lo, hi = wilson(STATED_CAMPAIGN["hsiao"][3][0], DRAWS)
    same("half-width of the balanced Wilson interval, points", pct((hi - lo) / 2), "0.307")
    same("draws for a half-width of 1 point at the balanced share", draws_normal(p, 0.01), 9453)
    same("draws for a half-width of 1 point at p = 1/2", draws_normal(0.5, 0.01), 9604)
    same("draws for 1 point with 95% certainty, Chebyshev", draws_chebyshev(0.01, 0.05), 50000)
    same("draws for 1 point with 95% certainty, Hoeffding", draws_hoeffding(0.01, 0.05), 18445)
    same("Chebyshev half-width at 100000 draws, points", pct(half_chebyshev(DRAWS, 0.05)), "0.707")
    same("Hoeffding half-width at 100000 draws, points", pct(half_hoeffding(DRAWS, 0.05)), "0.429")
    same("normal half-width at p = 1/2 and 100000 draws, points", pct(Z * sqrt(0.25 / DRAWS)), "0.310")

    h, d = Code("hsiao"), Code("hsiao_dictionary")
    only_h = only_d = 0
    for positions in combinations(range(N), 3):
        a, b = outcome(h, positions) == MISCORRECTED, outcome(d, positions) == MISCORRECTED
        only_h += a and not b
        only_d += b and not a
    p1, p2 = STATED_EXACT["hsiao"][0] / TRIPLES, STATED_EXACT["hsiao_dictionary"][0] / TRIPLES
    diff = p2 - p1
    same("triples miscorrected by the balanced matrix alone", only_h, 5583)
    same("triples miscorrected by the dictionary matrix alone", only_d, 6047)
    same("difference of the two shares, points", pct(diff), "0.778")
    same("draws to resolve it, separate patterns for each matrix",
         ceil_whole(Z * Z * (p1 * (1 - p1) + p2 * (1 - p2)) / diff ** 2), 31167)
    same("draws to resolve it, the same patterns for both",
         ceil_whole(Z * Z * ((only_h + only_d) / TRIPLES - diff ** 2) / diff ** 2), 12373)
    c_h, c_d = paired_counts(run["cases"])
    dh = (c_d - c_h) / DRAWS
    half = Z * sqrt(((c_d + c_h) / DRAWS - dh * dh) / DRAWS)
    same("campaign's paired difference, points", pct(dh), "0.762")
    same("its 95% interval, points", f"{pct(dh - half)} to {pct(dh + half)}", "0.488 to 1.036")


def check_bursts() -> None:
    section("Triples of adjacent cells, 70 windows in a 72-cell word")
    for kind in ("natural", "stride5", "stride7", "stride11", "stride13", "stride29", "searched"):
        got = tuple(adjacent_outcomes(Code(name), layout(kind, name))[MISCORRECTED] for name in CODES)
        same(f"miscorrected windows, {kind} layout, three matrices", " ".join(map(str, got)),
             " ".join(map(str, STATED_ADJACENT[kind])))
    for name in ("hsiao", "hsiao_dictionary"):
        _, steps = search_layout(columns(name))
        same(f"{label(name)}: steps the layout search took", steps, STATED_SEARCH_STEPS[name])
    require("Extended Hamming: no search runs, its searched layout repeats stride 11",
            layout("searched", "hamming") == layout("stride11"))
    ham = Code("hamming")
    require("Extended Hamming: three low positions always miscorrect",
            all(outcome(ham, t) == MISCORRECTED for t in combinations(range(64), 3)))
    require("Extended Hamming: no layout tried flags more than 24 windows",
            all(adjacent_outcomes(ham, layout(k))[FLAGGED] <= 24
                for k in ("natural", "stride5", "stride7", "stride11", "stride13", "stride29")))
    for name in CODES:
        same(f"{label(name)}: expected miscorrected windows, random layout",
             f"{70 * STATED_EXACT[name][0] / TRIPLES:.2f}", STATED_EXPECTED_WINDOWS[name])


def check_designer() -> None:
    section(f"A designer's number: {WORDS} words at {RATE:g} upsets per word-hour")
    same("chance three strikes land on three different cells", f"{DISTINCT:.4f}", "0.9587")
    for name in ("hsiao", "hamming"):
        share = STATED_EXACT[name][0] / TRIPLES
        for (when, hours), stated in zip(SCRUBS, STATED_SILENT[name]):
            same(f"{label(name)}: silent corruptions a year, scrub every {when}",
                 f"{silent_per_year(share, hours):.2e}", stated)
    same("ratio of the two, any interval", f"{STATED_EXACT['hamming'][0] / STATED_EXACT['hsiao'][0]:.3f}",
         "1.350")
    same("factor a 30-day to 1-day scrub change buys, balanced",
         f"{silent_per_year(STATED_EXACT['hsiao'][0] / TRIPLES, 720) / silent_per_year(STATED_EXACT['hsiao'][0] / TRIPLES, 24):.1f}",
         "899.4")


# ---------------------------------------------------------------- the options

def option_runs(args) -> None:
    if args.stride is not None:
        from math import gcd
        if gcd(args.stride, N) != 1:
            sys.exit(f"--stride {args.stride} shares a factor with {N}, so it does not visit every position")
        section(f"Adjacent triples when cell k holds position {args.stride} k mod {N}")
        for name in CODES:
            c = adjacent_outcomes(Code(name), layout(f"stride{args.stride}"))
            print(f"  {label(name):<17} miscorrected {c[MISCORRECTED]:>3} of 70, flagged {c[FLAGGED]:>3}")
    weights = (args.weight,) if args.weight else (3, 4)
    exact = {}
    for name in CODES:
        code = Code(name)
        counts = syndrome_counts(code.cols, max(weights))
        exact[name] = {w: exact_split(code, counts, w) for w in weights}
    if args.weight:
        w = args.weight
        section(f"Every pattern of weight {w}, {comb(N, w)} of them, from the syndrome counts")
        for name in CODES:
            c = exact[name][w]
            print(f"  {label(name):<17} " + "  ".join(f"{o} {c[o]}" for o in OUTCOMES)
                  + f"  ({pct(c[tracked(w)] / comb(N, w))} percent {tracked(w)})")
    if args.stride is None or args.weight or args.seed != SEED or args.draws != DRAWS:
        section(f"Software campaign: seed {args.seed}, {args.draws} draws of each weight")
        run = campaign_lines(args.seed, args.draws, weights, exact)
        if 3 in weights:
            c_h, c_d = paired_counts(run["cases"])
            dh = (c_d - c_h) / args.draws
            half = Z * sqrt(max((c_d + c_h) / args.draws - dh * dh, 0) / args.draws)
            print(f"  dictionary minus balanced, weight 3: {pct(dh)} points, 95% interval"
                  f" {pct(dh - half)} to {pct(dh + half)}, {'excludes' if dh - half > 0 or dh + half < 0 else 'includes'} 0")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=SEED, help=f"seed of the software campaign (default {SEED})")
    parser.add_argument("--draws", type=int, default=DRAWS, help=f"draws of each weight (default {DRAWS})")
    parser.add_argument("--weight", type=int, choices=(5,), help="campaign and exact split for five-bit errors")
    parser.add_argument("--stride", type=int, help="adjacent triples when cell k holds position S k mod 72")
    args = parser.parse_args()
    if args.seed != SEED or args.draws != DRAWS or args.weight or args.stride is not None:
        option_runs(args)
        return
    if not DATA.exists():
        sys.exit(f"{DATA} is missing; run tb/miscorrection/run_miscorrection.py to measure it")
    data = json.loads(DATA.read_text())
    check_matrices(data)
    exact = check_enumeration()
    check_weight_distribution(exact)
    run = check_campaign(exact)
    check_hardware(data, run)
    check_rates(run)
    check_bursts()
    check_designer()
    print("\nAll chapter 31 numbers reproduced.")


if __name__ == "__main__":
    main()
