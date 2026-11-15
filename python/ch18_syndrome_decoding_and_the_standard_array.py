"""Reproduce every number in chapter 18 of Error Correction Coding, Syndrome Decoding and the Standard Array.

    python3 python/ch18_syndrome_decoding_and_the_standard_array.py

Each check prints the value computed here beside the value the chapter states,
and the script stops with a non-zero exit at the first disagreement. Python 3.10,
standard library only. The Octave twin in octave/ prints the same lines.

Two codes carry the chapter. C6 is the (6,3) code with generator rows 100110,
010011 and 001111. The (7,4) Hamming code is the one of Linear Codes, with generator
rows 1000110, 0100011, 0010101 and 0001111 and parity-check columns 110, 011, 101,
111, 100, 010 and 001.
"""

import sys
from math import gcd
from itertools import combinations, product


def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


def parse(text: str) -> list:
    return [tuple(int(b) for b in chunk) for chunk in text.split()]


def word(v) -> str:
    return "".join(str(b) for b in v)


def words(vs) -> str:
    return " ".join(word(v) for v in vs)


def f6(x: float) -> str:
    return f"{x:.6f}"


def xor(a, b) -> tuple:
    return tuple((x + y) % 2 for x, y in zip(a, b))


def wt(v) -> int:
    return sum(v)


def unit(n: int, j: int) -> tuple:
    return tuple(1 if i == j else 0 for i in range(n))


def all_words(n: int) -> list:
    """Every binary word of length n in increasing binary value."""
    return [tuple(v) for v in product((0, 1), repeat=n)]


def codewords(G: list) -> list:
    """The codewords u G for the messages u in increasing binary value."""
    n = len(G[0])
    out = []
    for u in product((0, 1), repeat=len(G)):
        c = tuple(0 for _ in range(n))
        for bit, row in zip(u, G):
            if bit:
                c = xor(c, row)
        out.append(c)
    return out


def syndrome(H: list, v) -> tuple:
    """s = v H^T, one bit per row of H."""
    return tuple(sum(h * x for h, x in zip(row, v)) % 2 for row in H)


def columns(H: list) -> list:
    return [tuple(row[j] for row in H) for j in range(len(H[0]))]


def standard_array(code: list, n: int):
    """Slepian's layout, built as the chapter's definition states. The first row is the
    code. Each later row is led by a word of least weight among those not yet written,
    ties going to the smaller binary value, and holds that leader XORed with each
    codeword in the order of the first row."""
    listed = set(code)
    rows, leaders = [list(code)], [code[0]]
    order = sorted(all_words(n), key=lambda v: (wt(v), v))
    for v in order:
        if v in listed:
            continue
        leader = v
        row = [xor(leader, c) for c in code]
        rows.append(row)
        leaders.append(leader)
        listed.update(row)
    return rows, leaders


def syndrome_table(H: list, leaders: list) -> dict:
    return {syndrome(H, e): e for e in leaders}


def decode(H: list, table: dict, r) -> tuple:
    return xor(r, table[syndrome(H, r)])


def likelihood(p: float, n: int, w: int) -> float:
    return p ** w * (1 - p) ** (n - w)


def leader_counts(leaders: list, n: int) -> list:
    counts = [0] * (n + 1)
    for e in leaders:
        counts[wt(e)] += 1
    return counts


def p_correct(counts: list, n: int, p: float) -> float:
    return sum(a * likelihood(p, n, i) for i, a in enumerate(counts))


def p_correct_by_enumeration(H: list, table: dict, code: list, n: int, p: float) -> float:
    """Sum the probabilities of the error patterns after which the decoder returns the
    codeword that was sent, with the all-zeros codeword sent."""
    zero = code[0]
    return sum(likelihood(p, n, wt(e)) for e in all_words(n) if decode(H, table, xor(zero, e)) == zero)


def erasure_solutions(H: list, r: list, erased: list) -> list:
    """All codewords agreeing with r off the erased positions, found by trying every
    filling of the erased positions and keeping those whose syndrome is zero."""
    out = []
    for fill in product((0, 1), repeat=len(erased)):
        v = list(r)
        for pos, bit in zip(erased, fill):
            v[pos] = bit
        if syndrome(H, v) == tuple(0 for _ in H):
            out.append(tuple(v))
    return out


def independent(cols: list) -> bool:
    """True when no nonempty subset of the columns adds to zero."""
    m = len(cols[0])
    for size in range(1, len(cols) + 1):
        for subset in combinations(cols, size):
            total = tuple(0 for _ in range(m))
            for c in subset:
                total = xor(total, c)
            if wt(total) == 0:
                return False
    return True


G6 = parse("100110 010011 001111")
H6 = parse("101100 111010 011001")
C6 = codewords(G6)
G7 = parse("1000110 0100011 0010101 0001111")
H7 = parse("1011100 1101010 0111001")
C7 = codewords(G7)


def decoder_netlist(H: list) -> list:
    """The AND-gate decoder as a list of gates (kind, output, inputs), in an order that evaluates."""
    m, n = len(H), len(H[0])
    gates = []
    for i in range(m):
        level = [f"r{j + 1}" for j in range(n) if H[i][j]]
        k = 0
        while len(level) > 1:
            nxt = []
            for a in range(0, len(level) - 1, 2):
                k += 1
                out = f"s{i + 1}" if len(level) == 2 else f"t{i + 1}{k}"
                gates.append(("XOR", out, (level[a], level[a + 1])))
                nxt.append(out)
            if len(level) % 2:
                nxt.append(level[-1])
            level = nxt
    for i in range(m):
        gates.append(("NOT", f"n{i + 1}", (f"s{i + 1}",)))
    for j, col in enumerate(columns(H)):
        gates.append(("AND", f"e{j + 1}", tuple(f"s{i + 1}" if col[i] else f"n{i + 1}" for i in range(m))))
    for j in range(n):
        gates.append(("XOR", f"c{j + 1}", (f"r{j + 1}", f"e{j + 1}")))
    return gates


def run_netlist(gates: list, r) -> tuple:
    """Evaluate the netlist on the received word r and return the corrected word and each signal's depth."""
    value = {f"r{j + 1}": b for j, b in enumerate(r)}
    depth = {name: 0 for name in value}
    for kind, out, ins in gates:
        bits = [value[x] for x in ins]
        value[out] = (bits[0] ^ bits[1]) if kind == "XOR" else (1 - bits[0]) if kind == "NOT" else int(all(bits))
        depth[out] = 1 + max(depth[x] for x in ins)
    return tuple(value[f"c{j + 1}"] for j in range(len(r))), depth


def check_c6() -> None:
    section("The (6,3) code C6")
    same("G", words(G6), "100110 010011 001111")
    same("H = [P^T | I]", words(H6), "101100 111010 011001")
    require("every row of G passes every check of H", all(wt(syndrome(H6, g)) == 0 for g in G6))
    same("codewords for messages 000 to 111", words(C6),
         "000000 001111 010011 011100 100110 101001 110101 111010")
    same("weights of the codewords", " ".join(str(wt(c)) for c in C6), "0 4 3 3 3 3 4 4")
    same("minimum distance of C6", min(wt(c) for c in C6[1:]), 3)


def check_opening() -> None:
    section("Which codeword was sent")
    c = C6[6]
    same("codeword for message 110", word(c), "110101")
    r = xor(c, unit(6, 2))
    same("received word after an error in position 3", word(r), "111101")
    dist = [wt(xor(r, x)) for x in C6]
    same("distances from 111101 to the eight codewords", " ".join(map(str, dist)), "5 3 4 2 4 2 1 3")
    best = min(dist)
    same("nearest codeword", word(C6[dist.index(best)]), "110101")
    same("codewords at the least distance", dist.count(best), 1)
    for w, stated in [(1, "0.059049"), (2, "0.006561"), (3, "0.000729"), (4, "0.000081"), (5, "0.000009")]:
        same(f"likelihood at p = 0.1 of a codeword at distance {w}", f6(likelihood(0.1, 6, w)), stated)
    same("likelihood ratio, distance 1 over distance 2", f6(likelihood(0.1, 6, 1) / likelihood(0.1, 6, 2)), "9.000000")
    same("total likelihood over the eight codewords at p = 0.1",
         f6(sum(likelihood(0.1, 6, d) for d in dist)), "0.073800")
    same("share of that total held by 110101",
         f6(likelihood(0.1, 6, 1) / sum(likelihood(0.1, 6, d) for d in dist)), "0.800122")


def check_ml_cases() -> None:
    section("Maximum likelihood when p is not below one half")
    r = parse("111101")[0]
    dist = [wt(xor(r, x)) for x in C6]
    same("likelihood of every codeword at p = 1/2", f6(likelihood(0.5, 6, 3)), "0.015625")
    require("at p = 1/2 all eight likelihoods are equal",
            len({likelihood(0.5, 6, d) for d in dist}) == 1)
    scores = [likelihood(0.9, 6, d) for d in dist]
    same("ML codeword at p = 0.9", word(C6[scores.index(max(scores))]), "000000")
    same("its distance from 111101", dist[scores.index(max(scores))], 5)
    same("its likelihood at p = 0.9", f6(max(scores)), "0.059049")
    flipped = tuple(1 - b for b in r)
    same("111101 with every bit flipped", word(flipped), "000010")
    fd = [wt(xor(flipped, x)) for x in C6]
    same("nearest codeword to the flipped word", word(C6[fd.index(min(fd))]), "000000")


def check_standard_array() -> None:
    section("The standard array of C6")
    rows, leaders = standard_array(C6, 6)
    stated = [
        "000000 001111 010011 011100 100110 101001 110101 111010",
        "000001 001110 010010 011101 100111 101000 110100 111011",
        "000010 001101 010001 011110 100100 101011 110111 111000",
        "000100 001011 010111 011000 100010 101101 110001 111110",
        "001000 000111 011011 010100 101110 100001 111101 110010",
        "010000 011111 000011 001100 110110 111001 100101 101010",
        "100000 101111 110011 111100 000110 001001 010101 011010",
        "000101 001010 010110 011001 100011 101100 110000 111111",
    ]
    same("number of rows", len(rows), 8)
    for i, (row, text) in enumerate(zip(rows, stated)):
        same(f"row {i + 1}", words(row), text)
    flat = [v for row in rows for v in row]
    same("words in the array", len(flat), 64)
    require("every word of length 6 appears exactly once", sorted(flat) == all_words(6))
    same("leader weights", " ".join(str(wt(e)) for e in leaders), "0 1 1 1 1 1 1 2")
    require("each leader has the least weight in its own row",
            all(wt(e) == min(wt(v) for v in row) for e, row in zip(leaders, rows)))
    last = rows[7]
    same("least-weight words of the last coset", words(sorted(v for v in last if wt(v) == 2)),
         "000101 001010 110000")
    r = parse("111101")[0]
    i = next(k for k, row in enumerate(rows) if r in row)
    same("row holding 111101", i + 1, 5)
    same("codeword at the head of its column", word(C6[rows[i].index(r)]), "110101")
    require("every row has 8 words and the rows are disjoint", all(len(set(row)) == 8 for row in rows)
            and len(set(flat)) == 64)


def check_leader_decoding() -> None:
    section("Decoding with the array")
    rows, leaders = standard_array(C6, 6)
    w = parse("000101")[0]
    dist = [wt(xor(w, x)) for x in C6]
    same("distances from 000101 to the codewords", " ".join(map(str, dist)), "2 2 3 3 3 3 2 6")
    same("codewords at distance 2 from 000101", words([C6[k] for k in range(8) if dist[k] == 2]),
         "000000 001111 110101")
    same("leader decoding of 000101", word(xor(w, leaders[7])), "000000")
    same("decision region of 000000", words(leaders),
         "000000 000001 000010 000100 001000 010000 100000 000101")
    p = 0.1
    total = sum(likelihood(p, 6, d) for d in dist)
    same("total likelihood of the codewords given 000101 at p = 0.1", f6(total), "0.022600")
    same("share of that total held by each tied codeword", f6(likelihood(p, 6, 2) / total), "0.290310")
    require("each tied codeword holds less than 1/3 at p = k/1000 for k = 1 to 499",
            all(likelihood(k / 1000, 6, 2) / sum(likelihood(k / 1000, 6, d) for d in dist) < 1 / 3
                for k in range(1, 500)))
    region = [xor(C6[0], e) for e in leaders]
    same("weights of the words in the decision region of 000000",
         " ".join(str(wt(v)) for v in region), "0 1 1 1 1 1 1 2")
    for r in all_words(6):
        out = xor(r, leaders[next(k for k, row in enumerate(rows) if r in row)])
        if out not in C6 or wt(xor(r, out)) != min(wt(xor(r, x)) for x in C6):
            require("leader decoding returns a nearest codeword for every word", False)
    require("leader decoding returns a nearest codeword for every word", True)


def check_syndromes_c6() -> None:
    section("Syndromes of C6")
    rows, leaders = standard_array(C6, 6)
    same("columns of H", words(columns(H6)), "110 011 111 100 010 001")
    same("syndrome of 111101", word(syndrome(H6, parse("111101")[0])), "111")
    same("syndromes of the leaders, row by row", words(syndrome(H6, e) for e in leaders),
         "000 001 010 100 111 011 110 101")
    require("all words of a row share one syndrome",
            all(len({syndrome(H6, v) for v in row}) == 1 for row in rows))
    require("different rows have different syndromes",
            len({syndrome(H6, row[0]) for row in rows}) == 8)
    require("s(x) = s(y) exactly when x + y is a codeword",
            all((syndrome(H6, x) == syndrome(H6, y)) == (xor(x, y) in C6)
                for x in all_words(6) for y in all_words(6)))
    table = syndrome_table(H6, leaders)
    stated = {"000": "000000", "001": "000001", "010": "000010", "011": "010000",
              "100": "000100", "101": "000101", "110": "100000", "111": "001000"}
    for s in all_words(3):
        same(f"C6 table, syndrome {word(s)}", word(table[s]), stated[word(s)])
    same("syndrome missing from the columns of H", word(sorted(set(all_words(3)) - set(columns(H6))
                                                               - {(0, 0, 0)})[0]), "101")


def check_hamming() -> None:
    section("The (7,4) Hamming code")
    same("P", words(row[4:] for row in G7), "110 011 101 111")
    same("G = [I | P]", words(G7), "1000110 0100011 0010101 0001111")
    same("H = [P^T | I]", words(H7), "1011100 1101010 0111001")
    require("every row of G passes every check of H", all(wt(syndrome(H7, g)) == 0 for g in G7))
    same("columns of H", words(columns(H7)), "110 011 101 111 100 010 001")
    require("the seven columns are the seven nonzero words of length 3",
            sorted(columns(H7)) == all_words(3)[1:])
    same("codewords", len(C7), 16)
    counts = [0] * 8
    for c in C7:
        counts[wt(c)] += 1
    same("codewords of weight 0 to 7", " ".join(map(str, counts)), "1 0 0 7 7 0 0 1")
    rows, leaders = standard_array(C7, 7)
    same("rows of its standard array", len(rows), 8)
    require("its standard array lists all 128 words once",
            sorted(v for row in rows for v in row) == all_words(7))
    table = syndrome_table(H7, leaders)
    stated = {"000": "0000000", "001": "0000001", "010": "0000010", "011": "0100000",
              "100": "0000100", "101": "0010000", "110": "1000000", "111": "0001000"}
    for s in all_words(3):
        same(f"(7,4) table, syndrome {word(s)}", word(table[s]), stated[word(s)])


def check_lookup() -> None:
    section("Decoding the (7,4) code by lookup")
    _, leaders = standard_array(C7, 7)
    table = syndrome_table(H7, leaders)
    c = C7[int("1001", 2)]
    same("codeword for message 1001", word(c), "1001001")
    same("it is the first row of G plus the fourth", word(xor(G7[0], G7[3])), "1001001")
    r = xor(c, unit(7, 5))
    same("received word, error in position 6", word(r), "1001011")
    same("its syndrome", word(syndrome(H7, r)), "010")
    same("leader read from the table", word(table[syndrome(H7, r)]), "0000010")
    same("decoded codeword", word(decode(H7, table, r)), "1001001")
    e = parse("0100100")[0]
    r2 = xor(c, e)
    same("received word, errors in positions 2 and 5", word(r2), "1101101")
    same("its syndrome", word(syndrome(H7, r2)), "111")
    same("column 2 plus column 5 of H", word(xor(columns(H7)[1], columns(H7)[4])), "111")
    same("column of H equal to that syndrome", columns(H7).index((1, 1, 1)) + 1, 4)
    same("decoded codeword", word(decode(H7, table, r2)), "1100101")
    same("positions where it differs from the sent word", wt(xor(decode(H7, table, r2), c)), 3)
    same("message bits delivered", word(decode(H7, table, r2)[:4]), "1100")
    cl = columns(H7)
    require("(7,4): the XOR of two distinct columns is always a third column",
            all(xor(cl[i], cl[j]) in cl and xor(cl[i], cl[j]) not in (cl[i], cl[j])
                for i, j in combinations(range(7), 2)))
    r4 = parse("100100")[0]
    _, l6 = standard_array(C6, 6)
    t6 = syndrome_table(H6, l6)
    same("C6: syndrome of 100100", word(syndrome(H6, r4)), "010")
    same("C6: decoding of 100100", word(decode(H6, t6, r4)), "100110")
    r5 = parse("001010")[0]
    same("C6: syndrome of 001010", word(syndrome(H6, r5)), "101")
    same("C6: decoding of 001010", word(decode(H6, t6, r5)), "001111")
    d5 = [wt(xor(r5, x)) for x in C6]
    same("C6: codewords at distance 2 from 001010", words([C6[k] for k in range(8) if d5[k] == 2]),
         "000000 001111 111010")


def check_sizes() -> None:
    section("Table sizes")
    same("entries in the (7,4) syndrome table", 2 ** 3, 8)
    same("bits in the (7,4) syndrome table", 2 ** 3 * 7, 56)
    same("bits in a ROM addressed by the received word", 2 ** 7 * 4, 512)
    same("bits in the whole standard array", 2 ** 7 * 7, 896)
    same("entries for 72 bits with 8 check bits", 2 ** 8, 256)
    same("entries in a table addressed by 72 bits", f"{2.0 ** 72:.2e}", "4.72e+21")


def check_hardware() -> None:
    section("Syndrome hardware")
    same("inputs of s1", " ".join(str(j + 1) for j in range(7) if H7[0][j]), "1 3 4 5")
    same("inputs of s2", " ".join(str(j + 1) for j in range(7) if H7[1][j]), "1 2 4 6")
    same("inputs of s3", " ".join(str(j + 1) for j in range(7) if H7[2][j]), "2 3 4 7")
    gates = sum(wt(row) - 1 for row in H7)
    same("two-input XOR gates in the syndrome trees", gates, 9)
    same("depth of each tree in gates", max((wt(row) - 1).bit_length() for row in H7), 2)
    _, leaders = standard_array(C7, 7)
    table = syndrome_table(H7, leaders)
    same("ROM contents, addresses 0 to 7", words(table[s] for s in all_words(3)),
         "0000000 0000001 0000010 0100000 0000100 0010000 1000000 0001000")
    r = parse("1001011")[0]
    same("tree inputs for r = 1001011, s1", " ".join(str(r[j]) for j in range(7) if H7[0][j]), "1 0 1 0")
    same("tree inputs for r = 1001011, s2", " ".join(str(r[j]) for j in range(7) if H7[1][j]), "1 0 1 1")
    same("tree inputs for r = 1001011, s3", " ".join(str(r[j]) for j in range(7) if H7[2][j]), "0 0 1 1")
    same("first-level XORs, s1", f"{r[0] ^ r[2]} {r[3] ^ r[4]}", "1 1")
    same("first-level XORs, s2", f"{r[0] ^ r[1]} {r[3] ^ r[5]}", "1 0")
    same("first-level XORs, s3", f"{r[1] ^ r[2]} {r[3] ^ r[6]}", "0 0")
    same("AND gate taps, columns 1 to 7", words(columns(H7)), "110 011 101 111 100 010 001")
    ok = True
    for s in all_words(3):
        for j, col in enumerate(columns(H7)):
            fires = all((s[i] == 1) == (col[i] == 1) for i in range(3))
            ok = ok and (fires == (table[s][j] == 1))
    require("seven AND gates on s or its complement reproduce the ROM", ok)
    net = decoder_netlist(H7)
    tree_xors = sum(1 for g in net if g[0] == "XOR" and not g[1].startswith("c"))
    same("netlist: XOR gates in the syndrome trees", tree_xors, 9)
    same("netlist: inverters for the complements", sum(1 for g in net if g[0] == "NOT"), 3)
    same("netlist: three-input AND gates", sum(1 for g in net if g[0] == "AND" and len(g[2]) == 3), 7)
    same("netlist: XOR gates that apply the correction", sum(1 for g in net if g[1].startswith("c")), 7)
    same("netlist: gates in the whole decoder, inverters aside", sum(1 for g in net if g[0] != "NOT"), 23)
    _, depth = run_netlist(net, (0,) * 7)
    same("netlist: depth of s1, s2, s3 in gates", f"{depth['s1']} {depth['s2']} {depth['s3']}", "2 2 2")
    require("netlist: the gates decode all 128 received words as the table does",
            all(run_netlist(net, v)[0] == decode(H7, table, v) for v in all_words(7)))
    same("netlist: output for r = 1001011", word(run_netlist(net, r)[0]), "1001001")


def check_erasures() -> None:
    section("Erasures")
    c = parse("1001001")[0]
    erased = [0, 5]
    known = [b if k not in erased else 0 for k, b in enumerate(c)]
    same("received word with positions 1 and 6 erased", "".join("?" if k in erased else str(b)
                                                                 for k, b in enumerate(c)), "?0010?1")
    same("known part u with zeros in the gaps", word(known), "0001001")
    same("syndrome of the known part, positions 1 and 6 erased", word(syndrome(H7, known)), "110")
    cols = columns(H7)
    same("columns 4 and 7 of H added", word(xor(cols[3], cols[6])), "110")
    same("columns 1 and 6 of H", words([cols[0], cols[5]]), "110 010")
    aug = [[cols[0][i], cols[5][i], syndrome(H7, known)[i]] for i in range(3)]
    same("augmented matrix [h1 h6 | s]", words(aug), "101 111 000")
    aug[1] = list(xor(aug[1], aug[0]))
    same("after R2 <- R2 + R1", words(aug), "101 010 000")
    same("codewords that fit", words(erasure_solutions(H7, c, erased)), "1001001")
    erased3 = [0, 1, 2]
    same("columns 1, 2 and 3 added together", word(xor(xor(cols[0], cols[1]), cols[2])), "000")
    require("1110000 is a codeword", parse("1110000")[0] in C7)
    k3 = [b if k not in erased3 else 0 for k, b in enumerate(c)]
    same("known part with positions 1, 2 and 3 erased", word(k3), "0001001")
    same("its syndrome", word(syndrome(H7, k3)), "110")
    sols = erasure_solutions(H7, c, erased3)
    same("codewords that fit", words(sols), "0111001 1001001")
    same("their erased bits", " ".join(word(v[:3]) for v in sols), "011 100")
    same("the two candidates differ by", word(xor(sols[0], sols[1])), "1110000")
    for n, H, code, label in [(7, H7, C7, "(7,4)"), (6, H6, C6, "C6")]:
        cl = columns(H)
        stated = {(7, 1): "7 of 7", (7, 2): "21 of 21", (7, 3): "28 of 35", (7, 4): "0 of 35",
                  (6, 1): "6 of 6", (6, 2): "15 of 15", (6, 3): "16 of 20", (6, 4): "0 of 15"}
        for e in range(1, 5):
            sets = list(combinations(range(n), e))
            good = sum(1 for E in sets if independent([cl[j] for j in E]))
            same(f"{label}: sets of {e} erasures with independent columns", f"{good} of {len(sets)}",
                 stated[(n, e)])
        dep = [E for E in combinations(range(n), 3) if not independent([cl[j] for j in E])]
        supports = [tuple(j for j in range(n) if v[j]) for v in code if wt(v) == 3]
        require(f"{label}: the dependent triples are the supports of the weight-3 codewords",
                sorted(dep) == sorted(supports))
        same(f"{label}: supports of the weight-3 codewords",
             " ".join("{" + ",".join(str(j + 1) for j in E) + "}" for E in sorted(supports)),
             "{1,2,3} {1,4,7} {1,5,6} {2,4,5} {2,6,7} {3,4,6} {3,5,7}" if n == 7 else
             "{1,3,6} {1,4,5} {2,3,4} {2,5,6}")
        ok = all(len(erasure_solutions(H, v, list(E))) == 2 for v in code for E in combinations(range(n), 4))
        require(f"{label}: every set of 4 erasures leaves exactly 2 candidates", ok)
    ok = all(erasure_solutions(H7, v, list(E)) == [v]
             for v in C7 for e in (1, 2) for E in combinations(range(7), e))
    require("(7,4): every codeword is recovered from every 1 or 2 erasures", ok)
    c2 = C7[int("0110", 2)]
    same("codeword for message 0110", word(c2), "0110110")
    k2 = [b if k not in (1, 6) else 0 for k, b in enumerate(c2)]
    same("received word 0?1011? with zeros in the gaps", word(k2), "0010110")
    same("syndrome of the known part of 0?1011?", word(syndrome(H7, k2)), "011")
    same("columns 3, 5 and 6 added", word(xor(xor(cols[2], cols[4]), cols[5])), "011")
    same("codewords that fit 0?1011?", words(erasure_solutions(H7, c2, [1, 6])), "0110110")
    same("columns 2, 3, 5 and 6 added", word(xor(xor(cols[1], cols[2]), xor(cols[4], cols[5]))), "000")


def check_failure() -> None:
    section("Decoding failure")
    _, l6 = standard_array(C6, 6)
    _, l7 = standard_array(C7, 7)
    a6, a7 = leader_counts(l6, 6), leader_counts(l7, 7)
    same("C6 leaders of weight 0, 1, 2", " ".join(map(str, a6[:3])), "1 6 1")
    same("(7,4) leaders of weight 0, 1", " ".join(map(str, a7[:2])), "1 7")
    stated = {0.01: ("0.998636", "0.970299", "0.997969", "0.960596"),
              0.05: ("0.969262", "0.857375", "0.955619", "0.814506"),
              0.1: ("0.892296", "0.729000", "0.850306", "0.656100")}
    for p, (s6, u3, s7, u4) in stated.items():
        same(f"p = {p}: C6", f6(p_correct(a6, 6, p)), s6)
        same(f"p = {p}: three bits sent bare", f6((1 - p) ** 3), u3)
        same(f"p = {p}: (7,4) Hamming code", f6(p_correct(a7, 7, p)), s7)
        same(f"p = {p}: four bits sent bare", f6((1 - p) ** 4), u4)
    t6, t7 = syndrome_table(H6, l6), syndrome_table(H7, l7)
    require("C6 sums agree with a run over all 64 error patterns",
            all(abs(p_correct(a6, 6, p) - p_correct_by_enumeration(H6, t6, C6, 6, p)) < 1e-12 for p in stated))
    require("(7,4) sums agree with a run over all 128 error patterns",
            all(abs(p_correct(a7, 7, p) - p_correct_by_enumeration(H7, t7, C7, 7, p)) < 1e-12 for p in stated))
    ok = True
    for code, H, table, n in [(C6, H6, t6, 6), (C7, H7, t7, 7)]:
        for c in code:
            for e in all_words(n):
                ok = ok and ((decode(H, table, xor(c, e)) == c) == (e in table.values()))
    require("the decoder is right exactly when the error is a leader", ok)
    same("(7,4) over four bare bits at p = 1/2", f6(p_correct(a7, 7, 0.5) / 0.5 ** 4), "1.000000")
    same("factored form (1-p)^6 (1+6p) at p = 0.1", f6(0.9 ** 6 * 1.6), "0.850306")
    same("p where the ratio (1-p)^2 (1+6p) peaks", f6(2 / 9), "0.222222")
    same("largest ratio", f6((7 / 9) ** 2 * (1 + 12 / 9)), "1.411523")
    num, den = 7 ** 2 * (9 + 12), 9 ** 2 * 9
    g = gcd(num, den)
    same("largest ratio before reduction", f"{num}/{den}", "1029/729")
    same("its factors (7/9)^2 and 1 + 12/9", f"{7 ** 2}/{9 ** 2} {9 + 12}/9", "49/81 21/9")
    same("largest ratio as a fraction in lowest terms", f"{num // g}/{den // g}", "343/243")
    same("root of 4 - 18p in lowest terms", f"{4 // gcd(4, 18)}/{18 // gcd(4, 18)}", "2/9")
    same("percent gain at the peak", f"{((7 / 9) ** 2 * (1 + 12 / 9) - 1) * 100:.0f}", "41")
    require("the (7,4) code beats four bare bits at p = k/1000 for k = 1 to 499",
            all(p_correct(a7, 7, k / 1000) > (1 - k / 1000) ** 4 for k in range(1, 500)))
    same("(7,4) block error at p = 0.01", f"{1 - p_correct(a7, 7, 0.01):.5f}", "0.00203")
    same("C6 at p = 1/2", f6(p_correct(a6, 6, 0.5)), "0.125000")
    same("(7,4) at p = 1/2", f6(p_correct(a7, 7, 0.5)), "0.062500")
    p = 0.1
    s101 = [e for e in all_words(6) if syndrome(H6, e) == (1, 0, 1)]
    flag = sum(likelihood(p, 6, wt(e)) for e in s101)
    bounded = a6[0] * likelihood(p, 6, 0) + a6[1] * likelihood(p, 6, 1)
    same("no bit flips at p = 0.1", f6(likelihood(p, 6, 0)), "0.531441")
    same("one of six bits flips at p = 0.1", f6(6 * likelihood(p, 6, 1)), "0.354294")
    same("words of coset 101 of weight 2, 3 and 6",
         " ".join(str(sum(1 for e in s101 if wt(e) == w)) for w in (2, 3, 6)), "3 4 1")
    same("the pattern of weight 6 at p = 0.1", f6(likelihood(p, 6, 6)), "0.000001")
    same("bounded decoder at p = 0.1: correct", f6(bounded), "0.885735")
    same("bounded decoder at p = 0.1: failure flagged", f6(flag), "0.022600")
    same("bounded decoder at p = 0.1: wrong codeword", f6(1 - bounded - flag), "0.091665")
    same("complete decoder at p = 0.1: wrong codeword", f6(1 - p_correct(a6, 6, p)), "0.107704")
    same("syndrome 101: right guesses", f6(likelihood(p, 6, 2)), "0.006561")
    same("syndrome 101: wrong guesses", f6(flag - likelihood(p, 6, 2)), "0.016039")
    same("C6 at p = 0.1: error pattern is a nonzero codeword",
         f6(sum(likelihood(p, 6, wt(c)) for c in C6[1:])), "0.003159")


def check_exercises() -> None:
    section("Worked examples and exercises")
    _, l6 = standard_array(C6, 6)
    t6 = syndrome_table(H6, l6)
    r = parse("011011")[0]
    same("C6: syndrome of 011011", word(syndrome(H6, r)), "111")
    same("C6: decoding of 011011", word(decode(H6, t6, r)), "010011")
    same("ratio p/(1-p) at p = 0.2", f6(0.2 / 0.8), "0.250000")
    same("ratio p/(1-p) at p = 0.25", f6(0.25 / 0.75), "0.333333")
    R3 = [(0, 0, 0), (1, 1, 1)]
    H3 = parse("110 101")
    rows, leaders = standard_array(R3, 3)
    same("repetition code: leaders", words(leaders), "000 001 010 100")
    same("repetition code: their syndromes", words(syndrome(H3, e) for e in leaders), "00 01 10 11")
    same("repetition code: P(correct) at p = 0.1", f6(p_correct(leader_counts(leaders, 3), 3, 0.1)), "0.972000")
    same("repetition code: its two terms at p = 0.1", f"{0.9 ** 3:.3f} {3 * 0.1 * 0.9 ** 2:.3f}", "0.729 0.243")
    _, l7 = standard_array(C7, 7)
    a7 = leader_counts(l7, 7)
    same("(7,4) P(correct) at p = 0.02", f6(p_correct(a7, 7, 0.02)), "0.992143")
    same("(7,4) block error at p = 0.02", f6(1 - p_correct(a7, 7, 0.02)), "0.007857")
    same("four bits sent bare at p = 0.02", f6(0.98 ** 4), "0.922368")
    same("1 + 6p at p = 0.02", f"{1 + 6 * 0.02:.2f}", "1.12")
    c = parse("101001")[0]
    known = [b if k not in (1, 2) else 0 for k, b in enumerate(c)]
    same("C6: syndrome of the known part of 1??001", word(syndrome(H6, known)), "111")
    same("C6: codewords that fit 1??001", words(erasure_solutions(H6, c, [1, 2])), "101001")
    cols = columns(H7)
    same("(7,4): columns 3, 4 and 5 added", word(xor(xor(cols[2], cols[3]), cols[4])), "110")
    require("(7,4): columns 3, 4, 5 are independent", independent([cols[2], cols[3], cols[4]]))
    same("(7,4): columns 2 and 4 added", word(xor(cols[1], cols[3])), "100")
    require("(7,4): columns 2, 4, 5 are dependent", not independent([cols[1], cols[3], cols[4]]))
    same("(7,4): codeword on positions 2, 4, 5", word(C7[int("0101", 2)]), "0101100")
    same("syndrome table bits, length 15 with 4 checks", 2 ** 4 * 15, 240)
    same("ROM bits addressed by 15 received bits", 2 ** 15 * 11, 360448)
    same("C6 P(correct) at p = 0.05 (worked example)", f6(p_correct(leader_counts(l6, 6), 6, 0.05)), "0.969262")
    same("C6 at p = 0.05: the leader of weight 0", f6(0.95 ** 6), "0.735092")
    same("C6 at p = 0.05: the six leaders of weight 1", f6(6 * 0.05 * 0.95 ** 5), "0.232134")
    same("C6 at p = 0.05: the leader of weight 2", f6(0.05 ** 2 * 0.95 ** 4), "0.002036")
    same("0.95^5 and 0.05^2", f"{0.95 ** 5:.6f} {0.05 ** 2:.4f}", "0.773781 0.0025")
    same("three bits sent bare at p = 0.05: lost", f6(1 - 0.95 ** 3), "0.142625")
    same("C6 at p = 0.05: lost", f6(1 - p_correct(leader_counts(l6, 6), 6, 0.05)), "0.030738")
    same("four bits sent bare at p = 0.02: lost", f6(1 - 0.98 ** 4), "0.077632")
    same("ratio of the two failure rates at p = 0.02", f"{(1 - 0.98 ** 4) / (1 - p_correct(a7, 7, 0.02)):.2f}",
         "9.88")
    same("ROM bits over table bits", f"{2 ** 15 * 11 / (2 ** 4 * 15):.1f}", "1501.9")
    same("four bits sent bare at p = 0.01: block error", f6(1 - 0.99 ** 4), "0.039404")
    require("(1-p)^2 (1+6p) - 1 = p (2p-1) (3p-4) at p = k/1000 for k = 1 to 999",
            all(abs((1 - k / 1000) ** 2 * (1 + 6 * k / 1000) - 1
                    - (k / 1000) * (2 * k / 1000 - 1) * (3 * k / 1000 - 4)) < 1e-12 for k in range(1, 1000)))


def main() -> None:
    print("Chapter 18: Syndrome Decoding and the Standard Array")
    check_c6()
    check_opening()
    check_ml_cases()
    check_standard_array()
    check_leader_decoding()
    check_syndromes_c6()
    check_hamming()
    check_lookup()
    check_sizes()
    check_hardware()
    check_erasures()
    check_failure()
    check_exercises()
    print("\nAll chapter 18 numbers reproduced.")


if __name__ == "__main__":
    main()
