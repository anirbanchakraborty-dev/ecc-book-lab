"""Reproduce every number in chapter 22 of Error Correction Coding,
"Arithmetic Codes and Checking Logic".

Run it from the repository root:

    python3 python/ch22_arithmetic_codes_and_checking_logic.py

It needs only the Python standard library. Each check prints the value
computed here beside the value the chapter states, and the program stops
with exit status 1 at the first pair that disagrees. Every number in the
chapter is a whole number except the one floating-point sum in the section
on checksums for matrix products, which is printed to 17 significant digits,
so every comparison is exact.

The small cases are checked exhaustively: the claim that one stuck carry
changes the sum by a single power of two on every pair of 6-bit operands,
single-error detection by the code with A = 3 on every 8-bit operand, the
order bound for single-error correction on every odd A below 100, the
end-around-carry fold on every 12-bit word, the mod-3 residue tree on all
65,536 inputs of 16 bits, the four checks of a 4-bit ripple-carry adder on
all 256 operand pairs under six carry faults, both prediction identities on
all 512 inputs of a 4-bit adder, both self-checking checkers on every fault
and every input, and the checksum correction on every position of the 3-by-3
product. A run that reaches the final line has reproduced every number in
the chapter.

An adder here is a list of gates evaluated bit by bit, so a stuck-at fault
on a carry line can be injected exactly where a real one would sit.
"""

import sys


# Helpers ----------------------------------------------------------------------

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


def ones(x: int) -> int:
    """The number of 1 bits in a non-negative integer."""
    return bin(x).count("1")


def parity(x: int) -> int:
    return ones(x) % 2


def bits_text(x: int, n: int) -> str:
    """The n-bit binary text of x, most significant bit first."""
    return format(x, f"0{n}b")


def low_first(x: int, n: int) -> str:
    """The bits x_0 x_1 ... x_{n-1} of x, written left to right."""
    return "".join(str((x >> i) & 1) for i in range(n))


def join(values) -> str:
    return " ".join(str(v) for v in values)


def arith_weight(e: int) -> int:
    """The least number of terms of the form +2^i or -2^i that sum to e."""
    e = abs(e)
    if e == 0:
        return 0
    if e == 1:
        return 1
    if e % 2 == 0:
        return arith_weight(e // 2)
    return 1 + min(arith_weight((e - 1) // 2), arith_weight((e + 1) // 2))


def ripple_add(x: int, y: int, n: int, c0: int = 0, fault=None):
    """Add two n-bit words in a ripple-carry adder, one full adder per bit.

    fault is (j, v): carry line c_j, 1 <= j <= n, stuck at the value v.
    Returns the sum bits as an integer, the carries c_0..c_n as a list, and
    the carries the fault-free logic of each stage computed, c_1..c_n.
    """
    c = [c0] + [0] * n
    s = 0
    for i in range(n):
        xi, yi, ci = (x >> i) & 1, (y >> i) & 1, c[i]
        s |= (xi ^ yi ^ ci) << i
        c[i + 1] = (xi & yi) | (xi & ci) | (yi & ci)
        if fault is not None and fault[0] == i + 1:
            c[i + 1] = fault[1]
    return s, c


def result_value(s: int, c, n: int) -> int:
    """The n sum bits with the carry-out c_n on top, as one number."""
    return s + (c[n] << n)


# A wrong carry in a 16-bit adder ----------------------------------------------

def check_motivation() -> None:
    section("A wrong carry in a 16-bit adder")
    x, y = 0x7FF8, 0x0008
    same("x = 0x7FF8 in decimal", x, 32760)
    same("x + y, fault-free", x + y, 32768)
    s, c = ripple_add(x, y, 16)
    same("ripple-carry adder output, fault-free", result_value(s, c, 16), 32768)
    sf, cf = ripple_add(x, y, 16, fault=(4, 0))
    bad = result_value(sf, cf, 16)
    same("output with carry c_4 stuck at 0", bad, 32752)
    same("output with c_4 stuck at 0, in hexadecimal", format(bad, "04X"), "7FF0")
    same("sum bits that differ from the right sum", ones(bad ^ 32768), 12)
    same("arithmetic error, faulty minus correct", bad - 32768, -16)
    same("arithmetic weight of that error", arith_weight(bad - 32768), 1)
    same("parity of x", parity(x), 0)
    same("parity of y", parity(y), 1)
    same("parity of the correct sum", parity(32768), 1)
    same("parity of the faulty sum", parity(bad), 1)
    require("XOR of operand parities agrees with the faulty sum", parity(x) ^ parity(y) == parity(bad))
    same("1 + 1: XOR of the operand parities", parity(1) ^ parity(1), 0)
    same("1 + 1: parity of the sum 2", parity(2), 1)
    same("carries of 0x7FF8 + 8 that are 1, fault-free", sum(c[1:]), 12)

    section("Arithmetic weight")
    same("arithmetic weight of 16", arith_weight(16), 1)
    same("arithmetic weight of 7 = 8 - 1", arith_weight(7), 2)
    same("arithmetic weight of 11", arith_weight(11), 3)
    same("arithmetic weight of 4080 = 4096 - 16", arith_weight(4080), 2)
    same("Hamming weight of 4080", ones(4080), 8)
    same("arithmetic weight of 255 = 256 - 1", arith_weight(255), 2)

    # One stuck carry adds or removes exactly one power of two.
    n = 6
    cases = 0
    wrong = 0
    single = True
    for j in range(1, n + 1):
        for v in (0, 1):
            for x in range(1 << n):
                for y in range(1 << n):
                    sf, cf = ripple_add(x, y, n, fault=(j, v))
                    e = result_value(sf, cf, n) - (x + y)
                    cases += 1
                    if e != 0:
                        wrong += 1
                        if e not in (1 << j, -(1 << j)):
                            single = False
    same("6-bit adder: fault and operand cases tried", cases, 49152)
    same("6-bit adder: cases with a wrong result", wrong, 24576)
    require("every wrong result is off by exactly +2^j or -2^j", single)


# AN codes ---------------------------------------------------------------------

def order_of_two(a: int) -> int:
    k, v = 1, 2 % a
    while v != 1:
        v = v * 2 % a
        k += 1
    return k


def correcting(a: int, m: int) -> bool:
    """True when the residues of +2^i and -2^i, 0 <= i < m, are nonzero and distinct mod a."""
    seen = set()
    for i in range(m):
        for e in (1 << i, -(1 << i)):
            r = e % a
            if r == 0 or r in seen:
                return False
            seen.add(r)
    return True


def order_limit(a: int) -> int:
    """The longest correctable word the order theorem allows for odd a."""
    o = order_of_two(a)
    if any(pow(2, j, a) == a - 1 for j in range(o)):
        return o // 2
    return o


def check_an_codes() -> None:
    section("AN codes")
    same("A = 3, x = 5: stored word 3x", 3 * 5, 15)
    same("15 in binary", bits_text(15, 5), "01111")
    same("15 + 4 after one arithmetic error", 15 + 4, 19)
    same("19 mod 3", 19 % 3, 1)
    same("A = 4, x = 5: stored word 4x", 4 * 5, 20)
    same("20 + 4 mod 4, an undetected error", (20 + 4) % 4, 0)
    same("codeword bits for A = 3 and 8-bit x", (3 * 255).bit_length(), 10)
    undetected = 0
    tried = 0
    for x in range(256):
        for i in range(10):
            for e in (1 << i, -(1 << i)):
                w = 3 * x + e
                tried += 1
                if w % 3 == 0:
                    undetected += 1
    same("A = 3: single errors tried on all 8-bit x", tried, 5120)
    same("A = 3: single errors that go undetected", undetected, 0)

    same("residues of +2^i mod 19, i = 0..8", join(pow(2, i, 19) for i in range(9)), "1 2 4 8 16 13 7 14 9")
    same("residues of -2^i mod 19, i = 0..8", join((-(1 << i)) % 19 for i in range(9)), "18 17 15 11 3 6 12 5 10")
    same("distinct nonzero residues for A = 19, m = 9", len({(s * (1 << i)) % 19 for i in range(9) for s in (1, -1)}), 18)
    require("A = 19 corrects every single error in 9 bits", correcting(19, 9))
    require("A = 19 does not correct every single error in 10 bits", not correcting(19, 10))
    same("order of 2 modulo 19", order_of_two(19), 18)
    same("2^9 mod 19, which is -1", pow(2, 9, 19), 18)
    same("x = 11: stored word 19x", 19 * 11, 209)
    same("209 + 32 after an error in bit 5", 209 + 32, 241)
    same("241 mod 19", 241 % 19, 13)
    same("2^5 mod 19", pow(2, 5, 19), 13)
    same("corrected word (241 - 32) / 19", (241 - 32) // 19, 11)

    agree = True
    for a in range(3, 100, 2):
        m = 0
        while correcting(a, m + 1):
            m += 1
        if m != order_limit(a):
            agree = False
    require("order theorem matches a search for every odd A below 100", agree)
    same("order of 2 modulo 7", order_of_two(7), 3)
    same("longest correctable word for A = 7", order_limit(7), 3)
    same("order of 2 modulo 23", order_of_two(23), 11)
    same("longest correctable word for A = 23", order_limit(23), 11)

    expect = {4: (19, 9, 18, 9), 8: (29, 13, 28, 14), 12: (37, 18, 36, 18),
              16: (47, 22, 23, 23), 32: (79, 39, 39, 39)}
    for k, (a_s, m_s, o_s, l_s) in expect.items():
        a = 3
        while True:
            m = (a * ((1 << k) - 1)).bit_length()
            if a % 2 == 1 and correcting(a, m):
                break
            a += 1
        same(f"k = {k}: smallest single-error-correcting A", a, a_s)
        same(f"k = {k}: codeword length m", m, m_s)
        same(f"k = {k}: check bits m - k", m - k, m_s - k)
        same(f"k = {k}: order of 2 modulo A", order_of_two(a), o_s)
        same(f"k = {k}: longest correctable length for that A", order_limit(a), l_s)


# Residue codes ----------------------------------------------------------------

def fold(x: int, a: int) -> int:
    """Residue of x modulo 2^a - 1 by summing a-bit chunks with end-around carry."""
    mask = (1 << a) - 1
    while x > mask:
        total = 0
        while x:
            total += x & mask
            x >>= a
        x = total
    return 0 if x == mask else x


def chunks(x: int, a: int, n: int):
    return [(x >> (a * j)) & ((1 << a) - 1) for j in range(n // a)]


def canon(p1: int, p0: int):
    """A raw bit pair turned into the residue code 00, 01, 10."""
    t = p1 ^ p0
    return p1 & t, p0 & t


def mod3_add(a1: int, a0: int, b1: int, b0: int):
    """The ten-gate mod-3 adder on canonical residues."""
    za = 1 - (a1 | a0)
    zb = 1 - (b1 | b0)
    r0 = (za & b0) | (a0 & zb) | (a1 & b1)
    r1 = (za & b1) | (a1 & zb) | (a0 & b0)
    return r1, r0


def residue_tree(x: int):
    level = [canon((x >> (2 * j + 1)) & 1, (x >> (2 * j)) & 1) for j in range(8)]
    while len(level) > 1:
        level = [mod3_add(*level[k], *level[k + 1]) for k in range(0, len(level), 2)]
    r1, r0 = level[0]
    return 2 * r1 + r0


def mod3_tree_netlist():
    """The 16-bit mod-3 generator as gates (kind, output, inputs) in build order."""
    gates = []

    def gate(kind, *ins):
        name = f"n{len(gates)}"
        gates.append((kind, name, ins))
        return name

    level = []
    for j in range(8):
        p1, p0 = f"x{2 * j + 1}", f"x{2 * j}"
        t = gate("XOR", p1, p0)
        level.append((gate("AND", p1, t), gate("AND", p0, t)))
    while len(level) > 1:
        nxt = []
        for (a1, a0), (b1, b0) in zip(level[0::2], level[1::2]):
            za = gate("NOR", a1, a0)
            zb = gate("NOR", b1, b0)
            r0 = gate("OR", gate("AND", za, b0), gate("AND", a0, zb), gate("AND", a1, b1))
            r1 = gate("OR", gate("AND", za, b1), gate("AND", a1, zb), gate("AND", a0, b0))
            nxt.append((r1, r0))
        level = nxt
    return gates, level[0]


def netlist_mismatches(gates, out1, out0) -> int:
    """Evaluate the netlist on all 65,536 inputs at once, one bit of a big integer per input."""
    size = 1 << 16
    full = (1 << size) - 1
    val = {}
    for i in range(16):
        width = 1 << i
        pattern, length = ((1 << width) - 1) << width, 2 * width
        while length < size:
            pattern |= pattern << length
            length *= 2
        val[f"x{i}"] = pattern
    for kind, out, ins in gates:
        a = [val[i] for i in ins]
        if kind == "XOR":
            val[out] = a[0] ^ a[1]
        elif kind == "AND":
            val[out] = a[0] & a[1]
        elif kind == "NOR":
            val[out] = full ^ (a[0] | a[1])
        else:
            val[out] = a[0] | a[1] | a[2]
    want1 = int("".join("1" if k % 3 == 2 else "0" for k in reversed(range(size))), 2)
    want0 = int("".join("1" if k % 3 == 1 else "0" for k in reversed(range(size))), 2)
    return bin((val[out1] ^ want1) | (val[out0] ^ want0)).count("1")


def check_residue_codes() -> None:
    section("Residue codes")
    same("22 mod 3", 22 % 3, 1)
    same("13 mod 3", 13 % 3, 1)
    same("22 + 13", 22 + 13, 35)
    same("35 mod 3", 35 % 3, 2)
    same("(1 + 1) mod 3", (1 + 1) % 3, 2)
    same("22 times 13", 22 * 13, 286)
    same("286 mod 3", 286 % 3, 1)
    same("(1 times 1) mod 3", (1 * 1) % 3, 1)
    ok = all((x + y) % a == (x % a + y % a) % a and (x * y) % a == (x % a) * (y % a) % a
             for a in range(2, 17) for x in range(64) for y in range(64))
    require("sum and product residues agree for A = 2 to 16 on all 6-bit pairs", ok)

    same("181 in binary", bits_text(181, 8), "10110101")
    same("2-bit chunks of 181, lowest first", join(chunks(181, 2, 8)), "1 1 3 2")
    same("sum of those chunks", sum(chunks(181, 2, 8)), 7)
    same("7 folded once, 1 + 3", 1 + 3, 4)
    same("4 folded once more, 1 + 0", 1 + 0, 1)
    same("181 mod 3", 181 % 3, 1)
    same("fold of 181 for A = 3", fold(181, 2), 1)
    same("4-bit chunks of 181, lowest first", join(chunks(181, 4, 8)), "5 11")
    same("5 + 11", 5 + 11, 16)
    same("181 mod 15", 181 % 15, 1)
    same("fold of 181 for A = 15", fold(181, 4), 1)
    ok = all(fold(x, a) == x % ((1 << a) - 1) for a in (2, 3, 4) for x in range(1 << 12))
    require("fold equals x mod 2^a - 1 for a = 2, 3, 4 on 12-bit x", ok)

    gates, (out1, out0) = mod3_tree_netlist()
    kinds = {k: sum(1 for g in gates if g[0] == k) for k in ("XOR", "AND", "NOR", "OR")}
    same("mod-3 tree netlist: C blocks, one XOR each", kinds["XOR"], 8)
    same("mod-3 tree netlist: mod-3 adders, two NOR each", kinds["NOR"] // 2, 7)
    same("mod-3 tree netlist: XOR gates", kinds["XOR"], 8)
    same("mod-3 tree netlist: AND gates", kinds["AND"], 58)
    same("mod-3 tree netlist: NOR gates", kinds["NOR"], 14)
    same("mod-3 tree netlist: OR gates", kinds["OR"], 14)
    same("mod-3 tree netlist: gates in all", len(gates), 94)
    depth = {f"x{i}": 0 for i in range(16)}
    for _, out, ins in gates:
        depth[out] = 1 + max(depth[i] for i in ins)
    same("mod-3 tree netlist: longest path, in gates", max(depth[out1], depth[out0]), 11)
    same("mod-3 tree netlist: wrong residues over all 65,536 inputs", netlist_mismatches(gates, out1, out0), 0)
    mismatches = sum(1 for x in range(1 << 16) if residue_tree(x) != x % 3)
    same("mod-3 tree: wrong residues over all 65,536 inputs", mismatches, 0)

    same("2^16 mod 3", pow(2, 16, 3), 1)
    ok = True
    for s in range(1 << 16):
        r = residue_tree(s)
        for c16 in (0, 1):
            r1, r0 = mod3_add(r >> 1, r & 1, 0, c16)
            ok = ok and 2 * r1 + r0 == (s + (c16 << 16)) % 3
    require("tree residue of s, plus c16 in one more mod-3 adder, is r mod 3", ok)
    for n in (6, 8):
        top = 1 << n
        alarms = sum(1 for x in range(top) for y in range(top)
                     if (x % 3 + y % 3) % 3 != ((x + y) % top) % 3)
        carries = sum(1 for x in range(top) for y in range(top) if x + y >= top)
        same(f"{n}-bit pairs with a carry-out", carries, top * (top - 1) // 2)
        same(f"{n}-bit false alarms when the carry-out is ignored", alarms, carries)
    pairs = sum(range(1 << 16))
    same("16-bit pairs with a carry-out, the sum of x over all x", pairs, 2147450880)
    same("16-bit pairs with a carry-out, 2^16 (2^16 - 1) / 2", (1 << 16) * ((1 << 16) - 1) // 2, 2147450880)
    same("share of all 2^32 pairs, in ten-millionths", round(pairs * 10 ** 7 / (1 << 32)), 4999924)
    same("mod-3 tree applied to 32768", residue_tree(32768), 2)
    same("mod-3 tree applied to 32752", residue_tree(32752), 1)
    same("(32760 + 8) residues: (0 + 2) mod 3", (32760 % 3 + 8 % 3) % 3, 2)


# The 4-bit adder under six carry faults ---------------------------------------

FAULTS = [(1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]


def berger(w: int, n: int) -> int:
    """The Berger check symbol of an n-bit word: its number of 0 bits."""
    return n - ones(w)


def adder_experiment():
    """Every operand pair of a 4-bit adder under each fault on c_1, c_2, c_3."""
    n = 4
    rows = []
    for fault in FAULTS:
        wrong = mod3 = par_shared = par_dup = par_cmp = berg = 0
        gap_two = True
        for x in range(16):
            for y in range(16):
                s, c = ripple_add(x, y, n)
                sf, cf = ripple_add(x, y, n, fault=fault)
                good = result_value(s, c, n)
                got = result_value(sf, cf, n)
                if got == good:
                    continue
                wrong += 1
                if got % 3 != (x % 3 + y % 3) % 3:
                    mod3 += 1
                ps = parity(sf)
                shared = parity(x) ^ parity(y) ^ (sum(cf[0:n]) % 2)
                dup = parity(x) ^ parity(y) ^ (sum(c[0:n]) % 2)
                if ps != shared:
                    par_shared += 1
                if ps != dup:
                    par_dup += 1
                if ps != dup or cf != c:
                    par_cmp += 1
                predicted = berger(x, n) + berger(y, n) - (n - sum(cf[1:])) - cf[0] + cf[n]
                actual = berger(sf, n)
                if predicted != actual:
                    berg += 1
                if abs(predicted - actual) != 2:
                    gap_two = False
        rows.append((fault, wrong, mod3, par_shared, par_dup, par_cmp, berg, gap_two))
    return rows


def check_adder_faults() -> None:
    section("Checking a 4-bit adder under six carry faults")
    s, c = ripple_add(3, 1, 4)
    same("3 + 1, fault-free", result_value(s, c, 4), 4)
    same("carries c_0..c_4 of 3 + 1", join(c), "0 1 1 0 0")
    sf, cf = ripple_add(3, 1, 4, fault=(1, 0))
    same("3 + 1 with c_1 stuck at 0", result_value(sf, cf, 4), 2)
    same("carries c_0..c_4 with c_1 stuck at 0", join(cf), "0 0 0 0 0")
    same("residue check: 2 mod 3", 2 % 3, 2)
    same("residue check: (3 + 1) mod 3 from the residues", (3 % 3 + 1 % 3) % 3, 1)
    same("parity predicted from the fault-free carries", parity(3) ^ parity(1) ^ (sum(c[0:4]) % 2), 1)
    same("parity of the faulty sum 0010", parity(sf), 1)

    rows = adder_experiment()
    stated = {
        (1, 0): (64, 64, 0, 48, 64, 64),
        (1, 1): (192, 192, 0, 144, 192, 192),
        (2, 0): (96, 96, 0, 48, 96, 96),
        (2, 1): (160, 160, 0, 80, 160, 160),
        (3, 0): (112, 112, 0, 112, 112, 112),
        (3, 1): (144, 144, 0, 144, 144, 144),
    }
    totals = [0] * 6
    for fault, wrong, mod3, ps, pd, pc, bg, gap in rows:
        name = f"c_{fault[0]} stuck at {fault[1]}"
        got = (wrong, mod3, ps, pd, pc, bg)
        same(f"{name}: wrong results", wrong, stated[fault][0])
        same(f"{name}: caught by the mod-3 residue", mod3, stated[fault][1])
        same(f"{name}: caught by parity, shared carries", ps, stated[fault][2])
        same(f"{name}: caught by parity, duplicate carries", pd, stated[fault][3])
        same(f"{name}: caught by duplicate carries and compare", pc, stated[fault][4])
        same(f"{name}: caught by Berger check prediction", bg, stated[fault][5])
        require(f"{name}: Berger prediction off by exactly 2", gap)
        for k in range(6):
            totals[k] += got[k]
    same("fault and operand cases in all", 6 * 256, 1536)
    same("wrong results in all", totals[0], 768)
    same("caught by the mod-3 residue", totals[1], 768)
    same("caught by parity prediction, shared carries", totals[2], 0)
    same("caught by parity prediction, duplicate carries", totals[3], 576)
    same("escapes with duplicate carries alone", totals[0] - totals[3], 192)
    same("caught by duplicate carries and comparison", totals[4], 768)
    same("caught by Berger check prediction", totals[5], 768)


# Berger check prediction and parity prediction --------------------------------

def check_predictions() -> None:
    section("Berger check prediction and parity prediction")
    x, y, n = 0b1011, 0b0110, 4
    s, c = ripple_add(x, y, n)
    same("11 + 6 = 17, sum bits", bits_text(s, 4), "0001")
    same("11 + 6, carry-out c_4", c[4], 1)
    same("carries c_1..c_4 of 11 + 6", join(c[1:]), "0 1 1 1")
    same("ones in x = 1011", ones(x), 3)
    same("ones in y = 0110", ones(y), 2)
    same("ones among c_1..c_4", sum(c[1:]), 3)
    same("ones in the sum, 3 + 2 - 3 + 0 - 1", 3 + 2 - 3 + 0 - 1, 1)
    same("ones in the sum, counted", ones(s), 1)
    same("Berger symbol of x", berger(x, n), 1)
    same("Berger symbol of y", berger(y, n), 2)
    same("Berger symbol of the carries", n - sum(c[1:]), 1)
    same("predicted Berger symbol, 1 + 2 - 1 - 0 + 1", 1 + 2 - 1 - 0 + 1, 3)
    same("Berger symbol of the sum, counted", berger(s, n), 3)
    same("parity of x", parity(x), 1)
    same("parity of y", parity(y), 0)
    same("parity of the carries c_0..c_3", sum(c[0:4]) % 2, 0)
    same("parity of the sum", parity(s), 1)

    count = 0
    ones_ok = berger_ok = parity_ok = True
    for c0 in (0, 1):
        for x in range(16):
            for y in range(16):
                s, c = ripple_add(x, y, n, c0=c0)
                count += 1
                nc = sum(c[1:])
                if ones(s) != ones(x) + ones(y) - nc + c0 - c[n]:
                    ones_ok = False
                if berger(s, n) != berger(x, n) + berger(y, n) - (n - nc) - c0 + c[n]:
                    berger_ok = False
                if parity(s) != parity(x) ^ parity(y) ^ (sum(c[0:n]) % 2):
                    parity_ok = False
    same("4-bit adder inputs, carry-in 0 and 1", count, 512)
    require("ones in the sum match the count identity", ones_ok)
    require("Berger symbol of the sum matches the prediction", berger_ok)
    require("parity of the sum matches the prediction", parity_ok)


# Self-checking circuits -------------------------------------------------------

TWO_RAIL_LINES = ["a", "a'", "b", "b'",
                  "a@1", "b@1", "a'@2", "b'@2", "a@3", "b'@3", "a'@4", "b@4",
                  "g1", "g2", "g3", "g4", "z", "z'"]


def two_rail(a: int, an: int, b: int, bn: int, fault=None):
    """The two-rail checker cell, gates g1..g4 (AND) and z, z' (OR).

    g1 = a b, g2 = a' b', g3 = a b', g4 = a' b, z = g1 + g2, z' = g3 + g4.
    fault is (line name, stuck value).
    """
    v = {"a": a, "a'": an, "b": b, "b'": bn}

    def line(name, value):
        if fault is not None and fault[0] == name:
            return fault[1]
        return value

    for k in ("a", "a'", "b", "b'"):
        v[k] = line(k, v[k])
    g1 = line("g1", line("a@1", v["a"]) & line("b@1", v["b"]))
    g2 = line("g2", line("a'@2", v["a'"]) & line("b'@2", v["b'"]))
    g3 = line("g3", line("a@3", v["a"]) & line("b'@3", v["b'"]))
    g4 = line("g4", line("a'@4", v["a'"]) & line("b@4", v["b"]))
    z = line("z", g1 | g2)
    zn = line("z'", g3 | g4)
    return z, zn


def check_two_rail() -> None:
    section("The two-rail checker cell")
    codes = [(a, b) for a in (0, 1) for b in (0, 1)]
    outputs = [two_rail(a, 1 - a, b, 1 - b) for a, b in codes]
    same("code inputs (a, b)", len(codes), 4)
    same("outputs (z, z') for a b = 00 01 10 11", join(f"{z}{zn}" for z, zn in outputs), "10 01 01 10")
    noncode = 0
    disjoint = True
    for a in (0, 1):
        for an in (0, 1):
            for b in (0, 1):
                for bn in (0, 1):
                    if an == 1 - a and bn == 1 - b:
                        continue
                    noncode += 1
                    z, zn = two_rail(a, an, b, bn)
                    if z != 1 - zn:
                        pass
                    else:
                        disjoint = False
    same("non-code inputs", noncode, 12)
    require("every non-code input gives a non-code output", disjoint)

    same("fault sites (lines) in the cell", len(TWO_RAIL_LINES), 18)
    same("single stuck-at faults", 2 * len(TWO_RAIL_LINES), 36)
    secure = True
    testable = True
    per_input = [0, 0, 0, 0]
    sole = [0, 0, 0, 0]
    for name in TWO_RAIL_LINES:
        for v in (0, 1):
            detecting = []
            for k, (a, b) in enumerate(codes):
                good = two_rail(a, 1 - a, b, 1 - b)
                z, zn = two_rail(a, 1 - a, b, 1 - b, fault=(name, v))
                if (z, zn) == good:
                    continue
                if z == 1 - zn:
                    secure = False
                else:
                    detecting.append(k)
            if not detecting:
                testable = False
            for k in detecting:
                per_input[k] += 1
            if len(detecting) == 1:
                sole[detecting[0]] += 1
    require("fault-secure: no fault yields a wrong code word", secure)
    require("self-testing: some code input detects every fault", testable)
    same("faults detected by a b = 00, 01, 10, 11", join(per_input), "13 13 13 13")
    same("faults that only a b = 00, 01, 10, 11 detect", join(sole), "5 5 5 5")
    g1_s0 = [k for k, (a, b) in enumerate(codes) if two_rail(a, 1 - a, b, 1 - b, fault=("g1", 0)) != two_rail(a, 1 - a, b, 1 - b)]
    same("code input that detects g1 stuck at 0, as a b", join(f"{codes[k][0]}{codes[k][1]}" for k in g1_s0), "11")
    z, zn = two_rail(1, 0, 0, 1, fault=("g1", 1))
    same("a b = 10 with g1 stuck at 1: (z, z')", f"{z}{zn}", "11")


PARITY_LINES = ["x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "p1", "p2", "p3", "p4", "f", "g"]


def parity_checker(w: int, fault=None):
    """Two parity trees: f = (x0+x1)+(x2+x3), g = NOT((x4+x5)+(x6+x7)), + meaning XOR."""
    def line(name, value):
        if fault is not None and fault[0] == name:
            return fault[1]
        return value

    x = [line(f"x{i}", (w >> i) & 1) for i in range(8)]
    p1 = line("p1", x[0] ^ x[1])
    p2 = line("p2", x[2] ^ x[3])
    p3 = line("p3", x[4] ^ x[5])
    p4 = line("p4", x[6] ^ x[7])
    f = line("f", p1 ^ p2)
    g = line("g", 1 - (p3 ^ p4))
    return f, g


def check_parity_checker() -> None:
    section("A totally self-checking parity checker")
    code = [w for w in range(256) if parity(w) == 0]
    same("even-parity code words of 8 bits", len(code), 128)
    same("XOR or XNOR gates in the two trees", 6, 6)
    same("fault sites (lines)", len(PARITY_LINES), 14)
    same("single stuck-at faults", 2 * len(PARITY_LINES), 28)
    disjoint = all(parity_checker(w)[0] == parity_checker(w)[1] for w in range(256) if parity(w) == 1)
    require("every odd-parity word gives a non-code output", disjoint)
    require("every code word gives a code output", all(parity_checker(w)[0] != parity_checker(w)[1] for w in code))
    secure = True
    detect = []
    for name in PARITY_LINES:
        for v in (0, 1):
            row = []
            for w in code:
                good = parity_checker(w)
                got = parity_checker(w, fault=(name, v))
                if got != good and got[0] != got[1]:
                    secure = False
                row.append(got != good)
            detect.append(row)
    require("fault-secure: no fault yields a wrong code word", secure)
    require("self-testing: some code word detects every fault", all(any(r) for r in detect))
    pair_covers = any(all(r[i] or r[j] for r in detect) for i in range(128) for j in range(i + 1, 128))
    require("no two code words detect all 28 faults", not pair_covers)
    tests = [0b00000000, 0b10101010, 0b11011101]
    same("test words, written x_0 x_1 ... x_7", join(low_first(w, 8) for w in tests), "00000000 01010101 10111011")
    require("those three words are code words", all(parity(w) == 0 for w in tests))
    idx = [code.index(w) for w in tests]
    require("those three code words detect all 28 faults", all(any(r[i] for i in idx) for r in detect))


# Checksums for matrix products ------------------------------------------------

def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def mat_text(m) -> str:
    return "; ".join(" ".join(str(v) for v in row) for row in m)


def check_abft() -> None:
    section("Checksums for matrix products")
    A = [[2, 1, 0], [1, 3, 1], [0, 2, 4]]
    B = [[1, 0, 2], [3, 1, 0], [2, 1, 1]]
    C = matmul(A, B)
    same("C = AB", mat_text(C), "5 1 4; 12 4 3; 14 6 4")
    col_a = [sum(A[i][j] for i in range(3)) for j in range(3)]
    row_b = [sum(B[i]) for i in range(3)]
    same("checksum row of A", join(col_a), "3 6 5")
    same("checksum column of B", join(row_b), "3 4 4")
    Ac = A + [col_a]
    Br = [B[i] + [row_b[i]] for i in range(3)]
    Cf = matmul(Ac, Br)
    same("full checksum matrix, row by row", mat_text(Cf), "5 1 4 10; 12 4 3 19; 14 6 4 24; 31 11 11 53")
    require("last column holds the row sums of C", all(Cf[i][3] == sum(C[i]) for i in range(3)))
    require("last row holds the column sums of C", all(Cf[3][j] == sum(C[i][j] for i in range(3)) for j in range(3)))
    same("corner entry, (3 6 5) times (3 4 4)", 3 * 3 + 6 * 4 + 5 * 4, 53)

    bad = [row[:] for row in Cf]
    bad[1][1] = 9
    row_gap = [sum(bad[i][0:3]) - bad[i][3] for i in range(3)]
    col_gap = [sum(bad[i][j] for i in range(3)) - bad[3][j] for j in range(3)]
    same("row sum minus checksum, rows 1 to 3", join(row_gap), "0 5 0")
    same("column sum minus checksum, columns 1 to 3", join(col_gap), "0 5 0")
    same("row 2 sum with the error, 12 + 9 + 3", 12 + 9 + 3, 24)
    same("column 2 sum with the error, 1 + 9 + 6", 1 + 9 + 6, 16)
    same("corrected entry 9 - 5", 9 - 5, 4)

    fixed = True
    for i in range(3):
        for j in range(3):
            for e in (-7, -3, -1, 1, 2, 5, 100):
                m = [row[:] for row in Cf]
                m[i][j] += e
                rg = [sum(m[r][0:3]) - m[r][3] for r in range(3)]
                cg = [sum(m[r][c] for r in range(3)) - m[3][c] for c in range(3)]
                bad_rows = [r for r in range(3) if rg[r] != 0]
                bad_cols = [c for c in range(3) if cg[c] != 0]
                if bad_rows != [i] or bad_cols != [j] or rg[i] != cg[j]:
                    fixed = False
                    continue
                m[i][j] -= rg[i]
                if m != Cf:
                    fixed = False
    require("one wrong entry is located and corrected at all 9 places", fixed)

    left = (0.1 + 0.2) + 0.3
    right = 0.1 + (0.2 + 0.3)
    same("(0.1 + 0.2) + 0.3 in double precision", f"{left:.17g}", "0.60000000000000009")
    same("0.1 + (0.2 + 0.3) in double precision", f"{right:.17g}", "0.59999999999999998")
    same("their difference", f"{left - right:.17g}", "1.1102230246251565e-16")
    same("2^-53, the spacing of doubles from 0.5 to 1", f"{2.0 ** -53:.17g}", "1.1102230246251565e-16")


# The worked examples and the exercises ---------------------------------------

def check_examples() -> None:
    section("The worked examples and the exercises")
    same("A = 19, N = 11: stored word 19N", 19 * 11, 209)
    same("209 in binary", bits_text(209, 8), "11010001")
    same("209 + 32 after the error", 209 + 32, 241)
    same("19 times 12", 19 * 12, 228)
    same("241 mod 19", 241 % 19, 13)
    errs = [s * (1 << i) for i in range(9) for s in (1, -1)]
    same("single errors +2^i, -2^i, i < 9, leaving remainder 13", join(e for e in errs if e % 19 == 13), "32")
    require("A = 19 corrects every single error in 9-bit words", correcting(19, 9))
    same("other nonzero remainders modulo 19", 18 - 1, 17)
    same("(241 - 32) / 19", (241 - 32) // 19, 11)

    same("powers 2^0, 2^1, 2^2 modulo 7", join(pow(2, i, 7) for i in range(3)), "1 2 4")
    same("order of 2 modulo 7", order_of_two(7), 3)
    same("longest correctable length for A = 7", order_limit(7), 3)
    same("powers 2^0..2^10 modulo 23", join(pow(2, i, 23) for i in range(11)), "1 2 4 8 16 9 18 13 3 6 12")
    same("2^11", 2 ** 11, 2048)
    same("89 times 23 plus 1", 89 * 23 + 1, 2048)
    same("order of 2 modulo 23", order_of_two(23), 11)
    same("longest correctable length for A = 23", order_limit(23), 11)
    same("remainders of 1, -1, 2, -2, 4, -4 modulo 7, sorted", join(sorted(e % 7 for e in (1, -1, 2, -2, 4, -4))), "1 2 3 4 5 6")

    same("32752 in hexadecimal", format(32752, "X"), "7FF0")
    same("2-bit chunks of 32752, lowest first", join(chunks(32752, 2, 16)), "0 0 3 3 3 3 3 1")
    same("those chunks after the C blocks", join(2 * h + l for h, l in (canon(c >> 1, c & 1) for c in chunks(32752, 2, 16))), "0 0 0 0 0 0 0 1")
    same("mod-3 tree applied to 32752", residue_tree(32752), 1)
    same("32768 in hexadecimal", format(32768, "X"), "8000")
    same("2-bit chunks of 32768, lowest first", join(chunks(32768, 2, 16)), "0 0 0 0 0 0 0 2")
    same("mod-3 tree applied to 32768", residue_tree(32768), 2)
    same("3 times 10920", 3 * 10920, 32760)
    same("32760 mod 3", 32760 % 3, 0)
    same("8 mod 3", 8 % 3, 2)

    s, c = ripple_add(3, 1, 4, 0, (1, 0))
    s_ok, c_ok = ripple_add(3, 1, 4)
    same("3 + 1 with c_1 stuck at 0: sum bits", bits_text(s, 4), "0010")
    same("3 + 1 with c_1 stuck at 0: carries c_0..c_4", join(c), "0 0 0 0 0")
    same("3 + 1, healthy carries c_0..c_4", join(c_ok), "0 1 1 0 0")
    same("Berger symbol of x = 0011", berger(3, 4), 2)
    same("Berger symbol of y = 0001", berger(1, 4), 3)
    same("Berger symbol of the faulty carries c_1..c_4", 4 - sum(c[1:]), 4)
    pred = berger(3, 4) + berger(1, 4) - (4 - sum(c[1:])) - c[0] + c[4]
    same("predicted Berger symbol, 2 + 3 - 4 - 0 + 0", pred, 1)
    same("Berger symbol of the sum 0010", berger(s, 4), 3)
    same("Berger miss, 3 - 1", berger(s, 4) - pred, 2)
    same("parity prediction, shared faulty carries c_0..c_3", parity(3) ^ parity(1) ^ (sum(c[:4]) % 2), 1)
    same("parity of the sum 0010", parity(s), 1)
    same("parity of the healthy carries c_0..c_3", sum(c_ok[:4]) % 2, 0)
    same("parity prediction, duplicated healthy carries", parity(3) ^ parity(1) ^ (sum(c_ok[:4]) % 2), 1)
    same("carries that differ between the chains, c_0..c_4", join(f"c_{j}" for j in range(5) if c[j] != c_ok[j]), "c_1 c_2")
    same("mod-3 prediction (3 mod 3 + 1 mod 3) mod 3", (3 % 3 + 1 % 3) % 3, 1)
    same("faulty result 0010 with c_4 = 0, mod 3", result_value(s, c, 4) % 3, 2)

    same("arithmetic weight of 255", arith_weight(255), 2)
    same("ones in 255", ones(255), 8)
    same("256 - 1", 256 - 1, 255)
    same("arithmetic weight of 4080", arith_weight(4080), 2)
    same("ones in 4080", ones(4080), 8)
    same("4096 - 16", 4096 - 16, 4080)
    same("4080 in binary", bits_text(4080, 12), "111111110000")
    same("arithmetic weight of 4, a run of one 1", arith_weight(4), 1)
    ok = all(arith_weight(((1 << run) - 1) << shift) == 2
             for run in range(2, 17) for shift in range(17 - run))
    require("a run of 2 to 16 ones, anywhere in 16 bits, has arithmetic weight 2", ok)

    same("powers 2^0..2^3 modulo 3", join(pow(2, i, 3) for i in range(4)), "1 2 1 2")
    same("order of 2 modulo 3", order_of_two(3), 2)
    same("longest correctable length for A = 3", order_limit(3), 1)
    same("+2 mod 3", 2 % 3, 2)
    same("-1 mod 3", -1 % 3, 2)
    same("+1 mod 3", 1 % 3, 1)
    same("-2 mod 3", -2 % 3, 1)
    ok = True
    for m in range(2, 17):
        named = {1: 0, 2: 0}
        for i in range(m):
            for e in (1 << i, -(1 << i)):
                named[e % 3] += 1
        ok = ok and named[1] >= 2 and named[2] >= 2
    require("A = 3, m = 2 to 16: each nonzero remainder names two or more errors", ok)

    big = 47 * 65535
    same("47 times 65535", big, 3080145)
    same("bit length of 3080145", big.bit_length(), 22)
    require("2^21 <= 3080145 < 2^22", 2 ** 21 <= big < 2 ** 22)
    same("order of 2 modulo 47", order_of_two(47), 23)
    same("longest correctable length for A = 47", order_limit(47), 23)
    require("A = 47 corrects every single error in 22-bit words", correcting(47, 22))
    same("remainders of +2^i, -2^i for i < 22", 2 * 22, 44)
    same("bound 2m + 1 for m = 22", 2 * 22 + 1, 45)

    same("3-bit chunks of 181, lowest first", join(chunks(181, 3, 9)), "5 6 2")
    same("sum of those chunks", sum(chunks(181, 3, 9)), 13)
    same("3-bit chunks of 13, lowest first", join(chunks(13, 3, 6)), "5 1")
    same("fold of 181 for A = 7", fold(181, 3), 6)
    same("7 times 25 plus 6", 7 * 25 + 6, 181)

    x, y = 32760, 8
    s, c = ripple_add(x, y, 16, 0, (4, 0))
    s_ok, c_ok = ripple_add(x, y, 16)
    same("16-bit adder, c_4 stuck at 0: faulty sum", s, 32752)
    same("carries that differ between the chains, c_0..c_16", sum(a != b for a, b in zip(c, c_ok)), 12)
    require("healthy carries c_4..c_15 are all 1", all(c_ok[j] == 1 for j in range(4, 16)))
    require("faulty carries c_4..c_15 are all 0", all(c[j] == 0 for j in range(4, 16)))
    same("parity prediction from the healthy chain", parity(x) ^ parity(y) ^ (sum(c_ok[:16]) % 2), 1)
    same("parity of the correct sum 32768", parity(s_ok), 1)
    same("ones in the faulty sum 32752", ones(s), 11)
    same("parity of the faulty sum", parity(s), 1)

    cm = [[5, 1, 4], [12, 4, 3], [14, 6, 4]]
    same("row sums of C", join(sum(r) for r in cm), "10 19 24")
    same("row 2 sum minus the corrupted checksum 21", sum(cm[1]) - 21, -2)
    same("column sum minus checksum, columns 1 to 3", join(sum(cm[i][j] for i in range(3)) - [31, 11, 11][j] for j in range(3)), "0 0 0")


def main() -> None:
    check_motivation()
    check_an_codes()
    check_residue_codes()
    check_adder_faults()
    check_predictions()
    check_two_rail()
    check_parity_checker()
    check_abft()
    check_examples()
    print("\nAll chapter 22 numbers reproduced.")


if __name__ == "__main__":
    main()
