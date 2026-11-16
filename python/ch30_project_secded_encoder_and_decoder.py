"""Reproduce every number in chapter 30 of Error Correction Coding,
"Project: A SEC-DED Encoder and Decoder".

Run it from the repository root, or from anywhere, since every path is found
from this file:

    python3 python/ch30_project_secded_encoder_and_decoder.py

It needs Python 3.10 or later and nothing beyond the standard library, and it
has an Octave twin, octave/ch30_project_secded_encoder_and_decoder.m, that
prints the same lines.

Each check prints the value computed here beside the value the chapter
states, and the program stops with exit status 1 at the first pair that
disagrees. The check matrices are rebuilt from the constructions of
ch20_secded_and_hsiao_codes.py and compared with the generated package
hdl/secded72/secded72_pkg.sv, column for column. The flag counts come from
every error pattern of weight at most two and from every syndrome, not from a
sample. The measured gate counts and depths are read from
data/ch30/secded72_gates.json, which synth/synth_secded72.py writes from a
real Yosys netlist, and every percentage the chapter quotes is computed here
from that file. The bench model in tb/secded72/secded72_model.py is checked
against the same constructions without printing, so a disagreement there
stops the program too.

A column of a check matrix with r rows is held as a whole number whose bit i
is the entry in row i, and a bit pattern is printed with row 0 first, as in
"SECDED and Hsiao Codes". Positions are counted from 0.
"""

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent
PACKAGE = LAB / "hdl" / "secded72" / "secded72_pkg.sv"
GATES = LAB / "data" / "ch30" / "secded72_gates.json"

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LAB / "tb" / "secded72"))
import ch20_secded_and_hsiao_codes as ch20  # noqa: E402
import secded72_model as model  # noqa: E402

RANDOM_WORDS = 1000


# Helpers ----------------------------------------------------------------------

def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if str(computed) != str(stated):
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def quietly(label: str, condition: bool) -> None:
    """A check that prints nothing unless it fails."""
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


def weight(v: int) -> int:
    return bin(v).count("1")


def bits(v: int, r: int) -> str:
    return "".join(str((v >> i) & 1) for i in range(r))


def joined(values) -> str:
    return ",".join(str(v) for v in values)


def levels(w: int) -> int:
    """Levels of two-input gates in a balanced tree over w signals."""
    d, reach = 0, 1
    while reach < w:
        reach *= 2
        d += 1
    return d


def pct(part: int, whole: int) -> str:
    return f"{100 * part / whole:.1f}"


def row_loads(cols: list, r: int) -> list:
    return [sum((c >> i) & 1 for c in cols) for i in range(r)]


def package_columns(prefix: str) -> list:
    """The column list a generated package holds, position 0 first.

    Every code is stored in 72 entries of 8 bits, so a shorter code ends in
    entries that are zero, and those are returned too.
    """
    text = PACKAGE.read_text()
    block = re.search(r"localparam logic \[[^\]]*\] " + prefix + r"_COLUMNS = \{(.*?)\};", text, re.S)
    words = re.findall(r"8'b([01]+)", block.group(1))
    return [int(w, 2) for w in reversed(words)]


def check_fan_in(cols: list, r: int, checks: list) -> list:
    """Data bits in each check equation, in the order of the check positions.

    Every syndrome value is a sum of check columns in exactly one way, since
    the check columns are independent, so a table over all 2^r subsets of them
    solves each data column for the checks it involves.
    """
    table = {}
    for subset in range(1 << r):
        v = 0
        for j in range(r):
            if (subset >> j) & 1:
                v ^= cols[checks[j]]
        table[v] = subset
    data = [p for p in range(len(cols)) if p not in checks]
    return [sum((table[cols[p]] >> j) & 1 for p in data) for j in range(r)]


def parity_ok(name: str, s: int, r: int) -> bool:
    """The test the parity form of the flag applies: odd weight, or the parity row set."""
    if name == "hamming":
        return (s >> (r - 1)) & 1 == 1
    return weight(s) % 2 == 1


# The codes --------------------------------------------------------------------

CODES = [
    # (name, package prefix, construction named in the label, columns, r, check positions)
    ("hsiao", "HSIAO", "hsiao_72_64()", ch20.hsiao_72_64(), 8, list(range(64, 72))),
    ("hamming", "HAMMING", "extended_hamming(72, 7)", ch20.extended_hamming(72, 7), 8,
     [0] + [1 << i for i in range(7)]),
    ("hsiao_dictionary", "HSIAO_DICT", "hsiao_72_64_dictionary()", ch20.hsiao_72_64_dictionary(), 8,
     list(range(64, 72))),
    ("hsiao_39_32", "HSIAO_39_32", "hsiao_39_32()", ch20.hsiao_39_32(), 7, list(range(32, 39))),
]

STATED = {
    "hsiao": dict(ones=216, rows="27,27,27,27,27,27,27,27", largest=27, levels=5, sep=208,
                  fan="26,26,26,26,26,26,26,26", enc=200, enc_levels=5,
                  matched=72, nomatch=183, parity=127, only=56,
                  hand0=855, hand1=862, depth="9,15,16,9",
                  syn_whole=169, syn_whole_depth="5,5,5,5,5,5,5,5", enc_whole=159, enc_whole_depth=5,
                  dec0="550,9,15,16", dec1="457,9,15,10", cells0="139,154,257,39"),
    "hamming": dict(ones=284, rows="36,36,36,32,32,32,8,72", largest=72, levels=7, sep=276,
                    fan="35,35,35,35,31,31,31,7", enc=232, enc_levels=6,
                    matched=72, nomatch=183, parity=127, only=56,
                    hand0=923, hand1=922, depth="11,17,18,10",
                    syn_whole=137, syn_whole_depth="6,6,6,5,5,5,3,7", enc_whole=175, enc_whole_depth=6,
                    dec0="478,11,17,18", dec1="385,11,17,10", cells0="129,142,207,48"),
    "hsiao_dictionary": dict(ones=216, rows="30,30,30,26,26,25,25,24", largest=30, levels=5, sep=208,
                             fan="29,29,29,25,25,24,24,23", enc=200, enc_levels=5,
                             matched=72, nomatch=183, parity=127, only=56,
                             hand0=855, hand1=862, depth="9,15,16,9",
                             syn_whole=164, syn_whole_depth="5,5,5,5,5,5,5,5", enc_whole=156, enc_whole_depth=5,
                             dec0="549,9,15,16", dec1="457,9,15,10", cells0="139,153,257,39"),
    "hsiao_39_32": dict(ones=103, rows="14,14,15,15,15,15,15", largest=15, levels=4, sep=96,
                        fan="13,13,14,14,14,14,14", enc=89, enc_levels=4,
                        matched=39, nomatch=88, parity=63, only=25,
                        hand0=407, hand1=413, depth="8,13,14,8",
                        syn_whole=80, syn_whole_depth="4,4,4,4,4,4,4", enc_whole=73, enc_whole_depth=4,
                        dec0="303,8,13,14", dec1="251,8,13,9", cells0="80,91,132,23"),
}

LABELS = {"hsiao": "Hsiao", "hamming": "Hamming", "hsiao_dictionary": "Dictionary",
          "hsiao_39_32": "Hsiao (39,32)"}


def main() -> None:
    print("Chapter 30, Project: A SEC-DED Encoder and Decoder")

    section("The check matrices, rebuilt from chapter 20 and read from the package")
    fan_ins = {}
    for name, prefix, construction, cols, r, checks in CODES:
        label = LABELS[name]
        st = STATED[name]
        code = model.build(name)
        quietly(f"{name}: the bench model's columns equal {construction}", code.columns == cols)
        quietly(f"{name}: the bench model's check positions", code.check_positions == checks)
        stored = package_columns(prefix)
        require(f"{label}: package columns equal {construction}",
                stored[:len(cols)] == cols and not any(stored[len(cols):]))
        loads = row_loads(cols, r)
        same(f"{label}: ones in the matrix", sum(loads), st["ones"])
        same(f"{label}: ones in each row", joined(loads), st["rows"])
        same(f"{label}: largest row", max(loads), st["largest"])
        same(f"{label}: XOR levels for the largest row", levels(max(loads)), st["levels"])
        same(f"{label}: syndrome XORs as separate trees", sum(w - 1 for w in loads), st["sep"])
        fan = check_fan_in(cols, r, checks)
        quietly(f"{name}: the bench model's check equations",
                [weight(m) for m in code.check_masks] == fan)
        fan_ins[name] = fan
        same(f"{label}: data bits in each check equation", joined(fan), st["fan"])
        same(f"{label}: encoder XORs as separate trees", sum(f - 1 for f in fan), st["enc"])
        same(f"{label}: encoder levels", levels(max(fan)), st["enc_levels"])

    section("Two flags, compared on every syndrome and every pattern of weight <= 2")
    for name, _, _, cols, r, _ in CODES:
        label = LABELS[name]
        st = STATED[name]
        n = len(cols)
        column_set = set(cols)
        if name == "hamming":
            require(f"{label}: every column has the parity row set", all(parity_ok(name, c, r) for c in cols))
        else:
            require(f"{label}: every column has odd weight", all(parity_ok(name, c, r) for c in cols))
        require(f"{label}: columns distinct and nonzero", len(column_set) == n and 0 not in column_set)
        matched = nomatch = parity = only = 0
        parity_alone = False
        for s in range(1 << r):
            m = s in column_set
            no_match = s != 0 and not m
            par = s != 0 and not parity_ok(name, s, r)
            matched += m
            nomatch += no_match
            parity += par
            only += no_match and not par
            parity_alone = parity_alone or (par and not no_match)
        same(f"{label}: syndromes that match a column", matched, st["matched"])
        same(f"{label}: syndromes that set the no-match flag", nomatch, st["nomatch"])
        same(f"{label}: syndromes that set the parity flag", parity, st["parity"])
        same(f"{label}: syndromes that set the no-match flag alone", only, st["only"])
        require(f"{label}: the parity flag never fires alone", not parity_alone)
        singles = sum(1 for c in cols if c in column_set)
        doubles = disagree = 0
        for i in range(n):
            for j in range(i + 1, n):
                s = cols[i] ^ cols[j]
                no_match = s != 0 and s not in column_set
                par = s != 0 and not parity_ok(name, s, r)
                doubles += no_match
                disagree += no_match != par
        same(f"{label}: single errors whose syndrome is a column", singles, n)
        same(f"{label}: double errors that set the no-match flag", doubles, n * (n - 1) // 2)
        same(f"{label}: patterns of weight <= 2 where the flags disagree", disagree, 0)

    section("What the bench runs")
    for name, _, _, cols, r, checks in CODES:
        label = LABELS[name]
        n = len(cols)
        k = n - r
        patterns = 1 + n + n * (n - 1) // 2
        words = 2 + k + RANDOM_WORDS
        quietly(f"{name}: the bench model's envelope", len(model.envelope_patterns(n)) == patterns)
        quietly(f"{name}: the bench model's words", len(model.bench_words(k)) == words)
        stated_patterns = {72: 2629, 39: 781}[n]
        stated_words = {64: 1066, 32: 1034}[k]
        stated_cases = {72: 2802514, 39: 807554}[n]
        same(f"{label}: error patterns of weight <= 2", patterns, stated_patterns)
        same(f"{label}: data words", words, stated_words)
        same(f"{label}: cases", words * patterns, stated_cases)
        if n == 72:
            same(f"{label}: cases of weight 0, 1 and 2",
                 joined([words, words * n, words * n * (n - 1) // 2]), "1066,76752,2724696")

    section("The decoder counted by hand, before any synthesis")
    for name, _, _, cols, r, _ in CODES:
        label = LABELS[name]
        st = STATED[name]
        n = len(cols)
        k = n - r
        loads = row_loads(cols, r)
        syndrome = sum(w - 1 for w in loads)
        comparators = n * (r - 1)
        base = syndrome + comparators + k + (n - 1) + 1
        if name == "hamming":
            # nonzero over the r - 1 rows other than the parity row
            flag0 = base + (r - 1)
            flag1 = base + (r - 2)
        else:
            flag0 = base + (r - 1)
            flag1 = base + (r - 1) + (r - 1)
        same(f"{label}: comparator AND gates", comparators, {72: 504, 39: 234}[n])
        same(f"{label}: gates, no-match flag", flag0, st["hand0"])
        same(f"{label}: gates, parity flag", flag1, st["hand1"])
        ls = levels(max(loads))
        data = ls + levels(r) + 1
        corrected = ls + levels(r) + levels(n)
        unc0 = corrected + 1
        if name == "hamming":
            low = levels(max(loads[:r - 1]))
            unc1 = max(low + levels(r - 1), levels(loads[r - 1])) + 1
        else:
            unc1 = ls + levels(r) + 1
        same(f"{label}: depth of data, corrected, both flags", joined([data, corrected, unc0, unc1]),
             st["depth"])

    section("Measured in the Yosys netlist, data/ch30/secded72_gates.json")
    measured = {c["code"]: c for c in json.loads(GATES.read_text())["codes"]}
    build = re.search(r"\(git sha1 ([0-9a-f]{7})", json.loads(GATES.read_text())["yosys"])
    same("the Yosys build that measured the data file", f"git {build.group(1)}", "git 5628c13")
    for name, _, _, cols, r, _ in CODES:
        label = LABELS[name]
        st = STATED[name]
        m = measured[name]
        loads = row_loads(cols, r)
        require(f"{label}: the data file's rows match the matrix", m["row_weights"] == loads)
        require(f"{label}: the data file's check equations match", m["check_fan_in"] == fan_ins[name])
        trees = m["syndrome_trees"]
        require(f"{label}: each separate tree has W - 1 gates, ceil(log2 W) levels",
                all(t["gates"] == w - 1 and t["xor"] == w - 1 and t["depth"] == levels(w)
                    for t, w in zip(trees, loads)))
        same(f"{label}: syndrome XORs, separate trees", sum(t["gates"] for t in trees), st["sep"])
        same(f"{label}: syndrome depth, separate trees", max(t["depth"] for t in trees), st["levels"])
        enc = m["encoder_trees"]
        require(f"{label}: each separate encoder tree has F - 1 gates",
                all(t["gates"] == f - 1 and t["depth"] == levels(f) for t, f in zip(enc, fan_ins[name])))
        same(f"{label}: encoder XORs, separate trees", sum(t["gates"] for t in enc), st["enc"])
        sw = m["syndrome_whole"]
        same(f"{label}: syndrome two-input gates, whole block", sw["two_input"], st["syn_whole"])
        same(f"{label}: syndrome depth per bit, whole block", joined(sw["depth"]), st["syn_whole_depth"])
        ew = m["encoder_whole"]
        same(f"{label}: encoder two-input gates, whole block", ew["two_input"], st["enc_whole"])
        same(f"{label}: encoder depth, whole block", ew["depth"], st["enc_whole_depth"])
        for flag, key in ((0, "dec0"), (1, "dec1")):
            d = m[f"decoder_whole_flag{flag}"]
            quietly(f"{name}: flag {flag} two-input count", d["two_input"] == d["total"] - d["cells"]["not"])
            same(f"{label}: decoder flag {flag}: gates, data, corrected, flag depth",
                 joined([d["two_input"], d["data_depth"], d["corrected_depth"], d["uncorrectable_depth"]]),
                 st[key])
        c0 = m["decoder_whole_flag0"]["cells"]
        same(f"{label}: decoder flag 0 cells and, or, xor, not",
             joined([c0["and"], c0["or"], c0["xor"], c0["not"]]), st["cells0"])

    section("Hsiao against extended Hamming, computed from the data file")
    h, e = measured["hsiao"], measured["hamming"]

    def tree_sum(code, key):
        return sum(t["gates"] for t in code[key])

    hs, es = tree_sum(h, "syndrome_trees"), tree_sum(e, "syndrome_trees")
    he, ee = tree_sum(h, "encoder_trees"), tree_sum(e, "encoder_trees")
    hd = max(t["depth"] for t in h["syndrome_trees"])
    ed = max(t["depth"] for t in e["syndrome_trees"])
    same("syndrome XORs saved by Hsiao, separate trees, percent", pct(es - hs, es), "24.6")
    same("encoder XORs saved by Hsiao, separate trees, percent", pct(ee - he, ee), "13.8")
    same("syndrome levels saved by Hsiao, percent", pct(ed - hd, ed), "28.6")
    hw, ew = h["syndrome_whole"]["two_input"], e["syndrome_whole"]["two_input"]
    same("extra syndrome gates in Hsiao, whole block, percent", pct(hw - ew, ew), "23.4")
    hx, ex = h["encoder_whole"]["two_input"], e["encoder_whole"]["two_input"]
    same("encoder gates saved by Hsiao, whole block, percent", pct(ex - hx, ex), "9.1")
    for flag, stated in ((0, "15.1"), (1, "18.7")):
        a = h[f"decoder_whole_flag{flag}"]["two_input"]
        b = e[f"decoder_whole_flag{flag}"]["two_input"]
        same(f"extra decoder gates in Hsiao, flag {flag}, percent", pct(a - b, b), stated)
    a, b = h["decoder_whole_flag0"]["data_depth"], e["decoder_whole_flag0"]["data_depth"]
    same("data path levels saved by Hsiao, percent", pct(b - a, b), "18.2")
    for name, hand, stated in (("hsiao", 855, "35.7"), ("hamming", 923, "48.2")):
        got = measured[name]["decoder_whole_flag0"]["two_input"]
        same(f"{LABELS[name]}: share of the hand count shared away, percent", pct(hand - got, hand), stated)

    section("Worked examples and exercises")
    hc = ch20.hsiao_72_64()
    same("Hsiao: column at position 0", bits(hc[0], 8), "11100000")
    same("Hsiao: column at position 1", bits(hc[1], 8), "11010000")
    same("Hsiao: syndrome of errors at positions 0 and 1", bits(hc[0] ^ hc[1], 8), "00110000")
    p = hc.index(ch20.col((3, 5, 6)))
    same("Hsiao: position of the column 00010110", p, 49)
    s = hc[0] ^ hc[p] ^ hc[70]
    same("Hsiao: syndrome of errors at positions 0, 49 and 70", bits(s, 8), "11110100")
    require("Hsiao: that syndrome has odd weight and is no column", weight(s) % 2 == 1 and s not in hc)
    for name, stated in (("hsiao", "00000000"), ("hamming", "11111111")):
        code = model.build(name)
        word = code.encode((1 << 64) - 1)
        checks = [(word >> q) & 1 for q in code.check_positions]
        same(f"{LABELS[name]}: check bits of the all-ones data word", "".join(map(str, checks)), stated)
    index_rows = [sum(1 for q in range(72) if (q >> b) & 1) for b in range(7)]
    same("position output: match lines feeding each index bit", joined(index_rows), "36,36,36,32,32,32,8")
    same("position output: OR gates", sum(w - 1 for w in index_rows), 205)
    same("position output: levels after the match lines", levels(max(index_rows)), 6)

    print("\nAll chapter 30 numbers reproduced.")


if __name__ == "__main__":
    main()
