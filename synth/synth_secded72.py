"""Synthesize the project chapter's SEC-DED encoder and decoder with Yosys and count gates.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
"SECDED and Hsiao Codes" predicts the gates of each code by hand: a syndrome bit
whose row holds W ones is a tree of W - 1 two-input XOR gates, ceil(log2 W) gates
deep, which makes 208 XOR gates for the Hsiao (72,64) syndrome and 276 for the
extended Hamming (72,64) one. This script asks Yosys what it builds from
hdl/secded72/ in two ways. Each time Yosys runs "abc -g AND,OR,XOR", which maps
onto two-input AND, OR and XOR gates plus the NOT gate Yosys always adds, with the
ABC script "strash; balance -x; map". In the words of ABC's own help, strash
transforms the logic into an AIG, a graph of two-input ANDs and inverters, balance
-x makes that graph well balanced, minimizing levels and balancing multi-input
XORs, and map performs standard-cell mapping, with area recovery on and area-only
mapping off. The library it maps onto is the stdcells.genlib Yosys writes, where
every cell, the inverter included, has a delay of 1, and the areas are 2 for NOT,
4 for AND and OR and 5 for XOR.

    separate trees   each syndrome bit, and each check bit of the encoder, is
                     synthesized on its own, with every other output removed, so
                     no gate can serve two trees. This is the hand count's model.
    whole block      the syndrome stage, the encoder and the decoder are each
                     synthesized whole, and Yosys and ABC may share any gate
                     between outputs, as they do in a real netlist.

Depth is counted in two-input gates on the longest path to each output, read from
the netlist Yosys writes. An inverter adds no level here and is left out of the
gate count, as in the hand count, although ABC's library charges it a unit of delay.

    python synth/synth_secded72.py                  every code, compared with the data file
    python synth/synth_secded72.py --code hsiao_dictionary
                                                    one code, printed and checked
    python synth/synth_secded72.py --write          rewrite data/ch30/secded72_gates.json

Every run checks the separate-tree counts against W - 1 and ceil(log2 W) for the
code's own rows and check bits, and exits 1 on any difference. The default run,
which covers every code, also compares each measurement with
data/ch30/secded72_gates.json and exits 1 when one differs, which a different
Yosys release can cause for the whole-block counts. Needs Python 3.10 and a yosys
binary on the PATH, and exits 2 without one.
"""

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
HDL = LAB / "hdl" / "secded72"
DATA = LAB / "data" / "ch30" / "secded72_gates.json"
FILES = ("secded72_pkg.sv", "secded72_syndrome.sv", "secded72_enc.sv", "secded72_dec.sv")

sys.path.insert(0, str(LAB / "tb" / "secded72"))
import secded72_model as model                       # noqa: E402

# the CODE parameter number of each code, as secded72_pkg.sv numbers them
CODE_NUMBER = {model.HSIAO: 0, model.HAMMING: 1, model.HSIAO_DICT: 2, model.HSIAO_39_32: 3}
# the cell types abc -g AND,OR,XOR may produce, by the name the data file uses
CELL_NAMES = {"$_AND_": "and", "$_OR_": "or", "$_XOR_": "xor", "$_NOT_": "not"}


def levels(w: int) -> int:
    """Depth of a balanced tree of two-input gates over w inputs."""
    return math.ceil(math.log2(w)) if w > 1 else 0


def yosys_netlist(top: str, params: dict, keep: str | None = None) -> dict:
    """Synthesize one module and return its netlist as Yosys writes it in JSON.

    keep names one output bit, such as "syndrome[3]"; every other output is then
    removed before mapping, so only that bit's cone of logic is built.
    """
    chparams = " ".join(f"chparam -set {name} {value} {top};" for name, value in params.items())
    isolate = ""
    if keep is not None:
        isolate = f"splitnets -ports; delete -output o:* o:{keep} %d; opt_clean; "
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "netlist.json"
        script = (
            f"read_verilog -sv {' '.join(FILES)}; {chparams} "
            f"synth -flatten -noabc -top {top}; {isolate}"
            "abc -g AND,OR,XOR -script +strash;balance,-x;map; opt_clean; delete t:$scopeinfo; "
            f"write_json {out}"
        )
        subprocess.run(["yosys", "-q", "-p", script], cwd=HDL, capture_output=True, text=True, check=True)
        return json.loads(out.read_text())["modules"][top]


def measure(netlist: dict) -> dict:
    """Cell counts by type, and the depth in two-input gates of every output bit."""
    driver = {}
    cells = {name: 0 for name in CELL_NAMES.values()}
    for cell in netlist["cells"].values():
        if cell["type"] not in CELL_NAMES:
            raise SystemExit(f"unexpected cell type {cell['type']} in the netlist")
        cells[CELL_NAMES[cell["type"]]] += 1
        for bit in cell["connections"]["Y"]:
            driver[bit] = cell
    memo = {}

    def depth(bit) -> int:
        if bit in memo:
            return memo[bit]
        cell = driver.get(bit)
        if cell is None:                             # an input bit or a constant
            memo[bit] = 0
            return 0
        inputs = [b for pin, bits in cell["connections"].items() if pin != "Y" for b in bits]
        d = max(depth(b) for b in inputs) + (0 if cell["type"] == "$_NOT_" else 1)
        memo[bit] = d
        return d

    outputs = {}
    for name, port in netlist["ports"].items():
        if port["direction"] == "output":
            outputs[name] = [depth(b) for b in port["bits"]]
    total = sum(cells.values())
    return {"cells": cells, "total": total, "two_input": total - cells["not"], "depth": outputs}


def separate_trees(top: str, code_number: int, port: str, positions: list[int]) -> list[dict]:
    """Synthesize the given output bits one at a time and measure each."""
    found = []
    for p in positions:
        m = measure(yosys_netlist(top, {"CODE": code_number}, keep=f"{port}[{p}]"))
        found.append({"bit": p, "gates": m["total"], "xor": m["cells"]["xor"],
                      "depth": max(m["depth"][f"{port}[{p}]"])})
    return found


def synthesize_code(name: str) -> dict:
    code = model.build(name)
    number = CODE_NUMBER[name]
    row_weights = [m.bit_count() for m in code.row_masks]
    fan_in = [m.bit_count() for m in code.check_masks]
    result = {"code": name, "n": code.n, "k": code.k, "r": code.r,
              "row_weights": row_weights, "check_fan_in": fan_in}
    result["syndrome_trees"] = separate_trees("secded72_syndrome", number, "syndrome", list(range(code.r)))
    result["encoder_trees"] = separate_trees("secded72_enc", number, "codeword", code.check_positions)
    syn = measure(yosys_netlist("secded72_syndrome", {"CODE": number}))
    result["syndrome_whole"] = {"cells": syn["cells"], "total": syn["total"], "two_input": syn["two_input"],
                                "depth": syn["depth"]["syndrome"]}
    enc = measure(yosys_netlist("secded72_enc", {"CODE": number}))
    result["encoder_whole"] = {"cells": enc["cells"], "total": enc["total"], "two_input": enc["two_input"],
                               "depth": max(enc["depth"]["codeword"])}
    for flag in (0, 1):
        dec = measure(yosys_netlist("secded72_dec", {"CODE": number, "FLAG": flag}))
        result[f"decoder_whole_flag{flag}"] = {
            "cells": dec["cells"], "total": dec["total"], "two_input": dec["two_input"],
            "data_depth": max(dec["depth"]["data"]),
            "corrected_depth": dec["depth"]["corrected"][0],
            "uncorrectable_depth": dec["depth"]["uncorrectable"][0],
        }
    return result


def check_trees(result: dict) -> list[str]:
    """Differences between the separate-tree counts and W - 1, ceil(log2 W)."""
    problems = []
    for label, trees, weights in (("syndrome", result["syndrome_trees"], result["row_weights"]),
                                  ("encoder", result["encoder_trees"], result["check_fan_in"])):
        for tree, w in zip(trees, weights):
            if tree["gates"] != w - 1 or tree["xor"] != w - 1 or tree["depth"] != levels(w):
                problems.append(f"{result['code']} {label} bit {tree['bit']}: {tree['gates']} gates,"
                                f" depth {tree['depth']}, expected {w - 1} and {levels(w)}")
    return problems


def report(result: dict) -> None:
    syn_trees = sum(t["gates"] for t in result["syndrome_trees"])
    enc_trees = sum(t["gates"] for t in result["encoder_trees"])
    print(f"{result['code']}: n = {result['n']}, k = {result['k']}, r = {result['r']}")
    print(f"  row weights {result['row_weights']}, check-bit fan-in {result['check_fan_in']}")
    print(f"  syndrome, separate trees  {syn_trees:>4} XOR    depth {max(t['depth'] for t in result['syndrome_trees'])}"
          f"   per row {[t['gates'] for t in result['syndrome_trees']]}")
    syn = result['syndrome_whole']
    print(f"  syndrome, whole block     {syn['two_input']:>4} gates  depth"
          f" {max(syn['depth'])}   per row {syn['depth']}   (+{syn['cells']['not']} NOT)   {syn['cells']}")
    print(f"  encoder, separate trees   {enc_trees:>4} XOR    depth {max(t['depth'] for t in result['encoder_trees'])}")
    enc = result['encoder_whole']
    print(f"  encoder, whole block      {enc['two_input']:>4} gates  depth"
          f" {enc['depth']}   (+{enc['cells']['not']} NOT)   {enc['cells']}")
    for flag in (0, 1):
        d = result[f"decoder_whole_flag{flag}"]
        print(f"  decoder FLAG={flag}, whole      {d['two_input']:>4} gates  data {d['data_depth']},"
              f" corrected {d['corrected_depth']}, uncorrectable {d['uncorrectable_depth']}"
              f"   (+{d['cells']['not']} NOT)   {d['cells']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--code", choices=model.CODE_NAMES, help="synthesize one code only")
    parser.add_argument("--write", action="store_true", help="rewrite the data file from this run")
    args = parser.parse_args()
    if shutil.which("yosys") is None:
        print("yosys is not on the PATH")
        return 2
    version = subprocess.run(["yosys", "-V"], capture_output=True, text=True).stdout.strip()
    print(version)
    print("Gates are two-input AND, OR and XOR cells. Inverters are counted apart, as NOT, and add no level.")
    print()
    names = [args.code] if args.code else list(model.CODE_NAMES)
    results = []
    problems = []
    for name in names:
        result = synthesize_code(name)
        report(result)
        print()
        problems += check_trees(result)
        results.append(result)
    if args.write:
        DATA.parent.mkdir(parents=True, exist_ok=True)
        DATA.write_text(json.dumps({"yosys": version, "codes": results}, indent=1) + "\n")
        print(f"wrote {DATA.relative_to(LAB)}")
    elif not args.code:
        stored = json.loads(DATA.read_text())["codes"] if DATA.exists() else []
        if stored != results:
            problems.append(f"this run differs from {DATA.relative_to(LAB)}; a different Yosys"
                            " release can map the whole blocks differently")
    if problems:
        for p in problems:
            print(f"FAIL {p}")
        return 1
    print("Every separate tree has W - 1 XOR gates and depth ceil(log2 W).")
    if not args.code and not args.write:
        print(f"Every measurement matches {DATA.relative_to(LAB)}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
