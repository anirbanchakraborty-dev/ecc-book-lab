"""Synthesize the chapter 12 parity generator and checker with Yosys and count gates.

Chapter 12 of Error Correction Coding is "Lab: A First Hardware Checker". It proves
that a tree of two-input XOR gates over W bits needs W - 1 gates and has a deepest
path of ceil(log2 W) gates. This script asks Yosys what it actually builds from
hdl/parity/ for codeword lengths n = 8, 16, 32 and 64, where the generator takes
K = n - 1 data bits and the checker takes all n bits of the received word.

Needs only Python 3.10 and a yosys binary on the PATH:

    python synth/synth_parity.py

For each block and width it runs

    read_verilog -sv; chparam -set K <K>; synth -flatten; delete t:$scopeinfo;
    stat; ltp -noff

then reads the gate count, the split between XOR and XNOR cells, and the length of
the longest path (the number of cells on it) from the report. It exits non-zero if
any count differs from the formula, or if the cells are anything other than
two-input XOR and XNOR gates.
"""

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

HDL = Path(__file__).resolve().parent.parent / "hdl" / "parity"
N_VALUES = (8, 16, 32, 64)

# block name, files it needs, the width its tree spans for a given n
BLOCKS = (
    ("parity_gen", ("xor_tree.sv", "parity_gen.sv"), lambda n: n - 1),
    ("parity_check", ("xor_tree.sv", "parity_check.sv"), lambda n: n),
)


def synthesize(top: str, files: tuple, k: int) -> dict:
    """Run Yosys on one block at one width and return its cell counts and depth."""
    script = (
        f"read_verilog -sv {' '.join(files)}; "
        f"chparam -set K {k} {top}; "
        f"synth -flatten -top {top}; "
        "delete t:$scopeinfo; "
        "tee -o /dev/stdout stat; "
        "tee -o /dev/stdout ltp -noff"
    )
    report = subprocess.run(
        ["yosys", "-q", "-p", script], cwd=HDL, capture_output=True, text=True, check=True
    ).stdout
    # Recent Yosys releases print each count before its name ("6 $_XNOR_", "6 cells"),
    # older ones after it ("$_XNOR_ 6", "Number of cells: 6"). Both layouts are read.
    pairs = re.findall(r"^\s+(\d+)\s+\$_(\w+)_\s*$", report, re.M)
    pairs += [(count, name) for name, count in re.findall(r"^\s+\$_(\w+)_\s+(\d+)\s*$", report, re.M)]
    cells = {name: int(count) for count, name in pairs}
    found = re.search(r"^\s+(\d+) cells\s*$|Number of cells:\s+(\d+)", report, re.M)
    total = int(found.group(1) or found.group(2))
    depth = int(re.search(r"length=(\d+)", report).group(1))
    return {"cells": cells, "total": total, "depth": depth}


def main() -> int:
    if shutil.which("yosys") is None:
        print("yosys is not on the PATH")
        return 2
    version = subprocess.run(["yosys", "-V"], capture_output=True, text=True).stdout.strip()
    print(version)
    print()
    print(f"  {'block':<14}{'n':>4}{'K':>4}{'W':>4}{'gates':>7}{'W-1':>5}{'XOR':>5}{'XNOR':>6}{'depth':>7}{'ceil(log2 W)':>14}")
    failures = 0
    for top, files, span in BLOCKS:
        for n in N_VALUES:
            k = n - 1
            w = span(n)
            result = synthesize(top, files, k)
            xor = result["cells"].get("XOR", 0)
            xnor = result["cells"].get("XNOR", 0)
            bound = math.ceil(math.log2(w))
            print(
                f"  {top:<14}{n:>4}{k:>4}{w:>4}{result['total']:>7}{w - 1:>5}"
                f"{xor:>5}{xnor:>6}{result['depth']:>7}{bound:>14}"
            )
            if set(result["cells"]) - {"XOR", "XNOR"}:
                print(f"    unexpected cell types {sorted(result['cells'])}")
                failures += 1
            if result["total"] != w - 1 or xor + xnor != w - 1:
                print(f"    gate count {result['total']} is not W - 1 = {w - 1}")
                failures += 1
            if result["depth"] != bound:
                print(f"    depth {result['depth']} is not ceil(log2 W) = {bound}")
                failures += 1
            if xnor % 2:
                print(f"    {xnor} XNOR gates is odd, so the output would be inverted")
                failures += 1
    print()
    if failures:
        print(f"{failures} check(s) failed.")
        return 1
    print("Every gate count is W - 1, every depth is ceil(log2 W), every XNOR count is even.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
