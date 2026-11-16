"""Build the (72,64) SEC-DED bench with Verilator, once per code, and run its cocotb sweep.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
Needs Verilator 5.036 or later on the PATH and cocotb 2.1 in the Python that runs
this:

    python tb/secded72/run_secded72.py                         Hsiao, then extended Hamming
    python tb/secded72/run_secded72.py --code hsiao_dictionary
    python tb/secded72/run_secded72.py --code hsiao --flag parity

It first runs check_secded72.py, which stops everything if a column list differs
from the one "SECDED and Hsiao Codes" builds or the package is stale. Then, for
each code, it builds hdl/secded72/secded72_bench.sv with that code and the chosen
form of the uncorrectable flag (nomatch, the default, or parity) and runs the
sweep of test_secded72.py: every error pattern of weight 2 or less on the
codeword of every bench word, 2629 x 1066 = 2,802,514 cases for a (72,64) code.
It prints the number of cases of each weight, the cases on which the hardware
differs from the model or from the guarantee, and the time each sweep took, and
exits non-zero when any case differs or a test fails.

Every output lands under tb/secded72/sim_build/, which git ignores, one folder per
build.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from cocotb_tools.runner import get_results, get_runner

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
HDL = LAB / "hdl" / "secded72"
SOURCES = [HDL / name for name in ("secded72_pkg.sv", "secded72_syndrome.sv", "secded72_enc.sv",
                                   "secded72_dec.sv", "secded72_bench.sv")]

sys.path.insert(0, str(HERE))
import secded72_model as model                       # noqa: E402

CODE_NUMBER = {model.HSIAO: 0, model.HAMMING: 1, model.HSIAO_DICT: 2, model.HSIAO_39_32: 3}


def run(name: str, flag: str) -> tuple:
    """Build and sweep one code; return (tests run, tests failed, counts, build seconds)."""
    build_dir = HERE / "sim_build" / f"{name}_{flag}"
    results_json = build_dir / "sweep.json"
    runner = get_runner("verilator")
    start = time.perf_counter()
    runner.build(
        sources=SOURCES,
        hdl_toplevel="secded72_bench",
        parameters={"CODE": CODE_NUMBER[name], "FLAG": model.FLAGS.index(flag)},
        build_args=["-Wall"],
        build_dir=build_dir,
        timescale=("1ns", "1ps"),
        always=True,
    )
    built = time.perf_counter() - start
    results_json.unlink(missing_ok=True)
    xml = runner.test(
        test_module="test_secded72",
        hdl_toplevel="secded72_bench",
        testcase=["sweep"],
        build_dir=build_dir,
        test_dir=build_dir,
        extra_env={"SECDED72_CODE": name, "SECDED72_FLAG": flag, "SECDED72_RESULTS": str(results_json)},
    )
    ran, failed = get_results(xml)
    counts = json.loads(results_json.read_text()) if results_json.exists() else None
    return ran, failed, counts, built


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--code", choices=model.CODE_NAMES, action="append",
                        help="code to sweep; give it more than once for several (default: hsiao and hamming)")
    parser.add_argument("--flag", choices=model.FLAGS, default="nomatch",
                        help="form of the uncorrectable flag (default: nomatch)")
    args = parser.parse_args()
    names = args.code or [model.HSIAO, model.HAMMING]

    check = subprocess.run([sys.executable, str(HERE / "check_secded72.py")], capture_output=True, text=True)
    print(check.stdout.strip().splitlines()[-1])
    if check.returncode != 0:
        print(check.stdout)
        print("FAIL check_secded72.py failed, so nothing was simulated")
        return 1

    rows = []
    trouble = []
    for name in names:
        ran, failed, counts, built = run(name, args.flag)
        rows.append((name, counts, built, ran, failed))
        if failed:
            trouble.append(f"{name}: {failed} of {ran} tests failed")
        if counts is None:
            trouble.append(f"{name}: the sweep wrote no counts")
        elif counts["model_mismatches"] or counts["spec_mismatches"]:
            trouble.append(f"{name}: {counts['model_mismatches']} cases differ from the model and"
                           f" {counts['spec_mismatches']} from the guarantee; first {counts['first']}")
    print()
    print(f"Sweep of every error pattern of weight <= 2, uncorrectable flag in its {args.flag} form")
    print(f"  {'code':<17} {'words':>6} {'patterns':>9} {'cases':>10} {'weight 0':>9} {'weight 1':>9}"
          f" {'weight 2':>10} {'vs model':>9} {'vs guarantee':>13} {'build s':>8} {'sweep s':>8}")
    for name, counts, built, ran, failed in rows:
        if counts is None:
            print(f"  {name:<17} no counts")
            continue
        w = counts["by_weight"]
        print(f"  {name:<17} {counts['words']:>6} {counts['patterns']:>9} {counts['cases']:>10,}"
              f" {w['0']:>9,} {w['1']:>9,} {w['2']:>10,} {counts['model_mismatches']:>9}"
              f" {counts['spec_mismatches']:>13} {built:>8.1f} {counts['seconds']:>8.1f}")
    if trouble:
        for line in trouble:
            print(f"FAIL {line}")
        return 1
    print("Every case matches the model and the guarantee.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
