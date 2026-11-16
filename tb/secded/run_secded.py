"""Build the SEC-DED bench with Verilator, once per mutant, and run its cocotb tests.

The chapter is "Verifying ECC Hardware" in Error Correction Coding. Needs Verilator
5.036 or later on the PATH and cocotb 2.1 in the Python that runs this. From the
repository root, with the virtual environment of "Lab: A First Hardware Checker"
active, so that python names its interpreter:

    python tb/secded/run_secded.py

It builds hdl/secded/secded_bench.sv five times, with MUTANT = 0 (the correct
decoder) and with each planted bug 1 to 4, and runs the exhaustive envelope sweep
of 23,552 cases on every build. The correct build also runs the sweep over all
2,097,152 cases. For each build it prints the number of wrong cases the hardware
showed beside the number two software models predict: tb/secded/secded_model.py
and python/ch28_verifying_ecc_hardware.py, written separately. It exits non-zero
when any of the three disagree, when the correct decoder fails a test, or when a
mutant goes uncaught.

Every output lands under tb/secded/sim_build/, which git ignores, one folder per
build.
"""

import json
import sys
from pathlib import Path

from cocotb_tools.runner import get_results, get_runner

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
HDL = LAB / "hdl" / "secded"
SOURCES = [HDL / name for name in ("secded_enc.sv", "secded_dec.sv", "secded_bench.sv")]

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LAB / "python"))
import secded_model                                  # noqa: E402
import ch28_verifying_ecc_hardware as chapter        # noqa: E402

BUILDS = (
    # (MUTANT, tests to run)
    (0, ["test_envelope", "test_every_pattern"]),
    (1, ["test_envelope"]),
    (2, ["test_envelope"]),
    (3, ["test_envelope"]),
    (4, ["test_envelope"]),
)


def run(mutant: int, tests: list) -> tuple:
    """Build and test one mutant; return (tests run, tests failed, measured counts)."""
    build_dir = HERE / "sim_build" / f"m{mutant}"
    results_json = build_dir / "envelope.json"
    runner = get_runner("verilator")
    runner.build(
        sources=SOURCES,
        hdl_toplevel="secded_bench",
        parameters={"MUTANT": mutant},
        build_args=["-Wall"],
        build_dir=build_dir,
        timescale=("1ns", "1ps"),
        always=True,
    )
    results_json.unlink(missing_ok=True)
    xml = runner.test(
        test_module="test_secded",
        hdl_toplevel="secded_bench",
        testcase=tests,
        build_dir=build_dir,
        test_dir=build_dir,
        extra_env={"MUTANT": str(mutant), "SECDED_RESULTS": str(results_json)},
    )
    ran, failed = get_results(xml)
    measured = json.loads(results_json.read_text()) if results_json.exists() else None
    return ran, failed, measured


def main() -> int:
    rows = []
    trouble = []
    for mutant, tests in BUILDS:
        ran, failed, measured = run(mutant, tests)
        model_wrong, model_silent = secded_model.count_wrong(mutant)
        chapter_wrong, chapter_silent = chapter.mutant_counts(mutant)
        hw_wrong = measured["wrong"] if measured else None
        hw_silent = measured["silent"] if measured else None
        rows.append((mutant, measured["cases"] if measured else None, model_wrong, chapter_wrong,
                     hw_wrong, model_silent, hw_silent, ran, failed))
        if failed:
            trouble.append(f"mutant {mutant}: {failed} of {ran} tests failed")
        if measured is None:
            trouble.append(f"mutant {mutant}: the envelope sweep wrote no counts")
        elif not (hw_wrong == model_wrong == chapter_wrong and hw_silent == model_silent == chapter_silent):
            trouble.append(f"mutant {mutant}: hardware and models disagree")
    print()
    print("Envelope sweep, 2^8 data words x 92 error patterns of weight <= 2")
    print(f"  {'MUTANT':<7} {'cases':>6} {'model m':>8} {'ch28 m':>7} {'hardware m':>11} "
          f"{'model silent':>13} {'hardware silent':>16} {'tests':>6}")
    for mutant, cases, mw, cw, hw, ms, hs, ran, failed in rows:
        print(f"  {mutant:<7} {cases!s:>6} {mw:>8} {cw:>7} {hw!s:>11} {ms:>13} {hs!s:>16} "
              f"{f'{ran - failed}/{ran}':>6}")
    if trouble:
        for line in trouble:
            print(f"FAIL {line}")
        return 1
    print("Every count agrees and the correct decoder passed every test.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
