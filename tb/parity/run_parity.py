"""Build the chapter 12 parity bench with Verilator and run its cocotb tests.

Chapter 12 of Error Correction Coding is "Lab: A First Hardware Checker". Needs
Verilator 5.036 or later on the PATH and cocotb 2.1 in the Python that runs this:

    python tb/parity/run_parity.py

It builds hdl/parity/parity_bench.sv three times:

    k7        K = 7 data bits, the (8, 7) code: the fault-injection, exhaustive
              and clean-word tests
    k7-wave   the same bench with tracing on, running only the twelve-word
              waveform test, so the waveform stays small enough to read
    k3        K = 3, the (4, 3) code: the exhaustive and clean-word tests, which
              show the same design working at a second width

Every output lands under tb/parity/sim_build/, which git ignores, one folder per
build. The waveform is a VCD file. Surfer opens it with the bench's five signals
already in view when given the command file beside this script:

    surfer tb/parity/sim_build/k7-wave/dump.vcd --command-file tb/parity/parity_bench.sucl
"""

from pathlib import Path

from cocotb_tools.runner import get_runner

HERE = Path(__file__).resolve().parent
HDL = HERE.parent.parent / "hdl" / "parity"
SOURCES = [HDL / name for name in ("xor_tree.sv", "parity_gen.sv", "parity_check.sv", "parity_bench.sv")]

BUILDS = (
    # (folder, K, tests to run, record a waveform)
    ("k7", 7, ["test_bsc_faults", "test_exhaustive", "test_clean_words"], False),
    ("k7-wave", 7, ["test_waveform"], True),
    ("k3", 3, ["test_exhaustive", "test_clean_words"], False),
)


def run(folder: str, k: int, tests: list, waves: bool) -> None:
    build_dir = HERE / "sim_build" / folder
    runner = get_runner("verilator")
    runner.build(
        sources=SOURCES,
        hdl_toplevel="parity_bench",
        parameters={"K": k},
        build_args=["-Wall"],
        build_dir=build_dir,
        waves=waves,
        timescale=("1ns", "1ps"),
        always=True,
    )
    runner.test(
        test_module="test_parity",
        hdl_toplevel="parity_bench",
        testcase=tests,
        build_dir=build_dir,
        test_dir=build_dir,
        waves=waves,
    )


if __name__ == "__main__":
    for folder, k, tests, waves in BUILDS:
        run(folder, k, tests, waves)
