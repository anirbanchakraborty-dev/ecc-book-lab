"""cocotb tests for the SEC-DED bench of "Verifying ECC Hardware" in Error Correction Coding.

The bench, hdl/secded/secded_bench.sv, chains the Hsiao (13,8) encoder, a channel
that flips the codeword bits an error pattern marks, and the decoder, built with
one value of its MUTANT parameter. These tests drive it under Verilator and judge
every output twice: against the specification inside the design envelope, and
against the software model in secded_model.py built with the same mutant. Run
them through the runner beside this file:

    python tb/secded/run_secded.py

Two tests:

    test_envelope       every data word against every error pattern of weight at
                        most 2, 2^8 * 92 = 23,552 cases, counting the cases the
                        decoder gets wrong and the wrong cases that are silent
    test_every_pattern  every data word against all 2^13 error patterns,
                        2,097,152 cases, run on the correct decoder only, checking
                        that the hardware matches the model on every input

The runner passes the mutant number in the environment variable MUTANT and a file
name in SECDED_RESULTS, where test_envelope writes its counts as JSON so the runner
can print them beside the model's.
"""

import json
import os
from pathlib import Path

import cocotb
from cocotb.triggers import Timer

from secded_model import K, N, check_case, count_wrong, decode, encode, envelope_patterns, spec

STEP_NS = 1                  # each case is held for 1 ns of simulated time
MUTANT = int(os.environ.get("MUTANT", "0"))


async def apply(dut, d: int, e: int) -> tuple:
    """Drive one case and return the decoder's (data, corrected, uncorrectable)."""
    dut.data.value = d
    dut.error.value = e
    await Timer(STEP_NS, unit="ns")
    return (dut.data_out.value.to_unsigned(), int(dut.corrected.value), int(dut.uncorrectable.value))


@cocotb.test()
async def test_envelope(dut):
    """All 23,552 cases inside the envelope: count wrong and silent, compare with the model."""
    wrong = silent = cases = 0
    first_wrong = None
    for d in range(2**K):
        codeword = encode(d)
        for e in envelope_patterns():
            got = await apply(dut, d, e)
            assert dut.codeword.value.to_unsigned() == codeword, f"codeword for data {d:#04x}"
            assert got == decode(codeword ^ e, MUTANT), (
                f"hardware {got} differs from the mutant-{MUTANT} model for d={d:#04x} e={e:#06x}")
            is_wrong = got != spec(d, e)
            wrong += is_wrong
            silent += is_wrong and got[0] != d and not got[2]
            if is_wrong and first_wrong is None:
                first_wrong = (d, e)
            cases += 1
    model_wrong, model_silent = count_wrong(MUTANT)
    dut._log.info(
        f"mutant={MUTANT} cases={cases} wrong={wrong} silent={silent} "
        f"model wrong={model_wrong} silent={model_silent} first wrong case={first_wrong}")
    results = os.environ.get("SECDED_RESULTS")
    if results:
        Path(results).write_text(json.dumps(
            {"mutant": MUTANT, "cases": cases, "wrong": wrong, "silent": silent}) + "\n")
    assert (wrong, silent) == (model_wrong, model_silent)
    if MUTANT == 0:
        assert wrong == 0, "the correct decoder broke a guarantee"
    else:
        assert wrong > 0, f"mutant {MUTANT} was never caught"


@cocotb.test(skip=MUTANT != 0)
async def test_every_pattern(dut):
    """The correct decoder against all 2^8 * 2^13 = 2,097,152 cases, hardware against model."""
    cases = beyond = 0
    for d in range(2**K):
        codeword = encode(d)
        for e in range(2**N):
            got = await apply(dut, d, e)
            assert got == decode(codeword ^ e), f"d={d:#04x} e={e:#06x}"
            if bin(e).count("1") <= 2:
                assert not check_case(d, e)[0]
            else:
                beyond += 1
            cases += 1
    dut._log.info(f"cases={cases} inside envelope={cases - beyond} beyond envelope={beyond} "
                  f"all match the model")
