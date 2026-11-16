"""cocotb sweep of the project chapter's (72,64) SEC-DED encoder and decoder.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
run_secded72.py builds hdl/secded72/secded72_bench.sv under Verilator for one code
and one form of the uncorrectable flag and runs the single test below, which
drives every error pattern of weight 0, 1 and 2 (1 + 72 + 2556 = 2629 of them for
a (72,64) code) onto the codeword of every bench word: the all-zero word, the
all-ones word, the 64 walking-one words and 1,000 random words drawn with a fixed
seed. For each case it checks two things:

    the hardware's codeword and its three decoder outputs equal what
    secded72_model.py computes for the same data word and error pattern
    the decoder outputs equal the guarantee: weight 0 returns the data with no
    flag, weight 1 returns the data with corrected = 1, and weight 2 raises
    uncorrectable and returns the received data bits unchanged

The code and flag come from the environment, set by the runner, which also names
the file the counts are written to:

    SECDED72_CODE       hsiao, hamming, hsiao_dictionary or hsiao_39_32
    SECDED72_FLAG       nomatch or parity
    SECDED72_RESULTS    path of the JSON file the counts go to

Needs cocotb 2.1 and Verilator, through run_secded72.py.
"""

import json
import os
import time

import cocotb
from cocotb.triggers import Timer

import secded72_model as model

CODE = model.build(os.environ.get("SECDED72_CODE", model.HSIAO))
FLAG = os.environ.get("SECDED72_FLAG", "nomatch")
RESULTS = os.environ.get("SECDED72_RESULTS")


@cocotb.test()
async def sweep(dut):
    """Every pattern of weight 2 or less on every bench word, against the model and the guarantee."""
    start = time.perf_counter()
    words = model.bench_words(CODE.k)
    patterns = model.envelope_patterns(CODE.n)
    by_weight = {0: 0, 1: 0, 2: 0}
    cases = 0
    model_mismatches = 0
    spec_mismatches = 0
    first = None
    for d in words:
        codeword = CODE.encode(d)
        dut.data.value = d
        for e in patterns:
            dut.error.value = e
            await Timer(1, unit="ns")
            got_codeword = dut.codeword.value.to_unsigned()
            got = (dut.data_out.value.to_unsigned(), int(dut.corrected.value), int(dut.uncorrectable.value))
            want = CODE.decode(codeword ^ e, FLAG)
            cases += 1
            by_weight[e.bit_count()] += 1
            if got_codeword != codeword or got != want:
                model_mismatches += 1
                if first is None:
                    first = f"model: d={d:#x} e={e:#x} hardware {got} codeword {got_codeword:#x}, model {want}"
            if got != model.spec(CODE, d, e):
                spec_mismatches += 1
                if first is None:
                    first = f"guarantee: d={d:#x} e={e:#x} hardware {got}, guarantee {model.spec(CODE, d, e)}"
    seconds = time.perf_counter() - start
    if RESULTS:
        with open(RESULTS, "w") as f:
            json.dump({"code": CODE.name, "flag": FLAG, "words": len(words), "patterns": len(patterns),
                       "cases": cases, "by_weight": by_weight, "model_mismatches": model_mismatches,
                       "spec_mismatches": spec_mismatches, "first": first, "seconds": seconds}, f)
    assert model_mismatches == 0 and spec_mismatches == 0, first
