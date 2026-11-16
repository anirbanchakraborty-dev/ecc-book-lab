"""cocotb campaign of chapter 31, "Lab: Measuring Miscorrection", in Error Correction Coding.

run_miscorrection.py builds hdl/secded72/secded72_bench.sv under Verilator for one
of three (72,64) codes and runs the single test below. The test draws the error
patterns of the chapter's campaign with draw() from
python/ch31_lab_measuring_miscorrection.py, through random.Random(seed), so they
are the patterns the chapter's script classifies in software: a number of
patterns of each chosen weight, each written onto the codeword of its own random
data word. It then writes a fresh random data word under every window of three
adjacent cells in three layouts of the word, the natural order, stride 11 and
the searched order, 70 windows each. Extended Hamming has no searched order of
its own, since no layout takes it below 46 miscorrected windows, so its third
layout repeats stride 11 and it sees 140 distinct windows.

For each case it drives the data word and the error pattern into the bench,
reads the decoder's three outputs and checks two things:

    the outputs equal what secded72_model.py's decoder returns for the same
    codeword and error pattern
    the data bits the outputs get wrong, and the corrected and uncorrectable
    flags, equal what the error pattern alone predicts, whatever the data word

Each case is sorted into one of four outcomes from the outputs alone: flagged
when uncorrectable is 1, correct when data_out equals the data word,
miscorrected when corrected is 1 and the data are still wrong, and undetected
when the data are wrong and neither flag is set.

The runner sets five variables in the environment:

    MISCORRECTION_CODE      hsiao, hsiao_dictionary or hamming
    MISCORRECTION_SEED      the campaign's seed
    MISCORRECTION_DRAWS     patterns of each weight
    MISCORRECTION_WEIGHTS   the weights, comma separated, such as 3,4
    MISCORRECTION_RESULTS   path of the JSON file the counts go to

Needs cocotb 2.1 and Verilator, through run_miscorrection.py.
"""

import json
import os
import time

import cocotb
from cocotb.triggers import Timer

import ch31_lab_measuring_miscorrection as lab
import secded72_model as model

NAME = os.environ.get("MISCORRECTION_CODE", model.HSIAO)
SEED = int(os.environ.get("MISCORRECTION_SEED", lab.SEED))
DRAWS = int(os.environ.get("MISCORRECTION_DRAWS", lab.DRAWS))
WEIGHTS = tuple(int(w) for w in os.environ.get("MISCORRECTION_WEIGHTS", "3,4").split(","))
RESULTS = os.environ.get("MISCORRECTION_RESULTS")

CODE = model.build(NAME)
LAB_CODE = lab.Code(NAME)
LETTER = {lab.CORRECT: "C", lab.FLAGGED: "F", lab.MISCORRECTED: "M", lab.UNDETECTED: "U"}


def prediction(positions) -> tuple:
    """(wrong data bits, corrected, uncorrectable) from the error pattern alone.

    The syndrome of the received word is the syndrome of the error pattern,
    because the codeword's own syndrome is zero, so the decoder's choice and
    the data bits it gets wrong depend on the pattern and not on the data.
    """
    e = sum(1 << p for p in positions)
    s = CODE.syndrome(e)
    j = CODE.position_of.get(s)
    if j is not None:
        e ^= 1 << j
    return CODE.data_bits(e), int(j is not None), int(s != 0 and j is None)


def classify(d: int, data_out: int, corrected: int, uncorrectable: int) -> str:
    """The outcome the decoder's outputs show."""
    if uncorrectable:
        return lab.FLAGGED
    if data_out == d:
        return lab.CORRECT
    return lab.MISCORRECTED if corrected else lab.UNDETECTED


class Campaign:
    def __init__(self, dut):
        self.dut = dut
        self.cases = 0
        self.model_mismatches = 0
        self.prediction_mismatches = 0
        self.first = None

    async def case(self, d: int, positions) -> str:
        """Drive one data word and one error pattern, check the outputs, return the outcome."""
        e = sum(1 << p for p in positions)
        self.dut.data.value = d
        self.dut.error.value = e
        await Timer(1, unit="ns")
        got = (self.dut.data_out.value.to_unsigned(), int(self.dut.corrected.value),
               int(self.dut.uncorrectable.value))
        self.cases += 1
        want = CODE.decode(CODE.encode(d) ^ e)
        if got != want:
            self.model_mismatches += 1
            self.first = self.first or f"model: d={d:#x} positions={positions} hardware {got}, model {want}"
        found = classify(d, *got)
        wrong, corrected, uncorrectable = prediction(positions)
        if (got[0] ^ d, got[1], got[2]) != (wrong, corrected, uncorrectable) or \
                found != lab.outcome(LAB_CODE, positions):
            self.prediction_mismatches += 1
            self.first = self.first or f"prediction: d={d:#x} positions={positions} hardware {got} {found}"
        return found


@cocotb.test()
async def campaign(dut):
    """The seeded random campaign and the adjacent windows, against the model and the prediction."""
    start = time.perf_counter()
    run = Campaign(dut)
    counts = {w: dict.fromkeys(lab.OUTCOMES, 0) for w in WEIGHTS}
    running = {w: [] for w in WEIGHTS}
    letters = {w: [] for w in WEIGHTS}
    for w, d, positions in lab.campaign_cases(SEED, DRAWS, WEIGHTS):
        found = await run.case(d, positions)
        counts[w][found] += 1
        letters[w].append(LETTER[found])
        if len(letters[w]) in lab.CHECKPOINTS:
            running[w].append(counts[w][lab.tracked(w)])

    adjacent = {}
    data = iter(lab.adjacent_data(SEED, 70 * 3))
    for kind in ("natural", "stride11", "searched"):
        order = lab.layout(kind, NAME)
        adjacent[kind] = dict.fromkeys(lab.OUTCOMES, 0)
        for k in range(lab.N - 2):
            adjacent[kind][await run.case(next(data), order[k:k + 3])] += 1

    seconds = time.perf_counter() - start
    if RESULTS:
        with open(RESULTS, "w") as f:
            json.dump({"code": NAME, "seed": SEED, "draws": DRAWS, "weights": list(WEIGHTS),
                       "columns": CODE.columns, "cases": run.cases,
                       "campaign": {f"weight_{w}": counts[w] for w in WEIGHTS},
                       "running": {f"weight_{w}": running[w] for w in WEIGHTS},
                       "outcomes": {f"weight_{w}": "".join(letters[w]) for w in WEIGHTS},
                       "adjacent": adjacent, "model_mismatches": run.model_mismatches,
                       "prediction_mismatches": run.prediction_mismatches, "first": run.first,
                       "seconds": seconds}, f)
    assert run.model_mismatches == 0 and run.prediction_mismatches == 0, run.first
