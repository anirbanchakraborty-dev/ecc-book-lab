"""cocotb tests for the parity bench of "Lab: A First Hardware Checker" in Error Correction Coding.

The chapter is "Lab: A First Hardware Checker". The bench, hdl/parity/parity_bench.sv,
chains a parity generator, a channel that flips the bits an error mask marks, and a
parity checker. These tests drive it under Verilator and compare every output with
a software model written in a few lines of Python below. Run them through the
runner beside this file:

    python tb/parity/run_parity.py

Four tests:

    test_bsc_faults      random data and binary symmetric channel faults at eight
                         crossover probabilities, counting what the checker flags
                         and what slips past it
    test_exhaustive      every data word against every error mask
    test_clean_words     every data word with no fault, which must never be flagged
    test_waveform        the first twelve words of the p = 0.2 run, alone, so the
                         waveform it records is short enough to read

Every random number comes from random.Random(seed).random(), the one method whose
sequence for a given seed the Python documentation guarantees across versions, so
python/ch12_lab_first_hardware_checker.py can replay the same draws with no
simulator and reproduce every count these tests log.
"""

import random
from math import sqrt

import cocotb
from cocotb.triggers import Timer

SEED = 12                    # restarts the channel generator for every crossover probability
SOURCE_OFFSET = 10**9        # the data generator is seeded with SEED + SOURCE_OFFSET
STEP_NS = 10                 # each test vector is held for 10 ns of simulated time
Z95 = 1.96                   # two-sided 95% point of the standard normal
P_VALUES = (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
WORDS_PER_P = 100_000        # words sent through the channel at each crossover probability
WAVE_P = 0.2                 # crossover probability of the waveform window
WAVE_WORDS = 12              # words in the waveform window


# The software model ---------------------------------------------------------------

def weight(x: int) -> int:
    """The number of 1 bits in x."""
    return bin(x).count("1")


def model_codeword(data: int, k: int) -> int:
    """The k data bits with the even-parity bit placed on top, at bit k."""
    return data | ((weight(data) & 1) << k)


def model_error(word: int) -> int:
    """1 when the word has odd weight, so it cannot be a codeword."""
    return weight(word) & 1


# Statistics -----------------------------------------------------------------------

def wilson(successes: int, trials: int, z: float = Z95) -> tuple:
    """Wilson score interval for a proportion, as "Project: A Channel Simulator" builds it."""
    phat = successes / trials
    center = phat + z * z / (2 * trials)
    spread = z * sqrt(phat * (1 - phat) / trials + z * z / (4 * trials * trials))
    scale = 1 + z * z / trials
    return (center - spread) / scale, (center + spread) / scale


def coverage_formula(p: float, n: int) -> float:
    """Fraction of damaged words that single parity flags on a BSC(p)."""
    return (1 - (1 - 2 * p) ** n) / (2 * (1 - (1 - p) ** n))


# Driving the bench ----------------------------------------------------------------

async def apply(dut, data: int, mask: int, k: int) -> int:
    """Drive one vector, check all three outputs against the model, return the flag."""
    dut.data.value = data
    dut.error_mask.value = mask
    await Timer(STEP_NS, unit="ns")
    codeword = model_codeword(data, k)
    received = codeword ^ mask
    got_codeword = dut.codeword.value.to_unsigned()
    got_received = dut.received.value.to_unsigned()
    got_error = int(dut.error.value)
    assert got_codeword == codeword, f"codeword {got_codeword:#x} != {codeword:#x} for data {data:#x}"
    assert got_received == received, f"received {got_received:#x} != {received:#x}"
    assert got_error == model_error(received), f"error flag {got_error} for word {received:#x}"
    return got_error


def bsc_stream(p: float, k: int):
    """Yield (data, mask) pairs: the data from one generator, the faults from another.

    Both generators restart for every call, so every crossover probability sees the
    same data words and the same uniform draws, and the replay script can rebuild
    any stretch of any run.
    """
    source = random.Random(SEED + SOURCE_OFFSET)
    channel = random.Random(SEED)
    while True:
        data = sum(int(source.random() < 0.5) << i for i in range(k))
        mask = sum(int(channel.random() < p) << i for i in range(k + 1))
        yield data, mask


@cocotb.test()
async def test_bsc_faults(dut):
    """Random data through a binary symmetric channel, eight crossover probabilities."""
    k = len(dut.data)
    n = k + 1
    for p in P_VALUES:
        stream = bsc_stream(p, k)
        damaged = flagged = slipped = 0
        for _ in range(WORDS_PER_P):
            data, mask = next(stream)
            flag = await apply(dut, data, mask, k)
            flagged += flag
            if mask:
                damaged += 1
                slipped += 1 - flag
        caught = damaged - slipped
        low, high = wilson(caught, damaged)
        formula = coverage_formula(p, n)
        verdict = "inside" if low <= formula <= high else "OUTSIDE"
        dut._log.info(
            f"p={p:<5} words={WORDS_PER_P} damaged={damaged} flagged={flagged} "
            f"slipped={slipped} coverage={caught / damaged:.5f} "
            f"wilson=[{low:.5f}, {high:.5f}] formula={formula:.5f} {verdict}"
        )


@cocotb.test()
async def test_exhaustive(dut):
    """Every data word against every one of the 2^n error masks."""
    k = len(dut.data)
    n = k + 1
    flagged = slipped = 0
    for data in range(2**k):
        for mask in range(2**n):
            flag = await apply(dut, data, mask, k)
            assert flag == weight(mask) % 2, f"mask {mask:#x} of weight {weight(mask)}"
            flagged += flag
            if mask and not flag:
                slipped += 1
    dut._log.info(
        f"k={k} n={n} vectors={2**k * 2**n} flagged={flagged} "
        f"slipped={slipped} per word: flagged={flagged // 2**k} slipped={slipped // 2**k}"
    )


@cocotb.test()
async def test_clean_words(dut):
    """No fault: every codeword must arrive unflagged and with even weight."""
    k = len(dut.data)
    for data in range(2**k):
        assert await apply(dut, data, 0, k) == 0
        assert weight(dut.codeword.value.to_unsigned()) % 2 == 0
    dut._log.info(f"k={k} clean words={2**k} all unflagged")


@cocotb.test()
async def test_waveform(dut):
    """The first twelve words of the p = 0.2 run, one every 10 ns, for the waveform."""
    k = len(dut.data)
    stream = bsc_stream(WAVE_P, k)
    for t in range(WAVE_WORDS):
        data, mask = next(stream)
        flag = await apply(dut, data, mask, k)
        dut._log.info(
            f"t={t * STEP_NS:>3} ns data={data:0{k}b} mask={mask:0{k + 1}b} "
            f"received={model_codeword(data, k) ^ mask:0{k + 1}b} error={flag}"
        )
