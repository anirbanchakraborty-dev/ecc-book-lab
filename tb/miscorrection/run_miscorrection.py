"""Build the (72,64) SEC-DED bench with Verilator for three codes and run chapter 31's campaign.

The chapter is "Lab: Measuring Miscorrection" in Error Correction Coding. Needs
Verilator 5.036 or later on the PATH and cocotb 2.1 in the Python that runs
this:

    python tb/miscorrection/run_miscorrection.py                    seed 31, 100000 draws of weights 3 and 4
    python tb/miscorrection/run_miscorrection.py --seed 7
    python tb/miscorrection/run_miscorrection.py --weights 5 --draws 20000

For each of the balanced Hsiao matrix, the dictionary Hsiao matrix and the
extended Hamming matrix it builds hdl/secded72/secded72_bench.sv as it stands,
with the uncorrectable flag in its nomatch form, and runs the campaign of
test_miscorrection.py. The bench takes a data word and an error pattern and
returns the decoder's outputs, which is all a miscorrection campaign needs, so
no wrapper is added around it.

The runner then repeats the campaign in software with
python/ch31_lab_measuring_miscorrection.py, from the same seed, and prints for
each code and weight the counts of the four outcomes, the rate of the outcome
the chapter tracks with its 95% Wilson interval beside the exact share, the
difference between the two Hsiao matrices measured on the same patterns, the
outcomes of the adjacent windows, and the time each build and each campaign
took. It exits non-zero when any count differs from the software campaign's,
when any case differs from the model or from the prediction, or when a test
fails.

At the recorded settings, seed 31 and 100000 draws of weights 3 and 4, it also
writes the counts to data/ch31/miscorrection_counts.json, which the chapter's
script checks. Every other output lands under tb/miscorrection/sim_build/, which
git ignores, one folder per code.
"""

import argparse
import json
import platform
import re
import subprocess
import sys
import time
from math import comb, sqrt
from pathlib import Path

import cocotb
from cocotb_tools.runner import get_results, get_runner

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
HDL = LAB / "hdl" / "secded72"
SOURCES = [HDL / name for name in ("secded72_pkg.sv", "secded72_syndrome.sv", "secded72_enc.sv",
                                   "secded72_dec.sv", "secded72_bench.sv")]
DATA = LAB / "data" / "ch31" / "miscorrection_counts.json"

for folder in (HERE, LAB / "tb" / "secded72", LAB / "python"):
    sys.path.insert(0, str(folder))
import ch31_lab_measuring_miscorrection as lab     # noqa: E402
import secded72_model as model                      # noqa: E402

CODE_NUMBER = {model.HSIAO: 0, model.HAMMING: 1, model.HSIAO_DICT: 2}
RECORDED = (lab.SEED, lab.DRAWS, (3, 4))


def run(name: str, seed: int, draws: int, weights) -> tuple:
    """Build and run one code; return (tests run, tests failed, results, build seconds, campaign seconds)."""
    build_dir = HERE / "sim_build" / name
    results_json = build_dir / "campaign.json"
    runner = get_runner("verilator")
    start = time.perf_counter()
    runner.build(
        sources=SOURCES,
        hdl_toplevel="secded72_bench",
        parameters={"CODE": CODE_NUMBER[name], "FLAG": model.FLAGS.index("nomatch")},
        build_args=["-Wall"],
        build_dir=build_dir,
        timescale=("1ns", "1ps"),
        always=True,
    )
    built = time.perf_counter() - start
    results_json.unlink(missing_ok=True)
    start = time.perf_counter()
    xml = runner.test(
        test_module="test_miscorrection",
        hdl_toplevel="secded72_bench",
        testcase=["campaign"],
        build_dir=build_dir,
        test_dir=build_dir,
        extra_env={"MISCORRECTION_CODE": name, "MISCORRECTION_SEED": str(seed),
                   "MISCORRECTION_DRAWS": str(draws), "MISCORRECTION_WEIGHTS": ",".join(map(str, weights)),
                   "MISCORRECTION_RESULTS": str(results_json)},
    )
    ran, failed = get_results(xml)
    results = json.loads(results_json.read_text()) if results_json.exists() else None
    return ran, failed, results, built, time.perf_counter() - start


def paired(outcomes: dict) -> tuple:
    """Weight-3 cases miscorrected by the balanced matrix alone, and by the dictionary matrix alone."""
    h, d = outcomes[model.HSIAO], outcomes[model.HSIAO_DICT]
    only_h = sum(a == "M" and b != "M" for a, b in zip(h, d))
    only_d = sum(b == "M" and a != "M" for a, b in zip(h, d))
    return only_h, only_d


def tool_versions() -> dict:
    verilator = subprocess.run(["verilator", "--version"], capture_output=True, text=True).stdout.strip()
    return {"python": platform.python_version(), "cocotb": cocotb.__version__, "verilator": verilator}


def write_data(seed: int, draws: int, weights, rows: dict, only: tuple) -> None:
    record = {"chapter": 31, "seed": seed, "draws": draws, "weights": list(weights),
              "adjacent_seed": seed + 1, "checkpoints": list(lab.CHECKPOINTS), "tools": tool_versions(),
              "codes": {}, "paired": {"hsiao_only": only[0], "dictionary_only": only[1]}}
    for name in lab.CODES:
        r = rows[name]
        record["codes"][name] = {
            "bench_code": CODE_NUMBER[name], "columns": r["columns"], "cases": r["cases"],
            "campaign": r["campaign"], "running": r["running"], "adjacent": r["adjacent"],
            "model_mismatches": r["model_mismatches"], "prediction_mismatches": r["prediction_mismatches"]}
    text = json.dumps(record, indent=2)
    text = re.sub(r"\[\s+([-\d,\s]+?)\s+\]", lambda m: "[" + " ".join(m.group(1).split()) + "]", text)
    DATA.parent.mkdir(parents=True, exist_ok=True)
    DATA.write_text(text + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=lab.SEED, help=f"seed of the campaign (default {lab.SEED})")
    parser.add_argument("--draws", type=int, default=lab.DRAWS, help=f"patterns of each weight (default {lab.DRAWS})")
    parser.add_argument("--weights", default="3,4", help="weights to draw, comma separated (default 3,4)")
    args = parser.parse_args()
    weights = tuple(int(w) for w in args.weights.split(","))
    if any(w < 1 or w > lab.N for w in weights):
        parser.error("each weight must lie between 1 and 72")
    wall = time.perf_counter()

    rows, trouble, times = {}, [], {}
    for name in lab.CODES:
        ran, failed, results, built, swept = run(name, args.seed, args.draws, weights)
        times[name] = (built, swept)
        if failed:
            trouble.append(f"{name}: {failed} of {ran} tests failed")
        if results is None:
            trouble.append(f"{name}: the campaign wrote no counts")
            continue
        rows[name] = results
        if results["model_mismatches"] or results["prediction_mismatches"]:
            trouble.append(f"{name}: {results['model_mismatches']} cases differ from the model and"
                           f" {results['prediction_mismatches']} from the prediction; first {results['first']}")

    cases = lab.campaign_cases(args.seed, args.draws, weights)
    exact = {}
    for name in lab.CODES:
        code = lab.Code(name)
        counts = lab.syndrome_counts(code.cols, max(weights))
        exact[name] = {w: lab.exact_split(code, counts, w) for w in weights}
        if name not in rows:
            continue
        soft = lab.tally(cases, code)
        for w in weights:
            if rows[name]["campaign"][f"weight_{w}"] != soft["counts"][w]:
                trouble.append(f"{name}, weight {w}: hardware counts {rows[name]['campaign'][f'weight_{w}']}"
                               f" differ from the software campaign's {soft['counts'][w]}")
        for kind, got in rows[name]["adjacent"].items():
            if got != lab.adjacent_outcomes(code, lab.layout(kind, name)):
                trouble.append(f"{name}: adjacent windows in the {kind} layout differ from the software's")

    print(f"Campaign: seed {args.seed}, {args.draws} draws of each weight, on random data,"
          f" through hdl/secded72 with the nomatch flag")
    print(f"  {'code':<17} {'weight':>6} {'correct':>8} {'flagged':>8} {'miscorrected':>13} {'undetected':>11}"
          f"   {'rate and 95% Wilson interval, percent':<40} {'exact':>7}")
    for name in lab.CODES:
        if name not in rows:
            continue
        for w in weights:
            c = rows[name]["campaign"][f"weight_{w}"]
            k = c[lab.tracked(w)]
            lo, hi = lab.wilson(k, args.draws)
            share = exact[name][w][lab.tracked(w)] / comb(lab.N, w)
            verdict = "inside" if lo <= share <= hi else "outside"
            print(f"  {name:<17} {w:>6} {c['correct']:>8} {c['flagged']:>8} {c['miscorrected']:>13}"
                  f" {c['undetected']:>11}   {lab.tracked(w):<12} {lab.pct(k / args.draws):>7}"
                  f"  [{lab.pct(lo)}, {lab.pct(hi)}]  {lab.pct(share):>7}  {verdict}")

    only = None
    if 3 in weights and model.HSIAO in rows and model.HSIAO_DICT in rows:
        only = paired({n: rows[n]["outcomes"]["weight_3"] for n in (model.HSIAO, model.HSIAO_DICT)})
        soft_only = lab.paired_counts(cases)
        if only != soft_only:
            trouble.append(f"paired counts {only} differ from the software's {soft_only}")
        dh = (only[1] - only[0]) / args.draws
        half = lab.Z * sqrt(max((only[0] + only[1]) / args.draws - dh * dh, 0) / args.draws)
        print(f"\nDictionary minus balanced on the same weight-3 patterns: {only[0]} miscorrected by the"
              f" balanced matrix alone, {only[1]} by the dictionary alone,")
        print(f"  difference {lab.pct(dh)} points, 95% interval {lab.pct(dh - half)} to {lab.pct(dh + half)}")

    print("\nAdjacent windows, 70 in each layout: miscorrected and flagged")
    for name in lab.CODES:
        if name in rows:
            a = rows[name]["adjacent"]
            print(f"  {name:<17} " + "   ".join(f"{kind} {a[kind]['miscorrected']:>2} {a[kind]['flagged']:>2}"
                                                for kind in a))

    print("\nCases driven and time taken, seconds")
    for name in lab.CODES:
        if name in rows:
            built, swept = times[name]
            print(f"  {name:<17} {rows[name]['cases']:>7} cases   build {built:6.1f}   campaign {swept:6.1f}"
                  f"   inside the test {rows[name]['seconds']:6.1f}")
    print(f"  whole run {time.perf_counter() - wall:.1f}")

    if trouble:
        for line in trouble:
            print(f"FAIL {line}")
        return 1
    if (args.seed, args.draws, weights) == RECORDED:
        write_data(args.seed, args.draws, weights, rows, only)
        print(f"wrote {DATA.relative_to(LAB)}")
    print("Every case matches the model, the prediction and the software campaign.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
