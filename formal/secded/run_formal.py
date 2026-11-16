"""Run the SymbiYosys tasks for the SEC-DED decoder and check every verdict.

The chapter is "Verifying ECC Hardware" in Error Correction Coding. Needs
SymbiYosys, Yosys, yosys-smtbmc and the z3 solver on the PATH, and Python 3.10 or
later with nothing beyond the standard library. From the repository root, with
the virtual environment of "Lab: A First Hardware Checker" active, so that
python names its interpreter:

    python formal/secded/run_formal.py
    python formal/secded/run_formal.py --any-counterexample

It runs the ten tasks of secded.sby beside this file and expects:

    prove_m0           PASS, the correct decoder meets every guarantee
    prove_m1 .. m4     FAIL, each with a counterexample trace; the (d, e) the
                       solver chose is read from the trace and handed to the
                       reference model in tb/secded/secded_model.py, which must
                       agree that the mutant gets that case wrong and must
                       predict the very outputs the trace records, and the pair
                       must also be the one the chapter quotes
    cover_m0           PASS, an input of each weight 0, 1 and 2 exists
    narrow_prove_m2    PASS, the off-by-one assumption wt(e) < 2 hides mutant 2
    narrow_cover_m2    FAIL, with the weight-2 cover the one left unreached
    empty_prove_m3     ERROR, the solver reports that no input satisfies the
                       contradictory assumptions
    empty_cover_m3     FAIL, with all three covers unreached

It exits with status 1 if any result differs. The pairs the chapter quotes are the
ones z3 returned on the recorded run, and another solver build may choose other
wrong cases. A run whose pairs differ fails and names both, and the option
--any-counterexample accepts any pair the model agrees is wrong. SymbiYosys writes
each task's work folder, secded_<task>/, beside this file, and git ignores them.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
WRAPPER = LAB / "hdl" / "secded" / "secded_formal.sv"

sys.path.insert(0, str(LAB / "tb" / "secded"))
import secded_model as model                         # noqa: E402

# The (d, e) pair each mutant's counterexample carries in the chapter, as z3 returned it.
CHAPTER_PAIRS = {1: (0b11111100, 0b0000000000011), 2: (0b01000000, 0b0001000000100),
                 3: (0b10111100, 0b0000000100000), 4: (0b10011001, 0b0000000000001)}

DONE = re.compile(r"DONE \((PASS|FAIL|ERROR|UNKNOWN|TIMEOUT), rc=(\d+)\)")
VAR = re.compile(r"\$var\s+\S+\s+(\d+)\s+(\S+)\s+(\S+)(?:\s+\[[^\]]*\])?\s+\$end")


def cover_lines() -> dict:
    """The source line of each cover statement in the wrapper, by the weight it asks for."""
    lines = {}
    for number, text in enumerate(WRAPPER.read_text().splitlines(), start=1):
        found = re.search(r"cover\s*\(\s*wt\s*==\s*4'd(\d)\s*\)", text)
        if found:
            lines[int(found.group(1))] = number
    return lines


def run_task(task: str) -> tuple:
    """Run one task; return (status, its log text)."""
    shutil.rmtree(HERE / f"secded_{task}", ignore_errors=True)
    done = subprocess.run(["sby", "-f", "secded.sby", task], cwd=HERE,
                          capture_output=True, text=True)
    log = done.stdout + done.stderr
    found = DONE.search(log)
    return (found.group(1) if found else "NO VERDICT"), log


def trace_values(vcd: Path) -> dict:
    """The top-level signals of a counterexample trace at time 0, as integers."""
    names = {}
    depth = 0
    values = {}
    at_zero = False
    for line in vcd.read_text().splitlines():
        line = line.strip()
        if line.startswith("$scope"):
            depth += 1
        elif line.startswith("$upscope"):
            depth -= 1
        elif line.startswith("$var") and depth == 1:
            found = VAR.match(line)
            if found:
                names[found.group(2)] = found.group(3)
        elif line.startswith("#"):
            if at_zero:
                break
            at_zero = line == "#0"
        elif at_zero and line:
            if line[0] in "bB":
                bits, ident = line[1:].split()
            else:
                bits, ident = line[0], line[1:]
            if ident in names:
                values[names[ident]] = int(bits, 2)
    return values


def main() -> int:
    any_pair = "--any-counterexample" in sys.argv[1:]
    covers = cover_lines()
    trouble = []

    def expect(task: str, want: str, note: str = "") -> str:
        status, log = run_task(task)
        verdict = "as expected" if status == want else f"EXPECTED {want}"
        print(f"  {task:<18} {status:<6} {verdict}{'  ' + note if note else ''}")
        if status != want:
            trouble.append(f"{task}: {status}, expected {want}")
        return log

    print("SymbiYosys tasks for the Hsiao (13,8) decoder, engine smtbmc z3")
    expect("prove_m0", "PASS", "every guarantee holds on every input the assumption allows")

    for mutant in (1, 2, 3, 4):
        task = f"prove_m{mutant}"
        expect(task, "FAIL")
        vcd = HERE / f"secded_{task}" / "engine_0" / "trace.vcd"
        if not vcd.exists():
            trouble.append(f"{task}: no counterexample trace at {vcd}")
            continue
        v = trace_values(vcd)
        d, e = v["d"], v["e"]
        got = (v["data_out"], v["corrected"], v["uncorrectable"])
        predicted = model.decode(model.encode(d) ^ e, mutant)
        wrong, silent = model.check_case(d, e, mutant)
        inside = model.weight(e) <= 2
        print(f"      counterexample d={d:08b} e={e:013b} (wt {model.weight(e)}) "
              f"out=({got[0]:08b},{got[1]},{got[2]}) spec=({model.spec(d, e)[0]:08b},"
              f"{model.spec(d, e)[1]},{model.spec(d, e)[2]}) model: wrong={int(wrong)} silent={int(silent)}")
        if not inside:
            trouble.append(f"{task}: counterexample lies outside the envelope")
        elif not wrong:
            trouble.append(f"{task}: the model says mutant {mutant} gets d={d:#04x} e={e:#06x} right")
        if got != predicted:
            trouble.append(f"{task}: trace outputs {got} differ from the model's {predicted}")
        quoted = CHAPTER_PAIRS[mutant]
        print(f"      the pair the chapter quotes: {'yes' if (d, e) == quoted else 'NO'}")
        if (d, e) != quoted and not any_pair:
            trouble.append(f"{task}: the solver chose d={d:08b} e={e:013b}, the chapter quotes "
                           f"d={quoted[0]:08b} e={quoted[1]:013b}, and --any-counterexample accepts any wrong case")

    expect("cover_m0", "PASS", "inputs of weight 0, 1 and 2 all reachable")

    print("Vacuity checks")
    expect("narrow_prove_m2", "PASS", "wt(e) < 2 hides mutant 2, the danger")
    log = expect("narrow_cover_m2", "FAIL")
    unreached = unreached_covers(log)
    print(f"      unreached cover lines {sorted(unreached)}, the weight-2 cover is line {covers[2]}")
    if unreached != {covers[2]}:
        trouble.append(f"narrow_cover_m2: unreached covers {sorted(unreached)}, expected only line {covers[2]}")
    log = expect("empty_prove_m3", "ERROR")
    unsat = "Assumptions are unsatisfiable" in log
    print(f"      solver reports the assumptions unsatisfiable: {'yes' if unsat else 'NO'}")
    if not unsat:
        trouble.append("empty_prove_m3: the log does not report unsatisfiable assumptions")
    log = expect("empty_cover_m3", "FAIL")
    unreached = unreached_covers(log)
    print(f"      unreached cover lines {sorted(unreached)}")
    if unreached != set(covers.values()):
        trouble.append(f"empty_cover_m3: unreached covers {sorted(unreached)}, expected all three")

    if trouble:
        for line in trouble:
            print(f"FAIL {line}")
        return 1
    print("Every task gave the expected verdict and every counterexample agrees with the model.")
    return 0


def unreached_covers(log: str) -> set:
    """Source lines of the cover statements the summary lists as unreached."""
    lines = set()
    in_list = False
    for line in log.splitlines():
        if "unreached cover statements:" in line:
            in_list = True
            continue
        if in_list:
            found = re.search(r"secded_formal\.sv:(\d+)\.", line)
            if "summary:" in line and found and "reached cover statement" not in line:
                lines.add(int(found.group(1)))
            else:
                in_list = False
    return lines


if __name__ == "__main__":
    sys.exit(main())
