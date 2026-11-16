"""Run the SymbiYosys tasks for the (72,64) SEC-DED encoder and decoder and check every verdict.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
Needs SymbiYosys, Yosys with its yosys-abc, yosys-smtbmc and the z3 solver on the
PATH, and Python 3.10 or later with nothing beyond the standard library:

    python3 formal/secded72/run_formal72.py

It runs the nine tasks of secded72.sby beside this file and expects:

    prove_hsiao            PASS, every guarantee holds for all 2^64 data words
                           and all 2629 error patterns of weight <= 2
    prove_hsiao_parity     PASS, the same with the parity form of the flag
    prove_hamming          PASS, the extended Hamming decoder, no-match form
    prove_hamming_parity   PASS, the same with the parity form of the flag
    cover_hsiao            PASS, an input of each weight 0, 1 and 2 exists
    cover_hamming          PASS, the same for the extended Hamming code
    narrow_cover_hsiao     FAIL, with the weight-2 cover the one left unreached
    prove_dictionary       PASS, the Hsiao code in dictionary order
    prove_39_32            PASS, the Hsiao (39,32) code

It prints each task's verdict and time and exits with status 1 when any result
differs. When a proof fails, ABC's pdr engine has found a counterexample, and
SymbiYosys turns it into a trace with yosys-witness; a missing Python module
there makes the task end in ERROR, so a task whose log shows the engine
returning FAIL is reported as FAIL. SymbiYosys writes each task's work folder,
secded72_<task>/, beside this file, and git ignores them.
"""

import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
WRAPPER = LAB / "hdl" / "secded72" / "secded72_formal.sv"

DONE = re.compile(r"DONE \((PASS|FAIL|ERROR|UNKNOWN|TIMEOUT), rc=(\d+)\)")
ENGINE_FAIL = re.compile(r"engine_\d+ \([^)]*\) returned FAIL")

TASKS = (
    # (task, expected verdict, what it shows)
    ("prove_hsiao", "PASS", "Hsiao, no-match flag: every guarantee for all 2^64 words"),
    ("prove_hsiao_parity", "PASS", "Hsiao, parity flag: the same guarantees"),
    ("prove_hamming", "PASS", "extended Hamming, no-match flag"),
    ("prove_hamming_parity", "PASS", "extended Hamming, parity flag"),
    ("cover_hsiao", "PASS", "inputs of weight 0, 1 and 2 all reachable"),
    ("cover_hamming", "PASS", "inputs of weight 0, 1 and 2 all reachable"),
    ("narrow_cover_hsiao", "FAIL", "the narrowed assumption leaves weight 2 unreachable"),
    ("prove_dictionary", "PASS", "Hsiao in dictionary order, for an exercise"),
    ("prove_39_32", "PASS", "Hsiao (39,32), for an exercise"),
)


def cover_lines() -> dict:
    """The source line of each cover statement in the wrapper, by the weight it asks for."""
    lines = {}
    for number, text in enumerate(WRAPPER.read_text().splitlines(), start=1):
        found = re.search(r"cover\s*\(\s*wt\s*==\s*2'd(\d)\s*\)", text)
        if found:
            lines[int(found.group(1))] = number
    return lines


def run_task(task: str) -> tuple:
    """Run one task; return (status, seconds, its log text)."""
    shutil.rmtree(HERE / f"secded72_{task}", ignore_errors=True)
    start = time.perf_counter()
    done = subprocess.run(["sby", "-f", "secded72.sby", task], cwd=HERE, capture_output=True, text=True)
    seconds = time.perf_counter() - start
    log = done.stdout + done.stderr
    found = DONE.search(log)
    status = found.group(1) if found else "NO VERDICT"
    if status == "ERROR" and ENGINE_FAIL.search(log):
        status = "FAIL"
    return status, seconds, log


def unreached_covers(log: str) -> set:
    """Source lines of the cover statements the summary lists as unreached."""
    lines = set()
    in_list = False
    for line in log.splitlines():
        if "unreached cover statements:" in line:
            in_list = True
            continue
        if in_list:
            found = re.search(r"secded72_formal\.sv:(\d+)\.", line)
            if "summary:" in line and found and "reached cover statement" not in line:
                lines.add(int(found.group(1)))
            else:
                in_list = False
    return lines


def main() -> int:
    if shutil.which("sby") is None:
        print("sby is not on the PATH")
        return 2
    covers = cover_lines()
    trouble = []
    total = 0.0
    print("SymbiYosys tasks for the (72,64) decoders: proofs on abc pdr, covers on smtbmc z3")
    for task, want, note in TASKS:
        status, seconds, log = run_task(task)
        total += seconds
        verdict = "as expected" if status == want else f"EXPECTED {want}"
        print(f"  {task:<22} {status:<6} {seconds:>6.1f} s  {verdict}  {note}")
        if status != want:
            trouble.append(f"{task}: {status}, expected {want}")
        if task == "narrow_cover_hsiao":
            unreached = unreached_covers(log)
            print(f"      unreached cover lines {sorted(unreached)}, the weight-2 cover is line {covers[2]}")
            if unreached != {covers[2]}:
                trouble.append(f"{task}: unreached covers {sorted(unreached)}, expected only line {covers[2]}")
    print(f"  all tasks {total:.1f} s")
    if trouble:
        for line in trouble:
            print(f"FAIL {line}")
        return 1
    print("Every task gave the expected verdict.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
