"""Check that the project chapter's matrices are the ones "SECDED and Hsiao Codes" builds.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
For each code of secded72_model.py this script checks, and prints a line for:

    the column list equals, column for column and in the same order, the list
    python/ch20_secded_and_hsiao_codes.py builds for the same code
    the columns are distinct and nonzero and no three of them sum to zero, so
    the minimum distance is at least 4: the code corrects every single error
    and detects every double error
    every bench word's codeword has syndrome 0, so the check masks the
    encoder uses really solve the check equations

and then that hdl/secded72/secded72_pkg.sv is exactly what gen_secded72_pkg.py
writes from the model today. It exits 1 at the first failure, which is the
point of the exercise that swaps two columns. run_secded72.py runs it before
any simulation. Needs Python 3.10 and nothing else.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LAB / "python"))
import secded72_model as model                       # noqa: E402
import ch20_secded_and_hsiao_codes as ch20           # noqa: E402

# the list "SECDED and Hsiao Codes" builds for each code
CH20 = {
    model.HSIAO: ch20.hsiao_72_64,
    model.HAMMING: lambda: ch20.extended_hamming(72, 7),
    model.HSIAO_DICT: ch20.hsiao_72_64_dictionary,
    model.HSIAO_39_32: ch20.hsiao_39_32,
}


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"FAIL {label}")
        sys.exit(1)


def main() -> int:
    for name in model.CODE_NAMES:
        code = model.build(name)
        print(f"{name} ({code.n},{code.k})")
        theirs = CH20[name]()
        first = next((j for j, (a, b) in enumerate(zip(code.columns, theirs)) if a != b), None)
        label = "columns equal the chapter 20 list, column for column"
        if first is not None:
            label += f" (first difference at position {first})"
        require(label, code.columns == theirs)
        require("columns distinct and nonzero",
                len(set(code.columns)) == code.n and 0 not in code.columns)
        require("no three columns sum to zero, so minimum distance at least 4",
                ch20.distance_at_least_4(code.columns))
        require("every bench codeword has syndrome 0",
                all(code.syndrome(code.encode(d)) == 0 for d in model.bench_words(code.k)))
    print("package")
    done = subprocess.run([sys.executable, str(HERE / "gen_secded72_pkg.py"), "--check"],
                          capture_output=True, text=True)
    print(f"  {done.stdout.strip()}")
    require("secded72_pkg.sv matches the model", done.returncode == 0)
    print("Every check matrix is the chapter 20 matrix, and the package matches the model.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
