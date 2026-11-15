"""Reproduce every number in front-matter chapter f3 of Error Correction Coding, "The Octave Route".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/chf3_ecc_octave_route.py
    python3 python/chf3_ecc_octave_route.py --handoff

The chapter teaches GNU Octave to a reader who knows the book's Python, so
most of its numbers come in pairs: what an Octave expression gives and what
the nearest Python expression gives. Where the two languages agree, as for
bitxor and the ^ operator, this script computes the Python side. Where they
differ, as for rounding, integer types and remainders, it computes a small
Python model of Octave's rule, written from Octave's documentation, so the
chapter's Octave numbers are reproduced here too. The Octave twin,
octave/chf3_ecc_octave_route.m, computes the same lines with Octave itself,
and tools/compare-routes.sh checks that the two print the same text, which
is what shows each model to be faithful.

Every value computed here is printed beside the value the chapter states.
The script stops with exit status 1 at the first pair that disagrees by more
than half a unit in the last digit the chapter prints, and at the first
count, list or string that differs at all.

With --handoff the script also runs GNU Octave as a separate program, the
way a reader moving between the two routes would. It writes three data words
to a JSON file, has Octave compute their parity bits and read them back, and
then measures how many of 100,000 doubles survive each crossing exactly:
through JSON and through text written with 17 significant digits, in each
direction. It also compares 100,000 draws of Octave's rand after
rand ("state", 1) with Python's random.Random(1), and it has Octave's
idivide divide 3,146 pairs of integers with its default rounding, with
"fix" and with "floor", and fix(a ./ b) on doubles, which this script checks
against exact truncation and against Python's //, and it counts the elements
Octave's parser makes of [ 1 - 1 ] and [ 1 -1 ]. The counts the chapter
states were measured with Octave 11.3.0. Another Octave version may count
differently, so for any version but 11.3.0 a different count is reported
as a note rather than as a failure. The parity bits must agree everywhere.
That mode exits 1 when Octave cannot be found, which is the OCTAVE variable
if it is set and octave on the PATH otherwise.
"""

import argparse
import json
import math
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


def half_unit(stated: str) -> float:
    text = stated.replace(",", "").rstrip("%")
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.partition(".")[2])
    return 0.5 * 10.0 ** (int(exponent or 0) - places)


def check(label: str, computed, stated: str) -> None:
    value = float(stated.replace(",", "").rstrip("%"))
    if stated.endswith("%"):
        computed = 100 * float(computed)
    print(f"  {label:<64} {float(computed):<16.9g} {stated}")
    if abs(float(computed) - value) > half_unit(stated) * (1 + 1e-9):
        print(f"MISMATCH in '{label}': computed {float(computed)!r}, the chapter states {stated}")
        sys.exit(1)


def same(label: str, computed, stated) -> None:
    print(f"  {label:<64} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def spaced(values) -> str:
    """A row of numbers as Octave's disp shows it, one space apart."""
    return " ".join(str(v) for v in values)


# Python models of Octave's rules ----------------------------------------------


def octave_range(first: int, last: int, step: int = 1) -> list:
    """Octave's first:step:last, which includes last when the steps land on it."""
    return list(range(first, last + (1 if step > 0 else -1), step))


def round_half_away(x: float) -> int:
    """Octave's round: to the nearest integer, a tie away from zero."""
    a = abs(x)
    whole = math.floor(a)
    r = whole + (1 if a - whole >= 0.5 else 0)
    return r if x >= 0 else -r


def to_int(x: float, lo: int, hi: int) -> int:
    """Octave's conversion to an integer type: round half away, then saturate."""
    return min(max(round_half_away(x), lo), hi)


def saturate(n: int, lo: int, hi: int) -> int:
    return min(max(n, lo), hi)


UINT8 = (0, 255)
INT8 = (-128, 127)
INT32 = (-(2**31), 2**31 - 1)


# Indexing and ranges -----------------------------------------------------------


def bits_lsb_first(x: int, n: int) -> str:
    return spaced((x >> i) & 1 for i in range(n))


def check_indexing() -> None:
    print("Indexing and ranges")
    v = [5, 6, 7]
    same("first element of [5 6 7]: v[0] in Python, v(1) in Octave", v[0], 5)
    same("last element: v[-1] in Python, v(end) in Octave", v[-1], 7)
    same("elements in the Octave range 0:6", len(octave_range(0, 6)), 7)
    same("elements in Python's range(7)", len(range(7)), 7)
    same("elements in the Octave range 1:2:9", len(octave_range(1, 9, 2)), 5)
    same("elements in the Octave range 1:2:10", len(octave_range(1, 10, 2)), 5)
    same("bitget(100, 1:7), the bits of 100 from the lowest", bits_lsb_first(100, 7), "0 0 1 0 0 1 1")
    same("bitget(100, 3), Python (100 >> 2) & 1", (100 >> 2) & 1, 1)


# Doubles ------------------------------------------------------------------------


def check_doubles() -> None:
    print("Doubles")
    check("eps, the gap from 1 to the next double", math.ulp(1.0), "2.22e-16")
    same("eps equals 2^-52", math.ulp(1.0) == 2.0**-52, True)
    same("flintmax, the integer 2^53", int(2.0**53), 9007199254740992)
    same("2^53 + 1 stored as a double", int(float(2**53 + 1)), 9007199254740992)
    same("0.1 + 0.2 == 0.3", 0.1 + 0.2 == 0.3, False)
    same("0.1 + 0.2 written with 17 significant digits", f"{0.1 + 0.2:.17g}", "0.30000000000000004")
    same("0.1 + 0.2 written with 16 significant digits", f"{0.1 + 0.2:.16g}", "0.3")


# Integer types ------------------------------------------------------------------


def check_integer_types() -> None:
    print("Integer types saturate")
    same("uint8(250) + 10", saturate(250 + 10, *UINT8), 255)
    same("uint8(5) - 10", saturate(5 - 10, *UINT8), 0)
    same("int8(-100) - 100", saturate(-100 - 100, *INT8), -128)
    same("int8(127) + 1", saturate(127 + 1, *INT8), 127)
    same("250 + 10 wrapped modulo 256, as C's uint8_t would", (250 + 10) % 256, 4)
    same("uint8(300)", to_int(300, *UINT8), 255)
    same("uint8(2.5), a tie rounds away from zero", to_int(2.5, *UINT8), 3)
    same("int8(-2.5)", to_int(-2.5, *INT8), -3)
    same("uint8(1) + 1.5, whose result is still a uint8", to_int(1 + 1.5, *UINT8), 3)
    same("uint8(200) * 2", saturate(200 * 2, *UINT8), 255)
    same("(uint8(250) + 10) - 10, saturating at each step", saturate(saturate(250 + 10, *UINT8) - 10, *UINT8), 245)
    same("uint8(250) + (10 - 10)", saturate(250 + (10 - 10), *UINT8), 250)


# Integer division ---------------------------------------------------------------


def check_integer_division() -> None:
    print("Integer division")
    same("uint8(7) / 2, where 3.5 rounds to 4", to_int(7 / 2, *UINT8), 4)
    same("uint8(5) / 2", to_int(5 / 2, *UINT8), 3)
    same("int32(-7) / int32(2)", to_int(-7 / 2, *INT32), -4)
    same("Python 7 // 2, floor(7 / 2) in Octave", 7 // 2, 3)
    same("Python -7 // 2, floor(-7 / 2) in Octave", -7 // 2, -4)
    same("idivide(int32(-7), int32(2), 'floor')", saturate(-7 // 2, *INT32), -4)
    same("int(-7 / 2), fix(-7 / 2) in Octave", int(-7 / 2), -3)


# Rounding -----------------------------------------------------------------------


def check_rounding() -> None:
    print("Rounding")
    same("Octave round(2.5), a tie away from zero", round_half_away(2.5), 3)
    same("Octave round(-2.5)", round_half_away(-2.5), -3)
    same("Octave round(0.5)", round_half_away(0.5), 1)
    same("Python round(2.5), roundb(2.5) in Octave, a tie to even", round(2.5), 2)
    same("Python round(-2.5), roundb(-2.5)", round(-2.5), -2)
    same("Python round(0.5), roundb(0.5)", round(0.5), 0)
    same("Python round(-0.5), roundb(-0.5)", round(-0.5), 0)
    same("fix(-3.5) in Octave, math.trunc(-3.5) in Python", math.trunc(-3.5), -3)
    same("floor(-3.5) in Octave, math.floor(-3.5) in Python", math.floor(-3.5), -4)


# Remainders ---------------------------------------------------------------------


def check_remainders() -> None:
    print("Remainders")
    same("mod(-7, 3) in Octave, -7 % 3 in Python", -7 % 3, 2)
    same("rem(-7, 3) in Octave, math.fmod(-7, 3) in Python", int(math.fmod(-7, 3)), -1)
    same("mod(7, -3), which takes the divisor's sign", 7 % -3, -2)
    same("rem(7, -3), which takes the dividend's sign", int(math.fmod(7, -3)), 1)
    pairs = [(a, b) for a in range(-20, 21) for b in range(-7, 8) if b != 0]
    floor_ok = sum(a == b * (a // b) + a % b for a, b in pairs)
    trunc_ok = sum(a == b * math.trunc(a / b) + int(math.fmod(a, b)) for a, b in pairs)
    same("pairs (a, b) tried, a from -20 to 20, b from -7 to 7, b != 0", len(pairs), 574)
    same("pairs with a = b*floor(a/b) + mod(a, b)", floor_ok, 574)
    same("pairs with a = b*fix(a/b) + rem(a, b)", trunc_ok, 574)
    differ = sum(a % b != int(math.fmod(a, b)) for a, b in pairs)
    same("pairs where mod(a, b) and rem(a, b) differ", differ, 180)
    same("mod(-6:6, 3), Python [a % 3 for a in range(-6, 7)]", spaced(a % 3 for a in range(-6, 7)),
         "0 1 2 0 1 2 0 1 2 0 1 2 0")
    same("rem(-6:6, 3), math.fmod(a, 3) for the same a", spaced(int(math.fmod(a, 3)) for a in range(-6, 7)),
         "0 -2 -1 0 -2 -1 0 1 2 0 1 2 0")


# Bit operations -----------------------------------------------------------------


def check_bits() -> None:
    print("Bit operations")
    same("bitxor(100, 29), 100 ^ 29 in Python", 100 ^ 29, 121)
    same("bitand(100, 29), 100 & 29", 100 & 29, 4)
    same("bitor(100, 29), 100 | 29", 100 | 29, 125)
    same("bitshift(100, -2), 100 >> 2", 100 >> 2, 25)
    same("bitshift(100, 3), 100 << 3", 100 << 3, 800)
    same("dec2bin(100), format(100, 'b')", format(100, "b"), "1100100")
    same("dec2bin(100, 8), format(100, '08b')", format(100, "08b"), "01100100")
    same("bin2dec('1100100'), int('1100100', 2)", int("1100100", 2), 100)
    same("ones in 100, (100).bit_count()", (100).bit_count(), 3)


# Vectors and matrices -----------------------------------------------------------


WORDS = [100, 29, 121]


def mod2_product(row: list, matrix: list) -> list:
    return [sum(r * matrix[i][j] for i, r in enumerate(row)) % 2 for j in range(len(matrix[0]))]


def check_matrices() -> None:
    print("Vectors and matrices")
    a, b = [1, 2, 3], [4, 5, 6]
    same("[1 2 3] .* [4 5 6], element by element", spaced(x * y for x, y in zip(a, b)), "4 10 18")
    same("[1 2 3] * [4; 5; 6], a row times a column", sum(x * y for x, y in zip(a, b)), 32)
    same("mod([1 1] * [1 0 1; 0 1 1], 2), a codeword", spaced(mod2_product([1, 1], [[1, 0, 1], [0, 1, 1]])), "1 1 0")
    same("ones in each of the words [100 29 121]", spaced(w.bit_count() for w in WORDS), "3 4 5")
    same("even-parity bit of each word", spaced(w.bit_count() % 2 for w in WORDS), "1 0 1")


# Random numbers -----------------------------------------------------------------


def check_random() -> None:
    print("Random numbers")
    rng = random.Random(1)
    first = rng.random()
    for _ in range(99_998):
        rng.random()
    last = rng.random()
    check("first draw of random.Random(1), rand after rand('state', 1)", first, "0.13436424411240122")
    check("draw 100,000 of the same stream", last, "0.6355707617477253")
    words = [(2**33 + 1) >> (32 * i) & 0xFFFFFFFF for i in range(2)]
    same("32-bit words of the seed 2^33 + 1, lowest first", spaced(words), "1 2")
    check("first draw of random.Random(2**33 + 1), rand('state', [1 2])", random.Random(2**33 + 1).random(), "0.5845141492495678")


# JSON between the routes --------------------------------------------------------


def check_json() -> None:
    print("JSON between the routes")
    text = json.dumps(WORDS, separators=(",", ":"))
    same("JSON text of [100 29 121], written compactly", text, "[100,29,121]")
    same("sum of the decoded array", sum(json.loads(text)), 250)
    same("rows and columns of the array as Octave's jsondecode returns it", f"{len(json.loads(text))} 1", "3 1")
    same("2^53 + 1 read from JSON into a double", int(float(json.loads("9007199254740993"))), 9007199254740992)
    same("2^53 + 1 sent as a JSON string, read back as text", json.loads('"9007199254740993"'), "9007199254740993")
    same("JSON true, a logical 1 in Octave", int(json.loads("true")), 1)


# Worked examples and exercises --------------------------------------------------


def check_examples_and_exercises() -> None:
    print("Worked examples and exercises")
    print("Example: a byte that saturates")
    same("uint8(200) + uint8(100)", saturate(200 + 100, *UINT8), 255)
    same("the true sum", 200 + 100, 300)
    same("the same sum wrapped modulo 256", (200 + 100) % 256, 44)
    print("Example: porting a floor division")
    same("Python -17 // 5, floor(-17 / 5) in Octave", -17 // 5, -4)
    same("Python -17 % 5, mod(-17, 5) in Octave", -17 % 5, 3)
    same("rem(-17, 5), math.fmod(-17, 5)", int(math.fmod(-17, 5)), -2)
    same("int32(-17) / int32(5), where -3.4 rounds to -3", to_int(-17 / 5, *INT32), -3)
    print("Example: the parity of three bytes")
    same("dec2bin([100 29 121]) rows, ones in each", spaced(format(w, "07b").count("1") for w in WORDS), "3 4 5")
    same("bitxor of the three words", WORDS[0] ^ WORDS[1] ^ WORDS[2], 0)
    print("Exercises")
    same("int8(-128) - 1", saturate(-128 - 1, *INT8), -128)
    same("int8(-128) * -1, whose true value 128 does not fit", saturate(-128 * -1, *INT8), 127)
    same("uint8(9) / 4, where 2.25 rounds to 2", to_int(9 / 4, *UINT8), 2)
    same("uint8(10) / 4, where 2.5 rounds to 3", to_int(10 / 4, *UINT8), 3)
    same("mod(-1, 8), Python -1 % 8", -1 % 8, 7)
    same("rem(-1, 8)", int(math.fmod(-1, 8)), -1)
    same("numel(10:-3:0), the range 10, 7, 4, 1", len(octave_range(10, 0, -3)), 4)
    same("round(-0.5) in Octave", round_half_away(-0.5), -1)
    same("bitshift(1, 10), 1 << 10", 1 << 10, 1024)
    same("dec2bin(5, 4), format(5, '04b')", format(5, "04b"), "0101")


# --handoff: run Octave as a separate program ------------------------------------


OCTAVE_HANDOFF = r"""
w = jsondecode (fileread ('words.json'));
parity = mod (sum (dec2bin (w) == '1', 2), 2);
exact = hex2num (strsplit (strtrim (fileread ('doubles.hex')), "\n")(:), 'double');
from_json = jsondecode (fileread ('doubles.json'));
from_text = str2double (strsplit (strtrim (fileread ('doubles.txt')), "\n"))(:);
rand ('state', 1);
draws = rand (100000, 1);
[A, B] = ndgrid (-60:60, [-13:-1, 1:13]);
q_default = double (idivide (int32 (A), int32 (B)));
q_floor = double (idivide (int32 (A), int32 (B), 'floor'));
q_fix = double (idivide (int32 (A), int32 (B), 'fix'));
q_trunc = fix (A ./ B);
doc_example = double (idivide (int8 ([-3, 3]), int8 (4)));
doc_fix = double (idivide (int8 ([-3, 3]), int8 (4), 'fix'));
out = struct ('version', version (), 'parity', parity', ...
              'json_in', sum (from_json == exact), 'text_in', sum (from_text == exact), ...
              'big', jsonencode (1e19), 'json_ulps', max (abs (from_json - exact) ./ eps (exact)), 'q_default', q_default(:)', 'q_floor', q_floor(:)', ...
              'q_fix', q_fix(:)', 'q_trunc', q_trunc(:)', 'doc_fix', doc_fix, ...
              'doc_example', doc_example, 'spaced_minus', numel ([ 1 - 1 ]), 'unary_minus', numel ([ 1 -1 ]));
fid = fopen ('back.json', 'w');  fputs (fid, jsonencode (exact));  fclose (fid);
fid = fopen ('back.txt', 'w');  fprintf (fid, '%.17g\n', exact);  fclose (fid);
fid = fopen ('draws.txt', 'w');  fprintf (fid, '%.17g\n', draws);  fclose (fid);
fid = fopen ('result.json', 'w');  fputs (fid, jsonencode (out));  fclose (fid);
"""


def handoff() -> int:
    octave = shutil.which(os.environ.get("OCTAVE", "octave"))
    print("\nHanding work to GNU Octave")
    if octave is None:
        print("  Octave was not found. Install it as the chapter describes, or set OCTAVE to its path.")
        return 1
    rng = random.Random(7)
    doubles = [rng.random() for _ in range(50_000)] + [rng.uniform(-1e6, 1e6) for _ in range(50_000)]
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "words.json").write_text(json.dumps(WORDS, separators=(",", ":")))
        (d / "doubles.json").write_text(json.dumps(doubles))
        (d / "doubles.txt").write_text("".join(f"{x:.17g}\n" for x in doubles))
        (d / "doubles.hex").write_text("".join(struct.pack(">d", x).hex() + "\n" for x in doubles))
        (d / "handoff.m").write_text(OCTAVE_HANDOFF)
        run = subprocess.run([octave, "--no-window-system", "--norc", "--quiet", "handoff.m"],
                             cwd=d, capture_output=True, text=True)
        if run.returncode != 0:
            print(f"  Octave exited with status {run.returncode}:")
            print("  " + (run.stdout + run.stderr).strip().replace("\n", "\n  "))
            return 1
        result = json.loads((d / "result.json").read_text())
        back_json = json.loads((d / "back.json").read_text())
        back_text = [float(t) for t in (d / "back.txt").read_text().split()]
        draws = [float(t) for t in (d / "draws.txt").read_text().split()]
    version = result["version"]
    print(f"  Octave {version} at {octave}")
    same("parity bits Octave computed for [100 29 121]", spaced(result["parity"]), "1 0 1")
    ref = random.Random(1)
    # ndgrid runs a fastest, so the pairs come out column by column: every a for the first b, and so on.
    pairs = [(a, b) for b in [*range(-13, 0), *range(1, 14)] for a in range(-60, 61)]
    truncated = [(abs(a) // abs(b)) * (1 if (a < 0) == (b < 0) else -1) for a, b in pairs]
    wrong = [a for (a, b), q, t in zip(pairs, result["q_default"], truncated) if q != t]
    counts = [
        ("idivide default against truncation, wrong of 3146 pairs", len(wrong), 1184),
        ("of those wrong, how many have a nonnegative dividend", sum(a >= 0 for a in wrong), 0),
        ("idivide 'floor' against Python's //, wrong of 3146 pairs",
         sum(q != a // b for (a, b), q in zip(pairs, result["q_floor"])), 0),
        ("idivide 'fix' against truncation, wrong of 3146 pairs",
         sum(q != t for q, t in zip(result["q_fix"], truncated)), 1184),
        ("fix(a ./ b) on doubles against truncation, wrong of 3146 pairs",
         sum(q != t for q, t in zip(result["q_trunc"], truncated)), 0),
        ("idivide(int8([-3 3]), int8(4)), documented as 0 0", spaced(int(q) for q in result["doc_example"]), "-1 0"),
        ("idivide(int8([-3 3]), int8(4), 'fix'), documented as 0 0", spaced(int(q) for q in result["doc_fix"]), "-1 0"),
        ("elements in [ 1 - 1 ], where - is subtraction", result["spaced_minus"], 1),
        ("elements in [ 1 -1 ], where - is a sign", result["unary_minus"], 2),
        ("doubles Python json.dumps wrote as their repr text", sum(json.dumps(x) == repr(x) for x in doubles), 100000),
        ("doubles exact after Python json.dumps, Octave jsondecode", result["json_in"], 89493),
        ("largest of those misses, in units in the last place", result["json_ulps"], 1),
        ("doubles exact after Python %.17g text, Octave str2double", result["text_in"], 100000),
        ("doubles exact after Octave jsonencode, Python json.loads", sum(a == b for a, b in zip(back_json, doubles)), 100000),
        ("doubles exact after Octave %.17g text, Python float", sum(a == b for a, b in zip(back_text, doubles)), 100000),
        ("draws of rand after rand('state', 1) equal to random.Random(1)", sum(x == ref.random() for x in draws), 100000),
        ("JSON text Octave's jsonencode writes for 1e19", result["big"], "9223372036854775807"),
    ]
    for label, computed, stated in counts:
        if version == "11.3.0" or computed == stated:
            same(label, computed, stated)
        else:
            print(f"  {label:<64} {str(computed):<16} {stated}")
            print(f"  NOTE: the chapter's value was measured with Octave 11.3.0, and Octave {version} differs here.")
    print("Hand-off checks pass.")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--handoff", action="store_true", help="also run GNU Octave and measure what crosses between the routes")
    args = ap.parse_args()
    check_indexing()
    check_doubles()
    check_integer_types()
    check_integer_division()
    check_rounding()
    check_remainders()
    check_bits()
    check_matrices()
    check_random()
    check_json()
    check_examples_and_exercises()
    print("\nAll chapter f3 numbers reproduced.")
    if args.handoff:
        sys.exit(handoff())


if __name__ == "__main__":
    main()
