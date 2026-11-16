"""Write hdl/secded72/secded72_pkg.sv, the check matrices of the project chapter, from the model.

The chapter is "Project: A SEC-DED Encoder and Decoder" in Error Correction Coding.
The package holds every table the encoder and decoder read: for each code its
sizes, the column of the check matrix at every codeword position, the same matrix
by rows, the codeword position of every data bit and every check bit, and the data
bits each check bit is the XOR of. Every table is computed by secded72_model.py, so
the hardware and the model cannot be built from two different matrices.

    python tb/secded72/gen_secded72_pkg.py            write the package
    python tb/secded72/gen_secded72_pkg.py --check    compare it with the file on disk
    python tb/secded72/gen_secded72_pkg.py --print hsiao_39_32
                                                      print one code's tables

The --check form exits 1 when the file on disk differs from what the model gives,
which is what check_secded72.py runs. Needs Python 3.10 and nothing else.
"""

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent.parent
PACKAGE = LAB / "hdl" / "secded72" / "secded72_pkg.sv"

sys.path.insert(0, str(HERE))
import secded72_model as model                       # noqa: E402

MAX_N, MAX_K, MAX_R = 72, 64, 8

# (package constant prefix, model name, one line on where the columns come from)
CODES = (
    ("HSIAO", model.HSIAO,
     "balanced Hsiao (72,64), hsiao_72_64() of python/ch20_secded_and_hsiao_codes.py"),
    ("HAMMING", model.HAMMING,
     "extended Hamming (72,64), extended_hamming(72, 7) of the same script"),
    ("HSIAO_DICT", model.HSIAO_DICT,
     "Hsiao (72,64) in dictionary order, hsiao_72_64_dictionary() there"),
    ("HSIAO_39_32", model.HSIAO_39_32,
     "Hsiao (39,32), hsiao_39_32() there"),
)


def packed(values: list, width: int, per_line: int, fmt) -> list[str]:
    """A concatenation, highest index first, of the values in fields of the given width."""
    fields = [fmt(v) for v in reversed(values)]
    lines = []
    for start in range(0, len(fields), per_line):
        chunk = ", ".join(fields[start:start + per_line])
        top = len(values) - 1 - start
        bottom = max(top - per_line + 1, 0)
        comma = "," if start + per_line < len(fields) else ""
        lines.append(f"        {chunk}{comma}  // {top} .. {bottom}" if top != bottom
                     else f"        {chunk}{comma}  // {top}")
    return lines


def table(name: str, width: int, slots: int, values: list, per_line: int, fmt, note: str) -> list[str]:
    padded = values + [0] * (slots - len(values))
    out = [f"    // {note}"]
    out.append(f"    localparam logic [{slots}*{width}-1:0] {name} = {{")
    out += packed(padded, width, per_line, fmt)
    out.append("    };")
    return out


def code_section(prefix: str, name: str, source: str) -> list[str]:
    code = model.build(name)
    out = [f"    // ---- {prefix}: {source}", ""]
    out.append(f"    localparam int {prefix}_N = {code.n};")
    out.append(f"    localparam int {prefix}_K = {code.k};")
    out.append(f"    localparam int {prefix}_R = {code.r};")
    parity = code.parity_row if code.parity_row is not None else -1
    out.append(f"    localparam int {prefix}_PARITY_ROW = {parity};")
    out.append("")
    out += table(f"{prefix}_COLUMNS", MAX_R, MAX_N, code.columns, 8,
                 lambda v: f"8'b{v:08b}",
                 "column of the check matrix at each codeword position, row 0 in bit 0")
    out.append("")
    out += table(f"{prefix}_ROWS", MAX_N, MAX_R, code.row_masks, 1,
                 lambda v: f"72'h{v:018x}",
                 "the same matrix by rows: bit j of row i is set when position j is in row i")
    out.append("")
    out += table(f"{prefix}_DATA_POS", 8, MAX_K, code.data_positions, 16,
                 lambda v: f"8'd{v}", "codeword position of each data bit")
    out.append("")
    out += table(f"{prefix}_CHECK_POS", 8, MAX_R, code.check_positions, 8,
                 lambda v: f"8'd{v}", "codeword position of each check bit")
    out.append("")
    out += table(f"{prefix}_CHECK_MASKS", MAX_K, MAX_R, code.check_masks, 1,
                 lambda v: f"64'h{v:016x}", "data bits whose XOR is each check bit")
    return out


def selector(kind: str, field: str) -> list[str]:
    """A function that returns one table, or one size, of the code a CODE number names."""
    out = [f"    function automatic {kind} {field.lower()}_of(input int code);"]
    out.append("        case (code)")
    for prefix, _, _ in CODES:
        out.append(f"            {prefix}: {field.lower()}_of = {prefix}_{field};")
    out.append(f"            default: {field.lower()}_of = HSIAO_{field};")
    out.append("        endcase")
    out.append("    endfunction")
    out.append("")
    return out


def render() -> str:
    out = [
        "// secded72_pkg.sv",
        "//",
        "// Check matrices of the SEC-DED encoder and decoder in hdl/secded72/. Written",
        "// by tb/secded72/gen_secded72_pkg.py from the column lists in",
        "// tb/secded72/secded72_model.py; edit the model and run the generator, never",
        "// this file, and tb/secded72/check_secded72.py fails while the two disagree.",
        "//",
        "// A module picks its code with a CODE parameter set to one of the numbers",
        "// below. Every table is laid out for the widest code, 72 positions of 8 rows,",
        "// and a smaller code leaves the slots above its own size at 0. A column is",
        "// written row 0 in the low bit, so 8'b00000111 has ones in rows 0, 1 and 2.",
        "//",
        "// From \"Project: A SEC-DED Encoder and Decoder\" in Error Correction Coding.",
        "// Synthesizable: constants and constant functions only.",
        "",
        "package secded72_pkg;",
        "",
    ]
    for number, (prefix, _, _) in enumerate(CODES):
        out.append(f"    localparam int {prefix} = {number};")
    out.append("")
    out.append(f"    localparam int MAX_N = {MAX_N};")
    out.append(f"    localparam int MAX_K = {MAX_K};")
    out.append(f"    localparam int MAX_R = {MAX_R};")
    out.append("")
    for prefix, name, source in CODES:
        out += code_section(prefix, name, source)
    out.append("    // ---- one table, or one size, of the code a CODE number names")
    out.append("")
    for field in ("N", "K", "R", "PARITY_ROW"):
        out += selector("int", field)
    out += selector("logic [MAX_N*MAX_R-1:0]", "COLUMNS")
    out += selector("logic [MAX_R*MAX_N-1:0]", "ROWS")
    out += selector("logic [MAX_K*8-1:0]", "DATA_POS")
    out += selector("logic [MAX_R*8-1:0]", "CHECK_POS")
    out += selector("logic [MAX_R*MAX_K-1:0]", "CHECK_MASKS")
    out.append("endpackage")
    return "\n".join(out) + "\n"


def show(name: str) -> None:
    code = model.build(name)
    print(f"{name}: n = {code.n}, k = {code.k}, r = {code.r}")
    print("columns, row 0 first:")
    for j, c in enumerate(code.columns):
        print(f"  position {j:>2}  {''.join(str((c >> i) & 1) for i in range(code.r))}")
    print("row weights:", ",".join(str(m.bit_count()) for m in code.row_masks))
    print("check-bit fan-in:", ",".join(str(m.bit_count()) for m in code.check_masks))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="compare the package on disk with the model")
    parser.add_argument("--print", dest="show", choices=model.CODE_NAMES, help="print one code's tables")
    args = parser.parse_args()
    if args.show:
        show(args.show)
        return 0
    text = render()
    if args.check:
        on_disk = PACKAGE.read_text() if PACKAGE.exists() else ""
        if on_disk != text:
            print(f"{PACKAGE.relative_to(LAB)} differs from what secded72_model.py gives;"
                  " run tb/secded72/gen_secded72_pkg.py to rewrite it")
            return 1
        print(f"{PACKAGE.relative_to(LAB)} matches the model")
        return 0
    PACKAGE.parent.mkdir(parents=True, exist_ok=True)
    PACKAGE.write_text(text)
    print(f"wrote {PACKAGE.relative_to(LAB)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
