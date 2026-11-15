"""Reproduce every number in chapter 16 of Error Correction Coding, "Matrices and Elimination over Finite Fields".

Needs Python 3.10 or later and nothing beyond the standard library:

    python3 python/ch16_matrices_and_elimination_over_finite_fields.py

Each check prints the value computed here beside the value the chapter states,
and the program stops with exit status 1 at the first pair that disagrees.
Every number in this chapter is a whole number, a matrix over a small field, a
list of row operations or a fraction printed to six places, so every comparison
is exact. The small cases are checked exhaustively: the elimination theorem,
the rank theorem and the bound on row additions on all 4,096 binary matrices of
three rows and four columns, associativity on every triple of binary 2-by-2
matrices, and the count of invertible matrices on every binary square matrix up
to 4 by 4. A run that reaches the final line has reproduced every number in the
chapter.

A matrix is a list of rows, and each row a list of residues modulo a prime p.
The module also holds the elimination, solving and counting functions that
ch16_charts.py draws, so the figures and this check share one source.
"""

import sys
from itertools import product


# Helpers ----------------------------------------------------------------------

def same(label: str, computed, stated) -> None:
    print(f"  {label:<62} {str(computed):<16} {stated}")
    if computed != stated:
        print(f"MISMATCH in '{label}': computed {computed!r}, the chapter states {stated!r}")
        sys.exit(1)


def require(label: str, condition: bool) -> None:
    print(f"  {label:<62} {'holds' if condition else 'FAILS'}")
    if not condition:
        print(f"MISMATCH: '{label}' does not hold")
        sys.exit(1)


def section(title: str) -> None:
    print(f"\n{title}")


# Matrices over GF(p) -------------------------------------------------------------

def parse(text: str) -> list:
    """Rows written as digit strings separated by spaces, such as "1101 1011 0111"."""
    return [[int(ch) for ch in row] for row in text.split()]


def show(M: list) -> str:
    return " ".join("".join(str(v) for v in row) for row in M)


def word(v: list) -> str:
    return "".join(str(b) for b in v)


def identity(n: int) -> list:
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def mat_vec(A: list, x: list, p: int = 2) -> list:
    return [sum(a * b for a, b in zip(row, x)) % p for row in A]


def mat_mul(A: list, B: list, p: int = 2) -> list:
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p for j in range(len(B[0]))]
            for i in range(len(A))]


def augment(A: list, B: list) -> list:
    return [ra + rb for ra, rb in zip(A, B)]


def vectors(n: int, p: int = 2) -> list:
    return [list(v) for v in product(range(p), repeat=n)]


def xor_rows(rows: list, n: int) -> list:
    out = [0] * n
    for row in rows:
        out = [(a + b) % 2 for a, b in zip(out, row)]
    return out


# Row operations and elimination ------------------------------------------------

def apply_step(M: list, step: tuple, p: int = 2) -> list:
    R = [row[:] for row in M]
    if step[0] == "swap":
        _, i, j = step
        R[i], R[j] = R[j], R[i]
    elif step[0] == "scale":
        _, i, c = step
        R[i] = [v * c % p for v in R[i]]
    else:
        _, i, f, j = step
        R[i] = [(a - f * b) % p for a, b in zip(R[i], R[j])]
    return R


def undo_step(step: tuple, p: int = 2) -> tuple:
    """The row operation of the same kind that reverses step."""
    if step[0] == "swap":
        return step
    if step[0] == "scale":
        return ("scale", step[1], pow(step[2], -1, p))
    _, i, f, j = step
    return ("add", i, (p - f) % p, j)


def step_text(step: tuple, p: int = 2) -> str:
    if step[0] == "swap":
        return f"R{step[1] + 1} <-> R{step[2] + 1}"
    if step[0] == "scale":
        return f"R{step[1] + 1} <- {step[2]} R{step[1] + 1}"
    _, i, f, j = step
    if p == 2:
        return f"R{i + 1} <- R{i + 1} + R{j + 1}"
    coef = "" if f == 1 else f"{f} "
    return f"R{i + 1} <- R{i + 1} - {coef}R{j + 1}"


def steps_text(steps: list, p: int = 2) -> str:
    return "; ".join(step_text(s, p) for s in steps)


def rref(A: list, p: int = 2):
    """Gauss-Jordan elimination over GF(p), column by column.

    Returns the reduced row echelon form, the pivot columns (counted from 0), every
    row operation in the order it was made, and one stage per pivot column: the
    column, the operations that cleared it and the matrix after them.
    """
    R = [row[:] for row in A]
    m, n = len(R), len(R[0])
    steps, pivots, stages = [], [], []
    r = 0
    for c in range(n):
        if r == m:
            break
        k = next((i for i in range(r, m) if R[i][c] % p), None)
        if k is None:
            continue
        done = []
        if k != r:
            done.append(("swap", r, k))
            R = apply_step(R, done[-1], p)
        if R[r][c] != 1:
            done.append(("scale", r, pow(R[r][c], -1, p)))
            R = apply_step(R, done[-1], p)
        for i in range(m):
            if i != r and R[i][c]:
                done.append(("add", i, R[i][c], r))
                R = apply_step(R, done[-1], p)
        pivots.append(c)
        steps += done
        stages.append((c, done, [row[:] for row in R]))
        r += 1
    return R, pivots, steps, stages


def is_rref(R: list) -> bool:
    lead_cols = []
    seen_zero = False
    for row in R:
        nz = [c for c, v in enumerate(row) if v]
        if not nz:
            seen_zero = True
            continue
        if seen_zero or row[nz[0]] != 1:
            return False
        lead_cols.append(nz[0])
    if lead_cols != sorted(set(lead_cols)):
        return False
    return all(sum(1 for row in R if row[c]) == 1 for c in lead_cols)


def rank(A: list, p: int = 2) -> int:
    return len(rref(A, p)[1])


def null_basis(A: list, p: int = 2) -> list:
    """One null-space vector per free column, by back substitution."""
    R, pivots, _, _ = rref(A, p)
    n = len(A[0])
    basis = []
    for f in (c for c in range(n) if c not in pivots):
        v = [0] * n
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = (-R[i][f]) % p
        basis.append(v)
    return basis


def null_space(A: list, p: int = 2) -> list:
    return [v for v in vectors(len(A[0]), p) if not any(mat_vec(A, v, p))]


def span(rows: list, p: int = 2) -> set:
    n = len(rows[0])
    out = set()
    for coeffs in product(range(p), repeat=len(rows)):
        out.add(tuple(sum(c * row[j] for c, row in zip(coeffs, rows)) % p for j in range(n)))
    return out


def solve(A: list, b: list, p: int = 2):
    """A solution with every free variable set to 0, or None when there is none."""
    n = len(A[0])
    R, pivots, _, _ = rref(augment(A, [[v] for v in b]), p)
    if n in pivots:
        return None
    x = [0] * n
    for i, c in enumerate(pivots):
        x[c] = R[i][n]
    return x


def solutions(A: list, b: list, p: int = 2) -> list:
    return sorted(word(x) for x in vectors(len(A[0]), p) if mat_vec(A, x, p) == b)


def inverse(A: list, p: int = 2):
    n = len(A)
    R = rref(augment(A, identity(n)), p)[0]
    if [row[:n] for row in R] != identity(n):
        return None
    return [row[n:] for row in R]


def systematic(G: list, p: int = 2):
    """Row-reduce G, then move the pivot columns to the front: the form [I | P]."""
    R, pivots, _, _ = rref(G, p)
    order = pivots + [c for c in range(len(G[0])) if c not in pivots]
    return [[R[i][c] for c in order] for i in range(len(pivots))], order


def coset_table(A: list) -> list:
    """For each right side b in counting order, the sorted solutions of A x = b."""
    return [(word(b), solutions(A, b)) for b in vectors(len(A))]


# Counting invertible matrices ----------------------------------------------------

def rank_of_bit_rows(rows: list) -> int:
    lead = {}
    r = 0
    for v in rows:
        while v:
            h = v.bit_length() - 1
            if h in lead:
                v ^= lead[h]
            else:
                lead[h] = v
                r += 1
                break
    return r


def rank_counts(n: int) -> list:
    """How many n-by-n binary matrices have each rank 0..n, by trying all of them."""
    counts = [0] * (n + 1)
    mask = (1 << n) - 1
    for code in range(1 << (n * n)):
        counts[rank_of_bit_rows([(code >> (n * i)) & mask for i in range(n)])] += 1
    return counts


def invertible_count(n: int, q: int = 2) -> int:
    c = 1
    for i in range(n):
        c *= q ** n - q ** i
    return c


def invertible_fraction(n: int) -> float:
    f = 1.0
    for j in range(1, n + 1):
        f *= 1 - 2.0 ** -j
    return f


def f6(x: float) -> str:
    return f"{x:.6f}"


# A network of two-input XOR gates ---------------------------------------------------

# Signals 0 to 3 are the inputs x1 to x4, and gate g writes signal 4 + g, so the
# first gate forms x3 XOR x4 once and the second and third outputs both reuse it.
SHARED = [(2, 3), (0, 1), (5, 3), (0, 4), (1, 4)]
SHARED_OUT = [6, 7, 8]


def run_network(gates: list, outs: list, x: list) -> list:
    s = list(x)
    for a, b in gates:
        s.append(s[a] ^ s[b])
    return [s[o] for o in outs]


# The chapter's matrices ------------------------------------------------------------

A = parse("1101 1011 0111")             # 3 by 4, the running example
X = [1, 0, 1, 1]
E2 = parse("11 01")
F2 = parse("10 11")
B = parse("011011 101100 110111")       # 3 by 6, rank 2
G = parse("0110111 1110001 1101010 1100101")   # 4 by 7, rank 4
M = parse("110 011 111")                # 3 by 3, invertible
U = parse("110 011 101")                # 3 by 3, singular
M5 = parse("213 131 412")               # over GF(5)
B5 = [3, 0, 2]

# Worked examples and exercises
C = parse("10110 11001 01100")
CB = [1, 0, 1]
Q = parse("1100 1110 0111 0011")
D = parse("110101 011110 101100")
T3 = parse("211 120 112")
TB = [0, 1, 2]
AE = parse("10110 01101")
XE = [1, 0, 1, 1, 0]
K = parse("11010 01101 10111")
J2 = parse("11 10")
M5E = parse("23 11")
B5E = [1, 2]


# The checks ----------------------------------------------------------------------

def check_matrix_vector() -> None:
    section("Matrix times vector as parity sums")
    same("A", show(A), "1101 1011 0111")
    same("A x for x = 1011, one parity sum per row", word(mat_vec(A, X)), "010")
    cols = [[row[j] for row in A] for j in range(4)]
    same("XOR of the columns of A that x selects (1, 3, 4)",
         word(xor_rows([cols[j] for j in range(4) if X[j]], 3)), "010")
    same("row weights of A", " ".join(str(sum(r)) for r in A), "3 3 3")
    same("two-input XOR gates, one chain per row", sum(sum(r) - 1 for r in A), 6)

    require("the 5-gate network with x3 XOR x4 shared gives A x on all 16 inputs",
            all(run_network(SHARED, SHARED_OUT, x) == mat_vec(A, x) for x in vectors(4)))
    same("gates in the shared network", len(SHARED), 5)


def check_products() -> None:
    section("Multiplying two matrices")
    same("E F", show(mat_mul(E2, F2)), "01 11")
    same("F E", show(mat_mul(F2, E2)), "11 10")
    require("E F differs from F E", mat_mul(E2, F2) != mat_mul(F2, E2))
    mats = [[list(v[:2]), list(v[2:])] for v in product(range(2), repeat=4)]
    require("(P Q) R = P (Q R) for all 4,096 triples of binary 2-by-2 matrices",
            all(mat_mul(mat_mul(P, Qm), R) == mat_mul(P, mat_mul(Qm, R))
                for P in mats for Qm in mats for R in mats))
    same("number of triples checked", len(mats) ** 3, 4096)
    require("I A = A and A I = A for the 3-by-4 A",
            mat_mul(identity(3), A) == A and mat_mul(A, identity(4)) == A)
    same("message 1011 times G, the XOR of rows 1, 3 and 4 of G",
         word(mat_mul([[1, 0, 1, 1]], G)[0]), "0111000")


def check_row_operations() -> None:
    section("Three row operations")
    ops2 = [("swap", i, j) for i in range(3) for j in range(3) if i != j]
    ops2 += [("add", i, 1, j) for i in range(3) for j in range(3) if i != j]
    require("every swap and row addition on A is undone by repeating it",
            all(apply_step(apply_step(A, s), s) == A for s in ops2))
    ops5 = [("swap", i, j) for i in range(3) for j in range(3) if i != j]
    ops5 += [("scale", i, c) for i in range(3) for c in range(1, 5)]
    ops5 += [("add", i, f, j) for i in range(3) for j in range(3) if i != j for f in range(1, 5)]
    require("every row operation on the GF(5) matrix is undone by its reverse",
            all(apply_step(apply_step(M5, s, 5), undo_step(s, 5), 5) == M5 for s in ops5))
    same("reverse of R1 <- 2 R1 over GF(5)", step_text(undo_step(("scale", 0, 2), 5), 5),
         "R1 <- 3 R1")
    require("each operation is left multiplication by the operated identity",
            all(mat_mul(apply_step(identity(3), s), A) == apply_step(A, s) for s in ops2)
            and all(mat_mul(apply_step(identity(3), s, 5), M5, 5) == apply_step(M5, s, 5)
                    for s in ops5))
    same("elementary matrix of R2 <- R2 + R1 on three rows",
         show(apply_step(identity(3), ("add", 1, 1, 0))), "100 110 001")
    for name, Mx in (("A", A), ("B", B), ("G", G)):
        R = rref(Mx)[0]
        require(f"elimination keeps the null space of {name}", null_space(Mx) == null_space(R))
        require(f"elimination keeps the row space of {name}", span(Mx) == span(R))


def check_elimination_of_a() -> None:
    section("Elimination of the 3-by-4 matrix A")
    R, pivots, steps, stages = rref(A)
    same("operations, in order", steps_text(steps),
         "R2 <- R2 + R1; R1 <- R1 + R2; R3 <- R3 + R2; R1 <- R1 + R3")
    same("after column 1", show(stages[0][2]), "1101 0110 0111")
    same("after column 2", show(stages[1][2]), "1011 0110 0001")
    same("after column 4, the RREF", show(R), "1010 0110 0001")
    same("pivot columns", " ".join(str(c + 1) for c in pivots), "1 2 4")
    same("free column", " ".join(str(c + 1) for c in range(4) if c not in pivots), "3")
    require("the result is in reduced row echelon form", is_rref(R))
    m34 = [[list(v[0:4]), list(v[4:8]), list(v[8:12])] for v in product(range(2), repeat=12)]
    require("elimination ends in RREF on all 4,096 binary 3-by-4 matrices",
            all(is_rref(rref(Mx)[0]) for Mx in m34))


def check_rank() -> None:
    section("Rank and the rank theorem")
    R, pivots, steps, stages = rref(B)
    same("B", show(B), "011011 101100 110111")
    same("operations on B", steps_text(steps), "R1 <-> R2; R3 <- R3 + R1; R3 <- R3 + R2")
    same("B after column 1", show(stages[0][2]), "101100 011011 011011")
    same("RREF of B", show(R), "101100 011011 000000")
    table = [("A", A, 4, 3, 1), ("B", B, 6, 2, 4), ("G", G, 7, 4, 3), ("M", M, 3, 3, 0),
             ("C", C, 5, 3, 2)]
    for name, Mx, n, r, nul in table:
        same(f"{name}: columns, rank, nullity",
             f"{len(Mx[0])} {rank(Mx)} {len(null_basis(Mx))}", f"{n} {r} {nul}")
    m34 = [[list(v[0:4]), list(v[4:8]), list(v[8:12])] for v in product(range(2), repeat=12)]
    require("rank + nullity = 4 and the null space has 2^nullity words, all 3-by-4",
            all(len(null_space(Mx)) == 2 ** (4 - rank(Mx)) for Mx in m34))
    require("the row space has 2^rank words, all 3-by-4",
            all(len(span(Mx)) == 2 ** rank(Mx) for Mx in m34))


def check_null_space() -> None:
    section("The null space by back substitution")
    basis = null_basis(B)
    same("null-space basis of B, one word per free column 3, 4, 5, 6",
         " ".join(word(v) for v in basis), "111000 100100 010010 010001")
    require("B v = 0 for each basis word", all(not any(mat_vec(B, v)) for v in basis))
    same("words in the null space of B, by trying all 64", len(null_space(B)), 16)
    same("their count from the rank theorem, 2^(6 - 2)", 2 ** (6 - rank(B)), 16)
    require("the span of the basis is the whole null space",
            span(basis) == {tuple(v) for v in null_space(B)})
    same("null-space basis of A", " ".join(word(v) for v in null_basis(A)), "1110")


def check_solving() -> None:
    section("Solving a linear system")
    R = rref(augment(B, [[1], [0], [0]]))[0]
    same("RREF of [B | 100]", show(R), "1011000 0110110 0000001")
    same("solution of B x = 100", str(solve(B, [1, 0, 0])), "None")
    same("solutions of B x = 100, by trying all 64", len(solutions(B, [1, 0, 0])), 0)
    R = rref(augment(B, [[1], [1], [0]]))[0]
    same("RREF of [B | 110]", show(R), "1011001 0110111 0000000")
    same("particular solution of B x = 110", word(solve(B, [1, 1, 0])), "110000")
    same("solutions of B x = 110, by trying all 64", len(solutions(B, [1, 1, 0])), 16)
    same("right sides b with B x = b solvable", sum(1 for b in vectors(3) if solutions(B, b)), 4)
    same("solution of M x = 110", word(solve(M, [1, 1, 0])), "101")
    same("solutions of M x = 110, by trying all 8", " ".join(solutions(M, [1, 1, 0])), "101")
    table = coset_table(A)
    same("A x = b: solutions for each b, b = 000 to 111",
         "; ".join(f"{b}: {' '.join(s)}" for b, s in table),
         "000: 0000 1110; 001: 0111 1001; 010: 0101 1011; 011: 0010 1100; "
         "100: 0011 1101; 101: 0100 1010; 110: 0110 1000; 111: 0001 1111")
    require("each set is a particular solution plus {0000, 1110}",
            all(len(s) == 2 and int(s[0], 2) ^ int(s[1], 2) == 0b1110 for _, s in table))


def check_inverse() -> None:
    section("Inverting a square matrix")
    R, pivots, steps, stages = rref(augment(M, identity(3)))
    same("operations on [M | I]", steps_text(steps),
         "R3 <- R3 + R1; R1 <- R1 + R2; R1 <- R1 + R3; R2 <- R2 + R3")
    same("[M | I] after column 1", show(stages[0][2]), "110100 011010 001101")
    same("[M | I] after column 2", show(stages[1][2]), "101110 011010 001101")
    same("[M | I] after column 3", show(stages[2][2]), "100011 010111 001101")
    Minv = inverse(M)
    same("inverse of M", show(Minv), "011 111 101")
    require("M times its inverse is I, both ways",
            mat_mul(M, Minv) == identity(3) and mat_mul(Minv, M) == identity(3))
    same("RREF of U", show(rref(U)[0]), "101 011 000")
    same("rank of U", rank(U), 2)
    same("inverse of U", str(inverse(U)), "None")
    same("solution of M x = 110 as the inverse times 110", word(mat_vec(Minv, [1, 1, 0])), "101")


def check_counting() -> None:
    section("How many square matrices are invertible")
    for n, stated in ((1, 1), (2, 6), (3, 168), (4, 20160)):
        same(f"invertible {n}-by-{n} binary matrices, by trying all", rank_counts(n)[n], stated)
    stated = [1, 6, 168, 20160, 9999360, 20158709760, 163849992929280]
    for n in range(1, 8):
        same(f"product formula, n = {n}", invertible_count(n), stated[n - 1])
    fracs = ["0.500000", "0.375000", "0.328125", "0.307617", "0.298004", "0.293348", "0.291056"]
    for n in range(1, 8):
        same(f"fraction invertible, n = {n}", f6(invertible_fraction(n)), fracs[n - 1])
        require(f"product of (1 - 2^-j) equals count / 2^(n^2), n = {n}",
                abs(invertible_fraction(n) - invertible_count(n) / 2 ** (n * n)) < 1e-15)
    same("fraction invertible, n = 10", f6(invertible_fraction(10)), "0.289070")
    same("the limit of the products, to six places", f6(invertible_fraction(200)), "0.288788")
    require("the fraction exceeds 1/4 for every n from 1 to 60",
            all(invertible_fraction(n) > 0.25 for n in range(1, 61)))
    same("3-by-3 binary matrices of rank 0, 1, 2, 3", " ".join(map(str, rank_counts(3))),
         "1 49 294 168")
    same("4-by-4 binary matrices of rank 0, 1, 2, 3, 4", " ".join(map(str, rank_counts(4))),
         "1 225 7350 37800 20160")
    same("total 4-by-4 binary matrices", sum(rank_counts(4)), 65536)


def check_systematic() -> None:
    section("Systematic form")
    R, pivots, steps, stages = rref(G)
    same("G", show(G), "0110111 1110001 1101010 1100101")
    same("operations on G", steps_text(steps),
         "R1 <-> R2; R3 <- R3 + R1; R4 <- R4 + R1; R1 <- R1 + R2; "
         "R2 <- R2 + R3; R4 <- R4 + R3; R2 <- R2 + R4; R3 <- R3 + R4")
    same("G after column 1", show(stages[0][2]), "1110001 0110111 0011011 0010100")
    same("G after column 2", show(stages[1][2]), "1000110 0110111 0011011 0010100")
    same("G after column 3", show(stages[2][2]), "1000110 0101100 0011011 0001111")
    same("G after column 4, the form [I | P]", show(R), "1000110 0100011 0010100 0001111")
    same("operations on G, of which row additions", sum(1 for s in steps if s[0] == "add"), 7)
    same("P", show([row[4:] for row in R]), "110 011 100 111")
    S, order = systematic(G)
    require("no column swap is needed for G", order == list(range(7)))
    require("the leading 4-by-4 block of G is invertible", inverse([row[:4] for row in G]) is not None)
    require("G and [I | P] have the same row space", span(G) == span(S))
    same("codewords in the row space", len(span(G)), 16)
    same("codeword of message 1011 from [I | P]", word(mat_mul([[1, 0, 1, 1]], S)[0]), "1011101")
    same("message 1011 times P", word(mat_mul([[1, 0, 1, 1]], [row[4:] for row in S])[0]), "101")
    S, order = systematic(A)
    same("column order that puts the pivots of A first", " ".join(str(c + 1) for c in order),
         "1 2 4 3")
    same("A after the column swap, the form [I | P]", show(S), "1001 0101 0010")
    require("the leading 3-by-3 block of A is singular", inverse([row[:3] for row in A]) is None)
    require("A and [I | P] have different row spaces", span(A) != span(S))
    require("swapping positions 3 and 4 of every word of span(A) gives span([I | P])",
            {(v[0], v[1], v[3], v[2]) for v in span(A)} == span(S))


def check_gf5() -> None:
    section("Elimination over GF(5)")
    same("inverses of 1, 2, 3, 4 modulo 5", " ".join(str(pow(a, -1, 5)) for a in range(1, 5)),
         "1 3 2 4")
    R, pivots, steps, stages = rref(augment(M5, [[v] for v in B5]), 5)
    same("operations on [M5 | b]", steps_text(steps, 5),
         "R1 <- 3 R1; R2 <- R2 - R1; R3 <- R3 - 4 R1; R2 <-> R3; R2 <- 4 R2; "
         "R1 <- R1 - 3 R2; R3 <- 3 R3; R1 <- R1 - 2 R3; R2 <- R2 - 4 R3")
    same("after column 1", show(stages[0][2]), "1344 0021 0411")
    same("after column 2", show(stages[1][2]), "1022 0144 0021")
    same("after column 3", show(R), "1001 0102 0013")
    x = solve(M5, B5, 5)
    same("solution", word(x), "123")
    same("M5 x", word(mat_vec(M5, x, 5)), "302")
    same("solutions over GF(5), by trying all 125", " ".join(solutions(M5, B5, 5)), "123")


def check_hardware() -> None:
    section("Elimination in hardware")
    same("R2 + R1 for A, bit by bit", word([(a + b) % 2 for a, b in zip(A[1], A[0])]), "0110")
    adds_a = sum(1 for s in rref(A)[2] if s[0] == "add")
    same("row additions in the elimination of A", adds_a, 4)
    same("bit XORs in the elimination of A", adds_a * 4, 16)
    adds_g = sum(1 for s in rref(G)[2] if s[0] == "add")
    same("bit XORs in the elimination of G", adds_g * 7, 49)
    same("bound min(m, n) (m - 1) n on bit XORs for 4-by-7", min(4, 7) * 3 * 7, 84)
    m34 = [[list(v[0:4]), list(v[4:8]), list(v[8:12])] for v in product(range(2), repeat=12)]
    counts = [sum(1 for s in rref(Mx)[2] if s[0] == "add") for Mx in m34]
    worst = max(counts)
    require("row additions never exceed min(m, n)(m - 1) = 6, all 3-by-4", worst <= 6)
    same("most row additions over all 3-by-4 binary matrices", worst, 6)
    same("3-by-4 binary matrices that use all 6 additions", counts.count(6), 15)


def check_examples() -> None:
    section("Worked examples")
    R, pivots, steps, _ = rref(augment(C, [[v] for v in CB]))
    same("C", show(C), "10110 11001 01100")
    same("operations on [C | 101]", steps_text(steps),
         "R2 <- R2 + R1; R3 <- R3 + R2; R1 <- R1 + R3; R2 <- R2 + R3")
    same("RREF of [C | 101]", show(R), "101011 011001 000110")
    same("particular solution", word(solve(C, CB)), "11000")
    same("null-space basis of C", " ".join(word(v) for v in null_basis(C)), "11100 10011")
    same("all solutions of C x = 101", " ".join(solutions(C, CB)), "00100 01011 10111 11000")
    R, pivots, steps, _ = rref(augment(Q, identity(4)))
    same("operations on [Q | I]", steps_text(steps),
         "R2 <- R2 + R1; R2 <-> R3; R1 <- R1 + R2; R1 <- R1 + R3; R2 <- R2 + R3; "
         "R4 <- R4 + R3; R1 <- R1 + R4; R2 <- R2 + R4")
    Qinv = inverse(Q)
    same("inverse of Q", show(Qinv), "1011 0011 1100 1101")
    require("Q times its inverse is I", mat_mul(Q, Qinv) == identity(4))
    R, pivots, steps, _ = rref(D)
    same("operations on D", steps_text(steps),
         "R3 <- R3 + R1; R1 <- R1 + R2; R3 <- R3 + R2; R2 <- R2 + R3")
    same("RREF of D", show(R), "101011 011001 000111")
    S, order = systematic(D)
    same("column order for D", " ".join(str(c + 1) for c in order), "1 2 4 3 5 6")
    same("D in the form [I | P] after the swap", show(S), "100111 010101 001011")
    R, pivots, steps, stages = rref(augment(T3, [[v] for v in TB]), 3)
    same("operations on [T | 012] over GF(3)", steps_text(steps, 3),
         "R1 <- 2 R1; R2 <- R2 - R1; R3 <- R3 - R1; R2 <-> R3; R2 <- 2 R2; "
         "R1 <- R1 - 2 R2; R1 <- R1 - 2 R3")
    same("[T | 012] after column 1", show(stages[0][2]), "1220 0011 0202")
    same("RREF of [T | 012]", show(R), "1002 0101 0011")
    same("solution over GF(3)", word(solve(T3, TB, 3)), "211")
    same("T x", word(mat_vec(T3, [2, 1, 1], 3)), "012")


def check_exercises() -> None:
    section("Exercises")
    same("parity sums of 10110 01101 on x = 10110", word(mat_vec(AE, XE)), "11")
    same("XOR gates for those parity sums", sum(sum(r) - 1 for r in AE), 4)
    same("K", show(K), "11010 01101 10111")
    same("rank and nullity of K", f"{rank(K)} {len(null_basis(K))}", "2 3")
    same("words in the null space of K", len(null_space(K)), 8)
    s = ("add", 1, 2, 0)    # R2 <- R2 - 2 R1, which is R2 <- R2 + R1 over GF(3)
    T = apply_step(apply_step(apply_step(T3, s, 3), s, 3), s, 3)
    require("over GF(3), adding R1 to R2 three times restores the matrix", T == T3)
    require("over GF(3), adding R1 to R2 twice does not", apply_step(apply_step(T3, s, 3), s, 3) != T3)
    same("inverse of 11 10", show(inverse(J2)), "01 11")
    inv3 = sum(1 for v in product(range(3), repeat=4) if rank([list(v[:2]), list(v[2:])], 3) == 2)
    same("invertible 2-by-2 matrices over GF(3), by trying all 81", inv3, 48)
    same("(3^2 - 1)(3^2 - 3)", invertible_count(2, 3), 48)
    same("solution of 2x + 3y = 1, x + y = 2 over GF(5)", word(solve(M5E, B5E, 5)), "02")
    S, order = systematic(augment(M, identity(3)))
    same("systematic form of [M | I]", show(S), "100011 010111 001101")


def main() -> None:
    check_matrix_vector()
    check_products()
    check_row_operations()
    check_elimination_of_a()
    check_rank()
    check_null_space()
    check_solving()
    check_inverse()
    check_counting()
    check_systematic()
    check_gf5()
    check_hardware()
    check_examples()
    check_exercises()
    print("\nAll chapter 16 numbers reproduced.")


if __name__ == "__main__":
    main()
