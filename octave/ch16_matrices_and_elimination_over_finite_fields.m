% Reproduce every number in chapter 16 of Error Correction Coding,
% "Matrices and Elimination over Finite Fields", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch16_matrices_and_elimination_over_finite_fields.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch16_matrices_and_elimination_over_finite_fields.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in this chapter is a whole number, a matrix over a
% small field, a list of row operations or a fraction printed to six places,
% so every comparison is exact. The small cases are checked exhaustively: the
% elimination theorem, the rank theorem and the bound on row additions on all
% 4,096 binary matrices of three rows and four columns, associativity on every
% triple of binary 2-by-2 matrices, and the count of invertible matrices on
% every binary square matrix up to 4 by 4. A run that reaches the final line
% has reproduced every number in the chapter.
%
% A matrix is an ordinary array of residues modulo a prime p, and a row
% operation is a row of four numbers: [1 i j 0] swaps rows i and j, [2 i c 0]
% multiplies row i by c, and [3 i f j] subtracts f times row j from row i.
% Octave's own rref, rank and inverse work over the real numbers, so the
% functions here carry other names. The largest whole number printed is below
% 2^53, where a double still holds every integer exactly.

1;


% Helpers ----------------------------------------------------------------------

function s = as_text (x)
  % The text Python's str() prints for a str, a bool or an int.
  if (ischar (x))
    s = x;
  elseif (islogical (x))
    if (x)
      s = 'True';
    else
      s = 'False';
    end
  else
    s = sprintf ('%d', x);
  end
end

function s = as_repr (x)
  if (ischar (x))
    s = ['''' x ''''];
  else
    s = as_text (x);
  end
end

function same (label, computed, stated)
  c = as_text (computed);
  s = as_text (stated);
  printf ('  %-62s %-16s %s\n', label, c, s);
  if (~strcmp (c, s))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, as_repr (computed), as_repr (stated));
    exit (1);
  end
end

function require (label, condition)
  if (condition)
    verdict = 'holds';
  else
    verdict = 'FAILS';
  end
  printf ('  %-62s %s\n', label, verdict);
  if (~condition)
    printf ('MISMATCH: ''%s'' does not hold\n', label);
    exit (1);
  end
end

function section (title)
  printf ('\n%s\n', title);
end


% Matrices over GF(p) -------------------------------------------------------------

function M = parse_rows (text)
  % Rows written as digit strings separated by spaces, such as '1101 1011 0111'.
  parts = strsplit (text, ' ');
  M = zeros (numel (parts), numel (parts{1}));
  for i = 1:numel (parts)
    M(i, :) = parts{i} - '0';
  end
end

function s = show (M)
  parts = cell (1, size (M, 1));
  for i = 1:size (M, 1)
    parts{i} = sprintf ('%d', M(i, :));
  end
  s = strjoin (parts, ' ');
end

function s = word (v)
  s = sprintf ('%d', v);
end

function s = join_numbers (v)
  parts = arrayfun (@(t) sprintf ('%d', t), v, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function V = all_vectors (n, p)
  % Every vector of length n over GF(p), one per row, in counting order.
  N = p ^ n;
  V = zeros (N, n);
  codes = (0:N-1)';
  for j = n:-1:1
    V(:, j) = mod (codes, p);
    codes = floor (codes / p);
  end
end

function M = matrix_from_code (code, m, n)
  % The binary m-by-n matrix whose m*n bits, row after row, spell code.
  bits = mod (floor (code ./ 2 .^ (m*n-1:-1:0)), 2);
  M = reshape (bits, n, m)';
end

function y = mat_vec (A, x, p)
  y = mod (A * x(:), p)';
end

function y = inv_mod (a, p)
  y = find (mod (a * (1:p-1), p) == 1, 1);
end


% Row operations and elimination ------------------------------------------------

function R = apply_step (M, st, p)
  R = M;
  switch (st(1))
    case 1
      R([st(2) st(3)], :) = M([st(3) st(2)], :);
    case 2
      R(st(2), :) = mod (M(st(2), :) * st(3), p);
    otherwise
      R(st(2), :) = mod (M(st(2), :) - st(3) * M(st(4), :), p);
  end
end

function u = undo_step (st, p)
  % The row operation of the same kind that reverses st.
  u = st;
  if (st(1) == 2)
    u(3) = inv_mod (st(3), p);
  elseif (st(1) == 3)
    u(3) = mod (p - st(3), p);
  end
end

function s = step_text (st, p)
  switch (st(1))
    case 1
      s = sprintf ('R%d <-> R%d', st(2), st(3));
    case 2
      s = sprintf ('R%d <- %d R%d', st(2), st(3), st(2));
    otherwise
      if (p == 2)
        s = sprintf ('R%d <- R%d + R%d', st(2), st(2), st(4));
      elseif (st(3) == 1)
        s = sprintf ('R%d <- R%d - R%d', st(2), st(2), st(4));
      else
        s = sprintf ('R%d <- R%d - %d R%d', st(2), st(2), st(3), st(4));
      end
  end
end

function s = steps_text (steps, p)
  parts = cell (1, size (steps, 1));
  for t = 1:size (steps, 1)
    parts{t} = step_text (steps(t, :), p);
  end
  s = strjoin (parts, '; ');
end

function [R, pivots, steps, stages] = rref_p (A, p)
  % Gauss-Jordan elimination over GF(p), column by column. It returns the
  % reduced row echelon form, the pivot columns, every row operation in the
  % order it was made, and one stage per pivot column: the column, the
  % operations that cleared it and the matrix after them.
  R = A;
  [m, n] = size (R);
  steps = zeros (0, 4);
  pivots = zeros (1, 0);
  stages = {};
  r = 1;
  for c = 1:n
    if (r > m)
      break;
    end
    k = find (mod (R(r:m, c), p) ~= 0, 1);
    if (isempty (k))
      continue;
    end
    k = k + r - 1;
    done = zeros (0, 4);
    if (k ~= r)
      done(end+1, :) = [1 r k 0];
      R = apply_step (R, done(end, :), p);
    end
    if (R(r, c) ~= 1)
      done(end+1, :) = [2 r inv_mod(R(r, c), p) 0];
      R = apply_step (R, done(end, :), p);
    end
    for i = 1:m
      if (i ~= r && R(i, c) ~= 0)
        done(end+1, :) = [3 i R(i, c) r];
        R = apply_step (R, done(end, :), p);
      end
    end
    pivots(end+1) = c;
    steps = [steps; done];
    if (nargout > 3)
      stages{end+1} = {c, done, R};
    end
    r = r + 1;
  end
end

function ok = is_rref (R)
  ok = false;
  lead = zeros (1, 0);
  seen_zero = false;
  for i = 1:size (R, 1)
    nz = find (R(i, :) ~= 0);
    if (isempty (nz))
      seen_zero = true;
      continue;
    end
    if (seen_zero || R(i, nz(1)) ~= 1)
      return;
    end
    lead(end+1) = nz(1);
  end
  if (any (diff (lead) <= 0))
    return;
  end
  ok = all (sum (R(:, lead) ~= 0, 1) == 1);
end

function r = gf_rank (A, p)
  [~, pivots] = rref_p (A, p);
  r = numel (pivots);
end

function adds = row_additions (A, p)
  [~, ~, steps] = rref_p (A, p);
  adds = sum (steps(:, 1) == 3);
end

function Nb = null_basis (A, p)
  % One null-space vector per free column, by back substitution.
  [R, pivots] = rref_p (A, p);
  n = size (A, 2);
  free = setdiff (1:n, pivots);
  Nb = zeros (numel (free), n);
  for t = 1:numel (free)
    v = zeros (1, n);
    v(free(t)) = 1;
    for i = 1:numel (pivots)
      v(pivots(i)) = mod (-R(i, free(t)), p);
    end
    Nb(t, :) = v;
  end
end

function s = words_text (V)
  parts = cell (1, size (V, 1));
  for i = 1:size (V, 1)
    parts{i} = word (V(i, :));
  end
  s = strjoin (parts, ' ');
end

function Z = null_space (A, p)
  V = all_vectors (size (A, 2), p);
  Z = V(all (mod (A * V', p) == 0, 1), :);
end

function W = span_words (M, p)
  % The row space of M, one word per row, sorted and without repeats.
  C = all_vectors (size (M, 1), p);
  W = unique (mod (C * M, p), 'rows');
end

function x = gf_solve (A, b, p)
  % A solution with every free variable set to 0, or [] when there is none.
  n = size (A, 2);
  [R, pivots] = rref_p ([A, b(:)], p);
  x = [];
  if (any (pivots == n + 1))
    return;
  end
  x = zeros (1, n);
  for i = 1:numel (pivots)
    x(pivots(i)) = R(i, n + 1);
  end
end

function s = solve_text (A, b, p)
  x = gf_solve (A, b, p);
  if (isempty (x))
    s = 'None';
  else
    s = word (x);
  end
end

function S = solutions (A, b, p)
  % Every solution of A x = b, as words in counting order.
  V = all_vectors (size (A, 2), p);
  S = V(all (mod (A * V', p) == b(:), 1), :);
end

function X = gf_inverse (A, p)
  n = size (A, 1);
  R = rref_p ([A, eye(n)], p);
  X = [];
  if (isequal (R(:, 1:n), eye (n)))
    X = R(:, n+1:end);
  end
end

function [S, order] = systematic (G, p)
  % Row-reduce G, then move the pivot columns to the front: the form [I | P].
  [R, pivots] = rref_p (G, p);
  order = [pivots, setdiff(1:size (G, 2), pivots)];
  S = R(1:numel (pivots), order);
end

function y = run_network (gates, outs, x)
  % Signals 1 to 4 are the inputs, and gate g writes signal 4 + g.
  s = x(:)';
  for g = 1:size (gates, 1)
    s(end+1) = bitxor (s(gates(g, 1)), s(gates(g, 2)));
  end
  y = s(outs);
end


% Counting invertible matrices ----------------------------------------------------

function counts = rank_counts (n)
  % How many n-by-n binary matrices have each rank 0..n, by trying all of
  % them at once: each row is an n-bit number, and a basis of leading bits is
  % built for every matrix in parallel.
  N = 2 ^ (n * n);
  codes = (0:N-1)';
  mask = 2 ^ n - 1;
  lead = zeros (N, n);
  r = zeros (N, 1);
  for i = 0:n-1
    v = bitand (floor (codes / 2 ^ (n * i)), mask);
    for h = n-1:-1:0
      has = bitand (v, 2 ^ h) > 0;
      ex = lead(:, h + 1) > 0;
      a = has & ex;
      v(a) = bitxor (v(a), lead(a, h + 1));
      b = has & ~ex;
      lead(b, h + 1) = v(b);
      r(b) = r(b) + 1;
      v(b) = 0;
    end
  end
  counts = accumarray (r + 1, 1, [n + 1, 1])';
end

function c = invertible_count (n, q)
  c = 1;
  for i = 0:n-1
    c = c * (q ^ n - q ^ i);
  end
end

function f = invertible_fraction (n)
  f = 1.0;
  for j = 1:n
    f = f * (1 - 2 ^ -j);
  end
end

function s = f6 (x)
  s = sprintf ('%.6f', x);
end

function ok = all_3by4 (test)
  % Whether test holds on every binary 3-by-4 matrix.
  ok = true;
  for code = 0:4095
    if (~test (matrix_from_code (code, 3, 4)))
      ok = false;
      return;
    end
  end
end


% The chapter's matrices ------------------------------------------------------------

A = parse_rows ('1101 1011 0111');
X = [1 0 1 1];
E2 = parse_rows ('11 01');
F2 = parse_rows ('10 11');
B = parse_rows ('011011 101100 110111');
G = parse_rows ('0110111 1110001 1101010 1100101');
M = parse_rows ('110 011 111');
U = parse_rows ('110 011 101');
M5 = parse_rows ('213 131 412');
B5 = [3 0 2];
SHARED = [3 4; 1 2; 6 4; 1 5; 2 5];
SHARED_OUT = [7 8 9];

C = parse_rows ('10110 11001 01100');
CB = [1 0 1];
Q = parse_rows ('1100 1110 0111 0011');
D = parse_rows ('110101 011110 101100');
T3 = parse_rows ('211 120 112');
TB = [0 1 2];
AE = parse_rows ('10110 01101');
XE = [1 0 1 1 0];
K = parse_rows ('11010 01101 10111');
J2 = parse_rows ('11 10');
M5E = parse_rows ('23 11');
B5E = [1 2];


% The checks ----------------------------------------------------------------------

section ('Matrix times vector as parity sums');
same ('A', show (A), '1101 1011 0111');
same ('A x for x = 1011, one parity sum per row', word (mat_vec (A, X, 2)), '010');
same ('XOR of the columns of A that x selects (1, 3, 4)', word (mod (sum (A(:, X == 1), 2), 2)'), '010');
same ('row weights of A', join_numbers (sum (A, 2)'), '3 3 3');
same ('two-input XOR gates, one chain per row', sum (sum (A, 2) - 1), 6);
V4 = all_vectors (4, 2);
ok = true;
for t = 1:16
  ok = ok && isequal (run_network (SHARED, SHARED_OUT, V4(t, :)), mat_vec (A, V4(t, :), 2));
end
require ('the 5-gate network with x3 XOR x4 shared gives A x on all 16 inputs', ok);
same ('gates in the shared network', size (SHARED, 1), 5);

section ('Multiplying two matrices');
same ('E F', show (mod (E2 * F2, 2)), '01 11');
same ('F E', show (mod (F2 * E2, 2)), '11 10');
require ('E F differs from F E', ~isequal (mod (E2 * F2, 2), mod (F2 * E2, 2)));
mats = cell (1, 16);
for t = 0:15
  mats{t + 1} = matrix_from_code (t, 2, 2);
end
ok = true;
for a = 1:16
  for b = 1:16
    PQ = mod (mats{a} * mats{b}, 2);
    for c = 1:16
      ok = ok && isequal (mod (PQ * mats{c}, 2), mod (mats{a} * mod (mats{b} * mats{c}, 2), 2));
    end
  end
end
require ('(P Q) R = P (Q R) for all 4,096 triples of binary 2-by-2 matrices', ok);
same ('number of triples checked', numel (mats) ^ 3, 4096);
require ('I A = A and A I = A for the 3-by-4 A', isequal (mod (eye (3) * A, 2), A) && isequal (mod (A * eye (4), 2), A));
same ('message 1011 times G, the XOR of rows 1, 3 and 4 of G', word (mod ([1 0 1 1] * G, 2)), '0111000');

section ('Three row operations');
ops2 = zeros (0, 4);
ops5 = zeros (0, 4);
for i = 1:3
  for j = 1:3
    if (i ~= j)
      ops2(end+1, :) = [1 i j 0];
      ops2(end+1, :) = [3 i 1 j];
      ops5(end+1, :) = [1 i j 0];
      for f = 1:4
        ops5(end+1, :) = [3 i f j];
      end
    end
  end
  for c = 1:4
    ops5(end+1, :) = [2 i c 0];
  end
end
ok = true;
for t = 1:size (ops2, 1)
  ok = ok && isequal (apply_step (apply_step (A, ops2(t, :), 2), ops2(t, :), 2), A);
end
require ('every swap and row addition on A is undone by repeating it', ok);
ok = true;
for t = 1:size (ops5, 1)
  ok = ok && isequal (apply_step (apply_step (M5, ops5(t, :), 5), undo_step (ops5(t, :), 5), 5), M5);
end
require ('every row operation on the GF(5) matrix is undone by its reverse', ok);
same ('reverse of R1 <- 2 R1 over GF(5)', step_text (undo_step ([2 1 2 0], 5), 5), 'R1 <- 3 R1');
ok = true;
for t = 1:size (ops2, 1)
  ok = ok && isequal (mod (apply_step (eye (3), ops2(t, :), 2) * A, 2), apply_step (A, ops2(t, :), 2));
end
for t = 1:size (ops5, 1)
  ok = ok && isequal (mod (apply_step (eye (3), ops5(t, :), 5) * M5, 5), apply_step (M5, ops5(t, :), 5));
end
require ('each operation is left multiplication by the operated identity', ok);
same ('elementary matrix of R2 <- R2 + R1 on three rows', show (apply_step (eye (3), [3 2 1 1], 2)), '100 110 001');
names = {'A', 'B', 'G'};
mlist = {A, B, G};
for t = 1:3
  R = rref_p (mlist{t}, 2);
  require (sprintf ('elimination keeps the null space of %s', names{t}), isequal (null_space (mlist{t}, 2), null_space (R, 2)));
  require (sprintf ('elimination keeps the row space of %s', names{t}), isequal (span_words (mlist{t}, 2), span_words (R, 2)));
end

section ('Elimination of the 3-by-4 matrix A');
[R, pivots, steps, stages] = rref_p (A, 2);
same ('operations, in order', steps_text (steps, 2), 'R2 <- R2 + R1; R1 <- R1 + R2; R3 <- R3 + R2; R1 <- R1 + R3');
same ('after column 1', show (stages{1}{3}), '1101 0110 0111');
same ('after column 2', show (stages{2}{3}), '1011 0110 0001');
same ('after column 4, the RREF', show (R), '1010 0110 0001');
same ('pivot columns', join_numbers (pivots), '1 2 4');
same ('free column', join_numbers (setdiff (1:4, pivots)), '3');
require ('the result is in reduced row echelon form', is_rref (R));
require ('elimination ends in RREF on all 4,096 binary 3-by-4 matrices', all_3by4 (@(Mx) is_rref (rref_p (Mx, 2))));

section ('Rank and the rank theorem');
[R, pivots, steps, stages] = rref_p (B, 2);
same ('B', show (B), '011011 101100 110111');
same ('operations on B', steps_text (steps, 2), 'R1 <-> R2; R3 <- R3 + R1; R3 <- R3 + R2');
same ('B after column 1', show (stages{1}{3}), '101100 011011 011011');
same ('RREF of B', show (R), '101100 011011 000000');
names = {'A', 'B', 'G', 'M', 'C'};
mlist = {A, B, G, M, C};
stated = {'4 3 1', '6 2 4', '7 4 3', '3 3 0', '5 3 2'};
for t = 1:5
  Mx = mlist{t};
  same (sprintf ('%s: columns, rank, nullity', names{t}), sprintf ('%d %d %d', size (Mx, 2), gf_rank (Mx, 2), size (null_basis (Mx, 2), 1)), stated{t});
end
require ('rank + nullity = 4 and the null space has 2^nullity words, all 3-by-4', all_3by4 (@(Mx) size (null_space (Mx, 2), 1) == 2 ^ (4 - gf_rank (Mx, 2))));
require ('the row space has 2^rank words, all 3-by-4', all_3by4 (@(Mx) size (span_words (Mx, 2), 1) == 2 ^ gf_rank (Mx, 2)));

section ('The null space by back substitution');
basis = null_basis (B, 2);
same ('null-space basis of B, one word per free column 3, 4, 5, 6', words_text (basis), '111000 100100 010010 010001');
require ('B v = 0 for each basis word', all (all (mod (B * basis', 2) == 0)));
same ('words in the null space of B, by trying all 64', size (null_space (B, 2), 1), 16);
same ('their count from the rank theorem, 2^(6 - 2)', 2 ^ (6 - gf_rank (B, 2)), 16);
require ('the span of the basis is the whole null space', isequal (span_words (basis, 2), null_space (B, 2)));
same ('null-space basis of A', words_text (null_basis (A, 2)), '1110');

section ('Solving a linear system');
same ('RREF of [B | 100]', show (rref_p ([B, [1; 0; 0]], 2)), '1011000 0110110 0000001');
same ('solution of B x = 100', solve_text (B, [1 0 0], 2), 'None');
same ('solutions of B x = 100, by trying all 64', size (solutions (B, [1 0 0], 2), 1), 0);
same ('RREF of [B | 110]', show (rref_p ([B, [1; 1; 0]], 2)), '1011001 0110111 0000000');
same ('particular solution of B x = 110', solve_text (B, [1 1 0], 2), '110000');
same ('solutions of B x = 110, by trying all 64', size (solutions (B, [1 1 0], 2), 1), 16);
V3 = all_vectors (3, 2);
count = 0;
for t = 1:8
  count = count + (size (solutions (B, V3(t, :), 2), 1) > 0);
end
same ('right sides b with B x = b solvable', count, 4);
same ('solution of M x = 110', solve_text (M, [1 1 0], 2), '101');
same ('solutions of M x = 110, by trying all 8', words_text (solutions (M, [1 1 0], 2)), '101');
parts = cell (1, 8);
ok = true;
for t = 1:8
  S = solutions (A, V3(t, :), 2);
  parts{t} = sprintf ('%s: %s', word (V3(t, :)), words_text (S));
  ok = ok && size (S, 1) == 2 && isequal (mod (S(1, :) + S(2, :), 2), [1 1 1 0]);
end
same ('A x = b: solutions for each b, b = 000 to 111', strjoin (parts, '; '), ['000: 0000 1110; 001: 0111 1001; 010: 0101 1011; 011: 0010 1100; ' '100: 0011 1101; 101: 0100 1010; 110: 0110 1000; 111: 0001 1111']);
require ('each set is a particular solution plus {0000, 1110}', ok);

section ('Inverting a square matrix');
[R, pivots, steps, stages] = rref_p ([M, eye(3)], 2);
same ('operations on [M | I]', steps_text (steps, 2), 'R3 <- R3 + R1; R1 <- R1 + R2; R1 <- R1 + R3; R2 <- R2 + R3');
same ('[M | I] after column 1', show (stages{1}{3}), '110100 011010 001101');
same ('[M | I] after column 2', show (stages{2}{3}), '101110 011010 001101');
same ('[M | I] after column 3', show (stages{3}{3}), '100011 010111 001101');
Minv = gf_inverse (M, 2);
same ('inverse of M', show (Minv), '011 111 101');
require ('M times its inverse is I, both ways', isequal (mod (M * Minv, 2), eye (3)) && isequal (mod (Minv * M, 2), eye (3)));
same ('RREF of U', show (rref_p (U, 2)), '101 011 000');
same ('rank of U', gf_rank (U, 2), 2);
if (isempty (gf_inverse (U, 2)))
  ntext = 'None';
else
  ntext = show (gf_inverse (U, 2));
end
same ('inverse of U', ntext, 'None');
same ('solution of M x = 110 as the inverse times 110', word (mat_vec (Minv, [1 1 0], 2)), '101');

section ('How many square matrices are invertible');
stated = [1 6 168 20160];
for n = 1:4
  counts = rank_counts (n);
  same (sprintf ('invertible %d-by-%d binary matrices, by trying all', n, n), counts(n + 1), stated(n));
end
stated = [1 6 168 20160 9999360 20158709760 163849992929280];
for n = 1:7
  same (sprintf ('product formula, n = %d', n), invertible_count (n, 2), stated(n));
end
fracs = {'0.500000', '0.375000', '0.328125', '0.307617', '0.298004', '0.293348', '0.291056'};
for n = 1:7
  same (sprintf ('fraction invertible, n = %d', n), f6 (invertible_fraction (n)), fracs{n});
  require (sprintf ('product of (1 - 2^-j) equals count / 2^(n^2), n = %d', n), abs (invertible_fraction (n) - invertible_count (n, 2) / 2 ^ (n * n)) < 1e-15);
end
same ('fraction invertible, n = 10', f6 (invertible_fraction (10)), '0.289070');
same ('the limit of the products, to six places', f6 (invertible_fraction (200)), '0.288788');
ok = true;
for n = 1:60
  ok = ok && invertible_fraction (n) > 0.25;
end
require ('the fraction exceeds 1/4 for every n from 1 to 60', ok);
same ('3-by-3 binary matrices of rank 0, 1, 2, 3', join_numbers (rank_counts (3)), '1 49 294 168');
counts4 = rank_counts (4);
same ('4-by-4 binary matrices of rank 0, 1, 2, 3, 4', join_numbers (counts4), '1 225 7350 37800 20160');
same ('total 4-by-4 binary matrices', sum (counts4), 65536);

section ('Systematic form');
[R, pivots, steps, stages] = rref_p (G, 2);
same ('G', show (G), '0110111 1110001 1101010 1100101');
same ('operations on G', steps_text (steps, 2), ['R1 <-> R2; R3 <- R3 + R1; R4 <- R4 + R1; R1 <- R1 + R2; ' 'R2 <- R2 + R3; R4 <- R4 + R3; R2 <- R2 + R4; R3 <- R3 + R4']);
same ('G after column 1', show (stages{1}{3}), '1110001 0110111 0011011 0010100');
same ('G after column 2', show (stages{2}{3}), '1000110 0110111 0011011 0010100');
same ('G after column 3', show (stages{3}{3}), '1000110 0101100 0011011 0001111');
same ('G after column 4, the form [I | P]', show (R), '1000110 0100011 0010100 0001111');
same ('operations on G, of which row additions', sum (steps(:, 1) == 3), 7);
same ('P', show (R(:, 5:7)), '110 011 100 111');
[S, order] = systematic (G, 2);
require ('no column swap is needed for G', isequal (order, 1:7));
require ('the leading 4-by-4 block of G is invertible', ~isempty (gf_inverse (G(:, 1:4), 2)));
require ('G and [I | P] have the same row space', isequal (span_words (G, 2), span_words (S, 2)));
same ('codewords in the row space', size (span_words (G, 2), 1), 16);
same ('codeword of message 1011 from [I | P]', word (mod ([1 0 1 1] * S, 2)), '1011101');
same ('message 1011 times P', word (mod ([1 0 1 1] * S(:, 5:7), 2)), '101');
[S, order] = systematic (A, 2);
same ('column order that puts the pivots of A first', join_numbers (order), '1 2 4 3');
same ('A after the column swap, the form [I | P]', show (S), '1001 0101 0010');
require ('the leading 3-by-3 block of A is singular', isempty (gf_inverse (A(:, 1:3), 2)));
require ('A and [I | P] have different row spaces', ~isequal (span_words (A, 2), span_words (S, 2)));
WA = span_words (A, 2);
require ('swapping positions 3 and 4 of every word of span(A) gives span([I | P])', isequal (unique (WA(:, [1 2 4 3]), 'rows'), span_words (S, 2)));

section ('Elimination over GF(5)');
same ('inverses of 1, 2, 3, 4 modulo 5', join_numbers (arrayfun (@(a) inv_mod (a, 5), 1:4)), '1 3 2 4');
[R, pivots, steps, stages] = rref_p ([M5, B5'], 5);
same ('operations on [M5 | b]', steps_text (steps, 5), ['R1 <- 3 R1; R2 <- R2 - R1; R3 <- R3 - 4 R1; R2 <-> R3; R2 <- 4 R2; ' 'R1 <- R1 - 3 R2; R3 <- 3 R3; R1 <- R1 - 2 R3; R2 <- R2 - 4 R3']);
same ('after column 1', show (stages{1}{3}), '1344 0021 0411');
same ('after column 2', show (stages{2}{3}), '1022 0144 0021');
same ('after column 3', show (R), '1001 0102 0013');
x = gf_solve (M5, B5, 5);
same ('solution', word (x), '123');
same ('M5 x', word (mat_vec (M5, x, 5)), '302');
same ('solutions over GF(5), by trying all 125', words_text (solutions (M5, B5, 5)), '123');

section ('Elimination in hardware');
same ('R2 + R1 for A, bit by bit', word (mod (A(2, :) + A(1, :), 2)), '0110');
adds_a = row_additions (A, 2);
same ('row additions in the elimination of A', adds_a, 4);
same ('bit XORs in the elimination of A', adds_a * 4, 16);
same ('bit XORs in the elimination of G', row_additions (G, 2) * 7, 49);
same ('bound min(m, n) (m - 1) n on bit XORs for 4-by-7', min (4, 7) * 3 * 7, 84);
worst = 0;
sixes = 0;
for code = 0:4095
  adds = row_additions (matrix_from_code (code, 3, 4), 2);
  worst = max (worst, adds);
  sixes = sixes + (adds == 6);
end
require ('row additions never exceed min(m, n)(m - 1) = 6, all 3-by-4', worst <= 6);
same ('most row additions over all 3-by-4 binary matrices', worst, 6);
same ('3-by-4 binary matrices that use all 6 additions', sixes, 15);

section ('Worked examples');
[R, pivots, steps] = rref_p ([C, CB'], 2);
same ('C', show (C), '10110 11001 01100');
same ('operations on [C | 101]', steps_text (steps, 2), 'R2 <- R2 + R1; R3 <- R3 + R2; R1 <- R1 + R3; R2 <- R2 + R3');
same ('RREF of [C | 101]', show (R), '101011 011001 000110');
same ('particular solution', solve_text (C, CB, 2), '11000');
same ('null-space basis of C', words_text (null_basis (C, 2)), '11100 10011');
same ('all solutions of C x = 101', words_text (solutions (C, CB, 2)), '00100 01011 10111 11000');
[R, pivots, steps] = rref_p ([Q, eye(4)], 2);
same ('operations on [Q | I]', steps_text (steps, 2), ['R2 <- R2 + R1; R2 <-> R3; R1 <- R1 + R2; R1 <- R1 + R3; R2 <- R2 + R3; ' 'R4 <- R4 + R3; R1 <- R1 + R4; R2 <- R2 + R4']);
Qinv = gf_inverse (Q, 2);
same ('inverse of Q', show (Qinv), '1011 0011 1100 1101');
require ('Q times its inverse is I', isequal (mod (Q * Qinv, 2), eye (4)));
[R, pivots, steps] = rref_p (D, 2);
same ('operations on D', steps_text (steps, 2), 'R3 <- R3 + R1; R1 <- R1 + R2; R3 <- R3 + R2; R2 <- R2 + R3');
same ('RREF of D', show (R), '101011 011001 000111');
[S, order] = systematic (D, 2);
same ('column order for D', join_numbers (order), '1 2 4 3 5 6');
same ('D in the form [I | P] after the swap', show (S), '100111 010101 001011');
[R, pivots, steps, stages] = rref_p ([T3, TB'], 3);
same ('operations on [T | 012] over GF(3)', steps_text (steps, 3), ['R1 <- 2 R1; R2 <- R2 - R1; R3 <- R3 - R1; R2 <-> R3; R2 <- 2 R2; ' 'R1 <- R1 - 2 R2; R1 <- R1 - 2 R3']);
same ('[T | 012] after column 1', show (stages{1}{3}), '1220 0011 0202');
same ('RREF of [T | 012]', show (R), '1002 0101 0011');
same ('solution over GF(3)', solve_text (T3, TB, 3), '211');
same ('T x', word (mat_vec (T3, [2 1 1], 3)), '012');

section ('Exercises');
same ('parity sums of 10110 01101 on x = 10110', word (mat_vec (AE, XE, 2)), '11');
same ('XOR gates for those parity sums', sum (sum (AE, 2) - 1), 4);
same ('K', show (K), '11010 01101 10111');
same ('rank and nullity of K', sprintf ('%d %d', gf_rank (K, 2), size (null_basis (K, 2), 1)), '2 3');
same ('words in the null space of K', size (null_space (K, 2), 1), 8);
s3 = [3 2 2 1];
T = apply_step (apply_step (apply_step (T3, s3, 3), s3, 3), s3, 3);
require ('over GF(3), adding R1 to R2 three times restores the matrix', isequal (T, T3));
require ('over GF(3), adding R1 to R2 twice does not', ~isequal (apply_step (apply_step (T3, s3, 3), s3, 3), T3));
same ('inverse of 11 10', show (gf_inverse (J2, 2)), '01 11');
V9 = all_vectors (4, 3);
inv3 = 0;
for t = 1:81
  inv3 = inv3 + (gf_rank ([V9(t, 1:2); V9(t, 3:4)], 3) == 2);
end
same ('invertible 2-by-2 matrices over GF(3), by trying all 81', inv3, 48);
same ('(3^2 - 1)(3^2 - 3)', invertible_count (2, 3), 48);
same ('solution of 2x + 3y = 1, x + y = 2 over GF(5)', solve_text (M5E, B5E, 5), '02');
[S, order] = systematic ([M, eye(3)], 2);
same ('systematic form of [M | I]', show (S), '100011 010111 001101');

printf ('\nAll chapter 16 numbers reproduced.\n');
