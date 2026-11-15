% Reproduce every number in chapter 24 of Error Correction Coding,
% "Modifying and Combining Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch24_modifying_and_combining_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch24_modifying_and_combining_codes.py and
% prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in this chapter is a whole number, a binary word, a
% binary matrix, a weight distribution or a rate printed to three places, so
% every comparison is exact. Every code is small enough to list in full, so
% each minimum distance is found by computing the weight of every codeword
% rather than taken from a theorem. The largest lists are the 65,536
% codewords of the (25,16) rectangular parity code and of the (49,16) product
% of the Hamming code with itself.
%
% A word is a row of 0s and 1s, and a code is the matrix of its distinct
% codewords, one per row, sorted the way Python sorts tuples. A list of
% messages is in the order Python's itertools.product gives, the first bit
% the most significant.

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


% Words, matrices and codes over GF(2) -------------------------------------------

function M = parse (text)
  % Rows written as digit strings separated by spaces, such as "1101 1011".
  rows = strsplit (text, ' ');
  M = zeros (numel (rows), numel (rows{1}));
  for i = 1:numel (rows)
    M(i, :) = rows{i} - '0';
  end
end

function s = word (v)
  s = char (v(:)' + '0');
end

function s = show (M)
  parts = cell (1, rows (M));
  for i = 1:rows (M)
    parts{i} = word (M(i, :));
  end
  s = strjoin (parts, ' ');
end

function s = join_numbers (v)
  s = strjoin (arrayfun (@(x) sprintf ('%d', x), v, 'UniformOutput', false), ' ');
end

function P = mat_mul (A, B)
  P = mod (A * B, 2);
end

function r = gf2_rank (A)
  % Rank over GF(2) by elimination.
  r = 0;
  for c = 1:columns (A)
    pivot = find (A(r + 1:end, c), 1);
    if (isempty (pivot))
      continue;
    end
    pivot = pivot + r;
    A([r + 1, pivot], :) = A([pivot, r + 1], :);
    for i = 1:rows (A)
      if (i ~= r + 1 && A(i, c))
        A(i, :) = xor (A(i, :), A(r + 1, :));
      end
    end
    r = r + 1;
    if (r == rows (A))
      break;
    end
  end
end

function M = messages (k)
  % All 2^k messages, in itertools.product order.
  M = dec2bin (0:2 ^ k - 1, k) - '0';
end

function S = span (G)
  % Every codeword m G, in message order.
  S = mod (messages (rows (G)) * G, 2);
end

function C = code (G)
  C = unique (span (G), 'rows');
end

function s = wdist (C)
  w = sum (C, 2);
  u = unique (w)';
  parts = arrayfun (@(x) sprintf ('%d:%d', x, sum (w == x)), u, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function d = dmin (C)
  w = sum (C, 2);
  d = min (w(w > 0));
end

function k = dimension (C)
  [~, e] = log2 (rows (C));
  k = e - 1;
end

function s = params (C)
  s = sprintf ('(%d,%d,%d)', columns (C), dimension (C), dmin (C));
end

function P = puncture (C, i)
  % Delete position i.
  P = unique (C(:, [1:i - 1, i + 1:end]), 'rows');
end

function S = shorten (C, i)
  % Keep the codewords that are 0 at position i, then delete it.
  K = C(C(:, i) == 0, :);
  S = unique (K(:, [1:i - 1, i + 1:end]), 'rows');
end

function E = extend (C)
  E = sortrows ([C, mod(sum (C, 2), 2)]);
end

function X = expurgate (C)
  X = sortrows (C(mod (sum (C, 2), 2) == 0, :));
end

function A = augment (C)
  A = unique ([C; 1 - C], 'rows');
end

function L = lengthen (C)
  L = extend (augment (C));
end

function D = direct_sum (C1, C2)
  D = zeros (rows (C1) * rows (C2), columns (C1) + columns (C2));
  r = 0;
  for i = 1:rows (C1)
    for j = 1:rows (C2)
      r = r + 1;
      D(r, :) = [C1(i, :), C2(j, :)];
    end
  end
  D = sortrows (D);
end

function D = u_u_plus_v (C1, C2)
  D = zeros (rows (C1) * rows (C2), 2 * columns (C1));
  r = 0;
  for i = 1:rows (C1)
    for j = 1:rows (C2)
      r = r + 1;
      D(r, :) = [C1(i, :), xor(C1(i, :), C2(j, :))];
    end
  end
  D = sortrows (D);
end

function X = product_array (G1, G2, M)
  % The n1-by-n2 array G1^T M G2: rows encoded by the row code C2, columns by C1.
  X = mat_mul (mat_mul (G1', M), G2);
end

function X = rows_then_columns (G1, G2, M)
  R = mat_mul (M, G2);
  X = mat_mul (R', G1)';
end

function X = columns_then_rows (G1, G2, M)
  X = mat_mul (mat_mul (G1', M), G2);
end

function v = flatten (A)
  % Read an array row by row.
  v = reshape (A', 1, []);
end

function P = product_code (G1, G2)
  % Every product array, read row by row, sorted.
  k1 = rows (G1);
  k2 = rows (G2);
  Ms = messages (k1 * k2);
  P = zeros (rows (Ms), columns (G1) * columns (G2));
  for t = 1:rows (Ms)
    M = reshape (Ms(t, :), k2, k1)';
    P(t, :) = flatten (product_array (G1, G2, M));
  end
  P = sortrows (P);
end

function P = product_generator (G1, G2)
  % Rows are the outer products of a row of G1 and a row of G2, read row by row.
  P = zeros (rows (G1) * rows (G2), columns (G1) * columns (G2));
  r = 0;
  for i = 1:rows (G1)
    for j = 1:rows (G2)
      r = r + 1;
      P(r, :) = kron (G1(i, :), G2(j, :));
    end
  end
end

function w = big_weights (G)
  % The weight of every codeword of a larger code.
  w = sum (mod (messages (rows (G)) * G, 2), 2);
end

function g = gates (checks)
  % XOR gates for a set of parity equations, one tree per distinct equation.
  keys = cellfun (@(e) mat2str (unique (e)), checks, 'UniformOutput', false);
  [~, idx] = unique (keys);
  g = 0;
  for i = idx(:)'
    g = g + max (numel (unique (checks{i})) - 1, 0);
  end
end

function p = check_positions (G)
  % Positions whose column of a systematic G is not a unit vector, so they hold checks.
  p = find (sum (G, 1) ~= 1);
end

function eqs = check_equations (G, positions, kept)
  % The message bits that each check position sums, read down its column of G.
  % Message bits outside kept are tied to 0 and drop out of every equation.
  if (nargin < 3)
    kept = 1:rows (G);
  end
  eqs = cell (1, numel (positions));
  for t = 1:numel (positions)
    e = find (G(:, positions(t)))';
    eqs{t} = e(ismember (e, kept));
  end
end

function s = rate (k, n)
  s = sprintf ('%.3f', k / n);
end

function t = has_row (C, v)
  t = ismember (v, C, 'rows');
end

function t = all_zero (A)
  t = all (A(:) == 0);
end


% Checks -------------------------------------------------------------------------

function check_hamming_code (G7, H7, C7)
  section ('The (7,4) Hamming code of Linear Codes');
  same ('generator G', show (G7), '1000110 0100011 0010101 0001111');
  same ('parity-check matrix H', show (H7), '1011100 1101010 0111001');
  require ('H G^T = 0', all_zero (mat_mul (H7, G7')));
  same ('number of codewords', rows (C7), 16);
  same ('weight distribution', wdist (C7), '0:1 3:7 4:7 7:1');
  same ('minimum distance', dmin (C7), 3);
  require ('the all-ones word 1111111 is a codeword', has_row (C7, ones (1, 7)));
  same ('weights of the rows of H', join_numbers (sum (H7, 2)'), '4 4 4');
  pairs = [1 2; 1 3; 2 3];
  shared = arrayfun (@(t) sum (H7(pairs(t, 1), :) & H7(pairs(t, 2), :)), 1:3);
  same ('1s shared by rows 1 and 2, 1 and 3, 2 and 3 of H', join_numbers (shared), '2 2 2');
  require ('H H^T = 0', all_zero (mat_mul (H7, H7')));
  same ('rank of H', gf2_rank (H7), 3);
end

function check_memory_word ()
  section ('A 64-bit memory word');
  same ('Hamming code with m = 6 checks, (n,k)', sprintf ('(%d,%d)', 2 ^ 6 - 1, 2 ^ 6 - 1 - 6), '(63,57)');
  same ('Hamming code with m = 7 checks, (n,k)', sprintf ('(%d,%d)', 2 ^ 7 - 1, 2 ^ 7 - 1 - 7), '(127,120)');
  same ('data positions to remove from 120', 120 - 64, 56);
  same ('shortened Hamming code, (n,k)', sprintf ('(%d,%d)', 127 - 56, 120 - 56), '(71,64)');
  same ('extended shortened code, (n,k)', sprintf ('(%d,%d)', 127 - 56 + 1, 120 - 56), '(72,64)');
  same ('extended Hamming code, then shortened, (n,k)', sprintf ('(%d,%d)', 128 - 56, 120 - 56), '(72,64)');
  same ('check bits of the (72,64) code', 72 - 64, 8);
  same ('C(72,2)', floor (72 * 71 / 2), 2556);
  same ('error patterns of weight at most 2 in 72 bits', 1 + 72 + floor (72 * 71 / 2), 2629);
  same ('syndromes of 8 check bits', 2 ^ 8, 256);
  require ('2629 patterns exceed 256 syndromes, so d is at most 4', 1 + 72 + floor (72 * 71 / 2) > 2 ^ 8);
  same ('rate of the (127,120) code', rate (120, 127), '0.945');
  same ('rate of the (72,64) code', rate (64, 72), '0.889');
end

function check_extending (G7, H7, C7, E8, G43, EVEN4)
  section ('Extending');
  GE = [G7, mod(sum (G7, 2), 2)];
  HE = [H7, zeros(3, 1); ones(1, 8)];
  same ('extended generator', show (GE), '10001101 01000111 00101011 00011110');
  same ('extended parity-check matrix', show (HE), '10111000 11010100 01110010 11111111');
  require ('extended H times extended G^T = 0', all_zero (mat_mul (HE, GE')));
  same ('rank of the extended parity-check matrix', gf2_rank (HE), 4);
  same ('extended code, weight distribution', wdist (E8), '0:1 4:14 8:1');
  same ('extended code, minimum distance', dmin (E8), 4);
  require ('every extended codeword has even weight', all (mod (sum (E8, 2), 2) == 0));
  require ('the extended code equals the span of the extended generator', isequal (code (GE), E8));
  same ('weight-3 codeword 1000110 extended', word (GE(1, :)), '10001101');
  same ('overall parity bit as message coefficients', word (GE(:, 8)), '1110');
  M = messages (4);
  S = span (G7);
  require ('c5 XOR c6 XOR c7 = m4 for every message', all (xor (xor (S(:, 5), S(:, 6)), S(:, 7)) == M(:, 4)));
  require ('overall parity equals m1 XOR m2 XOR m3 for every message', ...
           all (mod (sum (S, 2), 2) == mod (M(:, 1) + M(:, 2) + M(:, 3), 2)));
  same ('XOR gates in the overall parity tree', gates (check_equations (GE, 8)), 2);
  GEV = [G43, mod(sum (G43, 2), 2)];
  same ('even-weight (4,3) code, extended generator', show (GEV), '10010 01010 00110');
  EE = extend (EVEN4);
  require ('every appended bit of the even-weight code is 0', all (EE(:, 5) == 0));
  same ('even-weight code, distance before and after', sprintf ('%d %d', dmin (EVEN4), dmin (EE)), '2 2');
end

function check_puncturing (G7, C7, E8, D7)
  section ('Puncturing');
  require ('puncturing the extended code at position 8 gives the (7,4) code', isequal (puncture (E8, 8), C7));
  P6 = puncture (C7, 7);
  same ('(7,4) punctured at position 7, (n,k,d)', params (P6), '(6,4,2)');
  same ('its weight distribution', wdist (P6), '0:1 2:3 3:8 4:3 6:1');
  same ('weight-2 word from the row 0100011', word (G7(2, 1:6)), '010001');
  three = C7(sum (C7, 2) == 3, :);
  same ('weight-3 codewords covering each position', join_numbers (sum (three, 1)), '3 3 3 3 3 3 3');
  same ('code D, weight distribution', wdist (D7), '0:1 3:1 4:1 7:1');
  same ('code D, (n,k,d)', params (D7), '(7,2,3)');
  PD = puncture (D7, 7);
  same ('code D punctured at position 7, (n,k,d)', params (PD), '(6,2,3)');
  same ('its weight distribution', wdist (PD), '0:1 3:2 6:1');
end

function check_shortening (G7, H7, S6, S5)
  section ('Shortening');
  G6 = parse ('100011 010101 001111');
  H6 = H7(:, 2:end);
  same ('(6,3) generator', show (G6), '100011 010101 001111');
  require ('the (6,3) code is the span of its generator', isequal (code (G6), S6));
  same ('(6,3) parity-check matrix, column 1 of H deleted', show (H6), '011100 101010 111001');
  require ('(6,3) H times G^T = 0', all_zero (mat_mul (H6, G6')));
  require ('columns of the (6,3) H are distinct and nonzero', ...
           rows (unique (H6', 'rows')) == 6 && ~any (all (H6 == 0, 1)));
  same ('(6,3) codewords in message order', show (span (G6)), ...
        '000000 001111 010101 011010 100011 101100 110110 111001');
  same ('(6,3) weight distribution', wdist (S6), '0:1 3:4 4:3');
  same ('(6,3) code, (n,k,d)', params (S6), '(6,3,3)');
  G5 = parse ('10101 01111');
  H5 = H7(:, 3:end);
  require ('the (5,2) code is the span of its generator', isequal (code (G5), S5));
  same ('(5,2) parity-check matrix, columns 1 and 2 deleted', show (H5), '11100 01010 11001');
  require ('(5,2) H times G^T = 0', all_zero (mat_mul (H5, G5')));
  same ('(5,2) weight distribution', wdist (S5), '0:1 3:2 4:1');
  same ('(5,2) code, (n,k,d)', params (S5), '(5,2,3)');
  require ('(5,2) check bits 3 and 5 are equal', all (S5(:, 3) == S5(:, 5)));
  checks7 = check_positions (G7);
  same ('check positions of the (7,4) generator', join_numbers (checks7), '5 6 7');
  same ('encoder XOR gates, (7,4)', gates (check_equations (G7, checks7)), 6);
  same ('encoder XOR gates, (6,3) with m1 tied to 0', gates (check_equations (G7, checks7, [2 3 4])), 4);
  same ('encoder XOR gates, (5,2) with m1 and m2 tied to 0', gates (check_equations (G7, checks7, [3 4])), 1);
  eq6 = check_equations (G7, checks7, [2 3 4]);
  require ('(6,3) check c7 = c5 XOR m2', isequal (eq6{3}, setxor (eq6{1}, 2)));
  same ('encoder XOR gates, (6,3) with c7 built as c5 XOR m2', gates (eq6(1:2)) + 1, 3);
end

function check_expurgating (G7, C7, SIMPLEX)
  section ('Expurgating');
  X = expurgate (C7);
  same ('even-weight codewords of the (7,4) code', rows (X), 8);
  same ('expurgated code, weight distribution', wdist (X), '0:1 4:7');
  same ('expurgated code, (n,k,d)', params (X), '(7,3,4)');
  require ('the expurgated code equals the row space of H', isequal (X, SIMPLEX));
  M = messages (4);
  S = span (G7);
  p = mod (M(:, 1) + M(:, 2) + M(:, 3), 2);
  require ('a codeword has even weight exactly when m1 XOR m2 XOR m3 = 0', ...
           all ((mod (sum (S, 2), 2) == 0) == (p == 0)));
  same ('messages kept by expurgating', sum (p == 0), 8);
  odd = find (mod (sum (G7, 2), 2))';
  same ('message bits of the odd-weight rows of G', join_numbers (odd), '1 2 3');
  same ('XOR gates for m3 = m1 XOR m2', gates ({odd(1:end - 1)}), 1);
  ab = parse ('1000110 0100011');
  a = ab(1, :);
  b = ab(2, :);
  same ('w(a), w(b), w(a AND b), w(a XOR b)', ...
        sprintf ('%d %d %d %d', sum (a), sum (b), sum (a & b), sum (xor (a, b))), '3 3 1 4');
  W = messages (7);
  w = sum (W, 2);
  holds = true;
  for i = 1:rows (W)
    lhs = sum (xor (W, repmat (W(i, :), rows (W), 1)), 2);
    rhs = w(i) + w - 2 * (W * W(i, :)');
    holds = holds && all (lhs == rhs);
  end
  require ('w(a XOR b) = w(a) + w(b) - 2 w(a AND b) for all words of length 7', holds);
end

function check_augmenting (C7, SIMPLEX, S5)
  section ('Augmenting');
  wmax = max (sum (SIMPLEX, 2));
  same ('largest weight in the (7,3) code', wmax, 4);
  require ('1111111 is not in the (7,3) code', ~has_row (SIMPLEX, ones (1, 7)));
  A = augment (SIMPLEX);
  require ('augmenting the (7,3) code gives the (7,4) code', isequal (A, C7));
  same ('min(d, n - w_max) for the (7,3) code', min (dmin (SIMPLEX), 7 - wmax), 3);
  same ('augmented code, minimum distance', dmin (A), 3);
  wmax5 = max (sum (S5, 2));
  same ('largest weight in the (5,2) code', wmax5, 4);
  same ('min(d, n - w_max) for the (5,2) code', min (dmin (S5), 5 - wmax5), 1);
  same ('augmented (5,2) code, (n,k,d)', params (augment (S5)), '(5,3,1)');
  same ('augmented (5,2) code, weight distribution', wdist (augment (S5)), '0:1 1:1 2:2 3:2 4:1 5:1');
end

function check_lengthening (SIMPLEX, E8)
  section ('Lengthening');
  L = lengthen (SIMPLEX);
  same ('lengthened (7,3) code, (n,k,d)', params (L), '(8,4,4)');
  require ('lengthening the (7,3) code gives the extended (8,4) code', isequal (L, E8));
end

function check_direct_sum (EVEN4, REP4, REP3, D7, DIRECT8)
  section ('Direct sum');
  same ('even-weight (4,3) code, weight distribution', wdist (EVEN4), '0:1 2:6 4:1');
  same ('repetition (4,1) code, weight distribution', wdist (REP4), '0:1 4:1');
  same ('direct sum, (n,k,d)', params (DIRECT8), '(8,4,2)');
  same ('direct sum, weight distribution', wdist (DIRECT8), '0:1 2:6 4:2 6:6 8:1');
  poly = conv ([1 0 6 0 1], [1 0 0 0 1]);
  same ('coefficients of (1 + 6x^2 + x^4)(1 + x^4)', join_numbers (poly), '1 0 6 0 2 0 6 0 1');
  require ('code D is the direct sum of the (3,1) and (4,1) repetition codes', isequal (direct_sum (REP3, REP4), D7));
  same ('rate of the direct sum', rate (4, 8), '0.500');
end

function check_u_u_plus_v (EVEN4, REP4, E8, UUV8)
  section ('The (u | u+v) construction');
  G = parse ('10011001 01010101 00110011 00001111');
  same ('(u | u+v) generator', show (G), '10011001 01010101 00110011 00001111');
  require ('the generator spans the (u | u+v) code', isequal (code (G), UUV8));
  same ('(u | u+v) code, (n,k,d)', params (UUV8), '(8,4,4)');
  same ('min(2 d1, d2)', min (2 * dmin (EVEN4), dmin (REP4)), 4);
  same ('(u | u+v) weight distribution', wdist (UUV8), '0:1 4:14 8:1');
  require ('its weight distribution equals the extended code''s', strcmp (wdist (UUV8), wdist (E8)));
  same ('upper-bound word (u, u) with u = 1001', word (repmat ([1 0 0 1], 1, 2)), '10011001');
  same ('upper-bound word (0, v) with v = 1111', word ([zeros(1, 4), ones(1, 4)]), '00001111');
  U = EVEN4(any (EVEN4, 2), :);
  v_zero = sum ([U, U], 2);
  V = REP4(any (REP4, 2), :);
  v_nonzero = [];
  for i = 1:rows (EVEN4)
    for j = 1:rows (V)
      v_nonzero(end + 1) = sum ([EVEN4(i, :), xor(EVEN4(i, :), V(j, :))]);
    end
  end
  same ('least weight with v = 0, then with v nonzero', sprintf ('%d %d', min (v_zero), min (v_nonzero)), '4 4');
  order = [1 2 3 5 7 6 8 4];
  require ('positions read as 1 2 3 5 7 6 8 4 give the extended code', isequal (sortrows (UUV8(:, order)), E8));
  u = [0 1 1 0];
  v = [1 1 1 1];
  same ('encoding u = 0110, v = 1111', word ([u, xor(u, v)]), '01101001');
  require ('01101001 is a codeword', has_row (UUV8, [u, xor(u, v)]));
end

function check_product (G32, PRODUCT9)
  section ('Product codes');
  M = [1 0; 1 1];
  X = product_array (G32, G32, M);
  same ('(3,2) x (3,2) array for data 10 11', show (X), '101 110 011');
  Ms = messages (4);
  agree = true;
  even_lines = true;
  even3 = code (G32);
  for t = 1:rows (Ms)
    D = [Ms(t, 1:2); Ms(t, 3:4)];
    agree = agree && isequal (rows_then_columns (G32, G32, D), columns_then_rows (G32, G32, D));
    A = product_array (G32, G32, D);
    even_lines = even_lines && all (ismember (A, even3, 'rows')) && all (ismember (A', even3, 'rows'));
  end
  require ('rows first and columns first agree on all 16 data arrays', agree);
  require ('every row and every column of every array is even', even_lines);
  same ('(3,2) x (3,2) product, (n,k,d)', params (PRODUCT9), '(9,4,4)');
  same ('(3,2) x (3,2) product, weight distribution', wdist (PRODUCT9), '0:1 4:9 6:6');
  same ('outer product of 110 with 110', show ([1 1 0]' * [1 1 0]), '110 110 000');
  G54 = parse ('10001 01001 00101 00011');
  G25 = product_generator (G54, G54);
  w = big_weights (G25);
  same ('rectangular parity, (n,k) from the generator''s rank', ...
        sprintf ('(%d,%d)', columns (G25), gf2_rank (G25)), '(25,16)');
  require ('rectangular parity, codewords listed once each', sum (w == 0) == 1);
  same ('rectangular parity, minimum distance', min (w(w > 0)), 4);
  same ('rectangular parity, codewords of weight 4', sum (w == 4), 100);
  eq32 = check_equations (G32, check_positions (G32));
  same ('product encoder XOR gates, row and column encoders reused', ...
        rows (G32) * gates (eq32) + columns (G32) * gates (eq32), 5);
  G9 = product_generator (G32, G32);
  same ('check positions of the (9,4) product generator', join_numbers (check_positions (G9)), '3 6 7 8 9');
  same ('product encoder XOR gates, one tree per check', gates (check_equations (G9, check_positions (G9))), 7);
end

function check_rate_and_distance (CHART_CODES)
  section ('Rate and distance together');
  n = numel (CHART_CODES);
  rates = cell (1, n);
  dists = cell (1, n);
  sizes = cell (1, n);
  for i = 1:n
    C = CHART_CODES{i};
    rates{i} = rate (dimension (C), columns (C));
    dists{i} = sprintf ('%d', dmin (C));
    sizes{i} = sprintf ('(%d,%d)', columns (C), dimension (C));
  end
  same ('rates of the chart codes', strjoin (rates, ' '), ...
        '0.571 0.500 0.667 0.500 0.400 0.429 0.500 0.500 0.444');
  same ('distances of the chart codes', strjoin (dists, ' '), '3 4 2 3 3 4 2 4 4');
  same ('(n,k) of the chart codes', strjoin (sizes, ' '), ...
        '(7,4) (8,4) (6,4) (6,3) (5,2) (7,3) (8,4) (8,4) (9,4)');
end

function check_worked_examples (E8, SIMPLEX, G43, G32, S6)
  section ('Worked examples');
  require ('shortening the extended code at position 8 gives the (7,3) code', isequal (shorten (E8, 8), SIMPLEX));
  same ('extended code shortened at position 8, (n,k,d)', params (shorten (E8, 8)), '(7,3,4)');
  same ('extended code punctured at position 8, (n,k,d)', params (puncture (E8, 8)), '(7,4,3)');
  PR = product_code (G43, G32);
  same ('(4,3) x (3,2) product, (n,k,d)', params (PR), '(12,6,4)');
  A = product_array (G43, G32, [1 0; 1 1; 0 0]);
  same ('(4,3) x (3,2) array for data 10 11 00', show (A), '101 110 000 011');
  same ('weight of that array', sum (A(:)), 6);
  same ('rate of the (12,6) product', rate (6, 12), '0.500');
  require ('111111 is not in the (6,3) code', ~has_row (S6, ones (1, 6)));
  wmax6 = max (sum (S6, 2));
  same ('largest weight in the (6,3) code', wmax6, 4);
  same ('min(d, n - w_max) for the (6,3) code', min (dmin (S6), 6 - wmax6), 2);
  same ('augmented (6,3) code, (n,k,d)', params (augment (S6)), '(6,4,2)');
  same ('lengthened (6,3) code, (n,k,d)', params (lengthen (S6)), '(7,4,2)');
  c = parse ('110110');
  require ('110110 is a (6,3) codeword', has_row (S6, c));
  same ('111111 XOR 110110', word (xor (ones (1, 6), c)), '001001');
  same ('its weight', sum (xor (ones (1, 6), c)), 2);
end

function check_exercises (G7, C7, E8, SIMPLEX, S6, PRODUCT9)
  section ('Exercises');
  same ('extended (6,3) code, (n,k,d)', params (extend (S6)), '(7,3,4)');
  same ('extended (6,3) code, weight distribution', wdist (extend (S6)), '0:1 4:7');
  PC = puncture (PRODUCT9, 9);
  same ('(9,4) product punctured at the corner, (n,k,d)', params (PC), '(8,4,3)');
  same ('weight-4 product codewords through the corner', sum (sum (PRODUCT9, 2) == 4 & PRODUCT9(:, 9) == 1), 4);
  same ('extended code shortened at position 1, (n,k,d)', params (shorten (E8, 1)), '(7,3,4)');
  require ('extend then shorten equals shorten then extend at position 1', ...
           isequal (shorten (extend (C7), 1), extend (shorten (C7, 1))));
  dual = puncture (SIMPLEX, 1);
  require ('the (6,3) shortened code is the dual of the punctured (7,3) code', ...
           all_zero (mat_mul (S6, dual')) && dimension (S6) + dimension (dual) == 6);
  EVEN8 = code ([eye(7), ones(7, 1)]);
  same ('even-weight (8,7) code, (n,k,d)', params (EVEN8), '(8,7,2)');
  same ('direct sum of (8,7,2) and (8,4,4), (n,k,d)', params (direct_sum (EVEN8, E8)), '(16,11,2)');
  same ('(u | u+v) of (8,7,2) and (8,4,4), (n,k,d)', params (u_u_plus_v (EVEN8, E8)), '(16,11,4)');
  G49 = product_generator (G7, G7);
  w = big_weights (G49);
  same ('(7,4) x (7,4) product, (n,k) from the generator''s rank', ...
        sprintf ('(%d,%d)', columns (G49), gf2_rank (G49)), '(49,16)');
  require ('(7,4) x (7,4) product, codewords listed once each', sum (w == 0) == 1);
  same ('(7,4) x (7,4) product, minimum distance', min (w(w > 0)), 9);
  same ('(7,4) x (7,4) product, codewords of weight 9', sum (w == 9), 49);
  same ('rate of the (49,16) product', rate (gf2_rank (G49), columns (G49)), '0.327');
end


% The codes of the chapter ---------------------------------------------------------

G7 = parse ('1000110 0100011 0010101 0001111');
H7 = parse ('1011100 1101010 0111001');
C7 = code (G7);
E8 = extend (C7);
SIMPLEX = code (H7);
EVEN4 = code (parse ('1001 0101 0011'));
REP4 = code (parse ('1111'));
REP3 = code (parse ('111'));
G32 = parse ('101 011');
G43 = parse ('1001 0101 0011');
D7 = code (parse ('1110000 0001111'));
S6 = shorten (C7, 1);
S5 = shorten (S6, 1);
DIRECT8 = direct_sum (EVEN4, REP4);
UUV8 = u_u_plus_v (EVEN4, REP4);
PRODUCT9 = product_code (G32, G32);

% Every code the rate-against-distance chart plots, in the chart's order.
CHART_CODES = {C7, E8, puncture(C7, 7), S6, S5, expurgate(C7), DIRECT8, UUV8, PRODUCT9};

check_hamming_code (G7, H7, C7);
check_memory_word ();
check_extending (G7, H7, C7, E8, G43, EVEN4);
check_puncturing (G7, C7, E8, D7);
check_shortening (G7, H7, S6, S5);
check_expurgating (G7, C7, SIMPLEX);
check_augmenting (C7, SIMPLEX, S5);
check_lengthening (SIMPLEX, E8);
check_direct_sum (EVEN4, REP4, REP3, D7, DIRECT8);
check_u_u_plus_v (EVEN4, REP4, E8, UUV8);
check_product (G32, PRODUCT9);
check_rate_and_distance (CHART_CODES);
check_worked_examples (E8, SIMPLEX, G43, G32, S6);
check_exercises (G7, C7, E8, SIMPLEX, S6, PRODUCT9);
printf ('\nAll chapter 24 numbers reproduced.\n');
