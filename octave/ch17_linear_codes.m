% Reproduce every number in chapter 17 of Error Correction Coding,
% "Linear Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch17_linear_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch17_linear_codes.py and prints the same lines,
% which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in this chapter is a whole number, a binary word, a
% binary matrix or a rate printed to six places, so every comparison is
% exact. The general claims are also checked exhaustively on small cases:
% closure on every subset of the eight words of length 3, G H^T = 0 and
% C = null space of H for all 4,096 blocks P of four rows and three columns,
% the column view of H x^T on all 128 words of length 7, the half lemma for
% every length up to 6, and the automorphisms of the (7,4) Hamming code on all
% 5,040 orderings of its seven positions. A run that reaches the final line
% has reproduced every number in the chapter.
%
% A word is a row of bits, a matrix is an ordinary array of 0s and 1s, a code
% is a sorted cell array of bit strings, and every sum is taken modulo 2, so
% addition is XOR and multiplication is AND. Positions are numbered from 1 on
% the left, as in the chapter. Octave's own rank works over the real numbers,
% so the rank over GF(2) carries another name.

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


% Words and matrices over GF(2) --------------------------------------------------

function M = parse_rows (text)
  % Rows written as bit strings separated by spaces, such as '10110 01101'.
  parts = strsplit (text, ' ');
  M = zeros (numel (parts), numel (parts{1}));
  for i = 1:numel (parts)
    M(i, :) = parts{i} - '0';
  end
end

function s = word (v)
  s = sprintf ('%d', v);
end

function s = show (M)
  parts = cell (1, size (M, 1));
  for i = 1:size (M, 1)
    parts{i} = sprintf ('%d', M(i, :));
  end
  s = strjoin (parts, ' ');
end

function s = join_numbers (v)
  parts = arrayfun (@(t) sprintf ('%d', t), v, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function s = join_words (code)
  s = strjoin (code, ' ');
end

function V = all_words (n)
  % Every word of n bits, one per row, in counting order with the first bit most significant.
  V = zeros (2 ^ n, n);
  for i = 0:(2 ^ n - 1)
    V(i + 1, :) = bitand (bitshift (i, -((n - 1):-1:0)), 1);
  end
end

function c = encode (m, G)
  c = mod (m * G, 2);
end

function code = rows_to_code (M)
  code = cell (1, size (M, 1));
  for i = 1:size (M, 1)
    code{i} = word (M(i, :));
  end
end

function code = span_of (G)
  % Every codeword m G, sorted, as bit strings.
  code = unique (rows_to_code (mod (all_words (size (G, 1)) * G, 2)));
end

function code = null_space (H)
  % Every word x with H x^T = 0, in counting order, as bit strings.
  W = all_words (size (H, 2));
  keep = all (mod (W * H', 2) == 0, 2);
  code = rows_to_code (W(keep, :));
end

function s = syndrome (H, x)
  s = word (mod (H * x(:), 2)');
end

function r = rank2 (A)
  r = 0;
  for c = 1:size (A, 2)
    pivot = 0;
    for i = (r + 1):size (A, 1)
      if (A(i, c))
        pivot = i;
        break;
      end
    end
    if (pivot == 0)
      continue;
    end
    r = r + 1;
    A([r pivot], :) = A([pivot r], :);
    for i = 1:size (A, 1)
      if (i ~= r && A(i, c))
        A(i, :) = mod (A(i, :) + A(r, :), 2);
      end
    end
  end
end

function G = systematic (P)
  G = [eye(size (P, 1)) P];
end

function H = check_matrix (P)
  H = [P' eye(size (P, 2))];
end

function A = weights (code)
  % The weight distribution A_0, ..., A_n of a cell array of bit strings.
  n = numel (code{1});
  A = zeros (1, n + 1);
  for i = 1:numel (code)
    w = sum (code{i} == '1');
    A(w + 1) = A(w + 1) + 1;
  end
end

function w = min_weight (code)
  w = Inf;
  for i = 1:numel (code)
    t = sum (code{i} == '1');
    if (t > 0 && t < w)
      w = t;
    end
  end
end

function d = min_distance (code)
  d = Inf;
  for i = 1:numel (code)
    for j = (i + 1):numel (code)
      d = min (d, sum (code{i} ~= code{j}));
    end
  end
end

function count = zero_column_sets (H, w)
  % How many sets of w columns of H add to the zero column.
  sets = nchoosek (1:size (H, 2), w);
  count = 0;
  for i = 1:size (sets, 1)
    if (~any (mod (sum (H(:, sets(i, :)), 2), 2)))
      count = count + 1;
    end
  end
end

function w = least_zero_column_set (H)
  w = 0;
  for t = 1:size (H, 2)
    if (zero_column_sets (H, t))
      w = t;
      return;
    end
  end
end

function s = support (v)
  s = join_numbers (find (v));
end

function s = rate6 (k, n)
  s = sprintf ('%.6f', k / n);
end

function g = tree_gates (P)
  % Two-input XOR gates for each parity bit: one fewer than its column weight.
  g = max (sum (P, 1) - 1, 0);
end

function s = add_words (a, b)
  s = word (mod ((a - '0') + (b - '0'), 2));
end

function tf = in_code (x, code)
  tf = any (strcmp (x, code));
end

function tf = closed_under_xor (code)
  tf = true;
  for i = 1:numel (code)
    for j = 1:numel (code)
      tf = tf && in_code (add_words (code{i}, code{j}), code);
    end
  end
end

function s = column_sum (H, js)
  s = word (mod (sum (H(:, js), 2), 2)');
end


% The codes of the chapter ------------------------------------------------------------

C5 = {'00000', '01101', '10110', '11011'};
G5 = parse_rows ('10110 01101');
P5 = parse_rows ('110 101');
H5 = parse_rows ('11100 10010 01001');

P = parse_rows ('110 011 101 111');
G = systematic (P);
H = check_matrix (P);
HAMMING = span_of (G);

P16 = parse_rows ('110 011 100 111');
G16 = systematic (P16);
H16 = check_matrix (P16);
CODE16 = span_of (G16);

SIGMA = [5 6 1 7 3 2 4];

printf ('Chapter 17, Linear Codes: every number in the chapter, recomputed\n');

section ('A code closed under XOR');
same ('C5 contains the zero word', in_code ('00000', C5), true);
sums = {};
for i = 1:4
  for j = 1:4
    sums{end + 1} = add_words (C5{i}, C5{j});
  end
end
same ('ordered pairs of C5 words checked', numel (sums), 16);
same ('every XOR of two C5 words lies in C5', all (cellfun (@(s) in_code (s, C5), sums)), true);
W3 = rows_to_code (all_words (3));
wt3 = cellfun (@(s) sum (s == '1'), W3);
even3 = W3(mod (wt3, 2) == 0);
odd3 = W3(mod (wt3, 2) == 1);
same ('even-weight words of length 3', join_words (even3), '000 011 101 110');
same ('even-weight words closed under XOR', closed_under_xor (even3), true);
same ('odd-weight words of length 3', join_words (odd3), '001 010 100 111');
same ('001 xor 010', add_words ('001', '010'), '011');
same ('odd-weight words closed under XOR', in_code ('011', odd3), false);
sizes = [];
all_zero = true;
subsets = 0;
for mask = 1:255
  S = find (bitand (mask, 2 .^ (0:7))) - 1;
  subsets = subsets + 1;
  ok = true;
  for a = S
    for b = S
      ok = ok && any (S == bitxor (a, b));
    end
  end
  if (ok)
    sizes(end + 1) = numel (S);
    all_zero = all_zero && any (S == 0);
  end
end
same ('nonempty subsets of the 8 words of length 3', subsets, 255);
same ('of those, closed under XOR', numel (sizes), 16);
same ('every closed subset contains 000', all_zero, true);
same ('sizes of the closed subsets', join_numbers (sort (sizes)), '1 2 2 2 2 2 2 2 4 4 4 4 4 4 4 8');

section ('The generator matrix');
same ('G5', show (G5), '10110 01101');
msgs = {'00', '01', '10', '11'};
cws = {'00000', '01101', '10110', '11011'};
for i = 1:4
  same (sprintf ('message %s encodes to', msgs{i}), word (encode (msgs{i} - '0', G5)), cws{i});
end
same ('the four codewords of G5 are C5', isequal (span_of (G5), C5), true);
ok = true;
W2 = all_words (2);
for i = 1:4
  m = W2(i, :);
  ok = ok && isequal (encode (m, G5), [m(1) m(2) mod(m(1) + m(2), 2) m(1) m(2)]);
end
same ('c = (m1, m2, m1+m2, m1, m2) for every message', ok, true);
same ('rank of G5', rank2 (G5), 2);
same ('codewords 2^k for k = 2', numel (span_of (G5)), 4);
same ('rate of C5', rate6 (2, 5), '0.400000');
same ('transpose of G5', show (G5'), '10 01 11 10 01');

section ('The parity-check matrix');
same ('H5', show (H5), '11100 10010 01001');
same ('H5 c^T = 000 for every codeword of C5', all (cellfun (@(c) strcmp (syndrome (H5, c - '0'), '000'), C5)), true);
same ('G5 H5^T', show (mod (G5 * H5', 2)), '000 000');
same ('rank of H5', rank2 (H5), 3);
same ('words x of length 5 with H5 x^T = 0', numel (null_space (H5)), 4);
same ('null space of H5 equals C5', isequal (null_space (H5), C5), true);
same ('syndrome of 10111', syndrome (H5, '10111' - '0'), '001');
same ('syndrome of 10110', syndrome (H5, '10110' - '0'), '000');

section ('From the generator to the checks');
same ('P5', show (P5), '110 101');
same ('G5 = [I2 | P5]', isequal (systematic (P5), G5), true);
same ('P5^T', show (P5'), '11 10 01');
same ('[P5^T | I3] = H5', isequal (check_matrix (P5), H5), true);
same ('P of the (7,4) Hamming code', show (P), '110 011 101 111');
same ('G = [I4 | P]', show (G), '1000110 0100011 0010101 0001111');
same ('H = [P^T | I3]', show (H), '1011100 1101010 0111001');
same ('G H^T', show (mod (G * H', 2)), '000 000 000 000');
ok_zero = true;
ok_null = true;
for bits = 0:4095
  Q = reshape (bitand (bitshift (bits, -(0:11)), 1), 3, 4)';
  GQ = systematic (Q);
  HQ = check_matrix (Q);
  ok_zero = ok_zero && ~any (any (mod (GQ * HQ', 2)));
  ok_null = ok_null && isequal (span_of (GQ), null_space (HQ));
end
require ('G H^T = 0 for all 4096 blocks P of 4 rows, 3 columns', ok_zero);
require ('row space of G = null space of H for all 4096 blocks P', ok_null);
same ('rank of H', rank2 (H), 3);
same ('null space of H has 2^4 words', numel (null_space (H)), 16);
same ('null space of H equals the row space of G', isequal (null_space (H), HAMMING), true);

section ('Systematic encoding with the (7,4) code');
feeds = {'1 3 4', '1 2 4', '2 3 4'};
for j = 1:3
  same (sprintf ('message bits feeding p%d', j), support (P(:, j)'), feeds{j});
end
same ('message 1011 encodes to', word (encode ('1011' - '0', G)), '1011100');
same ('message 0110 encodes to', word (encode ('0110' - '0', G)), '0110110');
same ('message 1110 encodes to', word (encode ('1110' - '0', G)), '1110000');
W4 = all_words (4);
listed = rows_to_code (mod (W4 * G, 2));
ok = true;
for i = 1:16
  ok = ok && strcmp (listed{i}(1:4), word (W4(i, :)));
end
same ('first four bits of every codeword are its message', ok, true);
same ('the 16 codewords in message order', join_words (listed), ...
      ['0000000 0001111 0010101 0011010 0100011 0101100 0110110 0111001 ' ...
       '1000110 1001001 1010011 1011100 1100101 1101010 1110000 1111111']);
same ('their weights', join_numbers (cellfun (@(s) sum (s == '1'), listed)), '0 4 3 3 3 3 4 4 3 3 4 4 4 4 3 7');

section ('The encoder as XOR trees');
same ('column weights of P', join_numbers (sum (P, 1)), '3 3 3');
same ('two-input gates per parity bit', join_numbers (tree_gates (P)), '2 2 2');
same ('gates in the encoder, one tree per parity bit', sum (tree_gates (P)), 6);
same ('gate depth of a three-input tree', ceil (log2 (3)), 2);
same ('gates in the C5 encoder', sum (tree_gates (P5)), 1);
shared_ok = true;
for i = 1:16
  m = W4(i, :);
  t = mod (m(1) + m(4), 2);
  u = mod (m(2) + m(3), 2);
  shared = [mod(t + m(3), 2) mod(t + m(2), 2) mod(u + m(4), 2)];
  c = encode (m, G);
  shared_ok = shared_ok && isequal (shared, c(5:7));
end
require ('five shared gates give p1, p2, p3 for all 16 messages', shared_ok);

section ('Distance is the smallest weight');
same ('weight distribution of the (7,4) code', join_numbers (weights (HAMMING)), '1 0 0 7 7 0 0 1');
same ('nonzero weights checked', numel (HAMMING) - 1, 15);
same ('pairs of codewords', nchoosek (numel (HAMMING), 2), 120);
same ('minimum nonzero weight', min_weight (HAMMING), 3);
same ('minimum distance over all 120 pairs', min_distance (HAMMING), 3);
ok = true;
hw = sort (cellfun (@(s) sum (s == '1'), HAMMING));
for i = 1:16
  dists = sort (cellfun (@(s) sum (s ~= HAMMING{i}), HAMMING));
  ok = ok && isequal (dists, hw);
end
same ('distances from each codeword repeat the weights', ok, true);
same ('minimum weight of C5', min_weight (C5), 3);
same ('minimum distance of C5', min_distance (C5), 3);
same ('G of the code from Matrices and Elimination', show (G16), '1000110 0100011 0010100 0001111');
same ('weight distribution of that code', join_numbers (weights (CODE16)), '1 0 1 6 5 2 1 0');
same ('its minimum weight', min_weight (CODE16), 2);
same ('0010100 is one of its codewords', in_code ('0010100', CODE16), true);
same ('message 1011 encodes to', word (encode ('1011' - '0', G16)), '1011101');

section ('Reading the distance from the columns of H');
same ('columns of H', show (H'), '110 011 101 111 100 010 001');
cols = rows_to_code (H');
same ('every column nonzero', all (any (H, 1)), true);
same ('all seven columns different', numel (unique (cols)), 7);
same ('nonzero words of length 3', 2 ^ 3 - 1, 7);
same ('columns 1 + 5 + 6', column_sum (H, [1 5 6]), '000');
same ('1000110 is a codeword', in_code ('1000110', HAMMING), true);
column_view = true;
W7 = all_words (7);
for i = 1:128
  x = W7(i, :);
  column_view = column_view && strcmp (column_sum (H, find (x)), syndrome (H, x));
end
require ('H x^T = sum of the columns at the 1s of x, all 128 x', column_view);
zs = zeros (1, 7);
for w = 1:7
  zs(w) = zero_column_sets (H, w);
end
same ('sets of w columns adding to zero, w = 1 to 7', join_numbers (zs), '0 0 7 7 0 0 1');
A = weights (HAMMING);
same ('codewords of weight w, w = 1 to 7', join_numbers (A(2:end)), '0 0 7 7 0 0 1');
same ('fewest columns of H adding to zero', least_zero_column_set (H), 3);
same ('columns of H5', show (H5'), '110 101 100 010 001');
same ('columns 1 + 3 + 4 of H5', column_sum (H5, [1 3 4]), '000');
same ('fewest columns of H5 adding to zero', least_zero_column_set (H5), 3);
same ('columns of H for the Matrices and Elimination code', show (H16'), '110 011 100 111 100 010 001');
same ('its columns 3 and 5 are equal', isequal (H16(:, 3), H16(:, 5)), true);
same ('fewest of its columns adding to zero', least_zero_column_set (H16), 2);

section ('The Venn diagram of the (7,4) code');
circles = {'1 3 4 5', '1 2 4 6', '2 3 4 7'};
for i = 1:3
  same (sprintf ('positions inside circle %d', i), support (H(i, :)), circles{i});
end
held = cell (1, 7);
for j = 1:7
  held{j} = sprintf ('%d', find (H(:, j))');
end
same ('circles holding positions 1 to 7', strjoin (held, ' '), '12 23 13 123 1 2 3');
same ('circle parities of 1011100', syndrome (H, '1011100' - '0'), '000');
same ('every codeword has three even circles', all (cellfun (@(c) strcmp (syndrome (H, c - '0'), '000'), HAMMING)), true);
same ('words of length 7 with three even circles', numel (null_space (H)), 16);
same ('every position lies in some circle', all (any (H, 1)), true);

section ('The dual code');
dual5 = span_of (H5);
same ('C5 dual, the row space of H5', join_words (dual5), '00000 00111 01001 01110 10010 10101 11011 11100');
W5 = all_words (5);
C5M = cell2mat (cellfun (@(c) c - '0', C5', 'UniformOutput', false));
orth = rows_to_code (W5(all (mod (W5 * C5M', 2) == 0, 2), :));
same ('C5 dual = words orthogonal to all of C5', isequal (dual5, orth), true);
same ('dimension of C5 + dimension of its dual', sprintf ('%d + %d = %d', rank2 (G5), rank2 (H5), rank2 (G5) + rank2 (H5)), '2 + 3 = 5');
same ('weight distribution of C5 dual', join_numbers (weights (dual5)), '1 0 2 4 1 0');
simplex = span_of (H);
same ('dual of the (7,4) code, the row space of H', join_words (simplex), ...
      '0000000 0001111 0110110 0111001 1010011 1011100 1100101 1101010');
same ('its weight distribution', join_numbers (weights (simplex)), '1 0 0 0 7 0 0 0');
same ('H H^T', show (mod (H * H', 2)), '000 000 000');
same ('H G^T', show (mod (H * G', 2)), '0000 0000 0000');
same ('the dual lies inside the (7,4) code', all (cellfun (@(y) in_code (y, HAMMING), simplex)), true);
same ('null space of G equals the dual', isequal (null_space (G), simplex), true);
same ('dimension of the code + dimension of its dual', sprintf ('%d + %d = %d', rank2 (G), rank2 (H), rank2 (G) + rank2 (H)), '4 + 3 = 7');
same ('rate of the dual', rate6 (3, 7), '0.428571');
half = true;
for m = 1:6
  Wm = all_words (m);
  for i = 2:size (Wm, 1)
    half = half && sum (mod (Wm * Wm(i, :)', 2)) == 2 ^ (m - 1);
  end
end
W3b = all_words (3);
counts = zeros (1, 7);
for i = 2:8
  counts(i - 1) = sum (mod (W3b * W3b(i, :)', 2));
end
same ('h in GF(2)^3 with y . h = 1, each nonzero y', join_numbers (counts), '4 4 4 4 4 4 4');
require ('exactly half of GF(2)^m meets y . h = 1, m = 1 to 6', half);

section ('Equivalent codes');
same ('the permutation sigma', join_numbers (SIGMA), '5 6 1 7 3 2 4');
Hp = H(:, SIGMA);
same ('H'' with columns sigma(1), ..., sigma(7) of H', show (Hp), '1010101 0110011 0001111');
same ('column j of H'' read as binary, row 1 lowest', join_numbers ([1 2 4] * Hp), '1 2 3 4 5 6 7');
x = '1011100';
same ('1011100 permuted', x(SIGMA), '1010101');
same ('1010101 is in the (7,4) code', in_code ('1010101', HAMMING), false);
same ('message 1010 encodes to', word (encode ('1010' - '0', G)), '1010011');
Cp = sort (cellfun (@(c) c(SIGMA), HAMMING, 'UniformOutput', false));
same ('permuted code is the null space of H''', isequal (Cp, null_space (Hp)), true);
same ('weight distribution of the permuted code', join_numbers (weights (Cp)), '1 0 0 7 7 0 0 1');
same ('minimum distance of the permuted code', min_distance (Cp), 3);
CM = char (HAMMING);
orders = perms (1:7);
keys = cell (1, size (orders, 1));
fixed = 0;
for i = 1:size (orders, 1)
  img = sortrows (CM(:, orders(i, :)));
  keys{i} = reshape (img', 1, []);
  if (isequal (img, CM))
    fixed = fixed + 1;
  end
end
orderings = size (orders, 1);
same ('orderings of 7 positions', orderings, 5040);
same ('orderings that map the (7,4) code to itself', fixed, 168);
same ('different codes equivalent to it', numel (unique (keys)), 30);
same ('5040 / 168', floor (orderings / fixed), numel (unique (keys)));

section ('How many codewords, and at what rate');
table = {
  'repetition', '111', '3 1 3', 2, '0.333333';
  'even weight, length 3', '101 011', '3 2 2', 4, '0.666667';
  'C5', '10110 01101', '5 2 3', 4, '0.400000';
  '(6,3) code', '100011 010101 001110', '6 3 3', 8, '0.500000';
  '(7,4) Hamming', show(G), '7 4 3', 16, '0.571429';
  'its dual', show(H), '7 3 4', 8, '0.428571';
  'Matrices and Elimination code', show(G16), '7 4 2', 16, '0.571429'};
for i = 1:size (table, 1)
  Gt = parse_rows (table{i, 2});
  code = span_of (Gt);
  n = size (Gt, 2);
  k = rank2 (Gt);
  same (sprintf ('%s: n k d', table{i, 1}), sprintf ('%d %d %d', n, k, min_weight (code)), table{i, 3});
  same (sprintf ('%s: codewords', table{i, 1}), numel (code), table{i, 4});
  same (sprintf ('%s: rate', table{i, 1}), rate6 (k, n), table{i, 5});
end
same ('fraction of 7-bit words in the (7,4) code', sprintf ('1/%d', 2 ^ 7 / numel (HAMMING)), '1/8');
same ('(15,11) code: codewords', 2 ^ 11, 2048);
same ('(15,11) code: rate', rate6 (11, 15), '0.733333');
same ('(15,11) code: fraction of words that are codewords', sprintf ('1/%d', 2 ^ 15 / 2 ^ 11), '1/16');

section ('Worked examples');
P6 = parse_rows ('011 101 110');
G6 = systematic (P6);
H6 = check_matrix (P6);
code6 = span_of (G6);
same ('(6,3): G', show (G6), '100011 010101 001110');
same ('(6,3): H', show (H6), '011100 101010 110001');
same ('(6,3): columns of H', show (H6'), '011 101 110 100 010 001');
same ('(6,3): columns distinct and nonzero', numel (unique (rows_to_code (H6'))) == 6 && all (any (H6, 1)), true);
same ('(6,3): columns 1 + 5 + 6', column_sum (H6, [1 5 6]), '000');
same ('(6,3): 100011 is a codeword', in_code ('100011', code6), true);
same ('(6,3): weight distribution', join_numbers (weights (code6)), '1 0 0 4 3 0 0');
same ('(6,3): the eight codewords', join_words (code6), '000000 001110 010101 011011 100011 101101 110110 111000');
same ('(6,3): message 101 encodes to', word (encode ('101' - '0', G6)), '101101');
same ('(6,3): gates in the encoder', sum (tree_gates (P6)), 3);
G4 = systematic (parse_rows ('1 1 1'));
H4 = check_matrix (parse_rows ('1 1 1'));
same ('even weight, length 4: G', show (G4), '1001 0101 0011');
same ('even weight, length 4: H', show (H4), '1111');
same ('even weight, length 4: codewords', join_words (span_of (G4)), '0000 0011 0101 0110 1001 1010 1100 1111');
same ('even weight, length 4: minimum distance', min_weight (span_of (G4)), 2);
same ('even weight, length 4: dual', join_words (span_of (H4)), '0000 1111');
same ('even weight, length 4: rate', rate6 (3, 4), '0.750000');
Gn = parse_rows ('11011 01101');
same ('non-systematic G spans C5', isequal (span_of (Gn), C5), true);
same ('row 1 + row 2 of it', word (mod (Gn(1, :) + Gn(2, :), 2)), '10110');
c = word (encode ('0110' - '0', G));
same ('0110: codeword', c, '0110110');
same ('0110: weight of the codeword', sum (c == '1'), 4);
same ('0110: H c^T', syndrome (H, '0110110' - '0'), '000');
m = '0110' - '0';
t = mod (m(1) + m(4), 2);
u = mod (m(2) + m(3), 2);
same ('0110: t, p1, p2, u, p3 from the shared gates', join_numbers ([t mod(t + m(3), 2) mod(t + m(2), 2) u mod(u + m(4), 2)]), '0 1 1 0 0');

section ('Exercises');
S = {'0000', '1100', '0011', '1111', '1010'};
same ('1100 xor 1010', add_words ('1100', '1010'), '0110');
same ('0110 is in the five-word set', in_code ('0110', S), false);
P53 = parse_rows ('11 10 01');
G53 = systematic (P53);
H53 = check_matrix (P53);
same ('(5,3): G', show (G53), '10011 01010 00101');
same ('(5,3): H', show (H53), '11010 10101');
same ('(5,3): columns of H', show (H53'), '11 10 01 10 01');
same ('(5,3): columns 2 and 4 equal', isequal (H53(:, 2), H53(:, 4)), true);
same ('(5,3): 01010 is a codeword', in_code ('01010', span_of (G53)), true);
same ('(5,3): minimum distance', min_weight (span_of (G53)), 2);
same ('1110: H c^T', syndrome (H, '1110000' - '0'), '000');
HH = H * H';
same ('dot products of the rows of H, 11 12 13 22 23 33', join_numbers ([HH(1, 1) HH(1, 2) HH(1, 3) HH(2, 2) HH(2, 3) HH(3, 3)]), '4 2 2 4 2 4');
hwt = cellfun (@(s) sum (s == '1'), HAMMING);
even = HAMMING(mod (hwt, 2) == 0);
same ('even-weight codewords of the (7,4) code', numel (even), 8);
same ('they are closed under XOR', closed_under_xor (even), true);
same ('they are exactly the dual code', isequal (sort (even), span_of (H)), true);

printf ('\nAll chapter 17 numbers reproduced.\n');
