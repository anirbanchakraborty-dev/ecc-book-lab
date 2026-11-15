% Reproduce every number in chapter 20 of Error Correction Coding,
% "SECDED and Hsiao Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch20_secded_and_hsiao_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch20_secded_and_hsiao_codes.py and prints the
% same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in the chapter is a whole number, a bit pattern or
% a percentage printed to one decimal place, so every comparison is exact.
% The distance of every code is checked by testing every set of one, two and
% three columns of its check matrix, and the triple-error counts come from
% all 59,640 three-bit errors of the 72-bit codes, not from a sample.
%
% A column of a check matrix with r rows is held as a whole number whose bit
% i is the entry in row i, and a code is a row vector of such columns, data
% columns first and check columns last. A bit pattern is printed with row 0
% first. Positions are counted from 0, as in the chapter, so position p of a
% code is element p + 1 here. Octave's nchoosek lists combinations in the
% same order as Python's itertools.combinations, which keeps the columns of
% every code in the order the Python script builds them.

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

function w = weight (v)
  % The number of ones in each entry of v.
  w = zeros (size (v));
  while (any (v(:)))
    w = w + bitand (v, 1);
    v = bitshift (v, -1);
  end
end

function v = col (rows)
  % The column with ones in the given rows.
  v = sum (2 .^ rows);
end

function v = cols_of (C)
  % One column for each row of C, a list of row sets.
  v = sum (2 .^ C, 2)';
end

function s = bits (v, r)
  % A column or syndrome as text, row 0 first.
  s = char ('0' + bitand (bitshift (v, -(0:r - 1)), 1));
end

function s = bits_list (cols, r)
  parts = cell (1, numel (cols));
  for i = 1:numel (cols)
    parts{i} = bits (cols(i), r);
  end
  s = strjoin (parts, ' ');
end

function s = joined (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ',');
end

function loads = row_loads (cols, r)
  loads = zeros (1, r);
  for i = 0:r - 1
    loads(i + 1) = sum (bitand (bitshift (cols, -i), 1));
  end
end

function d = levels (w, fan_in)
  % Levels of fan_in-input gates in a balanced tree over w signals.
  if (nargin < 2)
    fan_in = 2;
  end
  d = 0;
  reach = 1;
  while (reach < w)
    reach = reach * fan_in;
    d = d + 1;
  end
end

function q = ceil_div (a, b)
  q = ceil (a / b);
end

function s = pct (part, whole)
  s = sprintf ('%.1f', 100 * part / whole);
end

function r = sec_r (k)
  % The least r with 2^r >= k + r + 1.
  r = 1;
  while (2 ^ r < k + r + 1)
    r = r + 1;
  end
end

function r = secded_r (k)
  % The least r with 2^(r-1) >= k + r.
  r = 1;
  while (2 ^ (r - 1) < k + r)
    r = r + 1;
  end
end


% Codes ------------------------------------------------------------------------

function h = hsiao_8_4 ()
  h = [cols_of(nchoosek (0:3, 3)), 2 .^ (0:3)];
end

function h = hsiao_13_8 ()
  data = cols_of (nchoosek (0:4, 3));
  data = data(~ismember (data, [col([0 1 2]), col([0 3 4])]));
  h = [data, 2 .^ (0:4)];
end

function h = hsiao_22_16 ()
  data = cols_of (nchoosek (0:5, 3));
  data = data(~ismember (data, [col([0 1 2]), col([0 3 4]), col([1 3 5]), col([2 4 5])]));
  h = [data, 2 .^ (0:5)];
end

function h = hsiao_39_32 ()
  data = cols_of (nchoosek (0:6, 3));
  data = data(~ismember (data, [col([0 1 2]), col([0 3 4]), col([1 5 6])]));
  h = [data, 2 .^ (0:6)];
end

function h = hsiao_72_64 ()
  % All 56 weight-3 columns, then the 8 cyclic runs of five rows.
  runs = zeros (1, 8);
  for i = 0:7
    runs(i + 1) = col (mod (i + (0:4), 8));
  end
  h = [cols_of(nchoosek (0:7, 3)), runs, 2 .^ (0:7)];
end

function h = hsiao_72_64_dictionary ()
  % All 56 weight-3 columns, then the first 8 weight-5 columns in dictionary order.
  C5 = nchoosek (0:7, 5);
  h = [cols_of(nchoosek (0:7, 3)), cols_of(C5(1:8, :)), 2 .^ (0:7)];
end

function h = hsiao_137_128 ()
  % All 84 weight-3 columns of length 9, then 44 weight-5 columns chosen
  % greedily. Each step takes the unused weight-5 word whose rows hold the
  % fewest ones so far, the first in dictionary order on a tie.
  data = cols_of (nchoosek (0:8, 3));
  C5 = nchoosek (0:8, 5);
  loads = (1 + nchoosek (8, 2)) * ones (1, 9);
  used = false (rows (C5), 1);
  order = zeros (1, 128 - numel (data));
  for step = 1:numel (order)
    score = sum (loads(C5 + 1), 2);
    score(used) = Inf;
    [~, i] = min (score);
    used(i) = true;
    order(step) = i;
    loads(C5(i, :) + 1) = loads(C5(i, :) + 1) + 1;
  end
  h = [data, cols_of(C5(order, :)), 2 .^ (0:8)];
end

function h = extended_hamming (n, rh)
  % Positions 0 to n-1 of the extended Hamming code with rh Hamming rows.
  % Rows 0 to rh-1 hold the binary form of the position number, and row rh
  % is the overall parity row of all ones. Position 0 is the parity bit.
  h = bitor (0:n - 1, 2 ^ rh);
end


% Properties -------------------------------------------------------------------

function ok = distance_at_least_4 (cols)
  ok = false;
  if (any (cols == 0) || numel (unique (cols)) ~= numel (cols))
    return;
  end
  P = nchoosek (1:numel (cols), 2);
  a = cols(P(:, 1));
  b = cols(P(:, 2));
  x = bitxor (a, b);
  ok = ~any (ismember (x, cols) & x ~= a & x ~= b);
end

function ok = has_weight_4_codeword (cols)
  % True when two disjoint pairs of columns have the same sum.
  P = nchoosek (1:numel (cols), 2);
  x = bitxor (cols(P(:, 1)), cols(P(:, 2)));
  [xs, order] = sort (x(:));
  P = P(order, :);
  ok = false;
  n = numel (xs);
  i = 1;
  while (i <= n)
    j = i;
    while (j < n && xs(j + 1) == xs(i))
      j = j + 1;
    end
    for a = i:j - 1
      for b = a + 1:j
        if (numel (unique ([P(a, :), P(b, :)])) == 4)
          ok = true;
          return;
        end
      end
    end
    i = j + 1;
  end
end

function s = syndrome (cols, positions)
  s = 0;
  for p = positions
    s = bitxor (s, cols(p + 1));
  end
end

function [x, m] = triple_sums (cols)
  % The syndrome of every three-bit error, and whether it equals a column.
  T = nchoosek (1:numel (cols), 3);
  x = bitxor (bitxor (cols(T(:, 1)), cols(T(:, 2))), cols(T(:, 3)));
  m = ismember (x, cols);
end

function [mis, det] = triple_outcomes (cols)
  [~, m] = triple_sums (cols);
  mis = sum (m);
  det = numel (m) - mis;
end

function table = triple_by_weight (cols, r)
  % Row w + 1 holds the (miscorrected, detected) counts at syndrome weight w.
  [x, m] = triple_sums (cols);
  w = weight (x);
  table = zeros (r + 1, 2);
  for t = 0:r
    table(t + 1, :) = [sum(m & w == t), sum(~m & w == t)];
  end
end

function c = check_bits (h, data, r)
  k = numel (data);
  c = zeros (1, r);
  for i = 0:r - 1
    c(i + 1) = mod (sum (data(bitand (bitshift (h(1:k), -i), 1) == 1)), 2);
  end
end

function g = decoder_gates (cols, r, k, parity_row)
  % Two-input gate counts of the full decoder. parity_row is -1 for a Hsiao
  % code, whose single-error signal is the parity of the whole syndrome, and
  % the index of the all-ones row for an extended Hamming code, whose
  % single-error signal is that row itself. Besides correcting and raising
  % the double-error flag, the decoder matches the syndrome against the
  % n - k check columns too, ORs all n match signals into M, and raises the
  % no-match flag for an odd syndrome with M = 0.
  loads = row_loads (cols, r);
  n = numel (cols);
  lv = zeros (1, r);
  for i = 1:r
    lv(i) = levels (loads(i));
  end
  g.syndrome = sum (loads - 1);
  g.syndrome_depth = max (lv);
  if (parity_row < 0)
    zero_inputs = r;
    g.parity = r - 1;
    parity_depth = g.syndrome_depth + levels (r);
    zero_depth = g.syndrome_depth + levels (r);
  else
    zero_inputs = r - 1;
    g.parity = 0;
    parity_depth = levels (loads(parity_row + 1));
    others = lv;
    others(parity_row + 1) = [];
    zero_depth = max (others) + levels (r - 1);
  end
  g.zero = zero_inputs - 1;
  g.match = k * (r - 1);
  g.correct = k;
  g.ded = 1;
  g.check_match = (n - k) * (r - 1);
  g.any_match = n - 1;
  g.no_match = 1;
  g.total = g.syndrome + g.zero + g.parity + g.match + g.correct + g.ded ...
            + g.check_match + g.any_match + g.no_match;
  g.inverters = r + 2;
  match_depth = g.syndrome_depth + levels (r);
  g.path = match_depth + 1;
  g.ded_depth = max (parity_depth, zero_depth) + 1;
  g.m_depth = match_depth + levels (n);
  g.no_match_depth = max (g.m_depth, parity_depth) + 1;
end

function t = decoder_actions (cols, r, parity_row, weight_of_error)
  % Run the full decoder's logic on every error of the given weight. The
  % fields count the errors that end with no flag and no flip, a flip of the
  % one bit in error, a flip of some other bit, the double-error flag, the
  % no-match flag, and more than one of these at once.
  E = nchoosek (1:numel (cols), weight_of_error);
  s = zeros (rows (E), 1);
  for c = 1:weight_of_error
    s = bitxor (s, cols(E(:, c))(:));
  end
  if (parity_row < 0)
    nonzero = s ~= 0;
    odd = mod (weight (s), 2) == 1;
  else
    nonzero = bitand (s, (2 ^ r - 1) - 2 ^ parity_row) ~= 0;
    odd = bitand (bitshift (s, -parity_row), 1) == 1;
  end
  ded = nonzero & ~odd;
  [hit, where] = ismember (s, cols);
  no_match = odd & ~hit;
  conflict = (ded + no_match + hit) > 1;
  own = false (size (hit));
  if (weight_of_error == 1)
    own = hit & where == E(:, 1);
  end
  ok = ~conflict;
  t.conflict = sum (conflict);
  t.double = sum (ok & ded);
  t.no_match = sum (ok & ~ded & no_match);
  t.silent = sum (ok & ~ded & ~no_match & ~hit);
  t.corrected = sum (ok & ~ded & ~no_match & own);
  t.miscorrected = sum (ok & ~ded & ~no_match & hit & ~own);
end

function ok = distance_is_4 (h)
  ok = distance_at_least_4 (h) && has_weight_4_codeword (h);
end


% Checks -----------------------------------------------------------------------

function check_bit_budget ()
  section ('How many check bits a 64-bit word needs');
  k = 64;
  same ('SEC test 2^r >= k + r + 1 at r = 6', sprintf ('%d >= %d', 2 ^ 6, k + 6 + 1), '64 >= 71');
  same ('SEC test at r = 7', sprintf ('%d >= %d', 2 ^ 7, k + 7 + 1), '128 >= 72');
  same ('least r for single correction, k = 64', sec_r (k), 7);
  same ('SECDED test 2^(r-1) >= k + r at r = 7', sprintf ('%d >= %d', 2 ^ 6, k + 7), '64 >= 71');
  same ('SECDED test at r = 8', sprintf ('%d >= %d', 2 ^ 7, k + 8), '128 >= 72');
  same ('least r for SECDED, k = 64', secded_r (k), 8);
  same ('codeword length 64 + 8', 64 + 8, 72);
  same ('bytes in a 72-bit word', floor (72 / 8), 9);
  same ('check-bit overhead 8/64, percent', pct (8, 64), '12.5');
  for row = [4 4 8; 8 5 13; 32 7 39; 128 9 137]'
    kk = row(1);
    r = secded_r (kk);
    same (sprintf ('least r for SECDED, k = %d', kk), r, row(2));
    same (sprintf ('length for k = %d', kk), kk + r, row(3));
  end
  same ('least r for single correction, k = 32', sec_r (32), 6);
  same ('SECDED test for k = 32 at r = 6', sprintf ('%d >= %d', 2 ^ 5, 32 + 6), '32 >= 38');
  same ('check-bit overhead 7/32, percent', pct (7, 32), '21.9');
  same ('longest SECDED length with r = 8', 2 ^ 7, 128);
  same ('most data bits with r = 8', 2 ^ 7 - 8, 120);
  same ('longest SECDED length with r = 7', 2 ^ 6, 64);
  same ('most data bits with r = 7', 2 ^ 6 - 7, 57);
  same ('SECDED test for k = 128 at r = 8', sprintf ('%d >= %d', 2 ^ 7, 128 + 8), '128 >= 136');
  same ('check-bit overhead 9/128, percent', pct (9, 128), '7.0');
  same ('least r for single correction, k = 16', sec_r (16), 5);
  same ('least r for SECDED, k = 16', secded_r (16), 6);
  same ('longest SECDED length with r = 5', 2 ^ 4, 16);
  same ('most data bits with r = 5', 2 ^ 4 - 5, 11);
  table = {4, 3, 4, 8, '100.0'; 8, 4, 5, 13, '62.5'; 16, 5, 6, 22, '37.5'; ...
           32, 6, 7, 39, '21.9'; 64, 7, 8, 72, '12.5'; 128, 8, 9, 137, '7.0'};
  for t = 1:rows (table)
    [kk, r1, r2, nn, over] = table{t, :};
    a = sec_r (kk);
    b = secded_r (kk);
    same (sprintf ('table row k = %d: SEC r, SECDED r, n, overhead', kk), ...
          sprintf ('%d,%d,%d,%s', a, b, kk + b, pct (b, kk)), sprintf ('%d,%d,%d,%s', r1, r2, nn, over));
  end
  ok = true;
  for kk = 1:4096
    ok = ok && secded_r (kk) == sec_r (kk) + 1;
  end
  require ('SECDED needs exactly one check bit more than SEC, k = 1 to 4096', ok);
end

function check_extended_hamming ()
  section ('SECDED from the extended Hamming code');
  full = extended_hamming (128, 7);
  require ('(128,120) extended Hamming check matrix has distance 4', distance_is_4 (full));
  h = extended_hamming (72, 7);
  same ('positions kept by the shortening', numel (h), 72);
  same ('positions deleted from the (128,120) code', 128 - 72, 56);
  j = 1:71;
  checks = j(weight (j) == 1);
  same ('Hamming check positions among 1..71', joined (checks), '1,2,4,8,16,32,64');
  same ('data positions among 1..71', 71 - numel (checks), 64);
  require ('shortened (72,64) check matrix has distance 4', distance_is_4 (h));
  loads = row_loads (h, 8);
  same ('ones in each row, rows 0 to 7', joined (loads), '36,36,36,32,32,32,8,72');
  same ('ones in the extended Hamming (72,64) check matrix', sum (loads), 284);
  % the decision table on (s', p)
  s = syndrome (h, 37);
  same ('syndrome of an error at position 37, rows 0 to 6', bits (bitand (s, 127), 7), '1010010');
  same ('overall parity of that error', bitshift (s, -7), 1);
  s = syndrome (h, [5 9]);
  same ('syndrome of errors at positions 5 and 9, rows 0 to 6', bits (bitand (s, 127), 7), '0011000');
  same ('overall parity of that pair', bitshift (s, -7), 0);
  s = syndrome (h, 0);
  same ('syndrome of an error in the parity bit, rows 0 to 6', bits (bitand (s, 127), 7), '0000000');
  % every single and double error falls in its row of the table
  ok = all (bitshift (h, -7) == 1);
  P = nchoosek (1:72, 2);
  s = bitxor (h(P(:, 1)), h(P(:, 2)));
  ok = ok && all (bitshift (s, -7) == 0 & bitand (s, 127) ~= 0);
  require ('every single error has p = 1, every double p = 0 and s'' nonzero', ok);
  h39 = extended_hamming (39, 6);
  require ('extended Hamming (39,32) has distance 4', distance_is_4 (h39));
  loads = row_loads (h39, 7);
  same ('ones in each row of extended Hamming (39,32)', joined (loads), '19,19,19,16,16,7,39');
  same ('ones in the extended Hamming (39,32) check matrix', sum (loads), 135);
  % replace the all-ones row by the sum of all eight rows: the same code, odd columns
  low = bitand (h, 127);
  odd_form = bitor (low, bitshift (mod (1 + weight (low), 2), 7));
  require ('odd-column form: every column has odd weight', all (mod (weight (odd_form), 2) == 1));
  require ('odd-column form still has distance 4', distance_is_4 (odd_form));
  same ('positions 0..71 whose index has even weight', sum (mod (weight (0:71), 2) == 0), 36);
  same ('ones in the odd-column form, 212 + 36', sum (row_loads (odd_form, 8)), 248);
end

function check_odd_weight ()
  section ('Odd-weight columns');
  for row = [4 8; 5 16; 6 32; 7 64; 8 128]'
    r = row(1);
    same (sprintf ('odd-weight words of length %d', r), sum (mod (weight (0:2 ^ r - 1), 2) == 1), row(2));
  end
  w8 = weight (0:255);
  by_weight = [sum(w8 == 1), sum(w8 == 3), sum(w8 == 5), sum(w8 == 7)];
  same ('8-bit odd words of weight 1, 3, 5, 7', joined (by_weight), '8,56,56,8');
  [A, B] = meshgrid (1:63, 1:63);
  ok = all (all (weight (bitxor (A, B)) == weight (A) + weight (B) - 2 * weight (bitand (A, B))));
  require ('wt(a+b) = wt(a) + wt(b) - 2 wt(a AND b), all 6-bit pairs', ok);
end

function check_hsiao_8_4 ()
  section ('The Hsiao (8,4) code');
  h = hsiao_8_4 ();
  same ('columns d0..d3 c0..c3', bits_list (h, 4), '1110 1101 1011 0111 1000 0100 0010 0001');
  require ('every column has odd weight', all (mod (weight (h), 2) == 1));
  require ('distance is exactly 4', distance_is_4 (h));
  same ('ones in each row', joined (row_loads (h, 4)), '4,4,4,4');
  same ('ones in the check matrix', sum (row_loads (h, 4)), 16);
  same ('columns used out of the 8 odd words', numel (unique (h)), 8);
  checks = check_bits (h, [1 0 1 1], 4);
  same ('check bits for data 1011', sprintf ('%d', checks), '0010');
  same ('syndrome of an error in d2', bits (syndrome (h, 2), 4), '1011');
  same ('syndrome of errors in d0 and c1', bits (syndrome (h, [0 5]), 4), '1010');
  P = nchoosek (1:8, 2);
  s = bitxor (h(P(:, 1)), h(P(:, 2)));
  same ('double errors', numel (s), 28);
  u = unique (s);
  same ('distinct double-error syndromes', numel (u), 7);
  require ('every double-error syndrome is even and nonzero', all (u ~= 0 & mod (weight (u), 2) == 0));
  counts = zeros (1, numel (u));
  for i = 1:numel (u)
    counts(i) = sum (s == u(i));
  end
  same ('pairs behind each double-error syndrome', joined (unique (counts)), '4');
  names = {'d0', 'd1', 'd2', 'd3', 'c0', 'c1', 'c2', 'c3'};
  stated = {'1111', 'd0c3 d1c2 d2c1 d3c0'; '1100', 'd0c2 d1c3 d2d3 c0c1'; ...
            '1010', 'd0c1 d1d3 d2c3 c0c2'; '1001', 'd0d3 d1c1 d2c2 c0c3'; ...
            '0110', 'd0c0 d1d2 d3c3 c1c2'; '0101', 'd0d2 d1c0 d3c2 c1c3'; ...
            '0011', 'd0d1 d2c0 d3c1 c2c3'};
  for t = 1:rows (stated)
    syn = stated{t, 1};
    found = {};
    for q = 1:rows (P)
      if (strcmp (bits (s(q), 4), syn))
        found{end + 1} = [names{P(q, 1)}, names{P(q, 2)}];
      end
    end
    same (sprintf ('pairs with syndrome %s', syn), strjoin (found, ' '), stated{t, 2});
  end
  same ('syndromes: 1 zero + 8 single + 7 double', 1 + 8 + numel (u), 16);
end

function check_hsiao_13_8 ()
  section ('The Hsiao (13,8) code');
  h = hsiao_13_8 ();
  same ('data columns d0..d7', bits_list (h(1:8), 5), '11010 11001 10110 10101 01110 01101 01011 00111');
  require ('every column has odd weight', all (mod (weight (h), 2) == 1));
  require ('distance is exactly 4', distance_is_4 (h));
  loads = row_loads (h, 5);
  same ('ones in each row', joined (loads), '5,6,6,6,6');
  same ('ones in the check matrix', sum (loads), 29);
  same ('data bits in each check', joined (loads - 1), '4,5,5,5,5');
  v = 0:31;
  odd_unused = v(mod (weight (v), 2) == 1 & ~ismember (v, h));
  same ('odd words that are not columns', bits_list (odd_unused, 5), '11100 10011 11111');
  data = [1 0 1 1 0 0 1 0];
  checks = check_bits (h, data, 5);
  same ('check bits for data 10110010', sprintf ('%d', checks), '10010');
  same ('weight of that codeword', sum (data) + sum (checks), 6);
  same ('syndrome of an error in d5', bits (syndrome (h, 5), 5), '01101');
  same ('syndrome of errors in d0 and c2', bits (syndrome (h, [0 10]), 5), '11110');
  s = syndrome (h, [0 1 2]);
  same ('syndrome of errors in d0, d1, d2', bits (s, 5), '10101');
  same ('column it matches', find (h == s, 1) - 1, 3);
  same ('syndrome of errors in d0, d7, c3', bits (syndrome (h, [0 7 11]), 5), '11111');
  same ('syndrome of errors in d0, c2, c3', bits (syndrome (h, [0 10 11]), 5), '11100');
  same ('syndrome of d0, c0, c1, c3', bits (syndrome (h, [0 8 9 11]), 5), '00000');
  P = nchoosek (1:13, 2);
  doubles = unique (bitxor (h(P(:, 1)), h(P(:, 2))));
  same ('double errors in 13 bits, C(13,2)', nchoosek (13, 2), 78);
  same ('distinct double-error syndromes', numel (doubles), 15);
  require ('every double-error syndrome is even and nonzero', all (doubles ~= 0 & mod (weight (doubles), 2) == 0));
  same ('syndromes: 1 zero + 13 single + 3 unused odd + 15 even', 1 + 13 + 3 + numel (doubles), 32);
end

function check_hsiao_72_64 ()
  section ('The Hsiao (72,64) and (39,32) codes');
  h = hsiao_72_64 ();
  same ('columns', numel (h), 72);
  require ('all 72 columns distinct and odd', numel (unique (h)) == 72 && all (mod (weight (h), 2) == 1));
  require ('distance is exactly 4', distance_is_4 (h));
  same ('data columns of weight 3 and 5', sprintf ('%d,%d', sum (weight (h) == 3), sum (weight (h) == 5)), '56,8');
  same ('weight-5 columns', bits_list (h(57:64), 8), ...
        '11111000 01111100 00111110 00011111 10001111 11000111 11100011 11110001');
  loads = row_loads (h, 8);
  same ('ones in each row', joined (loads), '27,27,27,27,27,27,27,27');
  same ('ones in the check matrix, 8 + 3*56 + 5*8', sum (loads), 216);
  same ('8 + 3*56 + 5*8', 8 + 3 * 56 + 5 * 8, 216);
  same ('ones in the data part of each row', joined (loads - 1), '26,26,26,26,26,26,26,26');
  same ('weight-3 columns that contain one given row, C(7,2)', nchoosek (7, 2), 21);
  same ('cyclic runs of five rows that contain row 0', sum (bitand (h(57:64), 1)), 5);
  same ('least possible largest row, ceil(216/8)', ceil_div (216, 8), 27);
  same ('syndrome of the four-bit error d0 c0 c1 c2', bits (syndrome (h, [0 64 65 66]), 8), '00000000');
  hd = hsiao_72_64_dictionary ();
  require ('dictionary-order choice has distance 4', distance_is_4 (hd));
  loads_d = row_loads (hd, 8);
  same ('ones in each row, dictionary order', joined (loads_d), '30,30,30,26,26,25,25,24');
  same ('ones in the dictionary-order matrix', sum (loads_d), 216);
  same ('two-input levels for 27 and 30 inputs', sprintf ('%d,%d', levels (27), levels (30)), '5,5');
  same ('three-input levels for 27 and 30 inputs', sprintf ('%d,%d', levels (27, 3), levels (30, 3)), '3,4');
  same ('two-input levels for 72 inputs', levels (72), 7);
  h39 = hsiao_39_32 ();
  require ('Hsiao (39,32) has distance exactly 4', distance_is_4 (h39));
  same ('weight-3 columns of length 7, C(7,3)', nchoosek (7, 3), 35);
  loads39 = row_loads (h39, 7);
  same ('ones in each row of Hsiao (39,32)', joined (loads39), '14,14,15,15,15,15,15');
  same ('ones in Hsiao (39,32)', sum (loads39), 103);
  same ('largest row: ext. Hamming, dictionary, balanced (72,64)', ...
        sprintf ('%d,%d,%d', max (row_loads (extended_hamming (72, 7), 8)), max (loads_d), max (loads)), '72,30,27');
  same ('largest row: ext. Hamming and Hsiao (39,32)', ...
        sprintf ('%d,%d', max (row_loads (extended_hamming (39, 6), 7)), max (loads39)), '39,15');
  same ('lower bound 3k + r for k = 32, r = 7', 3 * 32 + 7, 103);
  same ('least possible largest row, ceil(103/7)', ceil_div (103, 7), 15);
  same ('weight-3 columns of length 7 that contain a row, C(6,2)', nchoosek (6, 2), 15);
  h22 = hsiao_22_16 ();
  require ('Hsiao (22,16) has distance exactly 4', distance_is_4 (h22));
  same ('ones in each row of Hsiao (22,16)', joined (row_loads (h22, 6)), '9,9,9,9,9,9');
  same ('ones in Hsiao (22,16), 6 + 3*16', sum (row_loads (h22, 6)), 54);
  same ('weight-3 columns of length 6, C(6,3)', nchoosek (6, 3), 20);
  % a 128-bit word
  same ('weight-3 columns of length 9, C(9,3)', nchoosek (9, 3), 84);
  same ('weight-5 columns needed for k = 128', 128 - nchoosek (9, 3), 44);
  same ('ones in Hsiao (137,128), 9 + 3*84 + 5*44', 9 + 3 * 84 + 5 * 44, 481);
  same ('least possible largest row, ceil(481/9)', ceil_div (481, 9), 54);
  same ('two-input levels for 54 inputs', levels (54), 6);
  h137 = hsiao_137_128 ();
  require ('Hsiao (137,128) columns are 137 distinct odd words', ...
           numel (unique (h137)) == 137 && all (mod (weight (h137), 2) == 1));
  same ('ones in each row of a balanced Hsiao (137,128)', joined (row_loads (h137, 9)), ...
        '54,54,54,53,54,53,53,53,53');
  % systematic SECDED codes have no data column of weight 2
  ok = true;
  P = nchoosek (0:7, 2);
  for q = 1:rows (P)
    ok = ok && ~distance_at_least_4 ([2 .^ (0:7), col(P(q, :))]);
  end
  require ('a weight-2 column beside I_8 always breaks distance 4', ok);
end

function check_decoder ()
  section ('The decoder in gates');
  same ('gates in a tree over w inputs, w = 27', 27 - 1, 26);
  same ('levels of a tree over 27 inputs', levels (27), 5);
  hs = decoder_gates (hsiao_72_64 (), 8, 64, -1);
  eh = decoder_gates (extended_hamming (72, 7), 8, 64, 7);
  same ('Hsiao (72,64) syndrome XOR gates', hs.syndrome, 208);
  same ('Hsiao (72,64) syndrome levels', hs.syndrome_depth, 5);
  same ('Hsiao (72,64) zero detector OR gates', hs.zero, 7);
  same ('Hsiao (72,64) parity detector XOR gates', hs.parity, 7);
  same ('Hsiao (72,64) match AND gates, 64 x 7', hs.match, 448);
  same ('Hsiao (72,64) correction XOR gates', hs.correct, 64);
  same ('Hsiao (72,64) double-error AND gate', hs.ded, 1);
  same ('Hsiao (72,64) two-input gates in all', hs.total, 863);
  same ('Hsiao (72,64) inverters', hs.inverters, 10);
  same ('Hsiao (72,64) levels to corrected data', hs.path, 9);
  same ('Hsiao (72,64) levels to the double-error flag', hs.ded_depth, 9);
  same ('Hsiao (72,64) encoder XOR gates, 8 x 25', 8 * (26 - 1), 200);
  same ('extended Hamming (72,64) syndrome XOR gates', eh.syndrome, 276);
  same ('extended Hamming (72,64) syndrome levels', eh.syndrome_depth, 7);
  same ('extended Hamming (72,64) zero detector OR gates', eh.zero, 6);
  same ('extended Hamming (72,64) parity detector gates', eh.parity, 0);
  same ('extended Hamming (72,64) two-input gates in all', eh.total, 923);
  same ('extended Hamming (72,64) levels to corrected data', eh.path, 11);
  same ('extended Hamming (72,64) levels to the double-error flag', eh.ded_depth, 10);
  same ('syndrome gates saved by Hsiao (72,64)', eh.syndrome - hs.syndrome, 68);
  same ('gates saved by Hsiao (72,64) in all', eh.total - hs.total, 60);
  same ('gates Hsiao (72,64) spends that extended Hamming does not', 68 - (eh.total - hs.total), 8);
  hd = decoder_gates (hsiao_72_64_dictionary (), 8, 64, -1);
  same ('dictionary-order syndrome XOR gates and levels', sprintf ('%d,%d', hd.syndrome, hd.syndrome_depth), '208,5');
  h39 = decoder_gates (hsiao_39_32 (), 7, 32, -1);
  e39 = decoder_gates (extended_hamming (39, 6), 7, 32, 6);
  same ('Hsiao (39,32) syndrome XOR gates', h39.syndrome, 96);
  same ('Hsiao (39,32) syndrome levels', h39.syndrome_depth, 4);
  same ('Hsiao (39,32) zero detector OR gates', h39.zero, 6);
  same ('Hsiao (39,32) parity detector XOR gates', h39.parity, 6);
  same ('Hsiao (39,32) match AND gates, 32 x 6', h39.match, 192);
  same ('Hsiao (39,32) two-input gates in all', h39.total, 414);
  same ('Hsiao (39,32) levels to corrected data', h39.path, 8);
  same ('Hsiao (39,32) levels to the double-error flag', h39.ded_depth, 8);
  same ('extended Hamming (39,32) syndrome XOR gates', e39.syndrome, 128);
  same ('extended Hamming (39,32) syndrome levels', e39.syndrome_depth, 6);
  same ('extended Hamming (39,32) two-input gates in all', e39.total, 439);
  same ('extended Hamming (39,32) zero detector OR gates', e39.zero, 5);
  same ('extended Hamming (39,32) levels to corrected data', e39.path, 10);
  same ('extended Hamming (39,32) levels to the double-error flag', e39.ded_depth, 9);
  same ('Hsiao (72,64) check-column match AND gates, 8 x 7', hs.check_match, 56);
  same ('Hsiao (72,64) any-match OR gates, 72 - 1', hs.any_match, 71);
  same ('Hsiao (72,64) no-match AND gate', hs.no_match, 1);
  same ('gates the no-match flag adds, 56 + 71 + 1', hs.check_match + hs.any_match + hs.no_match, 128);
  same ('Hsiao (72,64) gates before the no-match flag', hs.total - 128, 735);
  same ('Hsiao (72,64) levels to the any-match signal M', hs.m_depth, 15);
  same ('Hsiao (72,64) levels to the no-match flag', hs.no_match_depth, 16);
  same ('extended Hamming (72,64) levels to M and to the no-match flag', ...
        sprintf ('%d,%d', eh.m_depth, eh.no_match_depth), '17,18');
  same ('Hsiao (39,32) gates the no-match flag adds, 42 + 38 + 1', ...
        h39.check_match + h39.any_match + h39.no_match, 81);
  same ('Hsiao (39,32) levels to the no-match flag', h39.no_match_depth, 14);
  same ('extended Hamming (39,32) levels to the no-match flag', e39.no_match_depth, 16);
  same ('Hsiao (39,32) gates before the no-match flag', h39.total - 81, 333);
  same ('check match, any-match OR and no-match AND of the four decoders', ...
        sprintf ('%d,%d,%d %d,%d,%d %d,%d,%d %d,%d,%d', ...
                 hs.check_match, hs.any_match, hs.no_match, eh.check_match, eh.any_match, eh.no_match, ...
                 h39.check_match, h39.any_match, h39.no_match, e39.check_match, e39.any_match, e39.no_match), ...
        '56,71,1 56,71,1 42,38,1 42,38,1');
  names = {'Hsiao (72,64)', 'extended Hamming (72,64)'};
  codes = {hsiao_72_64(), extended_hamming(72, 7)};
  prow = [-1, 7];
  for q = 1:2
    one = decoder_actions (codes{q}, 8, prow(q), 1);
    two = decoder_actions (codes{q}, 8, prow(q), 2);
    same ([names{q} ' decoder: single errors corrected'], one.corrected, 72);
    same ([names{q} ' decoder: double errors flagged, C(72,2)'], two.double, 2556);
    n1 = one.silent + one.corrected + one.miscorrected + one.double + one.no_match + one.conflict;
    n2 = two.silent + two.corrected + two.miscorrected + two.double + two.no_match + two.conflict;
    require ([names{q} ' decoder: no single or double error ends any other way'], ...
             n1 == one.corrected && n2 == two.double);
  end
end

function check_triples ()
  section ('Three errors');
  h = hsiao_72_64 ();
  T = nchoosek (1:numel (h), 3);
  a = h(T(:, 1));
  b = h(T(:, 2));
  c = h(T(:, 3));
  x = bitxor (bitxor (a, b), c);
  ok = all (mod (weight (x), 2) == 1 & x ~= a & x ~= b & x ~= c);
  require ('a triple''s syndrome is odd and is none of its three columns', ok);
  same ('three-bit errors in 72 bits, C(72,3)', nchoosek (72, 3), 59640);
  [mis, det] = triple_outcomes (h);
  same ('Hsiao (72,64) triples miscorrected', mis, 33568);
  same ('Hsiao (72,64) triples detected', det, 26072);
  same ('Hsiao (72,64) percent miscorrected', pct (mis, mis + det), '56.3');
  same ('72 columns out of 128 odd words, percent', sprintf ('%.2f', 100 * 72 / 128), '56.25');
  tw = triple_by_weight (h, 8);
  same ('Hsiao (72,64) triples by syndrome weight 1, 3, 5, 7', joined (sum (tw([2 4 6 8], :), 2)), ...
        '3952,26112,26040,3536');
  same ('Hsiao (72,64) weight-1 and weight-3 triples detected', tw(2, 2) + tw(4, 2), 0);
  same ('Hsiao (72,64) weight-5 triples miscorrected', tw(6, 1), 3504);
  same ('Hsiao (72,64) weight-5 triples detected', tw(6, 2), 22536);
  same ('Hsiao (72,64) weight-7 triples miscorrected', tw(8, 1), 0);
  eh = extended_hamming (72, 7);
  [mis_e, det_e] = triple_outcomes (eh);
  same ('extended Hamming (72,64) triples miscorrected', mis_e, 45304);
  same ('extended Hamming (72,64) triples detected', det_e, 14336);
  same ('extended Hamming (72,64) percent miscorrected', pct (mis_e, mis_e + det_e), '76.0');
  x = triple_sums (eh(1:64));
  inside = sum (ismember (x, eh));
  same ('triples inside positions 0..63, C(64,3)', nchoosek (64, 3), 41664);
  same ('triples inside positions 0..63 that miscorrect', inside, 41664);
  [mis_d, det_d] = triple_outcomes (hsiao_72_64_dictionary ());
  same ('dictionary-order Hsiao (72,64) triples miscorrected', mis_d, 34032);
  same ('dictionary-order Hsiao (72,64) triples detected', det_d, 25608);
  same ('dictionary-order Hsiao (72,64) percent miscorrected', pct (mis_d, mis_d + det_d), '57.1');
  same ('extended Hamming miscorrected triples touching 64..71', mis_e - 41664, 3640);
  same ('more miscorrected triples in dictionary order than balanced', mis_d - mis, 464);
  [mis39, det39] = triple_outcomes (hsiao_39_32 ());
  same ('three-bit errors in 39 bits, C(39,3)', nchoosek (39, 3), 9139);
  same ('Hsiao (39,32) triples miscorrected', mis39, 5452);
  same ('Hsiao (39,32) triples detected', det39, 3687);
  same ('Hsiao (39,32) percent miscorrected', pct (mis39, 9139), '59.7');
  [mis39e, det39e] = triple_outcomes (extended_hamming (39, 6));
  same ('extended Hamming (39,32) triples miscorrected', mis39e, 6332);
  same ('extended Hamming (39,32) triples detected', det39e, 2807);
  same ('extended Hamming (39,32) percent miscorrected', pct (mis39e, 9139), '69.3');
  [mis13, det13] = triple_outcomes (hsiao_13_8 ());
  same ('three-bit errors in 13 bits, C(13,3)', nchoosek (13, 3), 286);
  same ('Hsiao (13,8) triples miscorrected', mis13, 220);
  same ('Hsiao (13,8) triples detected', det13, 66);
  same ('Hsiao (13,8) percent miscorrected', pct (mis13, 286), '76.9');
  [mis8, det8] = triple_outcomes (hsiao_8_4 ());
  same ('Hsiao (8,4) triples miscorrected', mis8, 56);
  same ('Hsiao (8,4) triples detected', det8, 0);
  names = {'Hsiao (72,64)', 'dictionary-order Hsiao (72,64)', 'extended Hamming (72,64)', ...
           'Hsiao (39,32)', 'extended Hamming (39,32)', 'Hsiao (13,8)', 'Hsiao (8,4)'};
  codes = {h, hsiao_72_64_dictionary(), eh, hsiao_39_32(), extended_hamming(39, 6), ...
           hsiao_13_8(), hsiao_8_4()};
  r = [8, 8, 8, 7, 7, 5, 4];
  prow = [-1, -1, 7, -1, 6, -1, -1];
  stated = {'33568,26072', '34032,25608', '45304,14336', '5452,3687', '6332,2807', '220,66', '56,0'};
  for q = 1:numel (names)
    t = decoder_actions (codes{q}, r(q), prow(q), 3);
    same ([names{q} ' full decoder on triples: miscorrected, no-match flag'], ...
          sprintf ('%d,%d', t.miscorrected, t.no_match), stated{q});
    require ([names{q} ' full decoder: no triple passes unflagged or raises the double-error flag'], ...
             t.silent + t.double + t.conflict == 0);
  end
end


% Run --------------------------------------------------------------------------

check_bit_budget ();
check_extended_hamming ();
check_odd_weight ();
check_hsiao_8_4 ();
check_hsiao_13_8 ();
check_hsiao_72_64 ();
check_decoder ();
check_triples ();
printf ('\nAll chapter 20 numbers reproduced.\n');
