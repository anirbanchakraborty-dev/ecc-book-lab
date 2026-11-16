% Reproduce every number in chapter 26 of Error Correction Coding,
% "Reed-Muller Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch26_reed_muller_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch26_reed_muller_codes.py and prints the same
% lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in the chapter is a whole number, a binary word, a
% list of integers or a probability printed to six places, so every
% comparison is exact.
%
% A word is a row of 0s and 1s, and position j of a word is its element
% j + 1. The points of GF(2)^m are numbered j = 0, 1, ..., 2^m - 1, and
% variable x_i is bit i - 1 of j, so x_1 is the least significant bit. A
% code is a matrix with one codeword per row. Every code with at most 2,048
% codewords is listed whole, and so is RM(2,5) with 65,536. RM(3,5) and
% RM(4,5) get their weights from their duals by the MacWilliams identity, and
% the distances of RM(2,7), RM(3,6) and RM(1,10) come from the proved
% formula. Reed's majority-logic decoder
% runs here on a whole matrix of received words at once, and the 4,514,873
% error patterns of weight 0 to 7 on the all-zero codeword of RM(1,5) are
% built in blocks with nchoosek, one block for each first error position.

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

function s = word (v)
  s = sprintf ('%d', v);
end

function s = vec (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ' ');
end

function c = popcount (x)
  c = zeros (size (x));
  while (any (x(:) > 0))
    c = c + bitand (x, 1);
    x = bitshift (x, -1);
  end
end

function b = bits_of (x, n)
  % Row of the n low bits of x, least significant first.
  b = bitand (bitshift (x, -(0:(n - 1))), 1);
end

function W = all_words (n)
  % Every word of length n, one per row, row x + 1 holding the bits of x.
  x = (0:(2 ^ n - 1))';
  W = zeros (2 ^ n, n);
  for j = 0:(n - 1)
    W(:, j + 1) = bitand (bitshift (x, -j), 1);
  end
end

function P = combs (v, k)
  % Every k-subset of the vector v as a row, in the order itertools.combinations gives.
  if (k == 0)
    P = zeros (1, 0);
  elseif (numel (v) == k)
    P = v(:)';
  else
    P = nchoosek (v, k);
  end
end

function C = combos (m, k)
  C = combs (1:m, k);
end


% The codes --------------------------------------------------------------------

function M = monomials (r, m)
  % Every monomial of degree at most r, by degree, then in dictionary order.
  M = {};
  for d = 0:r
    C = combos (m, d);
    for i = 1:size (C, 1)
      M{end + 1} = C(i, :);
    end
  end
end

function s = mono_name (S)
  if (isempty (S))
    s = '1';
  else
    s = sprintf ('x%d', S);
  end
end

function t = truth_table (S, m)
  j = 0:(2 ^ m - 1);
  t = ones (1, 2 ^ m);
  for i = S
    t = t .* bitand (bitshift (j, -(i - 1)), 1);
  end
end

function G = generator (r, m)
  M = monomials (r, m);
  G = zeros (numel (M), 2 ^ m);
  for i = 1:numel (M)
    G(i, :) = truth_table (M{i}, m);
  end
end

function c = encode (coeffs, r, m)
  c = mod (coeffs(:)' * generator (r, m), 2);
end

function C = span_rows (G)
  k = size (G, 1);
  C = mod (all_words (k) * G, 2);
end

function C = rm_code (r, m)
  C = span_rows (generator (r, m));
end

function k = k_rm (r, m)
  k = 0;
  for i = 0:r
    k = k + nchoosek (m, i);
  end
end

function s = dist_text (C)
  w = sum (C, 2);
  u = unique (w)';
  parts = cell (1, numel (u));
  for i = 1:numel (u)
    parts{i} = sprintf ('%d:%d', u(i), sum (w == u(i)));
  end
  s = strjoin (parts, ' ');
end

function c = weight_count (C, w)
  c = sum (sum (C, 2) == w);
end

function d = min_weight (C)
  w = sum (C, 2);
  d = min (w(w > 0));
end

function e = same_set (A, B)
  e = isequal (size (A), size (B)) && isequal (sortrows (A), sortrows (B));
end

function st = mobius_stages (f, m)
  % One variable per stage: a(j) = a(j) + a(j - 2^(i-1)) where bit i-1 of j is 1.
  a = f;
  n = numel (f);
  st = zeros (m, n);
  j = 0:(n - 1);
  for i = 0:(m - 1)
    idx = find (bitand (bitshift (j, -i), 1) == 1);
    b = a;
    b(idx) = mod (a(idx) + a(idx - 2 ^ i), 2);
    a = b;
    st(i + 1, :) = a;
  end
end

function s = anf_terms (coeffs, m)
  masks = find (coeffs) - 1;
  if (isempty (masks))
    s = '0';
    return;
  end
  keys = zeros (numel (masks), m + 1);
  for q = 1:numel (masks)
    vars = find (bits_of (masks(q), m));
    keys(q, 1) = numel (vars);
    keys(q, 2:(numel (vars) + 1)) = vars;
  end
  [~, order] = sortrows (keys);
  parts = cell (1, numel (masks));
  for q = 1:numel (order)
    parts{q} = mono_name (find (bits_of (masks(order(q)), m)));
  end
  s = strjoin (parts, ' + ');
end


% Reed's majority-logic decoder ------------------------------------------------

function Q = subcube (S, b, m)
  % Positions whose variables outside S equal the bits of b, the variables in S running free.
  fixed = setdiff (1:m, S);
  j = 0:(2 ^ m - 1);
  keep = true (1, 2 ^ m);
  for t = 1:numel (fixed)
    keep = keep & (bitand (bitshift (j, -(fixed(t) - 1)), 1) == bitand (bitshift (b, -(t - 1)), 1));
  end
  Q = j(keep);
end

function [coeffs, trace] = reed_decode (y, r, m)
  % Coefficients in monomial order, or [] on a tied vote, and the trace of every vote.
  n = 2 ^ m;
  coeffs = zeros (1, k_rm (r, m));
  trace = {};
  for deg = r:-1:0
    part = zeros (1, n);
    C = combos (m, deg);
    base = 0;
    for i = 0:(deg - 1)
      base = base + nchoosek (m, i);
    end
    nb = 2 ^ (m - deg);
    for q = 1:size (C, 1)
      S = C(q, :);
      sums = zeros (1, nb);
      for b = 0:(nb - 1)
        sums(b + 1) = mod (sum (y(subcube (S, b, m) + 1)), 2);
      end
      cnt = sum (sums);
      trace{end + 1} = word (sums);
      if (2 * cnt == nb)
        coeffs = [];
        return;
      end
      a = 2 * cnt > nb;
      coeffs(base + q) = a;
      if (a)
        part = mod (part + truth_table (S, m), 2);
      end
    end
    y = mod (y + part, 2);
    trace{end + 1} = word (y);
  end
end

function plan = reed_plan (r, m)
  % For each degree from r down to 0, the coset positions, monomial and vote size of every coefficient.
  plan = {};
  for deg = r:-1:0
    C = combos (m, deg);
    level = {};
    nb = 2 ^ (m - deg);
    for q = 1:size (C, 1)
      S = C(q, :);
      Q = zeros (nb, 2 ^ deg);
      for b = 0:(nb - 1)
        Q(b + 1, :) = subcube (S, b, m) + 1;
      end
      level{end + 1} = struct ('Q', Q, 'mono', logical (truth_table (S, m)), 'nsums', nb);
    end
    plan{end + 1} = level;
  end
end

function [D, tie] = reed_batch (Y, plan)
  % The decoder on every row of Y at once. D holds the decoded codewords, tie marks a tied vote.
  Y = logical (Y);
  N = size (Y, 1);
  D = false (size (Y));
  tie = false (N, 1);
  for L = 1:numel (plan)
    level = plan{L};
    part = false (size (Y));
    for q = 1:numel (level)
      st = level{q};
      P = Y(:, st.Q(:, 1));
      for t = 2:size (st.Q, 2)
        P = xor (P, Y(:, st.Q(:, t)));
      end
      cnt = sum (P, 2);
      tie = tie | (2 * cnt == st.nsums);
      a = 2 * cnt > st.nsums;
      if (any (a))
        part(a, :) = xor (part(a, :), repmat (st.mono, sum (a), 1));
      end
    end
    Y = xor (Y, part);
    D = xor (D, part);
  end
end

function E = pattern_rows (P, n)
  N = size (P, 1);
  E = false (N, n);
  if (size (P, 2) > 0)
    E(sub2ind ([N n], repmat ((1:N)', 1, size (P, 2)), P)) = true;
  end
end

function [npat, tried, failures] = exhaustive_reed (r, m, every_codeword)
  % Patterns tried, decodings tried and failures, for every error pattern up to the radius.
  n = 2 ^ m;
  t = 2 ^ (m - r - 1) - 1;
  plan = reed_plan (r, m);
  if (every_codeword)
    W = logical (rm_code (r, m));
  else
    W = false (1, n);
  end
  nw = size (W, 1);
  npat = 0;
  tried = 0;
  failures = 0;
  for w = 0:t
    if (w <= 4)
      firsts = 0;
    else
      firsts = 1:(n - w + 1);
    end
    for a = firsts
      if (a == 0)
        P = combs (1:n, w);
      else
        R = combs ((a + 1):n, w - 1);
        P = [repmat(a, size (R, 1), 1), R];
      end
      E = pattern_rows (P, n);
      ne = size (E, 1);
      npat = npat + ne;
      if (nw == 1)
        Crep = repmat (W, ne, 1);
        Y = xor (E, Crep);
      else
        Crep = W(kron ((1:nw)', ones (ne, 1)), :);
        Y = xor (Crep, repmat (E, nw, 1));
      end
      [D, tie] = reed_batch (Y, plan);
      fail = tie | any (xor (D, Crep), 2);
      tried = tried + size (Y, 1);
      failures = failures + sum (fail);
    end
  end
end


% The fast Hadamard transform --------------------------------------------------

function st = fht_stages (Y)
  % Butterflies (a + b, a - b) on positions j and j + h, for h = 2^(m-1) down to 1.
  y = Y;
  n = numel (Y);
  st = zeros (round (log2 (n)), n);
  j = 0:(n - 1);
  h = n / 2;
  s = 0;
  while (h >= 1)
    lo = find (bitand (j, h) == 0);
    z = y;
    z(lo) = y(lo) + y(lo + h);
    z(lo + h) = y(lo) - y(lo + h);
    y = z;
    s = s + 1;
    st(s, :) = y;
    h = h / 2;
  end
end

function s = signs (y)
  s = 1 - 2 * y;
end

function H = hadamard (m)
  n = 2 ^ m;
  A = repmat ((0:(n - 1))', 1, n);
  J = repmat (0:(n - 1), n, 1);
  H = 1 - 2 * mod (popcount (bitand (A, J)), 2);
end

function T = direct_correlation (Y, m)
  T = (hadamard (m) * Y(:))';
end

function k = fht_additions (m)
  k = m * 2 ^ m;
end

function k = direct_additions (m)
  k = 2 ^ m * (2 ^ m - 1);
end

function [best, a0, T] = fht_decode (y, m)
  % Index of the largest |T|, a0 from its sign, and the transform. Ties go to the lowest index.
  st = fht_stages (signs (y));
  T = st(end, :);
  [~, i] = max (abs (T));
  best = i - 1;
  a0 = double (T(i) < 0);
end

function c = rm1_codeword (a0, a, m)
  j = 0:(2 ^ m - 1);
  c = mod (a0 + popcount (bitand (a, j)), 2);
end

function p = word_error (n, t, q)
  p = 0;
  for j = (t + 1):n
    p = p + nchoosek (n, j) * q ^ j * (1 - q) ^ (n - j);
  end
end

function s = positions_text (S, m)
  parts = cell (1, 4);
  for b = 0:3
    parts{b + 1} = sprintf ('%d', subcube (S, b, m));
  end
  s = strjoin (parts, ' ');
end


% The checks -------------------------------------------------------------------

function B = macwilliams (C, n)
  % Weight distribution of the dual of the code whose words are the rows of C, by the MacWilliams identity.
  w = sum (C, 2);
  A = zeros (1, n + 1);
  for i = 0:n
    A(i + 1) = sum (w == i);
  end
  P = zeros (n + 1, n + 1);
  for a = 0:n
    for b = 0:a
      P(a + 1, b + 1) = nchoosek (a, b);
    end
  end
  B = zeros (1, n + 1);
  for j = 0:n
    total = 0;
    for i = find (A) - 1
      K = 0;
      for l = 0:min (j, i)
        if (j - l <= n - i)
          K = K + (-1) ^ l * P(i + 1, l + 1) * P(n - i + 1, j - l + 1);
        end
      end
      total = total + A(i + 1) * K;
    end
    if (mod (total, size (C, 1)) ~= 0)
      error ('MacWilliams sum not divisible by the code size');
    end
    B(j + 1) = total / size (C, 1);
  end
end

% Hardware models ---------------------------------------------------------------

function net = serial_encoder_rm13 ()
  % The serial RM(1,3) encoder as a netlist: gates (kind, input, input, output), evaluated once per clock.
  net = {'AND', 'a1', 'x1', 'p1'; 'AND', 'a2', 'x2', 'p2'; 'AND', 'a3', 'x3', 'p3'; ...
         'XOR', 'a0', 'p1', 's1'; 'XOR', 'p2', 'p3', 's2'; 'XOR', 's1', 's2', 'c'};
end

function sig = run_netlist (net, sig)
  for g = 1:size (net, 1)
    if (strcmp (net{g, 1}, 'AND'))
      sig.(net{g, 4}) = bitand (sig.(net{g, 2}), sig.(net{g, 3}));
    else
      sig.(net{g, 4}) = bitxor (sig.(net{g, 2}), sig.(net{g, 3}));
    end
  end
end

function [out, carry, cells] = ripple (a, b, w, subtract)
  % A w-bit ripple-carry adder or subtracter of full-adder cells, applied to every element of the
  % arrays a and b of w-bit two's-complement patterns. The subtracter is the adder fed the inverted b
  % and a carry-in of 1. Returns the w-bit results, the carries out and the cells in one unit.
  if (subtract)
    b = bitxor (b, 2 ^ w - 1);
  end
  carry = double (subtract) * ones (size (a));
  out = zeros (size (a));
  for i = 0:(w - 1)
    x = bitand (bitshift (a, -i), 1);
    y = bitand (bitshift (b, -i), 1);
    out = out + bitxor (bitxor (x, y), carry) * 2 ^ i;
    carry = bitor (bitand (x, y), bitand (carry, bitxor (x, y)));
  end
  cells = w;
end

function p = to_bits (v, w)
  p = mod (v, 2 ^ w);
end

function v = from_bits (p, w)
  v = p - (p >= 2 ^ (w - 1)) * 2 ^ w;
end

function w = width_for (bound)
  % Fewest two's-complement bits that hold every whole number from -bound to bound.
  w = 1;
  while (2 ^ (w - 1) - 1 < bound)
    w = w + 1;
  end
end

function [T, widths, units, cells] = butterfly_array (Y, carry_as_new_bit)
  % The fast transform of every row of Y, built from fixed-width two's-complement adders and
  % subtracters. Each stage sign-extends its inputs to the width that holds the stage's range and
  % drops the top carry; with carry_as_new_bit it adds in the old width and keeps the carry as the
  % new top bit. Returns the correlations, the widths after the sign map and each stage, the adders
  % and subtracters, and the full-adder cells of one array.
  if (nargin < 2)
    carry_as_new_bit = false;
  end
  n = size (Y, 2);
  bound = 1;
  w = width_for (bound);
  vals = to_bits (Y, w);
  widths = w;
  units = 0;
  cells = 0;
  j = 0:(n - 1);
  h = n / 2;
  while (h >= 1)
    bound = 2 * bound;
    wn = width_for (bound);
    lo = find (bitand (j, h) == 0);
    a = vals(:, lo);
    b = vals(:, lo + h);
    z = vals;
    if (carry_as_new_bit)
      [s, cs, c1] = ripple (a, b, w, false);
      [d, cd, c2] = ripple (a, b, w, true);
      z(:, lo) = s + cs * 2 ^ w;
      z(:, lo + h) = d + cd * 2 ^ w;
    else
      a = to_bits (from_bits (a, w), wn);
      b = to_bits (from_bits (b, w), wn);
      [s, ~, c1] = ripple (a, b, wn, false);
      [d, ~, c2] = ripple (a, b, wn, true);
      z(:, lo) = s;
      z(:, lo + h) = d;
    end
    units = units + 2 * numel (lo);
    cells = cells + numel (lo) * (c1 + c2);
    vals = z;
    w = wn;
    widths(end + 1) = w;
    h = h / 2;
  end
  T = from_bits (vals, w);
end

function [best, a0, comparators, levels] = max_finder (T)
  % A tree of two-input comparators on |T| that carries the index and sign of the winner.
  % Ties keep the left input, so the lowest index wins.
  mag = abs (T);
  idx = 0:(numel (T) - 1);
  sgn = double (T < 0);
  comparators = 0;
  levels = 0;
  while (numel (mag) > 1)
    L = 1:2:numel (mag);
    R = L + 1;
    pick = R;
    keep = mag(L) >= mag(R);
    pick(keep) = L(keep);
    mag = mag(pick);
    idx = idx(pick);
    sgn = sgn(pick);
    comparators = comparators + numel (L);
    levels = levels + 1;
  end
  best = idx;
  a0 = sgn;
end

function check_truth_tables ()
  section ('Boolean functions as codewords, m = 3');
  same ('point j = 3 as (x1, x2, x3)', vec (bits_of (3, 3)), '1 1 0');
  same ('point j = 6 as (x1, x2, x3)', vec (bits_of (6, 3)), '0 1 1');
  names = {[], 1, 2, 3, [1 2], [1 3], [2 3], [1 2 3]};
  stated = {'11111111', '01010101', '00110011', '00001111', '00010001', '00000101', '00000011', '00000001'};
  for i = 1:numel (names)
    S = names{i};
    same (sprintf ('truth table of %s', mono_name (S)), word (truth_table (S, 3)), stated{i});
  end
  same ('Boolean functions of 3 variables', 2 ^ 8, 256);
  same ('monomials in 3 variables', numel (monomials (3, 3)), 8);
  same ('codeword of 1 + x2 + x3 in RM(1,3)', word (encode ([1 0 1 1], 1, 3)), '11000011');
  same ('weight of that codeword', sum (encode ([1 0 1 1], 1, 3)), 4);
  G = generator (1, 3);
  same ('RM(1,3) generator, rows 1, x1, x2, x3', strjoin ({word(G(1, :)), word(G(2, :)), word(G(3, :)), word(G(4, :))}, ' '), '11111111 01010101 00110011 00001111');
  same ('codewords of RM(1,3)', size (rm_code (1, 3), 1), 16);
  same ('distinct words spanned by the 8 monomials', size (unique (rm_code (3, 3), 'rows'), 1), 256);
end

function check_mobius ()
  section ('The Moebius transform');
  f = mod (truth_table ([1 2], 3) + truth_table (3, 3), 2);
  same ('truth table of x3 + x1x2', word (f), '00011110');
  st = mobius_stages (f, 3);
  same ('after the x1 stage', word (st(1, :)), '00011011');
  same ('after the x2 stage', word (st(2, :)), '00011001');
  same ('after the x3 stage, the coefficients', word (st(3, :)), '00011000');
  same ('polynomial read from the coefficients', anf_terms (st(3, :), 3), 'x3 + x1x2');
  ok_inv = true;
  ok_sum = true;
  W = all_words (8);
  for x = 0:255
    t = W(x + 1, :);
    s1 = mobius_stages (t, 3);
    a = s1(3, :);
    s2 = mobius_stages (a, 3);
    ok_inv = ok_inv && isequal (s2(3, :), t);
    for s = 0:7
      u = 0:7;
      sub = u(bitand (u, s) == u);
      ok_sum = ok_sum && (a(s + 1) == mod (sum (t(sub + 1)), 2));
    end
  end
  require ('transform applied twice returns all 256 tables', ok_inv);
  require ('coefficient a_S equals the sum of f over subsets of S', ok_sum);
  same ('XOR gates in the transform for m = 3', 3 * 4, 12);
end

function check_dimensions ()
  section ('Dimension and the special cases');
  cases = [1 3 4; 2 3 7; 1 4 5; 2 4 11; 1 5 6; 2 5 16; 3 5 26];
  for i = 1:size (cases, 1)
    same (sprintf ('k of RM(%d,%d)', cases(i, 1), cases(i, 2)), k_rm (cases(i, 1), cases(i, 2)), cases(i, 3));
  end
  for m = 1:5
    for r = 0:m
      if (k_rm (r, m) <= 11)
        k = k_rm (r, m);
        require (sprintf ('RM(%d,%d) has 2^k distinct codewords, k = %d', r, m, k), size (unique (rm_code (r, m), 'rows'), 1) == 2 ^ k);
      end
    end
  end
  stated = {'1 2', '1 3 4', '1 4 7 8', '1 5 11 15 16', '1 6 16 26 31 32', '1 7 22 42 57 63 64', '1 8 29 64 99 120 127 128'};
  for m = 1:7
    ks = zeros (1, m + 1);
    for r = 0:m
      ks(r + 1) = k_rm (r, m);
    end
    same (sprintf ('k for r = 0..%d, m = %d', m, m), vec (ks), stated{m});
  end
  for m = 1:5
    n = 2 ^ m;
    require (sprintf ('RM(0,%d) is the repetition code', m), same_set (rm_code (0, m), [zeros(1, n); ones(1, n)]));
  end
  W = all_words (8);
  require ('RM(2,3) is the even-weight code of length 8', same_set (rm_code (2, 3), W(mod (sum (W, 2), 2) == 0, :)));
  require ('RM(3,3) is every word of length 8', same_set (rm_code (3, 3), W));
  for m = 1:7
    require (sprintf ('k(m-1,%d) = 2^m - 1', m), k_rm (m - 1, m) == 2 ^ m - 1);
  end
end

function check_recursion ()
  section ('The (u | u + v) recursion');
  cases = [1 2; 1 3; 2 3; 1 4; 2 4; 1 5];
  for i = 1:size (cases, 1)
    r = cases(i, 1);
    m = cases(i, 2);
    U = rm_code (r, m - 1);
    V = rm_code (r - 1, m - 1);
    nu = size (U, 1);
    nv = size (V, 1);
    Ur = U(kron ((1:nu)', ones (nv, 1)), :);
    Vr = repmat (V, nu, 1);
    built = unique ([Ur, mod(Ur + Vr, 2)], 'rows');
    require (sprintf ('RM(%d,%d) = (u | u + v), u in RM(%d,%d), v in RM(%d,%d)', r, m, r, m - 1, r - 1, m - 1), same_set (built, rm_code (r, m)));
  end
  for m = 2:7
    ok = true;
    for r = 1:(m - 1)
      ok = ok && (k_rm (r, m) == k_rm (r, m - 1) + k_rm (r - 1, m - 1));
    end
    require (sprintf ('k(r,%d) = k(r,%d) + k(r-1,%d) for 1 <= r <= %d', m, m - 1, m - 1, m - 1), ok);
  end
  G13 = generator (1, 3);
  G03 = generator (0, 3);
  G = [G13, G13; zeros(1, 8), G03(1, :)];
  same ('rows of the recursive RM(1,4) generator', size (G, 1), 5);
  same ('its last row, x4', word (G(5, :)), '0000000011111111');
  require ('it spans RM(1,4)', same_set (unique (span_rows (G), 'rows'), rm_code (1, 4)));
  g = encode ([1 1 0 0], 1, 3);
  v = encode (1, 0, 3);
  c = [g, mod(g + v, 2)];
  same ('u = 1 + x1 on m = 3', word (g), '10101010');
  same ('codeword (u | u + v) with v = 1', word (c), '1010101001010101');
  same ('the same word from 1 + x1 + x4', word (encode ([1 1 0 0 1], 1, 4)), '1010101001010101');
  same ('its weight', sum (c), 8);
end

function check_distance ()
  section ('Minimum distance');
  checked = 0;
  for m = 1:7
    for r = 0:m
      if (k_rm (r, m) <= 11)
        require (sprintf ('d of RM(%d,%d) = 2^%d = %d', r, m, m - r, 2 ^ (m - r)), min_weight (rm_code (r, m)) == 2 ^ (m - r));
        checked = checked + 1;
      end
    end
  end
  same ('codes enumerated whole', checked, 18);
  % RM(2,5) is small enough to list; RM(3,5) and RM(4,5) get their weights from their duals
  dd = containers.Map ();
  listed = [1 5; 2 4; 1 7; 2 5];
  for i = 1:size (listed, 1)
    dd(sprintf ('%d,%d', listed(i, 1), listed(i, 2))) = min_weight (rm_code (listed(i, 1), listed(i, 2)));
  end
  same ('RM(2,5): words listed', size (rm_code (2, 5), 1), 65536);
  cases = {3, 5, '67108864 4'; 4, 5, '2147483648 2'};
  for i = 1:size (cases, 1)
    r = cases{i, 1};
    m = cases{i, 2};
    B = macwilliams (rm_code (m - r - 1, m), 2 ^ m);
    d = find (B(2:end), 1);
    dd(sprintf ('%d,%d', r, m)) = d;
    same (sprintf ('RM(%d,%d) by MacWilliams from RM(%d,%d): words, and d', r, m, m - r - 1, m), sprintf ('%d %d', sum (B), d), cases{i, 3});
  end
  cases = {1, 5, '32 6 16'; 2, 5, '32 16 8'; 3, 5, '32 26 4'; 4, 5, '32 31 2'; 2, 4, '16 11 4'; 1, 7, '128 8 64'};
  for i = 1:size (cases, 1)
    r = cases{i, 1};
    m = cases{i, 2};
    same (sprintf ('(n, k, d) of RM(%d,%d)', r, m), vec ([2 ^ m, k_rm(r, m), dd(sprintf ('%d,%d', r, m))]), cases{i, 3});
  end
  same ('(n, k, d) of RM(2,7), d from the formula', vec ([2 ^ 7, k_rm(2, 7), 2 ^ (7 - 2)]), '128 29 32');
  cases = [1 3 1; 1 4 3; 2 4 1; 1 5 7; 2 5 3];
  for i = 1:size (cases, 1)
    r = cases(i, 1);
    m = cases(i, 2);
    key = sprintf ('%d,%d', r, m);
    if (isKey (dd, key))
      d = dd(key);
    else
      d = min_weight (rm_code (r, m));
    end
    same (sprintf ('guaranteed radius of RM(%d,%d), (d - 1) / 2 rounded down', r, m), floor ((d - 1) / 2), cases(i, 3));
  end
  same ('guaranteed radius of RM(3,6), from the formula', 2 ^ (6 - 3 - 1) - 1, 3);
end

function check_first_order ()
  section ('First-order weights');
  for m = 1:7
    C = rm_code (1, m);
    w = sum (C, 2);
    if (m == 1)
      ok = sum (w == 0) == 1 && sum (w == 1) == 2 && sum (w == 2) == 1 && numel (w) == 4;
    else
      ok = sum (w == 0) == 1 && sum (w == 2 ^ (m - 1)) == 2 ^ (m + 1) - 2 && sum (w == 2 ^ m) == 1 && numel (w) == 2 ^ (m + 1);
    end
    require (sprintf ('RM(1,%d): 1 at 0, 2^(m+1) - 2 at 2^(m-1), 1 at 2^m', m), ok);
  end
  same ('weight distribution of RM(1,3)', dist_text (rm_code (1, 3)), '0:1 4:14 8:1');
  same ('weight distribution of RM(1,4)', dist_text (rm_code (1, 4)), '0:1 8:30 16:1');
  same ('weight distribution of RM(1,5)', dist_text (rm_code (1, 5)), '0:1 16:62 32:1');
  same ('words of weight 32 in RM(1,6)', weight_count (rm_code (1, 6), 32), 126);
  W = all_words (8);
  keep = false (256, 1);
  for x = 1:256
    c = W(x, :);
    s = 0;
    for j = 1:7
      if (c(j + 1))
        s = bitxor (s, j);
      end
    end
    keep(x) = (s == 0) && (c(1) == mod (sum (c(2:8)), 2));
  end
  ext = W(keep, :);
  same ('weight distribution of the extended (8,4) Hamming code', dist_text (ext), '0:1 4:14 8:1');
  require ('RM(1,3) is that code, position for position', same_set (ext, rm_code (1, 3)));
  G = generator (1, 3);
  ok = all (sum (G, 1) >= 1);
  for j = 0:7
    ok = ok && isequal (G(:, j + 1)', [1, bits_of(j, 3)]);
  end
  require ('column j of the RM(1,3) generator is (1, bits of j)', ok);
end

function ok = ext_ok (C, m)
  % Each row: the positions holding 1s XOR to zero, and the weight is even.
  ok = true (size (C, 1), 1);
  for x = 1:size (C, 1)
    s = 0;
    for j = 0:(2 ^ m - 1)
      if (C(x, j + 1))
        s = bitxor (s, j);
      end
    end
    ok(x) = (s == 0) && (mod (sum (C(x, :)), 2) == 0);
  end
end

function check_duality ()
  section ('Duality');
  for m = 1:6
    for r = 0:(m - 1)
      A = generator (r, m);
      B = generator (m - r - 1, m);
      ok = all (all (mod (A * B', 2) == 0));
      require (sprintf ('RM(%d,%d) is orthogonal to RM(%d,%d), k sum = %d', r, m, m - r - 1, m, k_rm (r, m) + k_rm (m - r - 1, m)), ok && (k_rm (r, m) + k_rm (m - r - 1, m) == 2 ^ m));
    end
  end
  for m = 1:5
    G = generator (m - 1, m);
    require (sprintf ('every monomial of degree below %d has even weight, m = %d', m, m), all (mod (sum (G, 2), 2) == 0));
  end
  same ('weight of x1x2 at m = 4', sum (truth_table ([1 2], 4)), 4);
  for m = [3 4]
    words = rm_code (m - 2, m);
    require (sprintf ('every word of RM(%d,%d) passes both extended Hamming checks', m - 2, m), all (ext_ok (words, m)));
    same (sprintf ('words of length %d passing both checks, all %d listed', 2 ^ m, 2 ^ (2 ^ m)), sum (ext_ok (all_words (2 ^ m), m)), size (words, 1));
  end
  same ('weight distribution of RM(2,4)', dist_text (rm_code (2, 4)), '0:1 4:140 6:448 8:870 10:448 12:140 16:1');
  % the 140 words of weight 4, counted by the recursion from u in RM(2,3) and v in RM(1,3)
  U = rm_code (2, 3);
  V = rm_code (1, 3);
  by_v = [0 0 0];
  inside = [];
  for iv = 1:size (V, 1)
    v = V(iv, :);
    wt = sum (U, 2) + sum (xor (U, repmat (v, size (U, 1), 1)), 2);
    k = sum (v) / 4 + 1;
    by_v(k) = by_v(k) + sum (wt == 4);
    if (sum (v) == 4)
      inside(end + 1) = sum (all (U <= repmat (v, size (U, 1), 1), 2));
    end
  end
  same ('weight-4 words of RM(2,4) with v = 0, with wt(v) = 4, with v all ones', vec (by_v), '28 112 0');
  same ('even-weight u inside the 1s of each v of weight 4', vec (unique (inside)), '8');
  same ('x1x2 at m = 3, then at m = 4', [word(truth_table ([1 2], 3)) ' ' word(truth_table ([1 2], 4))], '00010001 0001000100010001');
  Q = nchoosek (0:15, 4);
  xs = bitxor (bitxor (Q(:, 1), Q(:, 2)), bitxor (Q(:, 3), Q(:, 4)));
  same ('sets of 4 of the 16 points whose numbers XOR to 0', sum (xs == 0), 140);
  C = rm_code (2, 4);
  W4 = C(sum (C, 2) == 4, :);
  xo = zeros (size (W4, 1), 1);
  for i = 1:size (W4, 1)
    q = find (W4(i, :)) - 1;
    xo(i) = bitxor (bitxor (q(1), q(2)), bitxor (q(3), q(4)));
  end
  require ('the 140 weight-4 words of RM(2,4) are exactly those sets', size (W4, 1) == 140 && all (xo == 0));
  same ('dimension of the dual of RM(1,5)', 32 - k_rm (1, 5), 26);
end

function check_reed ()
  section ('Reed''s majority-logic decoder');
  ok = true;
  cases = 0;
  for m = 1:5
    for s = 0:m
      C = combos (m, s);
      M = monomials (s, m);
      for q = 1:size (C, 1)
        S = C(q, :);
        for b = 0:(2 ^ (m - s) - 1)
          Q = subcube (S, b, m);
          for i = 1:numel (M)
            T = M{i};
            tt = truth_table (T, m);
            ok = ok && ((mod (sum (tt(Q + 1)), 2) == 1) == isequal (T, S));
            cases = cases + 1;
          end
        end
      end
    end
  end
  require ('subcube sum of x_T is 1 exactly when T = S, |T| <= |S|', ok);
  same ('cases checked for m = 1 to 5', cases, 3880);
  same ('check sums per degree-1 monomial of RM(1,3)', 2 ^ 2, 4);
  same ('positions in the x1 check sums of RM(1,3)', positions_text (1, 3), '01 23 45 67');
  same ('positions in the x2 check sums', positions_text (2, 3), '02 13 46 57');
  same ('positions in the x3 check sums', positions_text (3, 3), '04 15 26 37');
  y = encode ([1 0 1 1], 1, 3);
  y(6) = 1 - y(6);
  same ('received word, error in position 5', word (y), '11000111');
  [coeffs, trace] = reed_decode (y, 1, 3);
  same ('check sums for x1', trace{1}, '0010');
  same ('check sums for x2', trace{2}, '1110');
  same ('check sums for x3', trace{3}, '1011');
  same ('word after removing x2 + x3', trace{4}, '11111011');
  same ('votes for a0, the word itself', trace{5}, '11111011');
  same ('decoded (a0, a1, a2, a3)', vec (coeffs), '1 0 1 1');
  % RM(2,4): f = 1 + x2 + x4 + x1x2 + x2x3 + x3x4, one error in position 10
  y = encode ([1 0 1 0 1 1 0 0 1 0 1], 2, 4);
  same ('RM(2,4) codeword of 1 + x2 + x4 + x1x2 + x2x3 + x3x4', word (y), '1101111000101110');
  y(11) = 1 - y(11);
  same ('RM(2,4) received, error in position 10', word (y), '1101111000001110');
  [coeffs, trace] = reed_decode (y, 2, 4);
  same ('RM(2,4) check sums, x1x2 x1x3 x1x4 x2x3 x2x4 x3x4', strjoin (trace(1:6), ' '), ...
        '1101 0001 0100 1101 1000 1101');
  same ('RM(2,4) degree-2 part, x1x2 + x2x3 + x3x4', ...
        word (mod (truth_table ([1 2], 4) + truth_table ([2 3], 4) + truth_table ([3 4], 4), 2)), '0001001000011101');
  same ('RM(2,4) word after removing the degree-2 part', trace{7}, '1100110000010011');
  same ('RM(2,4) check sums, x1 x2 x3 x4', strjoin (trace(8:11), ' '), ...
        '00000100 11110111 00000010 11011111');
  same ('RM(2,4) word after removing x2 + x4', trace{12}, '1111111111011111');
  same ('RM(2,4) left after removing a0 = 1, the error', trace{14}, '0000000000100000');
  same ('RM(2,4) decoded coefficients in monomial order', vec (coeffs), '1 0 1 0 1 1 0 0 1 0 1');
  W = all_words (8);
  [D, tie] = reed_batch (W, reed_plan (1, 3));
  agree = true;
  for x = 1:256
    c = reed_decode (W(x, :), 1, 3);
    if (isempty (c))
      agree = agree && tie(x);
    else
      agree = agree && ~tie(x) && isequal (double (D(x, :)), encode (c, 1, 3));
    end
  end
  require ('bit-mask decoder matches the step-by-step one, all 256 words', agree);
  cases = [1 3 9 144 0; 1 4 697 22304 0; 2 4 17 34816 0];
  for i = 1:size (cases, 1)
    r = cases(i, 1);
    m = cases(i, 2);
    [npat, tried, failures] = exhaustive_reed (r, m, true);
    same (sprintf ('RM(%d,%d): error patterns up to the radius', r, m), npat, cases(i, 3));
    same (sprintf ('RM(%d,%d): decodings, every codeword times every pattern', r, m), tried, cases(i, 4));
    same (sprintf ('RM(%d,%d): decodings that went wrong', r, m), failures, cases(i, 5));
  end
  same ('weight-7 patterns of length 32', nchoosek (32, 7), 3365856);
  [npat, tried, failures] = exhaustive_reed (1, 5, false);
  same ('RM(1,5): patterns of weight 0 to 7, all decoded', tried, 4514873);
  same ('RM(1,5): decodings that went wrong', failures, 0);
  plan = reed_plan (1, 5);
  words = logical (rm_code (1, 5));
  E = false (960, 32);
  row = 0;
  for step = [3 5 7]
    for w = 0:9
      for s = 0:31
        row = row + 1;
        E(row, mod (s + step * (0:(w - 1)), 32) + 1) = true;
      end
    end
  end
  E = unique (E, 'rows');
  ne = size (E, 1);
  same ('distinct spread patterns, strides 3 5 7, weights 0 to 9, starts 0 to 31', ne, 801);
  [D0, tie0] = reed_batch (E, plan);
  wrong0 = ~tie0 & any (D0, 2);
  same ('on the all-zero word: ties, wrong codewords, successes', vec ([sum(tie0), sum(wrong0), sum(~tie0 & ~wrong0)]), '96 96 609');
  nw = size (words, 1);
  Crep = words(kron ((1:nw)', ones (ne, 1)), :);
  Y = xor (Crep, repmat (E, nw, 1));
  [D, tie] = reed_batch (Y, plan);
  T0 = repmat (tie0, nw, 1);
  D0r = repmat (D0, nw, 1);
  ok = isequal (tie, T0) && isequal (D(~tie, :), xor (Crep(~tie, :), D0r(~tie, :)));
  require ('D(c + e) = c + D(e) on all 64 codewords, a tie on one side exactly when on the other', ok);
  z = zeros (1, 8);
  z(1) = 1;
  z(2) = 1;
  [coeffs, trace] = reed_decode (z, 1, 3);
  same ('errors at 0 and 1: check sums for x1', trace{1}, '0000');
  same ('errors at 0 and 1: check sums for x2', trace{2}, '1100');
  require ('errors at 0 and 1: the x2 vote ties and the decoder stops', isempty (coeffs) && numel (trace) == 2);
  same ('errors at 0 and 1: check sums for x3, had it gone on', word (xor (z(1:4), z(5:8))), '1100');
  C13 = rm_code (1, 3);
  same ('errors at 0 and 1: codewords at distance 2', sum (sum (xor (C13, repmat (z, size (C13, 1), 1)), 2) == 2), 4);
  P2 = nchoosek (1:8, 2);
  near = zeros (28, 1);
  held = zeros (28, 1);
  for i = 1:28
    e = zeros (1, 8);
    e(P2(i, :)) = 1;
    near(i) = sum (sum (xor (C13, repmat (e, size (C13, 1), 1)), 2) == 2);
    held(i) = sum (sum (C13, 2) == 4 & all (C13(:, P2(i, :)), 2));
  end
  same ('RM(1,3): codewords at distance 2 from each of the 28 double errors', vec (unique (near)), '4');
  same ('RM(1,3): weight-4 codewords holding each pair of positions', vec (unique (held)), '3');
  stz = fht_stages (signs (z));
  same ('errors at 0 and 1: correlations', vec (stz(end, :)), '4 0 -4 0 -4 0 -4 0');
end

function check_fht ()
  section ('The fast Hadamard transform');
  H2 = hadamard (1);
  H4 = hadamard (2);
  same ('H4, row by row', strjoin ({vec(H4(1, :)), vec(H4(2, :)), vec(H4(3, :)), vec(H4(4, :))}, ' | '), '1 1 1 1 | 1 -1 1 -1 | 1 1 -1 -1 | 1 -1 -1 1');
  require ('H4 = H2 (x) H2', isequal (kron (H2, H2), H4));
  require ('H8 = H2 (x) H4', isequal (kron (H2, H4), hadamard (3)));
  require ('H32 = H2 (x) H16', isequal (kron (H2, hadamard (4)), hadamard (5)));
  y = encode ([1 0 1 1], 1, 3);
  y(6) = 1 - y(6);
  Y = signs (y);
  same ('received word as signs', vec (Y), '-1 -1 1 1 1 -1 -1 -1');
  st = fht_stages (Y);
  same ('after stage 1, pairs (j, j + 4)', vec (st(1, :)), '0 -2 0 0 -2 0 2 2');
  same ('after stage 2, pairs (j, j + 2)', vec (st(2, :)), '0 -2 0 -2 0 2 -4 -2');
  same ('after stage 3, pairs (j, j + 1), the transform T', vec (st(3, :)), '-2 2 -2 2 2 -2 -6 -2');
  require ('T equals the direct correlation', isequal (st(3, :), direct_correlation (Y, 3)));
  [best, a0, T] = fht_decode (y, 3);
  same ('index of the largest |T|', best, 6);
  same ('its value', T(best + 1), -6);
  same ('decoded a0 from the sign', a0, 1);
  same ('agreements with 1 + x2 + x3, (8 + 6) / 2', (8 + abs (T(best + 1))) / 2, 7);
  W = all_words (8);
  ok = true;
  ok2 = true;
  for x = 1:256
    Yx = signs (W(x, :));
    stx = fht_stages (Yx);
    Tx = stx(3, :);
    ok = ok && isequal (Tx, direct_correlation (Yx, 3));
    for a = 0:7
      agreements = sum (W(x, :) == rm1_codeword (0, a, 3));
      ok2 = ok2 && (8 + Tx(a + 1) == 2 * agreements);
    end
  end
  require ('FHT equals direct correlation for all 256 words of length 8', ok);
  require ('agreements = (2^m + T) / 2 for every word and every a', ok2);
  cases = [2 8 12; 3 24 56; 4 64 240; 5 160 992; 10 10240 1047552];
  for i = 1:size (cases, 1)
    m = cases(i, 1);
    same (sprintf ('m = %d: additions, fast and direct', m), sprintf ('%d %d', fht_additions (m), direct_additions (m)), sprintf ('%d %d', cases(i, 2), cases(i, 3)));
  end
  same ('m = 10: ratio, direct over fast', sprintf ('%.1f', direct_additions (10) / fht_additions (10)), '102.3');
  same ('butterflies for m = 3', sum (arrayfun (@(h) sum (bitand (0:7, h) == 0), [4 2 1])), 12);
  c = rm1_codeword (1, 13, 5);
  same ('codeword with a0 = 1, (a1..a5) = 1 0 1 1 0', word (c), '10100101010110101010010101011010');
  errs = [2 7 11 16 21 26 30];
  y = c;
  y(errs + 1) = 1 - y(errs + 1);
  same ('received word with errors at 2 7 11 16 21 26 30', word (y), '10000100010010100010000101111000');
  [best, a0, T] = fht_decode (y, 5);
  same ('T over all 32 linear functions', vec (T), '10 -2 -2 2 -2 2 2 6 6 2 2 -2 2 -18 -2 -6 2 -2 -10 2 6 2 -6 -10 -2 2 -6 -2 -6 -2 6 -6');
  same ('index of the largest |T|', best, 13);
  same ('T there', T(14), -18);
  others = abs (T);
  others(14) = [];
  same ('largest |T| anywhere else', max (others), 10);
  same ('indices where |T| = 10', vec (find (abs (T) == 10) - 1), '0 18 23');
  same ('correct |T| with 7 errors, 32 - 14', 32 - 2 * 7, 18);
  same ('bound on every other |T|, 2 * 7', 2 * 7, 14);
  same ('decoded a0', a0, 1);
  c = rm1_codeword (0, 11, 4);
  y = c;
  y([3 8 14] + 1) = 1 - y([3 8 14] + 1);
  same ('RM(1,4) codeword x1 + x2 + x4', word (c), '0110011010011001');
  same ('received with errors at 3, 8, 14', word (y), '0111011000011011');
  st = fht_stages (signs (y));
  same ('RM(1,4) after stage 1', vec (st(1, :)), '2 0 0 -2 0 0 -2 0 0 -2 -2 0 2 -2 0 2');
  same ('RM(1,4) after stage 2', vec (st(2, :)), '2 0 -2 -2 2 0 2 -2 2 -4 -2 2 -2 0 -2 -2');
  same ('RM(1,4) after stage 3', vec (st(3, :)), '0 -2 4 2 4 -2 0 2 0 -2 4 -6 -4 -2 0 2');
  same ('RM(1,4) after stage 4, T', vec (st(4, :)), '-2 2 6 2 2 6 2 -2 -2 2 -2 10 -6 -2 2 -2');
  [best, a0, T] = fht_decode (y, 4);
  same ('RM(1,4): index of the largest |T|', best, 11);
  same ('RM(1,4): T there', T(12), 10);
  same ('RM(1,4): received signs', vec (signs (y)), '1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 -1');
  others = abs (T);
  others(12) = [];
  same ('RM(1,4): largest |T| anywhere else', max (others), 6);
  same ('RM(1,4): correct |T| with 3 errors, 16 - 6', 16 - 2 * 3, 10);
end

function check_mariner ()
  section ('The Mariner code');
  same ('(n, k, d) of RM(1,5)', vec ([32, k_rm(1, 5), 16]), '32 6 16');
  same ('rate 6/32', sprintf ('%.4f', 6 / 32), '0.1875');
  same ('word error, RM(1,5) to 7 errors, p = 0.05', sprintf ('%.6f', word_error (32, 7, 0.05)), '0.000139');
  same ('word error, 6 bits uncoded, p = 0.05', sprintf ('%.4f', 1 - 0.95 ^ 6), '0.2649');
end

function check_hardware ()
  section ('Hardware');
  cases = [1 3 12; 1 5 80; 2 4 56];
  for i = 1:size (cases, 1)
    r = cases(i, 1);
    m = cases(i, 2);
    G = generator (r, m);
    gates = sum (max (0, sum (G, 1) - 1));
    same (sprintf ('two-input XOR gates in the parallel RM(%d,%d) encoder', r, m), gates, cases(i, 3));
  end
  net = serial_encoder_rm13 ();
  ok = true;
  for msg = 0:15
    a = bitand (bitshift (msg, -(0:3)), 1);
    out = zeros (1, 8);
    for j = 0:7
      sig = struct ('a0', a(1), 'a1', a(2), 'a2', a(3), 'a3', a(4), ...
                    'x1', bitand (j, 1), 'x2', bitand (bitshift (j, -1), 1), 'x3', bitand (bitshift (j, -2), 1));
      sig = run_netlist (net, sig);
      out(j + 1) = sig.c;
    end
    ok = ok && isequal (out, encode (a, 1, 3));
  end
  require ('serial RM(1,3) encoder netlist gives every codeword, 16 messages times 8 clocks', ok);
  same ('serial RM(1,3) encoder: AND gates, XOR gates', sprintf ('%d %d', sum (strcmp (net(:, 1), 'AND')), sum (strcmp (net(:, 1), 'XOR'))), '3 3');
  cases = [1 3 12; 1 5 80; 2 4 72];
  for i = 1:size (cases, 1)
    r = cases(i, 1);
    m = cases(i, 2);
    C = combos (m, r);
    top = 0;
    for q = 1:size (C, 1)
      for b = 0:(2 ^ (m - r) - 1)
        top = top + numel (subcube (C(q, :), b, m)) - 1;
      end
    end
    same (sprintf ('two-input XOR gates in the top-degree check sums, RM(%d,%d)', r, m), top, cases(i, 3));
  end
  cases = {3, '3 4 1 8'; 5, '5 16 1 32'};
  for i = 1:size (cases, 1)
    m = cases{i, 1};
    plan = reed_plan (1, m);
    same (sprintf ('RM(1,%d): majority gates and inputs at degree 1, then at degree 0', m), ...
          sprintf ('%d %d %d %d', numel (plan{1}), plan{1}{1}.nsums, numel (plan{2}), plan{2}{1}.nsums), cases{i, 2});
  end
  for mc = [3 12; 5 80]'
    m = mc(1);
    same (sprintf ('m = %d: two-input XOR gates that remove the degree-1 part', m), sum (popcount (0:(2 ^ m - 1))), mc(2));
  end
  % the butterfly array of fixed-width two's-complement adders, checked against the exact transform
  [s00, c00] = ripple (1, 3, 2, false);
  same ('1 + (-1) in 2 bits: sum bits, carry, and the 3-bit word read as a number', ...
        sprintf ('%s %d %d', dec2bin (s00, 2), c00, from_bits (s00 + c00 * 4, 3)), '00 1 -4');
  c = rm1_codeword (1, 13, 5);
  errs = [2 7 11 16 21 26 30];
  c(errs + 1) = 1 - c(errs + 1);
  tests = {signs(all_words (8)), [signs(rm_code (1, 5)); signs(c)]};
  ms = [3 5];
  stated = {'2 3 4 5 24 96', '2 3 4 5 6 7 160 800'};
  nwrong = [256 65];
  for i = 1:2
    m = ms(i);
    Y = tests{i};
    H = hadamard (m);
    [T, widths, units, cells] = butterfly_array (Y);
    require (sprintf ('m = %d: fixed-width butterflies equal the exact transform on all %d test words', m, size (Y, 1)), isequal (T, Y * H'));
    same (sprintf ('m = %d: widths after the sign map and each stage, adders and subtracters, full-adder cells', m), vec ([widths units cells]), stated{i});
    Tc = butterfly_array (Y, true);
    same (sprintf ('m = %d: test words wrong when the carry is taken as the new bit', m), sum (any (Tc ~= Y * H', 2)), nwrong(i));
  end
  same ('m = 5: butterflies in all stages', sum (arrayfun (@(h) sum (bitand (0:31, h) == 0), [16 8 4 2 1])), 80);
  T = butterfly_array (signs (c));
  [best, a0, comparators, levels] = max_finder (T);
  same ('m = 5: max-finder on the seven-error word: index, a0, comparators, levels', vec ([best a0 comparators levels]), '13 1 31 5');
  [b2, a2] = fht_decode (c, 5);
  require ('the max-finder agrees with the decoder', best == b2 && a0 == a2);
  [~, ~, comparators, levels] = max_finder (tests{1}(1, :));
  same ('m = 3: comparators and levels in the max-finder', vec ([comparators levels]), '7 3');
  stated = {'8 4', '32 6'};
  for i = 1:2
    T = butterfly_array (tests{i});
    top = max (abs (T(:)));
    same (sprintf ('largest |T| for m = %d, and the bits in its magnitude', ms(i)), sprintf ('%d %d', top, numel (dec2bin (top))), stated{i});
  end
end

function check_exercises ()
  section ('Exercises');
  same ('Exercise: (n, k, d) of RM(2,5), d from all its words', vec ([32, k_rm(2, 5), min_weight(rm_code (2, 5))]), '32 16 8');
  same ('Exercise: (n, k, d) of RM(3,6), d from the formula', vec ([64, k_rm(3, 6), 2 ^ (6 - 3)]), '64 42 8');
  f = mod (truth_table ([], 3) + truth_table (2, 3) + truth_table ([1 3], 3), 2);
  same ('Exercise: truth table of 1 + x2 + x1x3', word (f), '11001001');
  maj = double (popcount (0:7) >= 2);
  same ('Exercise: truth table of the majority of 3 bits', word (maj), '00010111');
  st = mobius_stages (maj, 3);
  same ('Exercise: its polynomial', anf_terms (st(3, :), 3), 'x1x2 + x1x3 + x2x3');
  plan = reed_plan (2, 5);
  stated = {'10 8 4', '5 16 2', '1 32 1'};
  for L = 1:3
    st = plan{L}{1};
    same (sprintf ('Exercise: RM(2,5) degree %d: coefficients, check sums each, points per sum', 3 - L), ...
          sprintf ('%d %d %d', numel (plan{L}), st.nsums, size (st.Q, 2)), stated{L});
  end
  same ('Exercise: RM(2,5) radius, half the top-stage votes less 1', plan{1}{1}.nsums / 2 - 1, 3);
  same ('Exercise: words of weight 8 in RM(1,4)', weight_count (rm_code (1, 4), 8), 30);
  B = macwilliams (rm_code (1, 5), 32);
  same ('Exercise: (n, k, d) of the dual of RM(1,5), d by MacWilliams', vec ([32, k_rm(3, 5), find(B(2:end), 1)]), '32 26 4');
  same ('Exercise: (n, k, d) of RM(1,10), d from the formula', vec ([2 ^ 10, k_rm(1, 10), 2 ^ 9]), '1024 11 512');
  same ('Exercise: RM(1,10) radius, from the formula', 2 ^ 8 - 1, 255);
end


% Main -------------------------------------------------------------------------

check_truth_tables ();
check_mobius ();
check_dimensions ();
check_recursion ();
check_distance ();
check_first_order ();
check_duality ();
check_reed ();
check_fht ();
check_mariner ();
check_hardware ();
check_exercises ();
printf ('\nAll chapter 26 numbers reproduced.\n');
