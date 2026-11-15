% Rebuild the codes of chapter 21 of Error Correction Coding, "Beyond
% SECDED", and recompute its numbers, in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch21_beyond_secded.m
%
% It needs only core Octave 7 or later, with no pkg load and no Forge
% packages. It is the Octave twin of python/ch21_beyond_secded.py and prints
% the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number is a whole number, a bit pattern or a percentage,
% so every comparison is exact.
%
% An element of GF(2^m) is held as a whole number whose bit j is the
% coefficient of alpha^j, printed with the highest coefficient first. A
% column of a binary check matrix with r rows is a whole number whose highest
% bit is row 0, so the printed string is the column read from top to bottom.
% Positions are counted from 0, as in the chapter, so position p is element
% p + 1 here. The search for the adjacent-error codes runs depth first over
% the candidate columns in increasing order, exactly as the Python script
% does, so both routes find the same codes.

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

function w = popcount (v)
  % The number of ones in each entry of v.
  w = zeros (size (v));
  while (any (v(:)))
    w = w + bitand (v, 1);
    v = bitshift (v, -1);
  end
end

function s = bits (v, width)
  s = dec2bin (v, width);
end

function s = bits_list (v, width)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = bits (v(i), width);
  end
  s = strjoin (parts, ' ');
end

function s = num_list (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ' ');
end

function s = pct (part, whole)
  % A share in percent to two decimals, trailing zeros dropped.
  s = sprintf ('%.2f', 100 * part / whole);
  s = regexprep (s, '0+$', '');
  s = regexprep (s, '\.$', '');
end

function s = pct1 (part, whole)
  s = sprintf ('%.1f', 100 * part / whole);
end

function h = top_bit (x)
  [~, e] = log2 (x);
  h = e - 1;
end


% Binary polynomials and finite fields -----------------------------------------

function r = poly_rem (a, b)
  db = top_bit (b);
  while (a ~= 0 && top_bit (a) >= db)
    a = bitxor (a, bitshift (b, top_bit (a) - db));
  end
  r = a;
end

function ok = irreducible (p)
  % No factor of degree 1 to deg(p) / 2, found by trial division.
  m = top_bit (p);
  ok = true;
  for d = 2:(2 ^ (floor (m / 2) + 1) - 1)
    if (poly_rem (p, d) == 0)
      ok = false;
      return;
    end
  end
end

function r = poly_mulmod (a, b, poly, m)
  % Shift-and-add product of two field elements, reduced by poly.
  r = 0;
  while (b ~= 0)
    if (bitand (b, 1))
      r = bitxor (r, a);
    end
    b = bitshift (b, -1);
    a = bitshift (a, 1);
    if (bitand (bitshift (a, -m), 1))
      a = bitxor (a, poly);
    end
  end
end

function [ex, lg, wrap] = power_table (poly, m)
  q1 = 2 ^ m - 1;
  ex = zeros (1, q1);
  lg = -ones (1, 2 ^ m);
  v = 1;
  for i = 0:q1 - 1
    ex(i + 1) = v;
    lg(v + 1) = i;
    v = poly_mulmod (v, 2, poly, m);
  end
  wrap = v;
end

function c = gmul (a, b, ex, lg)
  if (a == 0 || b == 0)
    c = 0;
  else
    c = ex(mod (lg(a + 1) + lg(b + 1), numel (ex)) + 1);
  end
end

function c = gdiv (a, b, ex, lg)
  if (a == 0)
    c = 0;
  else
    c = ex(mod (lg(a + 1) - lg(b + 1), numel (ex)) + 1);
  end
end

function rank = rank_gf2 (cols)
  basis = zeros (1, 64);
  rank = 0;
  for c = cols
    x = c;
    while (x ~= 0)
      h = top_bit (x);
      if (basis(h + 1) ~= 0)
        x = bitxor (x, basis(h + 1));
      else
        basis(h + 1) = x;
        rank = rank + 1;
        break;
      end
    end
  end
end

function cols = bch_columns (m, ex, n)
  % Columns (beta, beta^3) for beta = alpha^0, ..., alpha^(n-1), beta in the top m rows.
  q1 = numel (ex);
  i = 0:n - 1;
  cols = bitor (bitshift (ex(mod (i, q1) + 1), m), ex(mod (3 * i, q1) + 1));
end

function b = ball (n, t)
  b = 0;
  for i = 0:t
    b = b + nchoosek (n, i);
  end
end

function r = least_check_bits (k, t, extended)
  r = 1;
  while (true)
    shift = double (extended);
    if (2 ^ (r - shift) >= ball (k + r - shift, t))
      return;
    end
    r = r + 1;
  end
end

function s = pair_sums (cols)
  P = nchoosek (1:numel (cols), 2);
  s = bitxor (cols(P(:, 1)), cols(P(:, 2)));
end

function n = distinct_low_syndromes (cols)
  % Distinct syndromes of the patterns of weight at most 2.
  n = numel (unique ([0, cols(:)', pair_sums(cols(:)')(:)']));
end

function T = mul_table (poly, m)
  q = 2 ^ m;
  T = zeros (q, q);
  for a = 0:q - 1
    for b = 0:q - 1
      T(a + 1, b + 1) = poly_mulmod (a, b, poly, m);
    end
  end
end


function cols = hamming_columns ()
  % The (71,64) Hamming code: the first 71 nonzero 7-bit words, all distinct.
  cols = 1:71;
end

function cols = hsiao_columns ()
  % A (72,64) Hsiao code: the 8 unit words, all 56 words of weight 3 and the
  % first 8 words of weight 5, all of odd weight and all distinct.
  v = 0:255;
  w = popcount (v);
  five = v(w == 5);
  cols = [v(w == 1), v(w == 3), five(1:8)];
end

function check_baselines ()
  section ('The Hamming SEC and Hsiao SEC-DED codes for 64 data bits, built');
  ham = hamming_columns ();
  require ('Hamming (71,64): 71 distinct nonzero columns', numel (unique (ham)) == 71 && all (ham ~= 0));
  same ('Hamming (71,64): rank of the 7 x 71 check matrix', rank_gf2 (ham), 7);
  hs = hsiao_columns ();
  require ('Hsiao (72,64): 72 distinct columns of odd weight', numel (unique (hs)) == 72 && all (mod (popcount (hs), 2) == 1));
  same ('Hsiao (72,64): rank of the 8 x 72 check matrix', rank_gf2 (hs), 8);
  d = pair_sums (hs);
  require ('Hsiao (72,64): every double has a nonzero syndrome no single has', all (d ~= 0) && ~any (ismember (d, hs)));
end


% The opening problem ----------------------------------------------------------

function check_opening ()
  section ('One strike, two neighboring bits of a (72,64) word');
  same ('double errors in a 72-bit word', nchoosek (72, 2), 2556);
  P = nchoosek (0:71, 2);
  same ('adjacent pairs in a 72-bit word', sum (P(:, 2) - P(:, 1) == 1), 71);
  same ('check bits of the (72,64) SEC-DED code', rank_gf2 (hsiao_columns ()), 8);
end


% DEC-TED ----------------------------------------------------------------------

function check_sphere_counts ()
  section ('Sphere-packing minimum for 64 data bits');
  same ('ball of radius 1 in 71 bits', ball (71, 1), 72);
  same ('ball of radius 1 in 70 bits', ball (70, 1), 71);
  same ('least check bits for SEC', least_check_bits (64, 1, false), 7);
  same ('least check bits for SEC-DED', least_check_bits (64, 1, true), 8);
  same ('ball of radius 2 in 75 bits', ball (75, 2), 2851);
  same ('ball of radius 2 in 76 bits', ball (76, 2), 2927);
  require ('2851 exceeds 2^11 = 2048, so 11 check bits fail', ball (75, 2) > 2 ^ 11);
  require ('2927 fits in 2^12 = 4096, so 12 check bits pass', ball (76, 2) <= 2 ^ 12);
  same ('least check bits for DEC', least_check_bits (64, 2, false), 12);
  same ('least check bits for DEC-TED', least_check_bits (64, 2, true), 13);
  same ('least check bits for DEC of 32 data bits', least_check_bits (32, 2, false), 10);
  same ('ball of radius 2 in 41 bits', ball (41, 2), 862);
  same ('ball of radius 2 in 42 bits', ball (42, 2), 904);
end

function check_gf16 (ex, lg, wrap, T)
  section ('GF(16) built from x^4 + x + 1');
  stated = {'0001', '0010', '0100', '1000', '0011', '0110', '1100', '1011', ...
            '0101', '1010', '0111', '1110', '1111', '1101', '1001'};
  same ('x^4 + x + 1 has no factor of degree 1 or 2', irreducible (19), true);
  for i = 0:14
    same (sprintf ('alpha^%d as c3 c2 c1 c0', i), bits (ex(i + 1), 4), stated{i + 1});
  end
  same ('alpha^15', bits (wrap, 4), '0001');
  same ('distinct nonzero powers of alpha', numel (unique (ex)), 15);
  ok = true;
  for a = 0:15
    for b = 0:15
      ok = ok && gmul (a, b, ex, lg) == T(a + 1, b + 1);
    end
  end
  require ('table product equals shift-and-add product, 256 pairs', ok);
  require ('every nonzero element has an inverse', all (any (T(2:16, :) == 1, 2)));
  require ('no two nonzero elements multiply to zero', all (all (T(2:16, 2:16) ~= 0)));
  [A, B, C] = ndgrid (0:15, 0:15, 0:15);
  AB = T(sub2ind ([16 16], A + 1, B + 1));
  BC = T(sub2ind ([16 16], B + 1, C + 1));
  left = T(sub2ind ([16 16], AB + 1, C + 1));
  right = T(sub2ind ([16 16], A + 1, BC + 1));
  require ('multiplication is associative on all 4,096 triples', all (left(:) == right(:)));
  AC = T(sub2ind ([16 16], A + 1, C + 1));
  left = T(sub2ind ([16 16], A + 1, bitxor (B, C) + 1));
  require ('multiplication distributes over addition, 4,096 triples', all (left(:) == bitxor (AB(:), AC(:))));
  same ('alpha^6 times alpha^13, as a power of alpha', lg(gmul (ex(7), ex(14), ex, lg) + 1), 4);
  same ('alpha^6 times alpha^13, shift and add', bits (poly_mulmod (ex(7), ex(14), 19, 4), 4), '0011');
  same ('inverse of alpha^5, as a power of alpha', lg(gdiv (1, ex(6), ex, lg) + 1), 10);
end

function check_dec_identities (T)
  section ('The algebra behind the DEC proof, checked in GF(16)');
  m = @(a, b) T(sub2ind ([16 16], a + 1, b + 1));
  cube = @(a) m (a, m (a, a));
  [A, B] = ndgrid (0:15, 0:15);
  A = A(:);
  B = B(:);
  S = bitxor (A, B);
  lhs = cube (S);
  rhs = bitxor (bitxor (cube (A), m (m (A, A), B)), bitxor (m (A, m (B, B)), cube (B)));
  require ('(a + b)^3 = a^3 + a^2 b + a b^2 + b^3, all 256 pairs', all (lhs == rhs));
  lhs = bitxor (cube (A), cube (B));
  rhs = m (S, bitxor (m (S, S), m (A, B)));
  require ('a^3 + b^3 = s(s^2 + ab) with s = a + b, all 256 pairs', all (lhs == rhs));
  lhs = bitxor (bitxor (cube (A), cube (B)), cube (S));
  rhs = m (m (A, B), S);
  require ('a + b + c = 0 gives a^3 + b^3 + c^3 = abc, all 256 pairs', all (lhs == rhs));
  ok = true;
  X = (0:15)';
  for s = 0:15
    for p = 0:15
      vals = bitxor (bitxor (m (X, X), m (s * ones (16, 1), X)), p);
      ok = ok && sum (vals == 0) <= 2;
    end
  end
  require ('x^2 + s x + p has at most two roots, all 256 pairs (s, p)', ok);
end

function dist = codeword_weights (cols)
  % Weight distribution of the code whose check matrix has these columns.
  n = numel (cols);
  V = 0:2 ^ n - 1;
  syn = zeros (size (V));
  for j = 1:n
    syn = bitxor (syn, bitand (bitshift (V, -(j - 1)), 1) * cols(j));
  end
  w = popcount (V(syn == 0));
  dist = accumarray (w(:) + 1, 1, [n + 1, 1])';
end

function check_bch_15_7 (ex, lg)
  section ('The (15,7) double-error-correcting code');
  cols = bch_columns (4, ex, 15);
  same ('column 0, read top to bottom', bits (cols(1), 8), '00010001');
  same ('column 1, read top to bottom', bits (cols(2), 8), '00101000');
  same ('column 6, read top to bottom', bits (cols(7), 8), '11001000');
  same ('column 14, read top to bottom', bits (cols(15), 8), '10011111');
  same ('distinct bottom halves among the 15 columns', numel (unique (bitand (cols, 15))), 5);
  same ('rank of the 8 x 15 check matrix', rank_gf2 (cols), 8);
  same ('data bits k = 15 - 8', 15 - rank_gf2 (cols), 7);
  dist = codeword_weights (cols);
  same ('codewords', sum (dist), 128);
  nz = find (dist(2:end)) ;
  same ('minimum nonzero weight', nz(1), 5);
  same ('weight distribution A0 to A15', num_list (dist), '1 0 0 0 0 18 30 15 15 30 18 0 0 0 0 1');
  same ('distinct syndromes of patterns of weight at most 2', distinct_low_syndromes (cols), 121);
  same ('patterns of weight at most 2 in 15 bits', ball (15, 2), 121);

  section ('Locating two errors at positions 3 and 10');
  s = bitxor (cols(4), cols(11));
  s1 = bitshift (s, -4);
  s3 = bitand (s, 15);
  same ('syndrome, top half S1', bits (s1, 4), '1111');
  same ('syndrome, bottom half S3', bits (s3, 4), '1011');
  same ('S1 as a power of alpha', lg(s1 + 1), 12);
  same ('S3 as a power of alpha', lg(s3 + 1), 7);
  same ('S3 / S1 as a power of alpha', lg(gdiv (s3, s1, ex, lg) + 1), 10);
  same ('S1^2 as a power of alpha', lg(gmul (s1, s1, ex, lg) + 1), 9);
  p = bitxor (gdiv (s3, s1, ex, lg), gmul (s1, s1, ex, lg));
  same ('ab = S3 / S1 + S1^2, as bits', bits (p, 4), '1101');
  same ('ab as a power of alpha', lg(p + 1), 13);
  roots_found = [];
  for i = 0:14
    x = ex(i + 1);
    if (bitxor (bitxor (gmul (x, x, ex, lg), gmul (s1, x, ex, lg)), p) == 0)
      roots_found(end + 1) = i;
    end
  end
  same ('positions whose alpha^i solves x^2 + S1 x + ab', num_list (roots_found), '3 10');
  % The division-free form, multiplied through by S1, that each position tests at once.
  s1sq = gmul (s1, s1, ex, lg);
  c0 = bitxor (s3, gmul (s1sq, s1, ex, lg));
  free_found = [];
  for i = 0:14
    x = ex(i + 1);
    if (bitxor (bitxor (gmul (s1, gmul (x, x, ex, lg), ex, lg), gmul (s1sq, x, ex, lg)), c0) == 0)
      free_found(end + 1) = i;
    end
  end
  same ('positions whose alpha^i solves S1 x^2 + S1^2 x + S3 + S1^3', num_list (free_found), '3 10');
end

function ext = dec_ted_columns ()
  % The (79,64) DEC-TED code: the first 78 columns (beta, beta^3) over GF(128),
  % each with a parity bit 1 appended, then the unit column of that parity bit.
  [ex, lg, wrap] = power_table (131, 7);
  ext = [bitor(bitshift (bch_columns (7, ex, 78), 1), 1), 1];
end

function check_bch_79_64 ()
  section ('A DEC-TED code for 64 data bits from GF(128)');
  [ex, lg, wrap] = power_table (131, 7);
  same ('x^7 + x + 1 has no factor of degree 1 to 3', irreducible (131), true);
  same ('distinct nonzero powers of alpha in GF(128)', numel (unique (ex)), 127);
  same ('alpha^127', wrap, 1);
  cols = bch_columns (7, ex, 127);
  same ('rank of the 14 x 127 check matrix', rank_gf2 (cols), 14);
  same ('data bits of the (127,k) code', 127 - rank_gf2 (cols), 113);
  short = cols(1:78);
  same ('rank of the first 78 columns', rank_gf2 (short), 14);
  same ('data bits of the shortened code', 78 - rank_gf2 (short), 64);
  same ('distinct syndromes, weight at most 2, length 78', distinct_low_syndromes (short), 3082);
  ext = [bitor(bitshift (short, 1), 1), 1];
  require ('every column of the extended matrix has parity bit 1', all (bitand (ext, 1)));
  same ('distinct syndromes, weight at most 2, length 79', distinct_low_syndromes (ext), 3161);
  same ('check bits of the DEC-TED code, the rank of its matrix', rank_gf2 (ext), 15);
  same ('the same matrix as dec_ted_columns() builds it', isequal (ext, dec_ted_columns ()), true);
  same ('overhead of the (79,64) code, percent', pct (15, 64), '23.44');
  same ('check bits above the sphere-packing minimum of 13', rank_gf2 (dec_ted_columns ()) - least_check_bits (64, 2, true), 2);
  same ('overhead of the (78,64) DEC code, percent', pct (14, 64), '21.88');
end


% DAEC -------------------------------------------------------------------------

function cols = daec_search (r, n)
  % Order n distinct odd-weight r-bit columns so that the n - 1 sums of
  % neighbors are distinct. Depth first, candidates in increasing order.
  % Each (r, n) is searched once and the result kept for later calls.
  persistent memo;
  if (isempty (memo))
    memo = struct ();
  end
  key = sprintf ('r%dn%d', r, n);
  if (isfield (memo, key))
    cols = memo.(key);
    return;
  end
  cols = daec_search_run (r, n);
  memo.(key) = cols;
end

function cols = daec_search_run (r, n)
  all_v = 1:2 ^ r - 1;
  cands = all_v(mod (popcount (all_v), 2) == 1);
  m = numel (cands);
  idx = -ones (1, n);
  used = false (1, m);
  sum_used = false (1, 2 ^ r);
  depth = 0;
  while (true)
    if (depth == n)
      cols = cands(idx + 1);
      return;
    end
    j = idx(depth + 1);
    if (j >= 0)
      used(j + 1) = false;
      if (depth > 0)
        sum_used(bitxor (cands(idx(depth) + 1), cands(j + 1)) + 1) = false;
      end
    end
    j = j + 1;
    while (j < m)
      if (~used(j + 1) && (depth == 0 || ~sum_used(bitxor (cands(idx(depth) + 1), cands(j + 1)) + 1)))
        break;
      end
      j = j + 1;
    end
    if (j == m)
      idx(depth + 1) = -1;
      depth = depth - 1;
      if (depth < 0)
        cols = [];
        return;
      end
      continue;
    end
    idx(depth + 1) = j;
    used(j + 1) = true;
    if (depth > 0)
      sum_used(bitxor (cands(idx(depth) + 1), cands(j + 1)) + 1) = true;
    end
    depth = depth + 1;
  end
end

function [ok1, ok2, mis, flag] = daec_outcomes (cols, r)
  % The decoder on every single and double error: singles corrected,
  % adjacent doubles corrected, nonadjacent doubles miscorrected and flagged.
  n = numel (cols);
  single = -ones (1, 2 ^ r);
  single(cols + 1) = 0:n - 1;
  adj_sums = bitxor (cols(1:n - 1), cols(2:n));
  adjacent = -ones (1, 2 ^ r);
  adjacent(adj_sums + 1) = 0:n - 2;
  ok1 = sum (single(cols + 1) == 0:n - 1);
  ok2 = sum (adjacent(adj_sums + 1) == 0:n - 2);
  P = nchoosek (1:n, 2);
  P = P(P(:, 2) >= P(:, 1) + 2, :);
  s = bitxor (cols(P(:, 1)), cols(P(:, 2)));
  hit = single(s + 1) >= 0 | adjacent(s + 1) >= 0;
  mis = sum (hit);
  flag = sum (~hit);
end

function b = daec_lower_bound (r, n)
  % Each nonzero even syndrome is the sum of at least (2n - 2^(r-1)) / 2 pairs.
  pairs = floor (max (0, 2 * n - 2 ^ (r - 1)) / 2);
  b = (n - 1) * max (0, pairs - 1);
end

function check_daec ()
  section ('SEC-DED-DAEC with odd-weight columns');
  v = 0:15;
  same ('odd-weight columns with 4 rows', sum (mod (popcount (v), 2) == 1), 8);
  even = v(v > 0 & mod (popcount (v), 2) == 0);
  same ('nonzero even-weight syndromes with 4 rows', numel (even), 7);
  total = 0;
  for e = even
    total = bitxor (total, e);
  end
  same ('sum of all seven nonzero even 4-bit words', bits (total, 4), '0000');
  same ('a search for an (8,4) code finds one', ~isempty (daec_search (4, 8)), false);
  c7 = daec_search (4, 7);
  same ('columns of the (7,3) code', bits_list (c7, 4), '0001 0010 0100 1000 0111 1110 1011');
  sums = bitxor (c7(1:6), c7(2:7));
  same ('neighbor sums of the (7,3) code', bits_list (sums, 4), '0011 0110 1100 1111 1001 0101');
  unused = even(~ismember (even, sums));
  same ('even syndrome no neighbor sum uses', bits_list (unused, 4), '1010');
  same ('rank of the (7,3) check matrix', rank_gf2 (c7), 4);
  parts = {};
  for i = 0:6
    for j = i + 2:6
      if (bitxor (c7(i + 1), c7(j + 1)) == 10)
        parts{end + 1} = sprintf ('%d,%d', i, j);
      end
    end
  end
  same ('nonadjacent doubles the (7,3) code flags', strjoin (parts, ' '), '0,6 1,3 2,5');
  same ('syndrome of a double error at 1 and 4', bits (bitxor (c7(2), c7(5)), 4), '0101');
  same ('syndrome of a double error at 2 and 3', bits (bitxor (c7(3), c7(4)), 4), '1100');
  parts = {};
  for i = 0:6
    for j = i + 2:6
      if (bitxor (c7(i + 1), c7(j + 1)) == 12)
        parts{end + 1} = sprintf ('%d,%d', i, j);
      end
    end
  end
  same ('nonadjacent pairs whose syndrome is also 1100', strjoin (parts, ' '), '1,5 4,6');
  same ('syndrome of a single error at 5', bits (c7(6), 4), '1110');

  stated = [4, 7, 6, 15, 12, 3, 12;
            5, 13, 12, 66, 51, 15, 48;
            7, 39, 38, 703, 490, 213, 228;
            8, 72, 71, 2485, 1465, 1020, 497];
  for row = 1:4
    r = stated(row, 1);
    n = stated(row, 2);
    cols = daec_search (r, n);
    label = sprintf ('(%d,%d)', n, n - r);
    require ([label ' search finds a code'], ~isempty (cols));
    require ([label ' every column has odd weight'], all (mod (popcount (cols), 2) == 1));
    same ([label ' rank of the check matrix'], rank_gf2 (cols), r);
    [ok1, ok2, mis, flag] = daec_outcomes (cols, r);
    same ([label ' single errors corrected'], ok1, n);
    same ([label ' adjacent doubles corrected'], ok2, stated(row, 3));
    same ([label ' nonadjacent doubles'], mis + flag, stated(row, 4));
    same ([label ' nonadjacent doubles miscorrected'], mis, stated(row, 5));
    same ([label ' nonadjacent doubles flagged'], flag, stated(row, 6));
    same ([label ' lower bound on miscorrected doubles'], daec_lower_bound (r, n), stated(row, 7));
  end
  cols72 = daec_search (8, 72);
  same ('(72,64) share of nonadjacent doubles miscorrected', pct1 (1465, 2485), '59.0');
  same ('(72,64) share at the lower bound', pct1 (497, 2485), '20.0');
  w = popcount (cols72);
  same ('(72,64) unit columns among its 72 columns', sum (w == 1), 8);
  same ('(72,64) columns of weight 3', sum (w == 3), 47);
  same ('(72,64) columns of weight 5', sum (w == 5), 17);
  adj = bitxor (cols72(1:71), cols72(2:72));
  counts = accumarray (pair_sums (cols72)(:) + 1, 1, [256, 1]);
  same ('(72,64) fewest pairs sharing a neighbor sum', min (counts(adj + 1)), 8);
  same ('(72,64) neighbor sums with that fewest count, of 71', sum (counts(adj + 1) == min (counts(adj + 1))), 4);
  same ('(72,64) neighbor sums with more pairs than that', sum (counts(adj + 1) > min (counts(adj + 1))), 67);
  same ('pairs on each even syndrome at least (144 - 128) / 2', floor ((2 * 72 - 128) / 2), 8);
end

function x = interleave_extra_check_bits ()
  % Lay four (72,64) words out bit by bit, cell p holding bit floor(p / 4) of
  % word mod(p, 4), and count the cells holding check bits beyond the 4 x 8 the
  % words carry anyway.
  p = 0:4 * 72 - 1;
  bit = floor (p / 4);
  x = sum (bit >= 64) - 4 * (72 - 64);
end

function check_interleave ()
  section ('Four-way bit interleaving of four 72-bit words');
  cells = 4 * 72;
  same ('cells in the interleaved row', cells, 288);
  same ('word holding cell 13', mod (13, 4), 1);
  same ('bit of that word in cell 13', floor (13 / 4), 3);
  worst = zeros (1, 9);
  for len = 1:9
    top = 0;
    for start = 0:cells - len
      hits = accumarray (mod ((start:start + len - 1)', 4) + 1, 1, [4, 1]);
      top = max (top, max (hits));
    end
    worst(len) = top;
  end
  same ('most errors in one word, bursts of 1 to 9 cells', num_list (worst), '1 1 1 1 2 2 2 2 3');
  L = 1:9;
  same ('longest burst SEC-DED words always correct', max (L(worst <= 1)), 4);
  same ('longest burst SEC-DED words always detect', max (L(worst <= 2)), 8);
  same ('extra check bits the interleave costs', interleave_extra_check_bits (), 0);
end


% Single-symbol correction -----------------------------------------------------

function S = s4ec_r9 ()
  % A binary (73,64) code found by a greedy search: one line per symbol, each
  % column a 9-bit string read from top to bottom. Lines 1 to 16 are the data
  % symbols, lines 17 and 18 the two 4-bit check symbols and line 19 the 1-bit
  % check symbol.
  S = { ...
    '111010000 011100011 111111100 101111110', ...
    '000100110 011000001 110010100 100001101', ...
    '011111100 011011001 110110101 100000110', ...
    '101101000 010100011 111110100 001010010', ...
    '010111110 001111001 111000110 001000001', ...
    '110010110 011010010 101101111 110100101', ...
    '110011100 010001001 011111011 001010110', ...
    '111111101 001010011 101000101 111100000', ...
    '110011110 001100111 010110110 110111100', ...
    '011111110 101100001 101111010 010100010', ...
    '110111011 100011101 100110100 100001001', ...
    '001000011 101110001 000010111 000110110', ...
    '110111010 100000111 010001010 101101001', ...
    '100011010 110110000 001101100 110001011', ...
    '101000110 011011110 010110001 010000011', ...
    '000011110 100100010 111100001 001110110', ...
    '000000001 000000010 000000100 000001000', ...
    '000010000 000100000 001000000 010000000', ...
    '100000000'};
end

function n = symbol_patterns (n_bits, b)
  % Nonzero errors inside one symbol when n_bits bits are cut into b-bit
  % symbols in order, the last one shorter when b does not divide n_bits.
  n = floor (n_bits / b) * (2 ^ b - 1) + 2 ^ mod (n_bits, b) - 1;
end

function r = least_binary_symbol_bits (k, b)
  r = 1;
  while (1 + symbol_patterns (k + r, b) > 2 ^ r)
    r = r + 1;
  end
end

function G = binary_symbol_groups ()
  S = s4ec_r9 ();
  G = cell (1, numel (S));
  for i = 1:numel (S)
    G{i} = cellfun (@bin2dec, strsplit (S{i}, ' '));
  end
end

function p = group_patterns (g)
  % Syndromes of every nonzero error inside one symbol, in order of the error value.
  p = zeros (1, 2 ^ numel (g) - 1);
  for e = 1:2 ^ numel (g) - 1
    v = 0;
    for j = 0:numel (g) - 1
      if (bitand (bitshift (e, -j), 1))
        v = bitxor (v, g(j + 1));
      end
    end
    p(e) = v;
  end
end

function [mis, flag, unseen] = double_symbol_outcomes (P)
  % Errors in two different symbols: how many are miscorrected, flagged or unseen.
  single = unique ([P{:}]);
  mis = 0;
  flag = 0;
  unseen = 0;
  for i = 1:numel (P)
    for j = i + 1:numel (P)
      a = P{i}(:);
      b = P{j}(:)';
      v = bitxor (repmat (a, 1, numel (b)), repmat (b, numel (a), 1));
      v = v(:);
      z = (v == 0);
      m = ismember (v, single) & ~z;
      unseen = unseen + sum (z);
      mis = mis + sum (m);
      flag = flag + sum (~z & ~m);
    end
  end
end

function out = symbol_binary_columns ()
  % The 76 binary columns of the (76,64) code over GF(16), four per symbol.
  [ex, lg, wrap] = power_table (19, 4);
  T = mul_table (19, 4);
  cols = zeros (19, 3);
  for i = 0:14
    cols(i + 1, :) = [1, ex(i + 1), gmul(ex(i + 1), ex(i + 1), ex, lg)];
  end
  cols(16, :) = [0, 1, 1];
  cols(17:19, :) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
  out = zeros (1, 76);
  t = 0;
  for i = 1:19
    for j = 0:3
      S = symbol_scale (2 ^ j, cols(i, :), T);
      t = t + 1;
      out(t) = S(1) * 256 + S(2) * 16 + S(3);
    end
  end
end

function S = symbol_scale (lam, col, T)
  S = T(sub2ind ([16 16], lam * ones (1, 3) + 1, col + 1));
end

function check_symbol_code (ex, lg, T)
  section ('Single 4-bit symbol correction for 64 data bits');
  same ('nonzero error patterns inside one 4-bit symbol', sum ((0:15) ~= 0), 15);
  same ('longest code with 2 check symbols, (16^2 - 1) / 15', (16 ^ 2 - 1) / 15, 17);
  same ('data bits that 17 symbols leave room for', (17 - 2) * 4, 60);
  same ('syndromes 18 symbols need, 1 + 15 x 18', 1 + 15 * 18, 271);
  require ('271 exceeds 16^2 = 256', 1 + 15 * 18 > 256);
  same ('longest code with 3 check symbols, (16^3 - 1) / 15', (16 ^ 3 - 1) / 15, 273);
  same ('syndromes 19 symbols need, 1 + 15 x 19', 1 + 15 * 19, 286);
  cols = zeros (19, 3);
  for i = 0:14
    cols(i + 1, :) = [1, ex(i + 1), gmul(ex(i + 1), ex(i + 1), ex, lg)];
  end
  cols(16, :) = [0, 1, 1];
  cols(17:19, :) = [1, 0, 0; 0, 1, 0; 0, 0, 1];
  same ('symbols in the code', rows (cols), 19);
  require ('no column is zero', all (any (cols, 2)));
  ok = true;
  for i = 1:19
    for j = 1:19
      if (i ~= j)
        for lam = 1:15
          ok = ok && ~isequal (symbol_scale (lam, cols(i, :), T), cols(j, :));
        end
      end
    end
  end
  require ('no column is a nonzero multiple of another, 171 pairs', ok);
  synd = zeros (1, 19 * 15);
  t = 0;
  for i = 1:19
    for v = 1:15
      S = symbol_scale (v, cols(i, :), T);
      t = t + 1;
      synd(t) = S(1) * 256 + S(2) * 16 + S(3);
    end
  end
  same ('distinct nonzero syndromes of single-symbol errors', numel (setdiff (unique (synd), 0)), 285);
  binary_cols = zeros (1, 76);
  t = 0;
  for i = 1:19
    for j = 0:3
      S = symbol_scale (2 ^ j, cols(i, :), T);
      t = t + 1;
      binary_cols(t) = S(1) * 256 + S(2) * 16 + S(3);
    end
  end
  same ('rank of the 12 x 76 binary check matrix', rank_gf2 (binary_cols), 12);
  same ('binary length and data bits', sprintf ('%d,%d', numel (binary_cols), numel (binary_cols) - 12), '76,64');
  same ('overhead of the (76,64) code, percent', pct (12, 64), '18.75');
  require ('the same matrix as symbol_binary_columns() builds it', isequal (binary_cols, symbol_binary_columns ()));

  section ('A binary code for one bad 4-bit symbol with fewer check bits');
  same ('8 check bits: syndromes needed, 1 + 18 x 15', 1 + symbol_patterns (72, 4), 271);
  require ('271 exceeds 2^8 = 256, so 8 check bits fail', 1 + symbol_patterns (72, 4) > 2 ^ 8);
  same ('9 check bits: syndromes needed, 1 + 18 x 15 + 1', 1 + symbol_patterns (73, 4), 272);
  require ('272 fits in 2^9 = 512, so 9 check bits pass the count', 1 + symbol_patterns (73, 4) <= 2 ^ 9);
  same ('least check bits for a binary single-symbol code', least_binary_symbol_bits (64, 4), 9);
  G = binary_symbol_groups ();
  widths = cellfun (@numel, G);
  same ('symbols of 4 bits and of 1 bit in the 73-bit word', sprintf ('%d %d', sum (widths == 4), sum (widths == 1)), '18 1');
  flat = [G{:}];
  require ('every column has 9 rows', all (flat < 2 ^ 9));
  same ('rank of the 9 x 73 binary check matrix', rank_gf2 (flat), 9);
  same ('binary length and data bits', sprintf ('%d,%d', numel (flat), numel (flat) - rank_gf2 (flat)), '73,64');
  require ('the last 9 columns are the 9 unit columns, so the word is systematic', isequal (sort (flat(65:73)), 2 .^ (0:8)));
  P73 = cellfun (@group_patterns, G, 'UniformOutput', false);
  every = [P73{:}];
  same ('single-symbol errors in the (73,64) code', numel (every), 271);
  same ('distinct nonzero syndromes they give', numel (setdiff (unique (every), 0)), 271);
  same ('overhead of the (73,64) code, percent', pct (9, 64), '14.06');
  same ('4-bit chips the 73-bit and 76-bit words fill', sprintf ('%d %d', ceil (73 / 4), ceil (76 / 4)), '19 19');
  same ('4-bit chips the 9 check bits fill', ceil (rank_gf2 (flat) / 4), 3);
  same ('(76,64) check bits above the binary minimum', rank_gf2 (symbol_binary_columns ()) - least_binary_symbol_bits (64, 4), 3);

  section ('Errors in two different symbols');
  [mis73, flag73, unseen73] = double_symbol_outcomes (P73);
  same ('(73,64): double-symbol errors', mis73 + flag73 + unseen73, 34695);
  same ('(73,64): miscorrected, flagged, unseen', sprintf ('%d %d %d', mis73, flag73, unseen73), '17553 17142 0');
  same ('(73,64): share miscorrected, percent', pct1 (mis73, 34695), '50.6');
  P76 = cell (1, 19);
  for i = 1:19
    P76{i} = synd((i - 1) * 15 + (1:15));
  end
  [mis76, flag76, unseen76] = double_symbol_outcomes (P76);
  same ('(76,64): double-symbol errors', mis76 + flag76 + unseen76, 38475);
  same ('(76,64): miscorrected, flagged, unseen', sprintf ('%d %d %d', mis76, flag76, unseen76), '405 38070 0');
  same ('(76,64): share miscorrected, percent', pct1 (mis76, 38475), '1.1');

  section ('Correcting the error value alpha^2 in data symbol 5');
  S = symbol_scale (ex(3), cols(6, :), T);
  same ('syndrome S1, S2, S3 as bits', [bits(S(1), 4) ' ' bits(S(2), 4) ' ' bits(S(3), 4)], '0100 1011 1111');
  same ('S1, S2, S3 as powers of alpha', sprintf ('%d %d %d', lg(S(1) + 1), lg(S(2) + 1), lg(S(3) + 1)), '2 7 12');
  same ('S2 / S1 as a power of alpha, the locator', lg(gdiv (S(2), S(1), ex, lg) + 1), 5);
  same ('S3 / S2 as a power of alpha, the same locator', lg(gdiv (S(3), S(2), ex, lg) + 1), 5);
  same ('bits of symbol 5', sprintf ('d%d to d%d', 4 * 5, 4 * 5 + 3), 'd20 to d23');

  section ('Byte symbols for 64 data bits');
  same ('longest code with 2 byte check symbols, (256^2 - 1) / 255', (256 ^ 2 - 1) / 255, 257);
  same ('syndromes 10 byte symbols need, 1 + 255 x 10', 1 + 255 * 10, 2551);
  same ('overhead of the (80,64) byte code, percent', pct (16, 64), '25');
end


% Unidirectional errors --------------------------------------------------------

function z = zeros_in (word, k)
  z = k - popcount (word);
end

function r = berger_bits (k)
  % ceil(log2(k + 1)), the length of k written in binary.
  r = 0;
  while (2 ^ r < k + 1)
    r = r + 1;
  end
end

function c = covering_pairs (code)
  % Ordered pairs of distinct codewords X, Y where every 1 of Y is a 1 of X.
  [X, Y] = meshgrid (code, code);
  c = sum (sum (X ~= Y & bitand (X, Y) == Y));
end

function d = least_covering_distance (code)
  [X, Y] = meshgrid (code, code);
  cover = X ~= Y & bitand (X, Y) == Y;
  if (any (cover(:)))
    d = min (popcount (bitxor (X(cover), Y(cover))));
  else
    d = [];
  end
end

function check_berger ()
  section ('The Berger code');
  same ('check bits for 8 data bits', berger_bits (8), 4);
  same ('check bits for 64 data bits', berger_bits (64), 7);
  same ('check bits for 32 data bits', berger_bits (32), 6);
  same ('check bits for 7 data bits', berger_bits (7), 3);
  same ('overhead for 64 data bits, percent', pct (7, 64), '10.94');
  same ('overhead for 32 data bits, percent', pct (6, 32), '18.75');
  d = bin2dec ('10110010');
  same ('zeros in 10110010', zeros_in (d, 8), 4);
  same ('check of 10110010', bits (zeros_in (d, 8), 4), '0100');
  same ('check recomputed after 10110010 becomes 10100010', bits (zeros_in (bin2dec ('10100010'), 8), 4), '0101');
  same ('check recomputed after 10110010 becomes 11110010', bits (zeros_in (bin2dec ('11110010'), 8), 4), '0011');
  w = 0:255;
  code = w * 16 + zeros_in (w, 8);
  same ('codeword pairs where one covers the other, k = 8', covering_pairs (code), 0);
  ones_code = w * 16 + popcount (w);
  same ('the same count when the check counts ones', covering_pairs (ones_code) > 0, true);
  x = bin2dec ('11011010') * 16 + bin2dec ('0101');
  y = bin2dec ('11010010') * 16 + bin2dec ('0100');
  require ('11011010 0101 and 11010010 0100 are both ones-count words', any (ones_code == x) && any (ones_code == y));
  same ('bits the 1-to-0 error flips', popcount (bitxor (x, y)), 2);
  same ('least check bits for all-unidirectional detection, k = 64', 7, berger_bits (64));
end

function c = f_bl2 (z, r)
  F_MAP = [3, 5, 6, 9, 10, 12];
  low = r - 4;
  v = mod (z, 6 * 2 ^ low);
  c = F_MAP(floor (v / 2 ^ low) + 1) * 2 ^ low + mod (v, 2 ^ low);
end

function t = detects_up_to (f, k)
  % Largest t such that every unidirectional error of weight at most t is
  % detected, from the zero counts alone. Empty means every weight is detected.
  F = f (0:k);
  [Zw, Zz] = meshgrid (0:k, 0:k);
  keep = Zw > Zz;
  cz = F(Zz(keep) + 1);
  cw = F(Zw(keep) + 1);
  j = Zw(keep) - Zz(keep);
  cover = bitand (cz, cw) == cw;
  if (any (cover))
    d = j(cover)' + popcount (bitxor (cz(cover), cw(cover)));
    t = min (d) - 1;
  else
    t = [];
  end
end

function check_bose_lin ()
  section ('Counting zeros modulo 2^r');
  for r = 2:6
    same (sprintf ('zeros mod 2^%d detects up to, k = 64', r), detects_up_to (@(z) mod (z, 2 ^ r), 64), r);
  end
  same ('zeros mod 2^7, k = 64, is the Berger code: every weight', isempty (detects_up_to (@(z) mod (z, 128), 64)), true);
  d = bin2dec ('11011010');
  same ('zeros in 11011010', zeros_in (d, 8), 3);
  same ('its 2-bit check', bits (mod (zeros_in (d, 8), 4), 2), '11');
  same ('zeros after 11011010 becomes 10011010', zeros_in (bin2dec ('10011010'), 8), 4);
  same ('the check that word needs', bits (mod (zeros_in (bin2dec ('10011010'), 8), 4), 2), '00');

  section ('Bose-Lin codes');
  f_bl4 = @(z) mod (z, 8) + 4;
  same ('(zeros mod 8) + 4 detects up to, k = 64', detects_up_to (f_bl4, 64), 6);
  same ('check of a word with 2 zeros', bits (f_bl4 (2), 4), '0110');
  same ('check of a word with 8 zeros', bits (f_bl4 (8), 4), '0100');
  stated1 = [11, 20, 37];
  stated2 = [11, 22, 43];
  for r = 5:7
    f1 = @(z) mod (z, 2 ^ (r - 1)) + 2 ^ (r - 2);
    f2 = @(z) f_bl2 (z, r);
    same (sprintf ('offset count, r = %d, detects up to', r), detects_up_to (f1, 64), stated1(r - 4));
    same (sprintf ('offset count, r = %d, 2^(r-2) + r - 2', r), 2 ^ (r - 2) + r - 2, stated1(r - 4));
    same (sprintf ('second method, r = %d, detects up to', r), detects_up_to (f2, 64), stated2(r - 4));
    same (sprintf ('second method, r = %d, 5 x 2^(r-4) + r - 4', r), 5 * 2 ^ (r - 4) + r - 4, stated2(r - 4));
  end
  same ('offset 2^(r-1) instead, r = 5, detects only', detects_up_to (@(z) mod (z, 16) + 16, 64), 4);
  same ('Berger, k = 64, detects every weight', isempty (detects_up_to (@(z) z, 64)), true);

  section ('Every codeword pair, k = 8, against the zero-count rule');
  names = {'zeros mod 4', 'zeros mod 8', '(zeros mod 8) + 4', 'Berger'};
  fs = {@(z) mod (z, 4), @(z) mod (z, 8), f_bl4, @(z) z};
  rs = [2, 3, 4, 4];
  w = 0:255;
  for t = 1:4
    f = fs{t};
    code = w * 2 ^ rs(t) + f (zeros_in (w, 8));
    brute = least_covering_distance (code);
    if (isempty (brute))
      brute_t = 'every weight';
    else
      brute_t = sprintf ('%d', brute - 1);
    end
    rule = detects_up_to (f, 8);
    if (isempty (rule))
      rule_t = 'every weight';
    else
      rule_t = sprintf ('%d', rule);
    end
    same ([names{t} ': all 65,536 pairs agree with the rule'], brute_t, rule_t);
  end
end


% Choosing a code --------------------------------------------------------------

function r = least_daec_bits (k)
  % Least r with 1 + n + (n - 1) <= 2^r for n = k + r: zero, singles, adjacent doubles.
  r = 1;
  while (2 * (k + r) > 2 ^ r)
    r = r + 1;
  end
end

function c = least_symbol_checks (k, b)
  % Least check symbols c with k / b + c symbols within (q^c - 1) / (q - 1), q = 2^b.
  q = 2 ^ b;
  c = 1;
  while (floor (k / b) + c > (q ^ c - 1) / (q - 1))
    c = c + 1;
  end
end

function check_choice_table ()
  section ('Choosing a code for 64 data bits');
  models = {'one random bit', 'one random bit, two detected', 'two random bits, three detected', ...
            'two adjacent bits', 'one 4-bit symbol', 'unidirectional, any weight', ...
            'unidirectional, at most 6 bits'};
  codes = {'Hamming SEC (71,64)', 'Hsiao SEC-DED (72,64)', 'DEC-TED (79,64)', ...
           'SEC-DED-DAEC (72,64)', 'S4EC (76,64)', 'Berger (71,64)', 'Bose-Lin (68,64)'};
  % Each row's check bits are the rank of a check matrix this file builds, or
  % the width of a check it computes, so none of them is typed in.
  built = {hamming_columns(), hsiao_columns(), dec_ted_columns(), daec_search(8, 72), symbol_binary_columns()};
  ranks = cellfun (@rank_gf2, built);
  bose_lin_width = floor (log2 (max (mod (0:64, 8) + 4))) + 1;
  rbits = [ranks, berger_bits(64), bose_lin_width];
  stated = {'7 10.94', '8 12.5', '15 23.44', '8 12.5', '12 18.75', '7 10.94', '4 6.25'};
  for i = 1:7
    same ([models{i} ': ' codes{i}], sprintf ('%d %s', rbits(i), pct (rbits(i), 64)), stated{i});
  end
  same ('sphere-packing minimum, one random bit', least_check_bits (64, 1, false), 7);
  same ('sphere-packing minimum, two random bits, three detected', least_check_bits (64, 2, true), 13);
  same ('Bose-Lin with 4 check bits detects up to', detects_up_to (@(z) mod (z, 8) + 4, 64), 6);
  same ('SEC-DED-DAEC counting minimum, 1 + n + (n - 1) <= 2^r', least_daec_bits (64), 8);
  same ('counting at r = 7: 1 + 71 + 70', 1 + 71 + 70, 142);
  same ('single 4-bit symbol: least check symbols', least_symbol_checks (64, 4), 3);
  same ('offset count, r = 4, is the Bose-Lin (zeros mod 8) + 4', detects_up_to (@(z) mod (z, 2 ^ (4 - 1)) + 2 ^ (4 - 2), 64), 6);
  same ('4-way interleave: extra check bits', interleave_extra_check_bits (), 0);
  lengths = cellfun (@numel, built);
  same ('lengths of the five built matrices', num_list (lengths), '71 72 79 72 76');
  same ('data bits each leaves, length minus rank', num_list (lengths - ranks), '64 64 64 64 64');
  same ('single 4-bit symbol, binary counting minimum', least_binary_symbol_bits (64, 4), 9);
  same ('single 4-bit symbol, whole GF(16) check symbols', 4 * least_symbol_checks (64, 4), 12);
end


[ex16, lg16, wrap16] = power_table (19, 4);
T16 = mul_table (19, 4);
check_opening ();
check_sphere_counts ();
check_baselines ();
check_gf16 (ex16, lg16, wrap16, T16);
check_dec_identities (T16);
check_bch_15_7 (ex16, lg16);
check_bch_79_64 ();
check_daec ();
check_interleave ();
check_symbol_code (ex16, lg16, T16);
check_berger ();
check_bose_lin ();
check_choice_table ();
printf ('\nAll chapter 21 numbers reproduced.\n');
