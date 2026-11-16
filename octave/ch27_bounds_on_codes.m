% Reproduce every number in chapter 27 of Error Correction Coding,
% "Bounds on Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch27_bounds_on_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch27_bounds_on_codes.py and prints the same
% lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. The chapter's numbers are whole numbers, ratios printed to a
% fixed number of places, or lists of whole numbers, so every comparison is
% exact once the value is printed in the chapter's own format.
%
% A binary word of length n is held as a whole number whose bit j is
% position j of the word. A linear code is listed in full by XOR-ing every
% subset of the rows of its generator matrix, and its minimum distance is
% read as its smallest nonzero weight. The cyclic codes of length 15 have
% generator matrices built from the shifts of their generator polynomials.

1;


% Helpers ----------------------------------------------------------------------

function s = as_text (x)
  if (ischar (x))
    s = x;
  else
    s = sprintf ('%d', x);
  end
end

function same (label, computed, stated)
  c = as_text (computed);
  s = as_text (stated);
  printf ('  %-62s %-16s %s\n', label, c, s);
  if (~strcmp (c, s))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, c, s);
    exit (1);
  end
end

function require (label, condition)
  if (condition)
    printf ('  %-62s %s\n', label, 'holds');
  else
    printf ('  %-62s %s\n', label, 'FAILS');
    printf ('MISMATCH: ''%s'' does not hold\n', label);
    exit (1);
  end
end

function section (title)
  printf ('\n%s\n', title);
end

function s = join_numbers (v, sep)
  if (nargin < 2)
    sep = ' ';
  end
  parts = arrayfun (@(x) sprintf ('%d', x), v, 'UniformOutput', false);
  s = strjoin (parts, sep);
end


% Words and codes --------------------------------------------------------------

function w = weight (x)
  w = zeros (size (x));
  y = x;
  while (any (y(:) > 0))
    w = w + bitand (y, 1);
    y = bitshift (y, -1);
  end
end

function words = span_rows (rows)
  words = 0;
  for row = rows
    words = [words, bitxor(words, row)];
  end
end

function m = min_weight (words)
  w = weight (words);
  m = min (w(words ~= 0));
end

function m = min_pair_distance (words)
  m = Inf;
  for i = 1:numel (words)
    for j = i + 1:numel (words)
      m = min (m, weight (bitxor (words(i), words(j))));
    end
  end
end

function rows = poly_rows (g, n)
  deg = floor (log2 (g));
  rows = g * 2 .^ (0:(n - deg - 1));
end

function x = parse_word (bits)
  x = sum (2 .^ (find (bits == '1') - 1));
end

function g = times_x_plus_1 (g)
  g = bitxor (g, 2 * g);
end

function rows = simplex_rows (m)
  n = 2 ^ m - 1;
  rows = zeros (1, m);
  for i = 0:m - 1
    for j = 0:n - 1
      if (bitand (bitshift (j + 1, -i), 1))
        rows(i + 1) = rows(i + 1) + 2 ^ j;
      end
    end
  end
end

function rows = even_rows ()
  rows = 2 .^ (0:13) + 2 ^ 14;
end

function rows = rows_15_2_10 ()
  rows = [sum(2 .^ (0:9)), sum(2 .^ (5:14))];
end


% The bounds -------------------------------------------------------------------

function v = volume (n, t)
  v = 0;
  for i = 0:min (t, n)
    v = v + nchoosek (n, i);
  end
end

function k = largest_k (m)
  k = 0;
  while (2 ^ (k + 1) <= m)
    k = k + 1;
  end
end

function m = hamming_bound (n, d)
  m = floor (2 ^ n / volume (n, floor ((d - 1) / 2)));
end

function k = hamming_k (n, d)
  k = largest_k (hamming_bound (n, d));
end

function k = singleton_k (n, d)
  k = n - d + 1;
end

function m = plotkin_simple (n, d)
  m = floor ((2 * d) / (2 * d - n));
end

function m = plotkin_sharp (n, d)
  m = 2 * floor (d / (2 * d - n));
end

function m = plotkin_best (n, d)
  m = plotkin_sharp (n, d);
  if (mod (d, 2) == 1 && 2 * (d + 1) > n + 1)
    m = min (m, plotkin_sharp (n + 1, d + 1));
  end
end

function s = griesmer_sum (k, d)
  s = sum (ceil (d ./ 2 .^ (0:k - 1)));
end

function k = griesmer_k (n, d)
  k = 1;
  while (griesmer_sum (k + 1, d) <= n)
    k = k + 1;
  end
end

function g = gilbert_value (n, d)
  g = 2 ^ n / volume (n, d - 1);
end

function k = varshamov_k (n, d)
  v = volume (n - 1, d - 2);
  k = n - 1;
  while (k > 1 && ~(v < 2 ^ (n - k)))
    k = k - 1;
  end
end

function sums = forbidden_sums (columns, d)
  % Every sum of at most d - 2 of the given columns, the empty sum 0 included.
  sums = 0;
  for iter = 1:d - 2
    grown = sums;
    for c = columns
      grown = [grown, bitxor(sums, c)];
    end
    sums = unique (grown);
  end
end

function [counts, chosen] = varshamov_greedy (r, d)
  % Choose r-bit columns one at a time, each the smallest that no sum of at
  % most d - 2 earlier columns equals, and record how many were forbidden.
  chosen = [];
  counts = [];
  while (true)
    forbidden = forbidden_sums (chosen, d);
    counts(end + 1) = numel (forbidden);
    free = setdiff (0:2 ^ r - 1, forbidden);
    if (isempty (free))
      return;
    end
    chosen(end + 1) = free(1);
  end
end

function words = columns_code (columns)
  % The code whose parity-check matrix has the given columns, every word listed.
  n = numel (columns);
  words = [];
  for x = 0:2 ^ n - 1
    s = 0;
    for j = 0:n - 1
      if (bitand (bitshift (x, -j), 1))
        s = bitxor (s, columns(j + 1));
      end
    end
    if (s == 0)
      words(end + 1) = x;
    end
  end
end

function [value, best_r] = elias_bound (n, d)
  value = Inf;
  best_r = 0;
  for r = 1:floor (n / 2)
    den = r * r - n * r + n * d / 2;
    if (den > 0)
      v = (n * d / 2) / den * (2 ^ n / volume (n, r));
      if (v < value)
        value = v;
        best_r = r;
      end
    end
  end
end

function f = elias_floor (n, d)
  f = floor (elias_bound (n, d));
end

function k = elias_k (n, d)
  k = largest_k (elias_floor (n, d));
end

function k = exhibited_k (d)
  table = [14 11 10 7 6 5 4 2 2 1 1 1 1 1];
  k = table(d - 1);
end


% Hardware ---------------------------------------------------------------------

function r = sec_bits (k)
  r = 1;
  while (2 ^ r < k + r + 1)
    r = r + 1;
  end
end

function r = dec_bits (k)
  r = 1;
  while (2 ^ r < 1 + (k + r) + nchoosek (k + r, 2))
    r = r + 1;
  end
end


% The checks -------------------------------------------------------------------

function check_easy_values ()
  section ('Easy values of A_2(n,d)');
  even = span_rows (even_rows ());
  same ('even-weight code of length 15, number of words', numel (even), 16384);
  same ('even-weight code of length 15, minimum distance', min_weight (even), 2);
  same ('2^(n-1) at n = 15', 2 ^ 14, 16384);
  rep = span_rows (2 ^ 15 - 1);
  same ('repetition code of length 15, number of words', numel (rep), 2);
  same ('repetition code of length 15, minimum distance', min_weight (rep), 15);
end

function c = code_5_2_3 ()
  c = [parse_word('00000'), parse_word('10101'), parse_word('01111'), parse_word('11010')];
end

function check_small_exact ()
  section ('Small exact values');
  c523 = code_5_2_3 ();
  same ('(5,3): simple Plotkin 2d/(2d-n)', plotkin_simple (5, 3), 6);
  same ('(5,3): sharpened Plotkin on A_2(6,4), 2 floor(4/2)', plotkin_sharp (6, 4), 4);
  same ('(5,3): words in the exhibited code', numel (c523), 4);
  same ('(5,3): its minimum distance', min_pair_distance (c523), 3);
  closed = true;
  for a = c523
    for b = c523
      closed = closed && any (c523 == bitxor (a, b));
    end
  end
  require ('(5,3): the exhibited code is closed under XOR', closed);
  same ('(6,3): sharpened Plotkin on A_2(7,4), 2 floor(4/1)', plotkin_sharp (7, 4), 8);
  c633 = span_rows ([parse_word('100011'), parse_word('010101'), parse_word('001110')]);
  same ('(6,3): words in the shortened Hamming code', numel (c633), 8);
  same ('(6,3): its minimum distance', min_weight (c633), 3);
  same ('(7,3): Hamming bound 2^7/(1+7)', hamming_bound (7, 3), 16);
  c743 = span_rows ([parse_word('1000110'), parse_word('0100101'), parse_word('0010011'), parse_word('0001111')]);
  same ('(7,3): words in the Hamming code', numel (c743), 16);
  same ('(7,3): its minimum distance', min_weight (c743), 3);
  same ('exact values A_2(5,3), A_2(6,3), A_2(7,3)', ...
        join_numbers ([numel(c523), numel(c633), numel(c743)]), '4 8 16');
end

function check_nonlinear ()
  section ('A nonlinear code beats every linear one at (8,3)');
  c = [0 7 25 30 42 53 75 84 97 108 114 127 140 147 166 169 176 194 197 216];
  same ('words in the listed (8,3) code', numel (unique (c)), 20);
  same ('pairs checked', nchoosek (20, 2), 190);
  same ('smallest distance over all pairs', min_pair_distance (c), 3);
  require ('every listed word has 8 bits', all (c >= 0 & c < 256));
  ext = bitshift (c, 1) + mod (weight (c), 2);
  same ('extended by parity to length 9: words, distance', ...
        join_numbers ([numel(unique (ext)), min_pair_distance(ext)]), '20 4');
  same ('Hamming bound for (8,3), floor(256/9)', hamming_bound (8, 3), 28);
  same ('[8,5,3]: 2^5 (1 + 8)', 2 ^ 5 * (1 + 8), 288);
  require ('[8,5,3]: 288 > 2^8 = 256, so no such code', 2 ^ 5 * 9 > 2 ^ 8);
  c844 = span_rows ([parse_word('10001101'), parse_word('01001011'), parse_word('00100111'), parse_word('00011110')]);
  same ('extended Hamming code: words', numel (c844), 16);
  same ('extended Hamming code: minimum distance', min_weight (c844), 4);
  same ('largest linear code at (8,3), words', 2 ^ hamming_k (8, 3), 16);
end

function check_hamming ()
  section ('Hamming bound at n = 15');
  v = arrayfun (@(t) volume (15, t), 0:7);
  same ('V_2(15,t) for t = 0..7', join_numbers (v), '1 16 121 576 1941 4944 9949 16384');
  same ('largest k, d = 3..12', join_numbers (arrayfun (@(d) hamming_k (15, d), 3:12)), ...
        '11 11 8 8 5 5 4 4 2 2');
  same ('largest k, d = 13..15', join_numbers (arrayfun (@(d) hamming_k (15, d), 13:15)), '1 1 1');
  same ('Hamming bound on M at d = 5', hamming_bound (15, 5), 270);
  same ('Hamming bound on M at d = 9', hamming_bound (15, 9), 16);
  same ('Hamming bound on M at d = 11', hamming_bound (15, 11), 6);
  same ('sharpened Plotkin on M at d = 11', plotkin_best (15, 11), 2);
end

function check_singleton ()
  section ('Singleton bound');
  same ('largest k at n = 15, d = 2..15', join_numbers (arrayfun (@(d) singleton_k (15, d), 2:15)), ...
        '14 13 12 11 10 9 8 7 6 5 4 3 2 1');
  require ('repetition code (15,1,15) meets k = n - d + 1', 1 == 15 - 15 + 1);
  require ('even-weight code (15,14,2) meets k = n - d + 1', 14 == 15 - 2 + 1);
  require ('Hamming (7,4,3) misses it: 4 < 7 - 3 + 1 = 5', 4 < 5);
  same ('[5,3,3]: 2^3 V_2(5,1) against 2^5', join_numbers ([2 ^ 3 * volume(5, 1), 2 ^ 5]), '48 32');
end

function check_plotkin ()
  section ('Plotkin bound');
  rows = code_5_2_3 ();
  cols = arrayfun (@(j) sum (bitand (bitshift (rows, -j), 1)), 0:4);
  same ('(5,3) code: ones per column', join_numbers (cols), '2 2 2 2 2');
  contrib = 2 .* cols .* (4 - cols);
  same ('(5,3) code: 2w(M-w) per column', join_numbers (contrib), '8 8 8 8 8');
  ordered = 0;
  for a = rows
    for b = rows
      ordered = ordered + weight (bitxor (a, b));
    end
  end
  same ('(5,3) code: sum of distances over ordered pairs', ordered, 40);
  same ('(5,3) code: column total', sum (contrib), 40);
  same ('lower bound M(M-1)d', 4 * 3 * 3, 36);
  same ('upper bound n M^2 / 2', 5 * 16 / 2, 40);
  same ('simple form at n = 15, d = 8, 9, 10', ...
        join_numbers (arrayfun (@(d) plotkin_simple (15, d), [8 9 10])), '16 6 4');
  same ('sharpened, d = 9 through A_2(16,10)', plotkin_sharp (16, 10), 4);
  same ('sharpened, d = 8 and d = 10', join_numbers (arrayfun (@(d) plotkin_sharp (15, d), [8 10])), '16 4');
  same ('sharpened bound at n = 15, d = 8..15', ...
        join_numbers (arrayfun (@(d) plotkin_best (15, d), 8:15)), '16 4 4 2 2 2 2 2');
  simplex = span_rows (simplex_rows (4));
  same ('simplex code: words', numel (simplex), 16);
  same ('simplex code: minimum distance', min_weight (simplex), 8);
  w = weight (simplex);
  require ('simplex code: every nonzero weight is 8', all (w(simplex ~= 0) == 8));
  for d = 2:8
    require (sprintf ('A_2(2d,d) <= 4d at d = %d: shortened bound times 2', d), ...
             2 * plotkin_sharp (2 * d - 1, d) == 4 * d);
  end
end

function check_plotkin_extension ()
  section ('Exercise: Plotkin for n >= 2d at n = 15');
  ext = zeros (1, 6);
  vals = zeros (1, 6);
  for d = 2:7
    if (mod (d, 2) == 0)
      dd = d;
      nn = 15;
    else
      dd = d + 1;
      nn = 16;
    end
    vals(d - 1) = dd * 2 ^ (nn - 2 * dd + 2);
    ext(d - 1) = largest_k (vals(d - 1));
  end
  same ('bound on words, d = 2..7', join_numbers (vals), '16384 4096 2048 384 192 32');
  same ('largest k from d 2^(n-2d+2), d = 2..7', join_numbers (ext), '14 12 11 8 7 5');
  same ('Griesmer largest k, d = 2..7', join_numbers (arrayfun (@(d) griesmer_k (15, d), 2:7)), ...
        '14 12 11 8 7 5');
end

function check_griesmer ()
  section ('Griesmer bound');
  params = [7 4 3; 8 4 4; 15 4 8; 16 5 8; 32 6 16; 15 5 7; 15 11 3; 15 7 5; 23 12 7; 24 12 8];
  sums = [7 8 15 16 32 15 14 14 22 23];
  meets = {};
  for i = 1:rows (params)
    n = params(i, 1);
    k = params(i, 2);
    d = params(i, 3);
    same (sprintf ('sum for (%d,%d,%d)', n, k, d), griesmer_sum (k, d), sums(i));
    if (griesmer_sum (k, d) == n)
      meets{end + 1} = sprintf ('(%d,%d,%d)', n, k, d);
    end
  end
  same ('codes that meet it with equality', strjoin (meets, ' '), ...
        '(7,4,3) (8,4,4) (15,4,8) (16,5,8) (32,6,16) (15,5,7)');
  same ('(15,5,7) sum written out', join_numbers (ceil (7 ./ 2 .^ (0:4)), '+'), '7+4+2+1+1');
  same ('largest k at n = 15, d = 2..10', join_numbers (arrayfun (@(d) griesmer_k (15, d), 2:10)), ...
        '14 12 11 8 7 5 4 2 2');
  ds = 2:15;
  gk = arrayfun (@(d) griesmer_k (15, d), ds);
  hk = arrayfun (@(d) hamming_k (15, d), ds);
  same ('n = 15: Griesmer below Hamming at d', join_numbers (ds(gk < hk)), '2 6 8 9 10 11 12');
  same ('n = 15: Griesmer equal to Hamming at d', join_numbers (ds(gk == hk)), '4 5 7 13 14 15');
  same ('n = 15: Griesmer above Hamming at d', join_numbers (ds(gk > hk)), '3');
  same ('(15,2,10): sum 10 + 5', griesmer_sum (2, 10), 15);
  same ('Golay (23,12,7): volume V_2(23,3)', volume (23, 3), 2048);
  require ('Golay (23,12,7) is perfect: 2^12 * 2048 = 2^23', 2 ^ 12 * volume (23, 3) == 2 ^ 23);
end

function check_gv ()
  section ('Gilbert-Varshamov bound');
  same ('V_2(7,2)', volume (7, 2), 29);
  same ('Gilbert at (7,3), 128/29', sprintf ('%.2f', gilbert_value (7, 3)), '4.41');
  same ('so at least this many words', ceil (128 / 29), 5);
  n = 1:19;
  same ('Varshamov at d = 3, n - k = 3: largest n with n < 8', max (n(n < 2 ^ 3)), 7);
  same ('V_2(14, d-2) for d = 2..15', join_numbers (arrayfun (@(d) volume (14, d - 2), 2:15)), ...
        '1 15 106 470 1471 3473 6476 9908 12911 14913 15914 16278 16369 16383');
  same ('Varshamov largest k at n = 15, d = 2..8', ...
        join_numbers (arrayfun (@(d) varshamov_k (15, d), 2:8)), '14 11 8 6 4 3 2');
  same ('Varshamov largest k at n = 15, d = 9..15', ...
        join_numbers (arrayfun (@(d) varshamov_k (15, d), 9:15)), '1 1 1 1 1 1 1');
  g = arrayfun (@(d) sprintf ('%.3f', log2 (gilbert_value (15, d))), 2:8, 'UniformOutput', false);
  same ('Gilbert log2 at n = 15, d = 2..8', strjoin (g, ' '), ...
        '11.000 8.081 5.830 4.077 2.729 1.720 1.000');
  [counts, cols] = varshamov_greedy (3, 3);
  same ('greedy columns at d = 3, r = 3', join_numbers (cols), '1 2 3 4 5 6 7');
  same ('forbidden count at steps 1..8', join_numbers (counts), '1 2 3 4 5 6 7 8');
  words = columns_code (cols);
  same ('greedy code at d = 3, r = 3: (n, k, d)', ...
        sprintf ('(%d,%d,%d)', numel (cols), largest_k (numel (words)), min_weight (words)), '(7,4,3)');
  n = 2:19;
  ok = arrayfun (@(m) volume (m - 1, 2) < 2 ^ 4, n);
  same ('Varshamov at d = 4, r = 4: largest n with V(n-1,2) < 16', max (n(ok)), 5);
  [~, cols] = varshamov_greedy (4, 4);
  same ('greedy columns at d = 4, r = 4', join_numbers (cols), '1 2 4 7 8 11 13 14');
  words = columns_code (cols);
  same ('greedy code at d = 4, r = 4: (n, k, d)', ...
        sprintf ('(%d,%d,%d)', numel (cols), largest_k (numel (words)), min_weight (words)), '(8,4,4)');
end

function check_elias ()
  section ('Elias bound at n = 15');
  same ('Elias check at (10,3), floor', elias_floor (10, 3), 232);
  [~, r10] = elias_bound (10, 3);
  same ('(10,3): best r and V_2(10,1)', join_numbers ([r10, volume(10, 1)]), '1 11');
  [~, r15] = elias_bound (15, 4);
  same ('(15,4): best r and r^2 - 15r + 30 there', join_numbers ([r15, 2 ^ 2 - 15 * 2 + 30]), '2 4');
  same ('floor of the Elias bound, d = 2..8', join_numbers (arrayfun (@(d) elias_floor (15, d), 2:8)), ...
        '30720 5421 2031 883 284 104 30');
  same ('floor of the Elias bound, d = 9..15', join_numbers (arrayfun (@(d) elias_floor (15, d), 9:15)), ...
        '11 7 6 5 4 4 3');
  r = zeros (1, 14);
  for d = 2:15
    [~, r(d - 1)] = elias_bound (15, d);
  end
  same ('best r, d = 2..15', join_numbers (r), '1 1 2 2 3 4 7 7 7 7 7 7 7 7');
  same ('largest k, d = 2..15', join_numbers (arrayfun (@(d) elias_k (15, d), 2:15)), ...
        '14 12 10 9 8 6 4 3 2 2 2 2 2 1');
  ok = true;
  for d = 2:15
    ok = ok && elias_floor (15, d) >= 2 ^ exhibited_k (d);
  end
  require ('Elias never below an exhibited code at n = 15', ok);
  require ('Elias at d = 4 rules out [15,11,4]: 2031 < 2048', elias_floor (15, 4) < 2 ^ 11);
end

function check_codes_15 ()
  section ('Real codes at n = 15, every codeword enumerated');
  g_bch2 = parse_word ('100010111');
  labels = {'even-weight code', 'Hamming code, x^4+x+1', 'even-weight subcode, (x+1)(x^4+x+1)', ...
            'degree-8 generator', '(x+1) times the degree-8 generator', 'degree-10 generator', ...
            'simplex code', 'two rows of weight 10', 'repetition code'};
  rowsets = {even_rows(), poly_rows(19, 15), poly_rows(53, 15), poly_rows(g_bch2, 15), ...
             poly_rows(times_x_plus_1 (g_bch2), 15), poly_rows(parse_word ('11101100101'), 15), ...
             simplex_rows(4), rows_15_2_10(), 2 ^ 15 - 1};
  stated = [14 2 16384; 11 3 2048; 10 4 1024; 7 5 128; 6 6 64; 5 7 32; 4 8 16; 2 10 4; 1 15 2];
  for i = 1:numel (labels)
    words = span_rows (rowsets{i});
    k = largest_k (numel (words));
    d = min_weight (words);
    same ([labels{i} ': (k, d)'], sprintf ('(%d,%d)', k, d), sprintf ('(%d,%d)', stated(i, 1), stated(i, 2)));
    require ([labels{i} ': 2^k distinct words'], numel (words) == 2 ^ k && 2 ^ k == stated(i, 3));
  end
  for d = 2:15
    k = exhibited_k (d);
    up = min ([singleton_k(15, d), hamming_k(15, d), griesmer_k(15, d), elias_k(15, d)]);
    require (sprintf ('exhibited k = %d at d = %d within every upper bound', k, d), ...
             k <= up && (2 * d <= 15 || 2 ^ k <= plotkin_best (15, d)));
    require (sprintf ('exhibited k = %d at d = %d at least Varshamov''s k', k, d), ...
             k >= varshamov_k (15, d));
  end
  same ('[15,11,4] punctured to [14,11,3]: 2^(14-11)', 2 ^ 3, 8);
  require ('8 < 1 + 14 = 15, so no [14,11,3] and no [15,11,4]', 2 ^ 3 < 15);
  same ('Hamming at d = 3: 2^11 * 16', 2 ^ 11 * 16, 2 ^ 15);
  same ('simplex: Griesmer sum', griesmer_sum (4, 8), 15);
  same ('simplex: Plotkin value at (15,8)', plotkin_sharp (15, 8), 16);
  same ('(15,5,7): 32768/576', sprintf ('%.1f', 32768 / 576), '56.9');
  same ('(15,5,7): largest k from Hamming', hamming_k (15, 7), 5);
  same ('(15,5,7): Griesmer sum', griesmer_sum (5, 7), 15);
  same ('Hamming largest k at d = 2', hamming_k (15, 2), 15);
  same ('Griesmer largest k, d = 11..15', join_numbers (arrayfun (@(d) griesmer_k (15, d), 11:15)), '1 1 1 1 1');
  same ('Plotkin largest k, d = 8..15', ...
        join_numbers (arrayfun (@(d) largest_k (plotkin_best (15, d)), 8:15)), '4 2 2 1 1 1 1 1');
  upper = zeros (1, 14);
  shown = zeros (1, 14);
  for d = 2:15
    u = min ([singleton_k(15, d), hamming_k(15, d), griesmer_k(15, d), elias_k(15, d)]);
    if (2 * d > 15)
      u = min (u, largest_k (plotkin_best (15, d)));
    end
    upper(d - 1) = u;
    shown(d - 1) = exhibited_k (d);
  end
  same ('smallest upper bound on k, d = 2..15', join_numbers (upper), '14 11 10 8 7 5 4 2 2 1 1 1 1 1');
  same ('exhibited k, d = 2..15', join_numbers (shown), '14 11 10 7 6 5 4 2 2 1 1 1 1 1');
  ds = 2:15;
  same ('distances where no bound applied directly meets a code', join_numbers (ds(upper > shown)), '5 6');
  a154 = floor (2 ^ 14 / volume (14, 1));
  same ('A_2(15,4) = A_2(14,3) <= floor(2^14 / V_2(14,1))', a154, 1092);
  require ('1092 < 2048, and below the Elias value 2031', ...
           a154 < 2 ^ 11 && a154 < elias_floor (15, 4));
end

function check_residual ()
  section ('Residual codes close d = 5 and d = 6 at n = 15');
  same ('[15,8,5] residual: length, dimension, distance', join_numbers ([15 - 5, 8 - 1, ceil(5 / 2)]), ...
        '10 7 3');
  same ('[10,7,3]: 2^7 * V_2(10,1) against 2^10', join_numbers ([2 ^ 7 * volume(10, 1), 2 ^ 10]), ...
        '1408 1024');
  same ('[15,7,6] residual: length, dimension, distance', join_numbers ([15 - 6, 7 - 1, ceil(6 / 2)]), ...
        '9 6 3');
  same ('[9,6,3]: 2^6 * V_2(9,1) against 2^9', join_numbers ([2 ^ 6 * volume(9, 1), 2 ^ 9]), ...
        '640 512');
  m = 1:10;
  require ('2^7 (1 + m) > 2^m for m = 1..10, so no [15,8,5]', all (2 ^ 7 * (1 + m) > 2 .^ m));
  m = 1:9;
  require ('2^6 (1 + m) > 2^m for m = 1..9, so no [15,7,6]', all (2 ^ 6 * (1 + m) > 2 .^ m));
  require ('2^m / (1 + m) increases for m = 1..10', all (2 .^ m .* (2 + m) < 2 .^ (m + 1) .* (1 + m)));
  same ('Griesmer sums for k = 8, d = 5 and k = 7, d = 6', ...
        join_numbers ([griesmer_sum(8, 5), griesmer_sum(7, 6)]), '15 15');
end

function check_hardware ()
  section ('Check bits for 64 data bits');
  widths = [4 8 16 32 64 128 256];
  same ('SEC check bits at k = 64', sec_bits (64), 7);
  same ('SECDED check bits at k = 64', sec_bits (64) + 1, 8);
  same ('DEC check bits at k = 64, sphere-packing floor', dec_bits (64), 12);
  same ('DEC at r = 12: 1 + 76 + C(76,2)', 1 + 76 + nchoosek (76, 2), 2927);
  same ('DEC-TED floor at k = 64: DEC floor plus 1', dec_bits (64) + 1, 13);
  same ('DEC check bits at k = 32', dec_bits (32), 10);
  same ('DEC at k = 32, r = 10: 1 + 42 + C(42,2) against 2^10', join_numbers ([1 + 42 + nchoosek(42, 2), 2 ^ 10]), '904 1024');
  same ('DEC at k = 32, r = 9: 1 + 41 + C(41,2) against 2^9', join_numbers ([1 + 41 + nchoosek(41, 2), 2 ^ 9]), '862 512');
  same ('SEC at k = 4..256', join_numbers (arrayfun (@sec_bits, widths)), '3 4 5 6 7 8 9');
  same ('SECDED at k = 4..256', join_numbers (arrayfun (@sec_bits, widths) + 1), '4 5 6 7 8 9 10');
  same ('DEC floor at k = 4..256', join_numbers (arrayfun (@dec_bits, widths)), '6 7 9 10 12 14 16');
  same ('SECDED with 7 checks, punctured to [70,64,3]: 2^6', 2 ^ 6, 64);
  require ('64 < 1 + 70, so 7 check bits cannot give SECDED', 2 ^ 6 < 71);
  same ('(8,3) nonlinear code: log2 20 to three places', sprintf ('%.3f', log2 (20)), '4.322');
end


function check_examples ()
  section ('Worked examples and exercises');
  same ('(5,3): Hamming bound floor(32/6)', hamming_bound (5, 3), 5);
  c = code_5_2_3 ();
  dists = [];
  for i = 1:numel (c)
    for j = i + 1:numel (c)
      dists(end + 1) = weight (bitxor (c(i), c(j)));
    end
  end
  same ('(5,3) code: the six pairwise distances', join_numbers (sort (dists)), '3 3 3 3 4 4');
  same ('Griesmer sum for k = 6, d = 7', griesmer_sum (6, 7), 16);
  same ('(15,11): 32768/4944 to two places', sprintf ('%.2f', 32768 / 4944), '6.63');
  same ('128 data bits: SEC, SECDED, DEC floor', ...
        join_numbers ([sec_bits(128), sec_bits(128) + 1, dec_bits(128)]), '8 9 14');
  same ('SEC at k = 128, r = 7: 2^7 against 128 + 7 + 1', join_numbers ([2 ^ 7, 128 + 7 + 1]), '128 136');
  same ('DEC at k = 128, r = 14: 1 + 142 + C(142,2) against 2^14', ...
        join_numbers ([1 + 142 + nchoosek(142, 2), 2 ^ 14]), '10154 16384');
  same ('DEC at k = 128, r = 13: 1 + 141 + C(141,2) against 2^13', ...
        join_numbers ([1 + 141 + nchoosek(141, 2), 2 ^ 13]), '10012 8192');
  same ('Hamming bound at (10,3), floor(1024/11)', hamming_bound (10, 3), 93);
  same ('Singleton largest k at (5,3)', singleton_k (5, 3), 3);
  same ('[15,3,10]: Griesmer sum 10 + 5 + 3', griesmer_sum (3, 10), 18);
  same ('[15,2,10]: Plotkin value at (15,10)', plotkin_best (15, 10), 4);
end

check_easy_values ();
check_small_exact ();
check_nonlinear ();
check_hamming ();
check_singleton ();
check_plotkin ();
check_plotkin_extension ();
check_griesmer ();
check_gv ();
check_elias ();
check_codes_15 ();
check_residual ();
check_hardware ();
check_examples ();
printf ('\nAll chapter 27 numbers reproduced.\n');
