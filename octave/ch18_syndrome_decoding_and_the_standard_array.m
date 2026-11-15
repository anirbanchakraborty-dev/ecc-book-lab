% Reproduce every number in chapter 18 of Error Correction Coding,
% "Syndrome Decoding and the Standard Array", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch18_syndrome_decoding_and_the_standard_array.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch18_syndrome_decoding_and_the_standard_array.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. A binary word is a row of zeros and ones, a list of words is a
% matrix with one word per row, and a probability is printed to six places.
% Two codes carry the chapter. C6 is the (6,3) code with generator rows
% 100110, 010011 and 001111. The (7,4) Hamming code is the one of Linear
% Codes, with generator rows 1000110, 0100011, 0010101 and 0001111 and
% parity-check columns 110, 011, 101, 111, 100, 010 and 001. The small cases are checked exhaustively: every word of length 6
% and 7 in the standard arrays, every pair of words for the syndrome theorem,
% every error pattern for the probability of correct decoding, and every set
% of erasures for the erasure counts.

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

function M = parse (txt)
  parts = strsplit (txt, ' ');
  M = zeros (numel (parts), numel (parts{1}));
  for i = 1:numel (parts)
    M(i, :) = parts{i} - '0';
  end
end

function s = word (v)
  s = char (v + '0');
end

function s = words (M)
  parts = cell (1, rows (M));
  for i = 1:rows (M)
    parts{i} = word (M(i, :));
  end
  s = strjoin (parts, ' ');
end

function s = nums (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ' ');
end

function s = f6 (x)
  s = sprintf ('%.6f', x);
end

function z = bxor (a, b)
  z = mod (a + b, 2);
end

function k = wt (v)
  k = sum (v, 2);
end

function u = unit (n, j)
  % Position j counted from 1.
  u = zeros (1, n);
  u(j) = 1;
end

function W = all_words (n)
  % Every binary word of length n in increasing binary value.
  W = dec2bin (0:2^n - 1, n) - '0';
end

function v = value (w)
  n = columns (w);
  v = w * (2 .^ (n - 1:-1:0))';
end

function C = codewords (G)
  % The codewords u G for the messages u in increasing binary value.
  C = mod (all_words (rows (G)) * G, 2);
end

function s = syndrome (H, v)
  % s = v H^T, one bit per row of H.
  s = mod (v * H', 2);
end

function L = cols_of (H)
  L = H';
end

function [R, leaders] = standard_array (code, n)
  % Slepian's layout, built as the chapter's definition states. The first row
  % is the code. Each later row is led by a word of least weight among those
  % not yet written, ties going to the smaller binary value, and holds that
  % leader XORed with each codeword in the order of the first row.
  listed = false (1, 2^n);
  listed(value (code) + 1) = true;
  R = {code};
  leaders = code(1, :);
  W = all_words (n);
  [~, idx] = sort (wt (W));
  order = W(idx, :);
  for i = 1:rows (order)
    v = order(i, :);
    if (listed(value (v) + 1))
      continue;
    end
    leader = v;
    row = bxor (repmat (leader, rows (code), 1), code);
    R{end + 1} = row;
    leaders(end + 1, :) = leader;
    listed(value (row) + 1) = true;
  end
end

function T = syndrome_table (H, leaders)
  % Row s + 1 holds the leader whose syndrome has binary value s.
  T = zeros (2^rows (H), columns (H));
  for i = 1:rows (leaders)
    T(value (syndrome (H, leaders(i, :))) + 1, :) = leaders(i, :);
  end
end

function c = decode (H, T, r)
  c = bxor (r, T(value (syndrome (H, r)) + 1, :));
end

function L = likelihood (p, n, w)
  L = p ^ w * (1 - p) ^ (n - w);
end

function counts = leader_counts (leaders, n)
  counts = zeros (1, n + 1);
  for i = 1:rows (leaders)
    counts(wt (leaders(i, :)) + 1) += 1;
  end
end

function t = p_correct (counts, n, p)
  t = 0;
  for i = 0:n
    t += counts(i + 1) * likelihood (p, n, i);
  end
end

function t = p_correct_by_enumeration (H, T, code, n, p)
  % Sum the probabilities of the error patterns after which the decoder returns
  % the codeword that was sent, with the all-zeros codeword sent.
  zero = code(1, :);
  W = all_words (n);
  t = 0;
  for i = 1:rows (W)
    e = W(i, :);
    if (isequal (decode (H, T, bxor (zero, e)), zero))
      t += likelihood (p, n, wt (e));
    end
  end
end

function out = erasure_solutions (H, r, erased)
  % All codewords agreeing with r off the erased positions, found by trying
  % every filling of the erased positions and keeping those whose syndrome is
  % zero. Positions are counted from 1.
  out = zeros (0, columns (r));
  F = all_words (numel (erased));
  for i = 1:rows (F)
    v = r;
    v(erased) = F(i, :);
    if (~any (syndrome (H, v)))
      out(end + 1, :) = v;
    end
  end
end

function ok = independent (C)
  % True when no nonempty subset of the rows of C adds to zero.
  k = rows (C);
  M = all_words (k);
  ok = true;
  for i = 2:rows (M)
    if (~any (mod (M(i, :) * C, 2)))
      ok = false;
      return;
    end
  end
end

function v = blank (c, erased)
  v = c;
  v(erased) = 0;
end

function k = row_holding (R, r)
  k = 0;
  for i = 1:numel (R)
    if (any (all (R{i} == r, 2)))
      k = i;
      return;
    end
  end
end

function k = index_in (M, r)
  k = find (all (M == r, 2), 1);
end

function g = decoder_netlist (H)
  % The AND-gate decoder as a list of gates, each with a kind, an output name
  % and its input names, in an order that evaluates.
  m = rows (H);
  n = columns (H);
  g = struct ('kind', {}, 'out', {}, 'ins', {});
  for i = 1:m
    level = arrayfun (@(j) sprintf ('r%d', j), find (H(i, :)), 'UniformOutput', false);
    k = 0;
    while (numel (level) > 1)
      nxt = {};
      for a = 1:2:numel (level) - 1
        k += 1;
        if (numel (level) == 2)
          out = sprintf ('s%d', i);
        else
          out = sprintf ('t%d%d', i, k);
        end
        g(end + 1) = struct ('kind', 'XOR', 'out', out, 'ins', {{level{a}, level{a + 1}}});
        nxt{end + 1} = out;
      end
      if (mod (numel (level), 2))
        nxt{end + 1} = level{end};
      end
      level = nxt;
    end
  end
  for i = 1:m
    g(end + 1) = struct ('kind', 'NOT', 'out', sprintf ('n%d', i), 'ins', {{sprintf('s%d', i)}});
  end
  for j = 1:n
    ins = cell (1, m);
    for i = 1:m
      if (H(i, j))
        ins{i} = sprintf ('s%d', i);
      else
        ins{i} = sprintf ('n%d', i);
      end
    end
    g(end + 1) = struct ('kind', 'AND', 'out', sprintf ('e%d', j), 'ins', {ins});
  end
  for j = 1:n
    g(end + 1) = struct ('kind', 'XOR', 'out', sprintf ('c%d', j), ...
                         'ins', {{sprintf('r%d', j), sprintf('e%d', j)}});
  end
end

function [c, depth] = run_netlist (g, r)
  % Evaluate the netlist on the received word r and return the corrected word
  % and each signal's depth in gates.
  val = containers.Map ();
  depth = containers.Map ();
  for j = 1:numel (r)
    val(sprintf ('r%d', j)) = r(j);
    depth(sprintf ('r%d', j)) = 0;
  end
  for q = 1:numel (g)
    bits = cellfun (@(x) val(x), g(q).ins);
    switch (g(q).kind)
      case 'XOR'
        v = xor (bits(1), bits(2));
      case 'NOT'
        v = 1 - bits(1);
      otherwise
        v = all (bits);
    end
    val(g(q).out) = double (v);
    depth(g(q).out) = 1 + max (cellfun (@(x) depth(x), g(q).ins));
  end
  c = arrayfun (@(j) val(sprintf ('c%d', j)), 1:numel (r));
end


% The two codes ----------------------------------------------------------------

G6 = parse ('100110 010011 001111');
H6 = parse ('101100 111010 011001');
C6 = codewords (G6);
G7 = parse ('1000110 0100011 0010101 0001111');
H7 = parse ('1011100 1101010 0111001');
C7 = codewords (G7);


% Checks -----------------------------------------------------------------------

printf ('Chapter 18: Syndrome Decoding and the Standard Array\n');

section ('The (6,3) code C6');
same ('G', words (G6), '100110 010011 001111');
same ('H = [P^T | I]', words (H6), '101100 111010 011001');
require ('every row of G passes every check of H', ~any (syndrome (H6, G6)(:)));
same ('codewords for messages 000 to 111', words (C6), ...
      '000000 001111 010011 011100 100110 101001 110101 111010');
same ('weights of the codewords', nums (wt (C6)), '0 4 3 3 3 3 4 4');
same ('minimum distance of C6', min (wt (C6(2:end, :))), 3);

section ('Which codeword was sent');
c = C6(7, :);
same ('codeword for message 110', word (c), '110101');
r = bxor (c, unit (6, 3));
same ('received word after an error in position 3', word (r), '111101');
dist = wt (bxor (repmat (r, 8, 1), C6))';
same ('distances from 111101 to the eight codewords', nums (dist), '5 3 4 2 4 2 1 3');
best = min (dist);
same ('nearest codeword', word (C6(find (dist == best, 1), :)), '110101');
same ('codewords at the least distance', sum (dist == best), 1);
stated = {'0.059049', '0.006561', '0.000729', '0.000081', '0.000009'};
for w = 1:5
  same (sprintf ('likelihood at p = 0.1 of a codeword at distance %d', w), f6 (likelihood (0.1, 6, w)), stated{w});
end
same ('likelihood ratio, distance 1 over distance 2', f6 (likelihood (0.1, 6, 1) / likelihood (0.1, 6, 2)), '9.000000');
total = 0;
for k = 1:8
  total += likelihood (0.1, 6, dist(k));
end
same ('total likelihood over the eight codewords at p = 0.1', f6 (total), '0.073800');
same ('share of that total held by 110101', f6 (likelihood (0.1, 6, 1) / total), '0.800122');

section ('Maximum likelihood when p is not below one half');
r = parse ('111101');
dist = wt (bxor (repmat (r, 8, 1), C6))';
same ('likelihood of every codeword at p = 1/2', f6 (likelihood (0.5, 6, 3)), '0.015625');
scores = arrayfun (@(d) likelihood (0.5, 6, d), dist);
require ('at p = 1/2 all eight likelihoods are equal', numel (unique (scores)) == 1);
scores = arrayfun (@(d) likelihood (0.9, 6, d), dist);
[top, at] = max (scores);
same ('ML codeword at p = 0.9', word (C6(at, :)), '000000');
same ('its distance from 111101', dist(at), 5);
same ('its likelihood at p = 0.9', f6 (top), '0.059049');
flipped = 1 - r;
same ('111101 with every bit flipped', word (flipped), '000010');
fd = wt (bxor (repmat (flipped, 8, 1), C6))';
same ('nearest codeword to the flipped word', word (C6(find (fd == min (fd), 1), :)), '000000');

section ('The standard array of C6');
[R6, L6] = standard_array (C6, 6);
stated = {
  '000000 001111 010011 011100 100110 101001 110101 111010'
  '000001 001110 010010 011101 100111 101000 110100 111011'
  '000010 001101 010001 011110 100100 101011 110111 111000'
  '000100 001011 010111 011000 100010 101101 110001 111110'
  '001000 000111 011011 010100 101110 100001 111101 110010'
  '010000 011111 000011 001100 110110 111001 100101 101010'
  '100000 101111 110011 111100 000110 001001 010101 011010'
  '000101 001010 010110 011001 100011 101100 110000 111111'
};
same ('number of rows', numel (R6), 8);
for i = 1:8
  same (sprintf ('row %d', i), words (R6{i}), stated{i});
end
flat = vertcat (R6{:});
same ('words in the array', rows (flat), 64);
require ('every word of length 6 appears exactly once', isequal (sortrows (flat), all_words (6)));
same ('leader weights', nums (wt (L6)), '0 1 1 1 1 1 1 2');
require ('each leader has the least weight in its own row', ...
         all (arrayfun (@(i) wt (L6(i, :)) == min (wt (R6{i})), 1:8)));
last = R6{8};
same ('least-weight words of the last coset', words (sortrows (last(wt (last) == 2, :))), ...
      '000101 001010 110000');
r = parse ('111101');
i = row_holding (R6, r);
same ('row holding 111101', i, 5);
same ('codeword at the head of its column', word (C6(index_in (R6{i}, r), :)), '110101');
require ('every row has 8 words and the rows are disjoint', ...
         all (cellfun (@(M) rows (unique (M, 'rows')), R6) == 8) && rows (unique (flat, 'rows')) == 64);

section ('Decoding with the array');
w0 = parse ('000101');
dist = wt (bxor (repmat (w0, 8, 1), C6))';
same ('distances from 000101 to the codewords', nums (dist), '2 2 3 3 3 3 2 6');
same ('codewords at distance 2 from 000101', words (C6(dist == 2, :)), '000000 001111 110101');
same ('leader decoding of 000101', word (bxor (w0, L6(8, :))), '000000');
same ('decision region of 000000', words (L6), ...
      '000000 000001 000010 000100 001000 010000 100000 000101');
p = 0.1;
total = sum (arrayfun (@(d) likelihood (p, 6, d), dist));
same ('total likelihood of the codewords given 000101 at p = 0.1', f6 (total), '0.022600');
same ('share of that total held by each tied codeword', f6 (likelihood (p, 6, 2) / total), '0.290310');
ok = true;
for k = 1:499
  q = k / 1000;
  ok = ok && likelihood (q, 6, 2) / sum (arrayfun (@(d) likelihood (q, 6, d), dist)) < 1 / 3;
end
require ('each tied codeword holds less than 1/3 at p = k/1000 for k = 1 to 499', ok);
region = bxor (repmat (C6(1, :), 8, 1), L6);
same ('weights of the words in the decision region of 000000', nums (wt (region)'), '0 1 1 1 1 1 1 2');
ok = true;
W = all_words (6);
for k = 1:rows (W)
  r = W(k, :);
  out = bxor (r, L6(row_holding (R6, r), :));
  d = wt (bxor (repmat (r, 8, 1), C6));
  if (isempty (index_in (C6, out)) || wt (bxor (r, out)) ~= min (d))
    ok = false;
  end
end
require ('leader decoding returns a nearest codeword for every word', ok);

section ('Syndromes of C6');
same ('columns of H', words (cols_of (H6)), '110 011 111 100 010 001');
same ('syndrome of 111101', word (syndrome (H6, parse ('111101'))), '111');
same ('syndromes of the leaders, row by row', words (syndrome (H6, L6)), ...
      '000 001 010 100 111 011 110 101');
require ('all words of a row share one syndrome', ...
         all (cellfun (@(M) rows (unique (syndrome (H6, M), 'rows')), R6) == 1));
heads = zeros (8, 3);
for k = 1:8
  heads(k, :) = syndrome (H6, R6{k}(1, :));
end
require ('different rows have different syndromes', rows (unique (heads, 'rows')) == 8);
W = all_words (6);
sv = value (syndrome (H6, W));
cv = value (C6);
ok = true;
for a = 1:64
  for b = 1:64
    ok = ok && ((sv(a) == sv(b)) == any (cv == value (bxor (W(a, :), W(b, :)))));
  end
end
require ('s(x) = s(y) exactly when x + y is a codeword', ok);
T6 = syndrome_table (H6, L6);
stated = {'000000', '000001', '000010', '010000', '000100', '000101', '100000', '001000'};
S3 = all_words (3);
for k = 1:8
  same (sprintf ('C6 table, syndrome %s', word (S3(k, :))), word (T6(k, :)), stated{k});
end
missing = setdiff (S3(2:end, :), cols_of (H6), 'rows');
same ('syndrome missing from the columns of H', word (missing(1, :)), '101');

section ('The (7,4) Hamming code');
same ('P', words (G7(:, 5:7)), '110 011 101 111');
same ('G = [I | P]', words (G7), '1000110 0100011 0010101 0001111');
same ('H = [P^T | I]', words (H7), '1011100 1101010 0111001');
require ('every row of G passes every check of H', ~any (syndrome (H7, G7)(:)));
same ('columns of H', words (cols_of (H7)), '110 011 101 111 100 010 001');
require ('the seven columns are the seven nonzero words of length 3', ...
         isequal (sortrows (cols_of (H7)), S3(2:end, :)));
same ('codewords', rows (C7), 16);
counts = zeros (1, 8);
for k = 1:16
  counts(wt (C7(k, :)) + 1) += 1;
end
same ('codewords of weight 0 to 7', nums (counts), '1 0 0 7 7 0 0 1');
[R7, L7] = standard_array (C7, 7);
same ('rows of its standard array', numel (R7), 8);
require ('its standard array lists all 128 words once', isequal (sortrows (vertcat (R7{:})), all_words (7)));
T7 = syndrome_table (H7, L7);
stated = {'0000000', '0000001', '0000010', '0100000', '0000100', '0010000', '1000000', '0001000'};
for k = 1:8
  same (sprintf ('(7,4) table, syndrome %s', word (S3(k, :))), word (T7(k, :)), stated{k});
end

section ('Decoding the (7,4) code by lookup');
c = C7(bin2dec ('1001') + 1, :);
same ('codeword for message 1001', word (c), '1001001');
same ('it is the first row of G plus the fourth', word (bxor (G7(1, :), G7(4, :))), '1001001');
r = bxor (c, unit (7, 6));
same ('received word, error in position 6', word (r), '1001011');
same ('its syndrome', word (syndrome (H7, r)), '010');
same ('leader read from the table', word (T7(value (syndrome (H7, r)) + 1, :)), '0000010');
same ('decoded codeword', word (decode (H7, T7, r)), '1001001');
e = parse ('0100100');
r2 = bxor (c, e);
same ('received word, errors in positions 2 and 5', word (r2), '1101101');
same ('its syndrome', word (syndrome (H7, r2)), '111');
same ('column 2 plus column 5 of H', word (bxor (H7(:, 2)', H7(:, 5)')), '111');
same ('column of H equal to that syndrome', index_in (cols_of (H7), [1 1 1]), 4);
same ('decoded codeword', word (decode (H7, T7, r2)), '1100101');
same ('positions where it differs from the sent word', wt (bxor (decode (H7, T7, r2), c)), 3);
d2 = decode (H7, T7, r2);
same ('message bits delivered', word (d2(1:4)), '1100');
cl = cols_of (H7);
pairs = nchoosek (1:7, 2);
ok = true;
for k = 1:rows (pairs)
  x = bxor (cl(pairs(k, 1), :), cl(pairs(k, 2), :));
  ok = ok && ~isempty (index_in (cl, x)) && ~isequal (x, cl(pairs(k, 1), :)) && ~isequal (x, cl(pairs(k, 2), :));
end
require ('(7,4): the XOR of two distinct columns is always a third column', ok);
r4 = parse ('100100');
same ('C6: syndrome of 100100', word (syndrome (H6, r4)), '010');
same ('C6: decoding of 100100', word (decode (H6, T6, r4)), '100110');
r5 = parse ('001010');
same ('C6: syndrome of 001010', word (syndrome (H6, r5)), '101');
same ('C6: decoding of 001010', word (decode (H6, T6, r5)), '001111');
d5 = wt (bxor (repmat (r5, 8, 1), C6));
same ('C6: codewords at distance 2 from 001010', words (C6(d5 == 2, :)), '000000 001111 111010');

section ('Table sizes');
same ('entries in the (7,4) syndrome table', 2^3, 8);
same ('bits in the (7,4) syndrome table', 2^3 * 7, 56);
same ('bits in a ROM addressed by the received word', 2^7 * 4, 512);
same ('bits in the whole standard array', 2^7 * 7, 896);
same ('entries for 72 bits with 8 check bits', 2^8, 256);
same ('entries in a table addressed by 72 bits', sprintf ('%.2e', 2^72), '4.72e+21');

section ('Syndrome hardware');
same ('inputs of s1', nums (find (H7(1, :))), '1 3 4 5');
same ('inputs of s2', nums (find (H7(2, :))), '1 2 4 6');
same ('inputs of s3', nums (find (H7(3, :))), '2 3 4 7');
gates = sum (wt (H7) - 1);
same ('two-input XOR gates in the syndrome trees', gates, 9);
same ('depth of each tree in gates', max (floor (log2 (wt (H7) - 1)) + 1), 2);
same ('ROM contents, addresses 0 to 7', words (T7), ...
      '0000000 0000001 0000010 0100000 0000100 0010000 1000000 0001000');
r = parse ('1001011');
same ('tree inputs for r = 1001011, s1', nums (r(find (H7(1, :)))), '1 0 1 0');
same ('tree inputs for r = 1001011, s2', nums (r(find (H7(2, :)))), '1 0 1 1');
same ('tree inputs for r = 1001011, s3', nums (r(find (H7(3, :)))), '0 0 1 1');
same ('first-level XORs, s1', sprintf ('%d %d', xor (r(1), r(3)), xor (r(4), r(5))), '1 1');
same ('first-level XORs, s2', sprintf ('%d %d', xor (r(1), r(2)), xor (r(4), r(6))), '1 0');
same ('first-level XORs, s3', sprintf ('%d %d', xor (r(2), r(3)), xor (r(4), r(7))), '0 0');
same ('AND gate taps, columns 1 to 7', words (cols_of (H7)), '110 011 101 111 100 010 001');
ok = true;
CH = cols_of (H7);
for a = 1:8
  for j = 1:7
    fires = isequal (S3(a, :), CH(j, :));
    ok = ok && (fires == (T7(a, j) == 1));
  end
end
require ('seven AND gates on s or its complement reproduce the ROM', ok);
net = decoder_netlist (H7);
kinds = {net.kind};
outs = {net.out};
corr = cellfun (@(o) o(1) == 'c', outs);
same ('netlist: XOR gates in the syndrome trees', sum (strcmp (kinds, 'XOR') & ~corr), 9);
same ('netlist: inverters for the complements', sum (strcmp (kinds, 'NOT')), 3);
same ('netlist: three-input AND gates', sum (strcmp (kinds, 'AND') & arrayfun (@(x) numel (x.ins) == 3, net)), 7);
same ('netlist: XOR gates that apply the correction', sum (corr), 7);
same ('netlist: gates in the whole decoder, inverters aside', sum (~strcmp (kinds, 'NOT')), 23);
[~, depth] = run_netlist (net, zeros (1, 7));
same ('netlist: depth of s1, s2, s3 in gates', sprintf ('%d %d %d', depth('s1'), depth('s2'), depth('s3')), '2 2 2');
W = all_words (7);
ok = true;
for k = 1:rows (W)
  ok = ok && isequal (run_netlist (net, W(k, :)), decode (H7, T7, W(k, :)));
end
require ('netlist: the gates decode all 128 received words as the table does', ok);
same ('netlist: output for r = 1001011', word (run_netlist (net, r)), '1001001');

section ('Erasures');
c = parse ('1001001');
known = blank (c, [1 6]);
shown = char ('0' + c);
shown([1 6]) = '?';
same ('received word with positions 1 and 6 erased', shown, '?0010?1');
same ('known part u with zeros in the gaps', word (known), '0001001');
same ('syndrome of the known part, positions 1 and 6 erased', word (syndrome (H7, known)), '110');
CH = cols_of (H7);
same ('columns 4 and 7 of H added', word (bxor (CH(4, :), CH(7, :))), '110');
same ('columns 1 and 6 of H', words (CH([1 6], :)), '110 010');
aug = [CH(1, :)', CH(6, :)', syndrome(H7, known)'];
same ('augmented matrix [h1 h6 | s]', words (aug), '101 111 000');
aug(2, :) = bxor (aug(2, :), aug(1, :));
same ('after R2 <- R2 + R1', words (aug), '101 010 000');
same ('codewords that fit', words (erasure_solutions (H7, c, [1 6])), '1001001');
same ('columns 1, 2 and 3 added together', word (mod (CH(1, :) + CH(2, :) + CH(3, :), 2)), '000');
require ('1110000 is a codeword', ~isempty (index_in (C7, parse ('1110000'))));
k3 = blank (c, [1 2 3]);
same ('known part with positions 1, 2 and 3 erased', word (k3), '0001001');
same ('its syndrome', word (syndrome (H7, k3)), '110');
sols = erasure_solutions (H7, c, [1 2 3]);
same ('codewords that fit', words (sols), '0111001 1001001');
same ('their erased bits', words (sols(:, 1:3)), '011 100');
same ('the two candidates differ by', word (bxor (sols(1, :), sols(2, :))), '1110000');
cases = {7, H7, C7, '(7,4)', {'7 of 7', '21 of 21', '28 of 35', '0 of 35'}, ...
         '{1,2,3} {1,4,7} {1,5,6} {2,4,5} {2,6,7} {3,4,6} {3,5,7}'; ...
         6, H6, C6, 'C6', {'6 of 6', '15 of 15', '16 of 20', '0 of 15'}, ...
         '{1,3,6} {1,4,5} {2,3,4} {2,5,6}'};
for q = 1:2
  n = cases{q, 1};
  Hq = cases{q, 2};
  code = cases{q, 3};
  label = cases{q, 4};
  cl = cols_of (Hq);
  st = cases{q, 5};
  for e = 1:4
    sets = nchoosek (1:n, e);
    good = 0;
    for k = 1:rows (sets)
      good += independent (cl(sets(k, :), :));
    end
    same (sprintf ('%s: sets of %d erasures with independent columns', label, e), ...
          sprintf ('%d of %d', good, rows (sets)), st{e});
  end
  sets = nchoosek (1:n, 3);
  dep = zeros (0, 3);
  for k = 1:rows (sets)
    if (~independent (cl(sets(k, :), :)))
      dep(end + 1, :) = sets(k, :);
    end
  end
  heavy = code(wt (code) == 3, :);
  supports = zeros (rows (heavy), 3);
  for k = 1:rows (heavy)
    supports(k, :) = find (heavy(k, :));
  end
  supports = sortrows (supports);
  require (sprintf ('%s: the dependent triples are the supports of the weight-3 codewords', label), ...
           isequal (sortrows (dep), supports));
  txt = strjoin (arrayfun (@(k) sprintf ('{%d,%d,%d}', supports(k, :)), 1:rows (supports), ...
                           'UniformOutput', false), ' ');
  same (sprintf ('%s: supports of the weight-3 codewords', label), txt, cases{q, 6});
  ok = true;
  sets = nchoosek (1:n, 4);
  for k = 1:rows (code)
    for m = 1:rows (sets)
      ok = ok && rows (erasure_solutions (Hq, code(k, :), sets(m, :))) == 2;
    end
  end
  require (sprintf ('%s: every set of 4 erasures leaves exactly 2 candidates', label), ok);
end
ok = true;
for k = 1:16
  cw = C7(k, :);
  for e = 1:2
    sets = nchoosek (1:7, e);
    for m = 1:rows (sets)
      sol = erasure_solutions (H7, cw, sets(m, :));
      ok = ok && isequal (sol, cw);
    end
  end
end
require ('(7,4): every codeword is recovered from every 1 or 2 erasures', ok);
c2 = C7(bin2dec ('0110') + 1, :);
same ('codeword for message 0110', word (c2), '0110110');
k2 = blank (c2, [2 7]);
same ('received word 0?1011? with zeros in the gaps', word (k2), '0010110');
same ('syndrome of the known part of 0?1011?', word (syndrome (H7, k2)), '011');
same ('columns 3, 5 and 6 added', word (mod (CH(3, :) + CH(5, :) + CH(6, :), 2)), '011');
same ('codewords that fit 0?1011?', words (erasure_solutions (H7, c2, [2 7])), '0110110');
same ('columns 2, 3, 5 and 6 added', word (mod (CH(2, :) + CH(3, :) + CH(5, :) + CH(6, :), 2)), '000');

section ('Decoding failure');
a6 = leader_counts (L6, 6);
a7 = leader_counts (L7, 7);
same ('C6 leaders of weight 0, 1, 2', nums (a6(1:3)), '1 6 1');
same ('(7,4) leaders of weight 0, 1', nums (a7(1:2)), '1 7');
ps = [0.01 0.05 0.1];
pnames = {'0.01', '0.05', '0.1'};
stated = {'0.998636', '0.970299', '0.997969', '0.960596'; ...
          '0.969262', '0.857375', '0.955619', '0.814506'; ...
          '0.892296', '0.729000', '0.850306', '0.656100'};
for k = 1:3
  p = ps(k);
  same (sprintf ('p = %s: C6', pnames{k}), f6 (p_correct (a6, 6, p)), stated{k, 1});
  same (sprintf ('p = %s: three bits sent bare', pnames{k}), f6 ((1 - p) ^ 3), stated{k, 2});
  same (sprintf ('p = %s: (7,4) Hamming code', pnames{k}), f6 (p_correct (a7, 7, p)), stated{k, 3});
  same (sprintf ('p = %s: four bits sent bare', pnames{k}), f6 ((1 - p) ^ 4), stated{k, 4});
end
ok = true;
for k = 1:3
  ok = ok && abs (p_correct (a6, 6, ps(k)) - p_correct_by_enumeration (H6, T6, C6, 6, ps(k))) < 1e-12;
end
require ('C6 sums agree with a run over all 64 error patterns', ok);
ok = true;
for k = 1:3
  ok = ok && abs (p_correct (a7, 7, ps(k)) - p_correct_by_enumeration (H7, T7, C7, 7, ps(k))) < 1e-12;
end
require ('(7,4) sums agree with a run over all 128 error patterns', ok);
ok = true;
codes = {C6, H6, T6, 6; C7, H7, T7, 7};
for q = 1:2
  code = codes{q, 1};
  Hq = codes{q, 2};
  Tq = codes{q, 3};
  W = all_words (codes{q, 4});
  for k = 1:rows (code)
    cw = code(k, :);
    for m = 1:rows (W)
      e = W(m, :);
      right = isequal (decode (Hq, Tq, bxor (cw, e)), cw);
      ok = ok && (right == ~isempty (index_in (Tq, e)));
    end
  end
end
require ('the decoder is right exactly when the error is a leader', ok);
same ('(7,4) over four bare bits at p = 1/2', f6 (p_correct (a7, 7, 0.5) / 0.5 ^ 4), '1.000000');
same ('factored form (1-p)^6 (1+6p) at p = 0.1', f6 (0.9 ^ 6 * 1.6), '0.850306');
same ('p where the ratio (1-p)^2 (1+6p) peaks', f6 (2 / 9), '0.222222');
same ('largest ratio', f6 ((7 / 9) ^ 2 * (1 + 12 / 9)), '1.411523');
num = 7^2 * (9 + 12);
den = 9^2 * 9;
g = gcd (num, den);
same ('largest ratio before reduction', sprintf ('%d/%d', num, den), '1029/729');
same ('its factors (7/9)^2 and 1 + 12/9', sprintf ('%d/%d %d/9', 7^2, 9^2, 9 + 12), '49/81 21/9');
same ('largest ratio as a fraction in lowest terms', sprintf ('%d/%d', num / g, den / g), '343/243');
same ('root of 4 - 18p in lowest terms', sprintf ('%d/%d', 4 / gcd (4, 18), 18 / gcd (4, 18)), '2/9');
same ('percent gain at the peak', sprintf ('%.0f', ((7 / 9) ^ 2 * (1 + 12 / 9) - 1) * 100), '41');
ok = true;
for k = 1:499
  ok = ok && p_correct (a7, 7, k / 1000) > (1 - k / 1000) ^ 4;
end
require ('the (7,4) code beats four bare bits at p = k/1000 for k = 1 to 499', ok);
same ('(7,4) block error at p = 0.01', sprintf ('%.5f', 1 - p_correct (a7, 7, 0.01)), '0.00203');
same ('C6 at p = 1/2', f6 (p_correct (a6, 6, 0.5)), '0.125000');
same ('(7,4) at p = 1/2', f6 (p_correct (a7, 7, 0.5)), '0.062500');
p = 0.1;
W = all_words (6);
flag = 0;
for k = 1:64
  if (isequal (syndrome (H6, W(k, :)), [1 0 1]))
    flag += likelihood (p, 6, wt (W(k, :)));
  end
end
bounded = a6(1) * likelihood (p, 6, 0) + a6(2) * likelihood (p, 6, 1);
same ('no bit flips at p = 0.1', f6 (likelihood (p, 6, 0)), '0.531441');
same ('one of six bits flips at p = 0.1', f6 (6 * likelihood (p, 6, 1)), '0.354294');
c101 = W(all (syndrome (H6, W) == [1 0 1], 2), :);
same ('words of coset 101 of weight 2, 3 and 6', nums (arrayfun (@(w) sum (wt (c101) == w), [2 3 6])), '3 4 1');
same ('the pattern of weight 6 at p = 0.1', f6 (likelihood (p, 6, 6)), '0.000001');
same ('bounded decoder at p = 0.1: correct', f6 (bounded), '0.885735');
same ('bounded decoder at p = 0.1: failure flagged', f6 (flag), '0.022600');
same ('bounded decoder at p = 0.1: wrong codeword', f6 (1 - bounded - flag), '0.091665');
same ('complete decoder at p = 0.1: wrong codeword', f6 (1 - p_correct (a6, 6, p)), '0.107704');
same ('syndrome 101: right guesses', f6 (likelihood (p, 6, 2)), '0.006561');
same ('syndrome 101: wrong guesses', f6 (flag - likelihood (p, 6, 2)), '0.016039');
t = 0;
for k = 2:8
  t += likelihood (p, 6, wt (C6(k, :)));
end
same ('C6 at p = 0.1: error pattern is a nonzero codeword', f6 (t), '0.003159');

section ('Worked examples and exercises');
r = parse ('011011');
same ('C6: syndrome of 011011', word (syndrome (H6, r)), '111');
same ('C6: decoding of 011011', word (decode (H6, T6, r)), '010011');
same ('ratio p/(1-p) at p = 0.2', f6 (0.2 / 0.8), '0.250000');
same ('ratio p/(1-p) at p = 0.25', f6 (0.25 / 0.75), '0.333333');
R3 = [0 0 0; 1 1 1];
H3 = parse ('110 101');
[~, L3] = standard_array (R3, 3);
same ('repetition code: leaders', words (L3), '000 001 010 100');
same ('repetition code: their syndromes', words (syndrome (H3, L3)), '00 01 10 11');
same ('repetition code: P(correct) at p = 0.1', f6 (p_correct (leader_counts (L3, 3), 3, 0.1)), '0.972000');
same ('repetition code: its two terms at p = 0.1', sprintf ('%.3f %.3f', 0.9 ^ 3, 3 * 0.1 * 0.9 ^ 2), '0.729 0.243');
same ('(7,4) P(correct) at p = 0.02', f6 (p_correct (a7, 7, 0.02)), '0.992143');
same ('(7,4) block error at p = 0.02', f6 (1 - p_correct (a7, 7, 0.02)), '0.007857');
same ('four bits sent bare at p = 0.02', f6 (0.98 ^ 4), '0.922368');
same ('1 + 6p at p = 0.02', sprintf ('%.2f', 1 + 6 * 0.02), '1.12');
c = parse ('101001');
same ('C6: syndrome of the known part of 1??001', word (syndrome (H6, blank (c, [2 3]))), '111');
same ('C6: codewords that fit 1??001', words (erasure_solutions (H6, c, [2 3])), '101001');
same ('(7,4): columns 3, 4 and 5 added', word (mod (CH(3, :) + CH(4, :) + CH(5, :), 2)), '110');
require ('(7,4): columns 3, 4, 5 are independent', independent (CH([3 4 5], :)));
same ('(7,4): columns 2 and 4 added', word (bxor (CH(2, :), CH(4, :))), '100');
require ('(7,4): columns 2, 4, 5 are dependent', ~independent (CH([2 4 5], :)));
same ('(7,4): codeword on positions 2, 4, 5', word (C7(bin2dec ('0101') + 1, :)), '0101100');
same ('syndrome table bits, length 15 with 4 checks', 2^4 * 15, 240);
same ('ROM bits addressed by 15 received bits', 2^15 * 11, 360448);
same ('C6 P(correct) at p = 0.05 (worked example)', f6 (p_correct (leader_counts (L6, 6), 6, 0.05)), '0.969262');
same ('C6 at p = 0.05: the leader of weight 0', f6 (0.95 ^ 6), '0.735092');
same ('C6 at p = 0.05: the six leaders of weight 1', f6 (6 * 0.05 * 0.95 ^ 5), '0.232134');
same ('C6 at p = 0.05: the leader of weight 2', f6 (0.05 ^ 2 * 0.95 ^ 4), '0.002036');
same ('0.95^5 and 0.05^2', sprintf ('%.6f %.4f', 0.95 ^ 5, 0.05 ^ 2), '0.773781 0.0025');
same ('three bits sent bare at p = 0.05: lost', f6 (1 - 0.95 ^ 3), '0.142625');
same ('C6 at p = 0.05: lost', f6 (1 - p_correct (leader_counts (L6, 6), 6, 0.05)), '0.030738');
same ('four bits sent bare at p = 0.02: lost', f6 (1 - 0.98 ^ 4), '0.077632');
same ('ratio of the two failure rates at p = 0.02', sprintf ('%.2f', (1 - 0.98 ^ 4) / (1 - p_correct (a7, 7, 0.02))), ...
      '9.88');
same ('ROM bits over table bits', sprintf ('%.1f', 2^15 * 11 / (2^4 * 15)), '1501.9');
same ('four bits sent bare at p = 0.01: block error', f6 (1 - 0.99 ^ 4), '0.039404');
k = (1:999) / 1000;
require ('(1-p)^2 (1+6p) - 1 = p (2p-1) (3p-4) at p = k/1000 for k = 1 to 999', ...
         all (abs ((1 - k) .^ 2 .* (1 + 6 * k) - 1 - k .* (2 * k - 1) .* (3 * k - 4)) < 1e-12));

printf ('\nAll chapter 18 numbers reproduced.\n');
