% Reproduce every number in chapter 19 of Error Correction Coding,
% "Hamming Codes and Perfect Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch19_hamming_codes_and_perfect_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch19_hamming_codes_and_perfect_codes.py and
% prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in the chapter is a whole number, a binary or
% ternary word, a matrix or a code rate printed to three places, so every
% comparison is exact. The small codes are checked exhaustively: every one of
% the 2,048 messages of the (15,11) code is encoded and every single error on
% every codeword is corrected, the 16 balls of radius one around the (7,4)
% codewords are shown to cover all 128 words of length 7 exactly once, every
% double error on the extended (8,4) and (16,11) codes is shown to be
% detected, and the gate-level decoder is simulated on every correctable
% received word of both binary codes.
%
% A word is a row of 0s and 1s (or 0, 1, 2 over GF(3)), and position j of a
% word is its element j. Column j of the parity-check matrix is j written in
% m bits, most significant bit in the top row. A list of codewords is a
% matrix with one codeword per row, in the order Python's itertools.product
% gives. The extended codes keep the overall parity bit in front, as
% position 0, which is element 1 here.

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

function s = show (rows)
  parts = cell (1, rows (1) * 0 + size (rows, 1));
  for i = 1:size (rows, 1)
    parts{i} = word (rows(i, :));
  end
  s = strjoin (parts, ' ');
end

function s = join_numbers (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ' ');
end

function v = parse_word (text)
  v = double (text) - double ('0');
end

function w = weight (v)
  w = sum (v ~= 0);
end

function s = rate_text (k, n)
  s = sprintf ('%.3f', k / n);
end

function b = bit_length (x)
  b = 0;
  while (x > 0)
    x = floor (x / 2);
    b = b + 1;
  end
end


% Binary Hamming codes in counting order ------------------------------------------

function H = columns_H (n, m)
  H = zeros (m, n);
  for j = 1:n
    for i = 1:m
      H(i, j) = bitget (j, m - i + 1);
    end
  end
end

function H = hamming_H (m)
  H = columns_H (2 ^ m - 1, m);
end

function s = syndrome (H, v)
  s = mod (H * v(:), 2)';
end

function j = index_of (bits)
  j = 0;
  for b = bits
    j = 2 * j + b;
  end
end

function p = check_positions (n)
  p = [];
  for j = 1:n
    if (bitand (j, j - 1) == 0)
      p(end + 1) = j;
    end
  end
end

function p = data_positions (n)
  p = [];
  for j = 1:n
    if (bitand (j, j - 1) ~= 0)
      p(end + 1) = j;
    end
  end
end

function c = encode_h (n, data)
  c = zeros (1, n);
  c(data_positions (n)) = data;
  for p = check_positions (n)
    idx = find (bitand (1:n, p) & ((1:n) ~= p));
    c(p) = mod (sum (c(idx)), 2);
  end
end

function d = extract_data (n, c)
  d = c(data_positions (n));
end

function out = decode_h (H, v)
  j = index_of (syndrome (H, v));
  out = v;
  if (j)
    out(j) = 1 - out(j);
  end
end

function out = flip (v, positions)
  out = v;
  for p = positions
    out(p) = 1 - out(p);
  end
end

function C = all_codewords (n)
  k = numel (data_positions (n));
  C = zeros (2 ^ k, n);
  for idx = 0:2 ^ k - 1
    data = zeros (1, k);
    for i = 1:k
      data(i) = bitget (idx, k - i + 1);
    end
    C(idx + 1, :) = encode_h (n, data);
  end
end

function v = volume (n, t, q)
  if (nargin < 3)
    q = 2;
  end
  v = 0;
  for i = 0:t
    v = v + nchoosek (n, i) * (q - 1) ^ i;
  end
end

function m = min_redundancy (k)
  m = 1;
  while (k + m + 1 > 2 ^ m)
    m = m + 1;
  end
end

function ok = decodes_all_single (H, C)
  n = size (C, 2);
  ok = true;
  for j = 1:n
    R = C;
    R(:, j) = 1 - R(:, j);
    S = mod (R * H', 2);
    idx = S * (2 .^ (size (H, 1) - 1:-1:0))';
    ok = ok && all (idx == j);
  end
end


% The ternary Hamming code of length 4 ------------------------------------------------

function s = syndrome_q (H, v, q)
  s = mod (H * v(:), q)';
end

function C = ternary_codewords (H3)
  C = [];
  for idx = 0:80
    v = [floor(idx / 27), mod(floor (idx / 9), 3), mod(floor (idx / 3), 3), mod(idx, 3)];
    if (all (syndrome_q (H3, v, 3) == 0))
      C(end + 1, :) = v;
    end
  end
end

function [pos, val, out] = ternary_decode (H3, v)
  s = syndrome_q (H3, v, 3);
  out = v;
  pos = 0;
  val = 0;
  if (all (s == 0))
    return;
  end
  for j = 1:4
    col = H3(:, j)';
    for a = [1 2]
      if (all (mod (a * col, 3) == s))
        out(j) = mod (out(j) - a, 3);
        pos = j;
        val = a;
        return;
      end
    end
  end
  error ('no column matches');
end

function cols = projective_columns (m, q)
  cols = [];
  for idx = 0:q ^ m - 1
    v = zeros (1, m);
    x = idx;
    for i = m:-1:1
      v(i) = mod (x, q);
      x = floor (x / q);
    end
    nz = v(v ~= 0);
    if (~isempty (nz) && nz(1) == 1)
      cols(end + 1, :) = v;
    end
  end
end


% Extended Hamming codes ----------------------------------------------------------

function x = extend_word (c)
  x = [mod(sum (c), 2), c];
end

function d = extended_decision (v, H)
  s = index_of (syndrome (H, v(2:end)));
  p = mod (sum (v), 2);
  if (s == 0 && p == 0)
    d = 'no error';
  elseif (p == 1 && s == 0)
    d = 'flip the parity bit';
  elseif (p == 1)
    d = sprintf ('flip position %d', s);
  else
    d = 'double error detected';
  end
end

function out = flip_ext (v, positions)
  out = v;
  for p = positions
    out(p + 1) = 1 - out(p + 1);
  end
end

function ok = names_all_single (H, E)
  % Every single flip on every row of E, parity bit included, is named.
  ok = true;
  for i = 1:size (E, 1)
    for j = 0:size (E, 2) - 1
      if (j == 0)
        expected = 'flip the parity bit';
      else
        expected = sprintf ('flip position %d', j);
      end
      ok = ok && strcmp (extended_decision (flip_ext (E(i, :), j), H), expected);
    end
  end
end

function ok = detects_all_double (H, E)
  % Every pair of positions flipped on every row of E lands on s ~= 0, p = 0.
  len = size (E, 2);
  m = size (H, 1);
  ok = true;
  for i = 0:len - 2
    for j = i + 1:len - 1
      R = E;
      R(:, i + 1) = 1 - R(:, i + 1);
      R(:, j + 1) = 1 - R(:, j + 1);
      S = mod (R(:, 2:end) * H', 2);
      p = mod (sum (R, 2), 2);
      ok = ok && all (any (S, 2) & p == 0);
    end
  end
end


% Gate-level decoder --------------------------------------------------------------

function [lines, out] = gate_decoder (H, v)
  [m, n] = size (H);
  s = mod (H * v(:), 2)';
  lines = zeros (1, n);
  for j = 1:n
    a = 1;
    for i = 1:m
      if (H(i, j))
        a = a * s(i);
      else
        a = a * (1 - s(i));
      end
    end
    lines(j) = a;
  end
  out = mod (v + lines, 2);
end

function [weights, gates, depth] = encoder_counts (n)
  cp = check_positions (n);
  dp = data_positions (n);
  weights = zeros (1, numel (cp));
  for i = 1:numel (cp)
    weights(i) = sum (bitand (dp, cp(i)) ~= 0);
  end
  gates = sum (weights - 1);
  depth = 0;
  for w = weights
    depth = max (depth, bit_length (w - 1));
  end
end

function c = decoder_counts (H)
  [m, n] = size (H);
  weights = sum (H, 2)';
  c.syndrome = sum (weights - 1);
  c.depth = 0;
  for w = weights
    c.depth = max (c.depth, bit_length (w - 1));
  end
  c.and_gates = n;
  c.and2 = n * (m - 1);
  c.not_gates = m;
  c.xor_gates = n;
end


% Checks -------------------------------------------------------------------------

function check_counting ()
  section ('How many check bits one error needs');
  same ('syndromes of 3 check bits', 2 ^ 3, 8);
  same ('outcomes to tell apart for n = 7: no error plus 7 positions', 7 + 1, 8);
  ks = [1 4 8 11 16 26 32 57 64 120];
  ms = zeros (size (ks));
  for i = 1:numel (ks)
    ms(i) = min_redundancy (ks(i));
  end
  same ('least m for k = 1 4 8 11 16 26 32 57 64 120', join_numbers (ms), '2 3 4 4 5 5 6 6 7 7');
  same ('the matching lengths n = k + m', join_numbers (ks + ms), '3 7 12 15 21 31 38 63 71 127');
  require ('k = 4 with m = 3 meets n + 1 = 2^m exactly', 4 + 3 + 1 == 2 ^ 3);
  require ('k = 8 with m = 3 fails, 12 > 8', 8 + 3 + 1 > 2 ^ 3);
end

function check_hamming_7_4 (H7, C7)
  section ('The (7,4) code in counting order');
  same ('H for m = 3', show (H7), '0001111 0110011 1010101');
  same ('check positions', join_numbers (check_positions (7)), '1 2 4');
  same ('data positions', join_numbers (data_positions (7)), '3 5 6 7');
  same ('codewords', size (C7, 1), 16);
  require ('every codeword has zero syndrome', all (all (mod (C7 * H7', 2) == 0)));
  w = sum (C7, 2);
  same ('minimum weight', min (w(w > 0)), 3);
  c = encode_h (7, parse_word ('1011'));
  same ('message 1011 encodes to', word (c), '0110011');
  rx = flip (c, 6);
  same ('error at position 6 gives', word (rx), '0110001');
  same ('its syndrome', word (syndrome (H7, rx)), '110');
  same ('which is the binary for', index_of (syndrome (H7, rx)), 6);
  require ('every single error on every codeword is corrected', decodes_all_single (H7, C7));
  all7 = dec2bin (0:127, 7) - '0';
  same ('cosets of the (7,4) code, distinct syndromes of the 128 words', size (unique (mod (all7 * H7', 2), 'rows'), 1), 8);
end

function check_distance (H7, C7, C15)
  section ('Dimension and distance three');
  for m = 2:5
    H = hamming_H (m);
    n = 2 ^ m - 1;
    cols = H';
    distinct = size (unique (cols, 'rows'), 1) == n;
    require (sprintf ('m = %d: no zero column and no two columns equal', m), all (any (cols, 2)) && distinct);
    require (sprintf ('m = %d: columns 1, 2 and 3 add to zero', m), all (mod (cols(1, :) + cols(2, :) + cols(3, :), 2) == 0));
    units = cols(2 .^ (0:m - 1), :);
    require (sprintf ('m = %d: the columns at 1, 2, 4, ... are the m unit vectors', m), isequal (sortrows (units), sortrows (eye (m))));
  end
  codes = {all_codewords(3), C7, C15};
  for m = 2:4
    code = codes{m - 1};
    n = 2 ^ m - 1;
    same (sprintf ('m = %d: log2 of the number of codewords, n - m', m), log2 (size (code, 1)), n - m);
    w = sum (code, 2);
    same (sprintf ('m = %d: minimum distance by enumeration', m), min (w(w > 0)), 3);
  end
  same ('weight-3 codewords of the (7,4) code', sum (sum (C7, 2) == 3), 7);
  same ('C(7,2) / 3', nchoosek (7, 2) / 3, 7);
  w3 = flip (zeros (1, 7), [1 2 3]);
  same ('the weight-3 word on columns 1, 2, 3', word (w3), '1110000');
  require ('and it is a codeword, with syndrome zero', ~any (syndrome (H7, w3)));
end

function check_hamming_15_11 (H15, C15, M15, CW15, RX15)
  section ('The (15,11) code by hand');
  same ('H for m = 4', show (H15), '000000011111111 000111100001111 011001100110011 101010101010101');
  same ('data positions', join_numbers (data_positions (15)), '3 5 6 7 9 10 11 12 13 14 15');
  w = encoder_counts (15);
  same ('data bits feeding c1 c2 c4 c8', join_numbers (w), '7 7 7 7');
  same ('message', word (M15), '10110011010');
  same ('check bits c1 c2 c4 c8', word (CW15([1 2 4 8])), '1101');
  same ('codeword', word (CW15), '111001110011010');
  same ('weight of the codeword', weight (CW15), 9);
  same ('received with position 11 flipped', word (RX15), '111001110001010');
  same ('syndrome', word (syndrome (H15, RX15)), '1011');
  same ('which names position', index_of (syndrome (H15, RX15)), 11);
  same ('decoded data', word (extract_data (15, decode_h (H15, RX15))), '10110011010');
  require ('all 2,048 codewords have zero syndrome', all (all (mod (C15 * H15', 2) == 0)));
  require ('all 30,720 single errors are corrected', decodes_all_single (H15, C15));
  same ('weight-3 codewords of the (15,11) code', sum (sum (C15, 2) == 3), 35);
end

function check_family (FAMILY)
  section ('The Hamming family');
  same ('lengths', join_numbers (FAMILY(:, 2)), '3 7 15 31 63 127');
  same ('dimensions', join_numbers (FAMILY(:, 3)), '1 4 11 26 57 120');
  parts = cell (1, size (FAMILY, 1));
  for i = 1:size (FAMILY, 1)
    parts{i} = rate_text (FAMILY(i, 3), FAMILY(i, 2));
  end
  same ('rates', strjoin (parts, ' '), '0.333 0.571 0.733 0.839 0.905 0.945');
  require ('2^k (n + 1) = 2^n for every member', all (2 .^ FAMILY(:, 3) .* (FAMILY(:, 2) + 1) == 2 .^ FAMILY(:, 2)));
  same ('the (3,1) code', show (all_codewords (3)), '000 111');
end

function check_sphere_packing (C7)
  section ('The sphere-packing bound');
  same ('words within distance 1 of a word of length 7', volume (7, 1), 8);
  same ('16 balls of 8', size (C7, 1) * volume (7, 1), 128);
  covered = [];
  for i = 1:size (C7, 1)
    c = C7(i, :);
    covered(end + 1) = c * (2 .^ (6:-1:0))';
    for j = 1:7
      covered(end + 1) = flip (c, j) * (2 .^ (6:-1:0))';
    end
  end
  require ('the 16 balls of radius 1 are disjoint', numel (unique (covered)) == numel (covered));
  require ('and cover all 128 words of length 7', numel (unique (covered)) == 128);
  same ('words within distance 1 of a word of length 15', volume (15, 1), 16);
  same ('2^11 times 16', 2 ^ 11 * 16, 2 ^ 15);
  same ('C(23,1) C(23,2) C(23,3)', sprintf ('%d %d %d', nchoosek (23, 1), nchoosek (23, 2), nchoosek (23, 3)), '23 253 1771');
  same ('1 + 23 + 253 + 1771', volume (23, 3), 2048);
  require ('2^12 times 2^11 = 2^23', 2 ^ 12 * volume (23, 3) == 2 ^ 23);
  ok = true;
  for n = 3:2:15
    ok = ok && (2 * volume (n, (n - 1) / 2) == 2 ^ n);
  end
  require ('repetition codes of length 3 to 15: 2 V(n, (n - 1)/2) = 2^n', ok);
  same ('V(5,2)', volume (5, 2), 16);
  same ('1 + 90 + 4005', volume (90, 2), 4096);
  require ('and 4096 = 2^12', volume (90, 2) == 2 ^ 12);
  same ('ternary Golay: 1 + 11*2 + 55*4', volume (11, 2, 3), 243);
  require ('3^6 times 243 = 3^11', 3 ^ 6 * volume (11, 2, 3) == 3 ^ 11);
end

function check_ternary (H3)
  section ('Hamming codes over GF(3)');
  same ('one-dimensional subspaces of GF(3)^2', (3 ^ 2 - 1) / (3 - 1), 4);
  same ('their representatives as columns', show (projective_columns (2, 3)), '01 10 11 12');
  same ('H', show (H3), '0111 1012');
  code = ternary_codewords (H3);
  same ('codewords', size (code, 1), 9);
  same ('the nine codewords', show (code), '0000 0111 0222 1012 1120 1201 2021 2102 2210');
  w = sum (code ~= 0, 2);
  require ('every nonzero codeword has weight 3', all (w(w > 0) == 3));
  same ('words within distance 1 of a ternary word of length 4', volume (4, 1, 3), 9);
  same ('9 balls of 9', 9 * volume (4, 1, 3), 3 ^ 4);
  rx = [1 2 1 1];
  same ('syndrome of 1211', word (syndrome_q (H3, rx, 3)), '11');
  [pos, val, out] = ternary_decode (H3, rx);
  same ('position, value, corrected word', sprintf ('%d %d %s', pos, val, word (out)), '3 1 1201');
  rx = [1 1 0 1];
  same ('syndrome of 1101', word (syndrome_q (H3, rx, 3)), '20');
  [pos, val, out] = ternary_decode (H3, rx);
  same ('position, value, corrected word', sprintf ('%d %d %s', pos, val, word (out)), '2 2 1201');
  ok = true;
  for i = 1:size (code, 1)
    for j = 1:4
      for a = [1 2]
        v = code(i, :);
        v(j) = mod (v(j) + a, 3);
        [~, ~, out] = ternary_decode (H3, v);
        ok = ok && isequal (out, code(i, :));
      end
    end
  end
  require ('all 72 single errors on the 9 codewords are corrected', ok);
  qm = [3 2; 3 3; 4 2; 5 2; 7 2];
  parts = cell (1, size (qm, 1));
  for i = 1:size (qm, 1)
    q = qm(i, 1);
    m = qm(i, 2);
    n = (q ^ m - 1) / (q - 1);
    parts{i} = sprintf ('(%d,%d)', n, n - m);
    require (sprintf ('q = %d, m = %d: q^k (1 + n(q - 1)) = q^n', q, m), q ^ (n - m) * volume (n, 1, q) == q ^ n);
  end
  same ('(n, k) for q, m = 3 2, 3 3, 4 2, 5 2, 7 2', strjoin (parts, ' '), '(4,2) (13,10) (5,3) (6,4) (8,6)');
  cols = projective_columns (3, 3);
  same ('columns of the ternary code with m = 3', size (cols, 1), 13);
  proportional = false;
  for i = 1:12
    for j = i + 1:13
      for a = [1 2]
        if (all (mod (a * cols(i, :), 3) == cols(j, :)))
          proportional = true;
        end
      end
    end
  end
  require ('no two of them are multiples of each other', ~proportional);
end

function check_extended (H7, H15, C7, C15, FAMILY, CW15)
  section ('The extended Hamming code');
  E8 = [mod(sum (C7, 2), 2), C7];
  require ('every extended (8,4) codeword has even weight', all (mod (sum (E8, 2), 2) == 0));
  w = sum (E8, 2);
  same ('minimum distance of the (8,4) code', min (w(w > 0)), 4);
  same ('weight-4 codewords of the (8,4) code', sum (w == 4), 14);
  E16 = [mod(sum (C15, 2), 2), C15];
  w = sum (E16, 2);
  same ('minimum distance of the (16,11) code', min (w(w > 0)), 4);
  parts = cell (1, size (FAMILY, 1));
  for i = 1:size (FAMILY, 1)
    parts{i} = sprintf ('(%d,%d)', FAMILY(i, 2) + 1, FAMILY(i, 3));
  end
  same ('extended lengths and dimensions', strjoin (parts, ' '), '(4,1) (8,4) (16,11) (32,26) (64,57) (128,120)');
  require ('(8,4): every single error is named, parity bit included', names_all_single (H7, E8));
  require ('(16,11): every single error is named, parity bit included', names_all_single (H15, E16));
  same ('double-error patterns of length 8', nchoosek (8, 2), 28);
  require ('(8,4): every double error on every codeword is detected', detects_all_double (H7, E8));
  same ('double-error patterns of length 16', nchoosek (16, 2), 120);
  require ('(16,11): every double error on every codeword is detected', detects_all_double (H15, E16));
  x = extend_word (CW15);
  same ('parity bit of 111001110011010', x(1), 1);
  same ('extended codeword c0 c1 ... c15', word (x), '1111001110011010');
  rx = flip_ext (x, [3 12]);
  same ('syndrome with positions 3 and 12 flipped', word (syndrome (H15, rx(2:end))), '1111');
  same ('overall parity of that word', mod (sum (rx), 2), 0);
  same ('decision', extended_decision (rx, H15), 'double error detected');
  wrong = decode_h (H15, rx(2:end));
  same ('positions where the plain decoder''s output differs from c', weight (mod (wrong + CW15, 2)), 3);
end

function check_hardware (H7, H15, C7, C15)
  section ('The encoder and decoder as gates');
  [w, g, d] = encoder_counts (7);
  same ('(7,4) encoder: data bits per check, gates, depth', sprintf ('%s %d %d', join_numbers (w), g, d), '3 3 3 6 2');
  [w, g, d] = encoder_counts (15);
  same ('(15,11) encoder: gates and depth', sprintf ('%d %d', g, d), '24 3');
  gates = [];
  depths = [];
  for m = 3:7
    [~, g, d] = encoder_counts (2 ^ m - 1);
    require (sprintf ('m = %d: gates equal m(2^(m - 1) - 2)', m), g == m * (2 ^ (m - 1) - 2));
    gates(end + 1) = g;
    depths(end + 1) = d;
  end
  same ('encoder gates for m = 3 to 7', join_numbers (gates), '6 24 70 180 434');
  same ('encoder depths for m = 3 to 7', join_numbers (depths), '2 3 4 5 6');
  Hs = {H7, H15};
  names = {'(7,4)', '(15,11)'};
  stated = {'9 2 7 14 3 7', '28 3 15 45 4 15'};
  for i = 1:2
    c = decoder_counts (Hs{i});
    same (sprintf ('%s decoder: syndrome XORs, depth, ANDs, AND2s, NOTs, XORs', names{i}), sprintf ('%d %d %d %d %d %d', c.syndrome, c.depth, c.and_gates, c.and2, c.not_gates, c.xor_gates), stated{i});
  end
  totals = [33 92];
  levels = [5 6];
  for i = 1:2
    c = decoder_counts (Hs{i});
    same (sprintf ('%s decoder: two-input gates in all', names{i}), c.syndrome + c.and2 + c.not_gates + c.xor_gates, totals(i));
    same (sprintf ('%s decoder: gate levels, tree + NOT + AND + XOR', names{i}), c.depth + 3, levels(i));
  end
  codes = {C7, C15};
  for i = 1:2
    H = Hs{i};
    code = codes{i};
    n = size (code, 2);
    ok = true;
    for idx = 1:size (code, 1)
      c = code(idx, :);
      [lines, out] = gate_decoder (H, c);
      ok = ok && sum (lines) == 0 && isequal (out, c);
      for j = 1:n
        [lines, out] = gate_decoder (H, flip (c, j));
        ok = ok && sum (lines) == 1 && lines(j) == 1 && isequal (out, c);
      end
    end
    require (sprintf ('gate-level decoder of length %d on every correctable word', n), ok);
  end
  same ('data-only decoder ANDs and XORs, (7,4) and (15,11)', sprintf ('%d %d', numel (data_positions (7)), numel (data_positions (15))), '4 11');
end

function check_examples (H7)
  section ('Worked examples');
  H12 = columns_H (12, 4);
  same ('H of the (12,8) code, columns 1 to 12', show (H12), '000000011111 000111100001 011001100110 101010101010');
  c = encode_h (12, parse_word ('11001010'));
  same ('byte 11001010 encodes to', word (c), '001110001010');
  require ('its syndrome is zero', all (syndrome (H12, c) == 0));
  code12 = all_codewords (12);
  w = sum (code12, 2);
  same ('(12,8) codewords and minimum distance', sprintf ('%d %d', size (code12, 1), min (w(w > 0))), '256 3');
  same ('(7,4): data read from 0110011', word (extract_data (7, decode_h (H7, parse_word ('0110001')))), '1011');
  same ('V(10,2)', volume (10, 2), 56);
  same ('floor(1024 / 56)', floor (1024 / 56), 18);
  same ('largest k with 2^k <= 18', max (find (2 .^ (0:10) <= 18)) - 1, 4);
  same ('five data bits in 8 bits: 32 times 9', 32 * volume (8, 1), 288);
  same ('five data bits in 9 bits: 32 times 10', 32 * volume (9, 1), 320);
  x = extend_word (encode_h (7, parse_word ('1011')));
  same ('extended codeword of 1011', word (x), '00110011');
  rx = flip_ext (x, [2 5]);
  same ('positions 2 and 5 flipped', word (rx), '00010111');
  same ('its syndrome and overall parity', sprintf ('%s %d', word (syndrome (H7, rx(2:end))), mod (sum (rx), 2)), '111 0');
  wrong = decode_h (H7, rx(2:end));
  same ('the (7,4) decoder alone outputs', word (wrong), '0010110');
  require ('which is a codeword', all (syndrome (H7, wrong) == 0));
  same ('carrying the data', word (extract_data (7, wrong)), '1110');
end

function check_exercises (H7, H15, CW15)
  section ('Exercises');
  same ('least m for k = 26 and k = 27', sprintf ('%d %d', min_redundancy (26), min_redundancy (27)), '5 6');
  same ('C(15,2) / 3', nchoosek (15, 2) / 3, 35);
  rx = flip (CW15, 6);
  same ('received word', word (rx), '111000110011010');
  same ('its syndrome', word (syndrome (H15, rx)), '0110');
  same ('H for m = 2', show (hamming_H (2)), '011 101');
  same ('C(90,2)', nchoosek (90, 2), 4005);
  n3 = (3 ^ 3 - 1) / (3 - 1);
  same ('ternary m = 3: n, k, 1 + 2n', sprintf ('%d %d %d', n3, n3 - 3, 1 + 2 * n3), '13 10 27');
  [~, g, d] = encoder_counts (31);
  same ('(31,26) encoder gates and depth', sprintf ('%d %d', g, d), '70 4');
  rx = parse_word ('01110011');
  same ('syndrome and parity of 01110011', sprintf ('%s %d', word (syndrome (H7, rx(2:end))), mod (sum (rx), 2)), '001 1');
  same ('decision', extended_decision (rx, H7), 'flip position 1');
end


% Main ---------------------------------------------------------------------------

H7 = hamming_H (3);
H15 = hamming_H (4);
C7 = all_codewords (7);
C15 = all_codewords (15);
FAMILY = zeros (6, 3);
for m = 2:7
  FAMILY(m - 1, :) = [m, 2 ^ m - 1, 2 ^ m - 1 - m];
end
M15 = parse_word ('10110011010');
CW15 = encode_h (15, M15);
RX15 = flip (CW15, 11);
H3 = [0 1 1 1; 1 0 1 2];

check_counting ();
check_hamming_7_4 (H7, C7);
check_distance (H7, C7, C15);
check_hamming_15_11 (H15, C15, M15, CW15, RX15);
check_family (FAMILY);
check_sphere_packing (C7);
check_ternary (H3);
check_extended (H7, H15, C7, C15, FAMILY, CW15);
check_hardware (H7, H15, C7, C15);
check_examples (H7);
check_exercises (H7, H15, CW15);
printf ('\nAll chapter 19 numbers reproduced.\n');
