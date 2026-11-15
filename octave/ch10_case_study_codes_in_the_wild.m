1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

% Reproduce every number in chapter 10 of Error Correction Coding, "Case Study: Codes in the Wild".
%
% The Octave twin of python/ch10_case_study_codes_in_the_wild.py. It needs only
% core GNU Octave, with no pkg load and no Forge packages. From the repository
% root:
%
%     octave --no-window-system --norc --quiet octave/ch10_case_study_codes_in_the_wild.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. The codes are built
% from their definitions rather than typed in. The Hamming, Reed-Muller, Hsiao and
% 5G (32,K) codes are enumerated codeword by codeword, the free distance of the
% CCSDS convolutional code is found by a shortest-path search over its trellis,
% the CCSDS Reed-Solomon generator is multiplied out over GF(256), the compact
% disc's burst bound is checked against a frame-by-frame model of its delay
% lines, the QR version 1 symbol is laid out module by module, and the 5G
% lifting sizes and code block parameters follow the rules of TS 38.212. A run
% that reaches the final line has reproduced every number in the chapter.
%
% tools/compare-routes.sh runs this file and the Python script side by side and
% requires the two to print the same lines. The plotting functions of
% ch10_charts.py stay on the Python route.
%
% Two conventions run through the file. A Python int is carried as int64, so
% that py_str prints it as 5 rather than 5.0, and every arithmetic step is done
% in doubles first, because int64 division in Octave rounds where Python's //
% floors. A Python dict or set is carried in a small struct built by py_dict or
% py_set below, since py_str has no form for either.


% Helpers ------------------------------------------------------------------------

function check (label, computed, stated)
  text = strrep (stated, ',', '');
  while (~isempty (text) && text(end) == '%')
    text = text(1:end - 1);
  end
  value = str2double (text);
  computed = double (computed);
  if (stated(end) == '%')
    computed = 100 * computed;
  end
  printf ('  %-62s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  printf ('  %-62s %-16s %s\n', label, value_text (computed, false), value_text (stated, false));
  if (~py_equal (computed, stated))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, ...
            value_text (computed, true), value_text (stated, true));
    exit (1);
  end
end

function near (label, a, b, rel)
  % Two computations of one quantity that must agree to a relative tolerance.
  printf ('  %-62s %-16s %s\n', label, py_g (a), py_g (b));
  if (abs (a - b) > rel * abs (b))
    printf ('MISMATCH in ''%s'': %s against %s\n', label, py_str (a), py_str (b));
    exit (1);
  end
end

function section (title)
  printf ('\n%s\n', title);
end

function d = py_dict (keys, values)
  % A Python dict of int keys and int values, kept as two parallel rows in key order.
  d = struct ('dict_keys', int64 (keys(:).'), 'dict_values', int64 (values(:).'));
end

function s = py_set (values)
  % A Python set of small ints. For the one-element sets this script prints, the
  % text is the same whatever order Python iterates in, and sorting fixes one.
  s = struct ('set_items', int64 (unique (values(:).')));
end

function t = tup (varargin)
  % A Python tuple of the values given, in the form py_str prints.
  t = struct ('tuple', {varargin});
end

function t = int_tuple (values)
  % A Python tuple of ints from a numeric row.
  t = struct ('tuple', {num2cell(int64 (values(:).'))});
end

function s = value_text (x, quoted)
  % str() of a value, or repr() when quoted is true, with dicts and sets added to py_str.
  if (isstruct (x) && isfield (x, 'dict_keys'))
    parts = arrayfun (@(k, v) sprintf ('%d: %d', k, v), x.dict_keys, x.dict_values, 'UniformOutput', false);
    s = ['{' strjoin(parts, ', ') '}'];
  elseif (isstruct (x) && isfield (x, 'set_items'))
    parts = arrayfun (@(v) sprintf ('%d', v), x.set_items, 'UniformOutput', false);
    s = ['{' strjoin(parts, ', ') '}'];
  else
    s = py_str (x, quoted);
  end
end

function eq = py_equal (a, b)
  % Python's ==, for the kinds of value this script compares. Numbers compare by
  % value, so the float 4.0 equals the int 4 as it does in Python.
  if (isstruct (a) || isstruct (b))
    if (~(isstruct (a) && isstruct (b)))
      eq = false;
    elseif (isfield (a, 'tuple') && isfield (b, 'tuple'))
      eq = numel (a.tuple) == numel (b.tuple) && all (cellfun (@py_equal, a.tuple, b.tuple));
    elseif (isfield (a, 'dict_keys') && isfield (b, 'dict_keys'))
      [ka, ia] = sort (a.dict_keys);
      [kb_, ib] = sort (b.dict_keys);
      eq = isequal (ka, kb_) && isequal (a.dict_values(ia), b.dict_values(ib));
    elseif (isfield (a, 'set_items') && isfield (b, 'set_items'))
      eq = isequal (a.set_items, b.set_items);
    else
      eq = false;
    end
  elseif (iscell (a) || iscell (b))
    eq = iscell (a) && iscell (b) && numel (a) == numel (b) && all (cellfun (@py_equal, a, b));
  elseif (ischar (a) || ischar (b))
    eq = ischar (a) && ischar (b) && strcmp (a, b);
  else
    eq = numel (a) == numel (b) && all (double (a(:)) == double (b(:)));
  end
end

function total = py_sum (terms)
  % Python's built-in sum() of floats. Since Python 3.12 it carries a Neumaier
  % compensation term beside the running total, where Octave's sum() adds plainly.
  total = 0;
  c = 0;
  for x = terms(:).'
    t = total + x;
    if (abs (total) >= abs (x))
      c += (total - t) + x;
    else
      c += (x - t) + total;
    end
    total = t;
  end
  if (c ~= 0 && isfinite (c))
    total += c;
  end
end

function c = comb_float (n, k)
  % math.comb(n, k) as a double. Every partial product is itself a binomial
  % coefficient, so the value is exact while it stays below 2^53. Beyond that it
  % is close but not exact, which only reaches tail terms far below the printed
  % digits. nchoosek would print a warning there instead.
  c = 1;
  for i = 1:k
    c = c * (n - k + i) / i;
  end
end

function q = q_func (x)
  % The Q-function, the probability that a standard Gaussian exceeds x.
  q = 0.5 * erfc (x / sqrt (2.0));
end

function s = tail_at_least (n, t, p)
  % The probability of t or more errors among n bits of BSC(p).
  terms = zeros (1, n - t + 1);
  for k = t:n
    terms(k - t + 1) = comb_float (n, k) * p ^ k * (1 - p) ^ (n - k);
  end
  s = py_sum (terms);
end

function w = weight (v)
  % The number of ones in the binary form of each entry of v, as bin(v).count("1").
  w = reshape (sum (dec2bin (v(:)) == '1', 2), size (v));
end


% Parity in early memories ----------------------------------------------------------

function q = even_count_probability (n, p)
  % The probability that n bits of BSC(p) suffer an even number of errors, zero included.
  q = (1 + (1 - 2 * p) ^ n) / 2;
end

function [alarm, silent] = byte_parity_word (p)
  % A 72-bit word stored as eight bytes, each with its own parity bit: (alarm, silent corruption).
  all_even = even_count_probability (9, p) ^ 8;
  alarm = 1 - all_even;
  silent = all_even - (1 - p) ^ 72;
end

function q = secded_uncorrected (p)
  % The probability that a 72-bit SEC-DED word holds two or more errors, so is not corrected.
  q = tail_at_least (72, 2, p);
end

function check_parity ()
  section ('Parity in Early Memories');
  % Octave's 0b literals make integer-class values, so double () turns each into
  % an ordinary number before any arithmetic.
  byte = double (0b10110010);
  same ('ones in the byte 10110010', int64 (weight (byte)), int64 (4));
  same ('even-parity bit for 10110010', int64 (mod (weight (byte), 2)), int64 (0));
  same ('ones after flipping one bit', int64 (weight (bitxor (byte, double (0b00001000)))), int64 (5));
  check ('byte parity rate 8/9', 8 / 9, '0.889');
  check ('check overhead 8/64', 8 / 64, '12.5%');
  same ('bits in a 64-bit word stored with byte parity', int64 (8 * 9), int64 (72));
  % The even-count formula against direct enumeration of every 9-bit pattern.
  p = 0.01;
  terms = [];
  for e = 0:511
    if (mod (weight (e), 2) == 0)
      terms(end + 1) = p ^ weight (e) * (1 - p) ^ (9 - weight (e));
    end
  end
  direct = py_sum (terms);
  near ('even-count formula against enumeration, n = 9', even_count_probability (9, p), direct, 1e-12);
end


% Hamming's code -------------------------------------------------------------------

function word = hamming_encode (d3, d5, d6, d7)
  % Hamming's arrangement: check bits at positions 1, 2 and 4, data at 3, 5, 6 and 7.
  p1 = bitxor (bitxor (d3, d5), d7);
  p2 = bitxor (bitxor (d3, d6), d7);
  p4 = bitxor (bitxor (d5, d6), d7);
  word = [p1, p2, d3, p4, d5, d6, d7];
end

function s = hamming_syndrome (word)
  % Positions count from 1, as Hamming numbered them, and that is also the
  % Octave index of each bit.
  s = 0;
  for position = 1:numel (word)
    if (word(position))
      s = bitxor (s, position);
    end
  end
end

function q = uncoded_block (p, k)
  q = 1 - (1 - p) ^ k;
end

function q = hamming74_block (p)
  q = 1 - ((1 - p) ^ 7 + 7 * p * (1 - p) ^ 6);
end

function q = repetition3_block (p)
  bit = 3 * p ^ 2 * (1 - p) + p ^ 3;
  q = 1 - (1 - bit) ^ 4;
end

function check_hamming ()
  section ('Hamming''s Code at Bell Labs');
  word = hamming_encode (1, 0, 1, 1);
  same ('codeword for data 1011', sprintf ('%d', word), '0110011');
  bad = word;
  % Python's bad[4] counts from 0, so it is position 5, which Octave indexes as 5.
  bad(5) = bitxor (bad(5), 1);
  same ('received word after flipping position 5', sprintf ('%d', bad), '0110111');
  same ('syndrome of the received word', dec2bin (hamming_syndrome (bad), 3), '101');
  % itertools.product((0, 1), repeat=4) runs through 0000, 0001, ..., 1111, which
  % is the binary form of 0 to 15 read from the most significant bit.
  all_zero = true;
  for n = 0:15
    d = dec2bin (n, 4) - '0';
    all_zero = all_zero && hamming_syndrome (hamming_encode (d(1), d(2), d(3), d(4))) == 0;
  end
  same ('syndrome of every codeword is zero', all_zero, true);
  for pos = 1:7
    flipped = word;
    flipped(pos) = bitxor (flipped(pos), 1);
    if (hamming_syndrome (flipped) ~= pos)
      printf ('MISMATCH: flipping position %d gives syndrome %d\n', pos, hamming_syndrome (flipped));
      exit (1);
    end
  end
  same ('every single flip names its own position', true, true);
  counts = zeros (1, 8);   % counts(w + 1) holds how many codewords have weight w
  for n = 0:15
    d = dec2bin (n, 4) - '0';
    w = sum (hamming_encode (d(1), d(2), d(3), d(4)));
    counts(w + 1) += 1;
  end
  ws = find (counts) - 1;
  same ('weight distribution of the (7,4) code', py_dict (ws, counts(ws + 1)), py_dict ([0 3 4 7], [1 7 7 1]));
  check ('rate of the (7,4) code', 4 / 7, '0.571');
  check ('rate of three-fold repetition', 1 / 3, '0.333');
  p = 0.01;
  check ('uncoded 4 bits at p = 0.01', uncoded_block (p, 4), '0.0394');
  check ('Hamming (7,4) block error at p = 0.01', hamming74_block (p), '0.00203');
  check ('three-fold repetition of 4 bits at p = 0.01', repetition3_block (p), '0.00119');
  check ('uncoded over Hamming at p = 0.01', uncoded_block (p, 4) / hamming74_block (p), '19.4');
  % The block error of the (7,4) code against enumeration of all 128 error patterns
  % with a syndrome decoder. The total is a plain running sum, as in the Python.
  total = 0.0;
  for e = 0:127
    pattern = bitand (bitshift (e, -(0:6)), 1);
    received = bitxor (word, pattern);
    s = hamming_syndrome (received);
    if (s)
      received(s) = bitxor (received(s), 1);
    end
    if (~isequal (received, word))
      total += p ^ sum (pattern) * (1 - p) ^ (7 - sum (pattern));
    end
  end
  near ('(7,4) block error against syndrome decoding of all 128 patterns', hamming74_block (p), total, 1e-12);
  same ('data bits that r = 3 check bits protect', int64 (2 ^ 3 - 1 - 3), int64 (4));
  same ('data bits that r = 7 check bits protect', int64 (2 ^ 7 - 1 - 7), int64 (120));
  r = 1:9;
  same ('fewest check bits for 64 data bits', int64 (min (r(2 .^ r - 1 - r >= 64))), int64 (7));
end


% The Mariner Reed-Muller code ------------------------------------------------------

function s = MARINER_FIXED ()
  % A script's variables are not visible inside its functions, so each
  % module-level constant of the Python becomes a function of no arguments.
  s = '10001101110101000010010110011111';
end

function word = rm15_codeword (message)
  % Bit x of the codeword is a0 + a1 x1 + ... + a5 x5 mod 2, with x1 the low bit of x.
  % The 32 positions x are handled at once as a row, where the Python loops over them.
  a0 = bitand (message, 1);
  x = 0:31;
  bit = repmat (a0, 1, 32);
  for i = 0:4
    if (bitand (bitshift (message, -(i + 1)), 1))
      bit = bitxor (bit, bitand (bitshift (x, -i), 1));
    end
  end
  % The largest word is 2^32 - 1, well inside the 2^53 a double holds exactly.
  word = sum (bit .* 2 .^ x);
end

function q = rm15_word_error (p)
  % A hard-decision decoder that corrects up to seven errors fails on eight or more.
  q = tail_at_least (32, 8, p);
end

function q = repetition5_word_error (p)
  bit = tail_at_least (5, 3, p);
  q = 1 - (1 - bit) ^ 6;
end

function check_mariner ()
  section ('The Mariner Reed-Muller Code');
  words = arrayfun (@rm15_codeword, 0:63);
  same ('distinct codewords', int64 (numel (unique (words))), int64 (64));
  ws = weight (words);
  keys = unique (ws);
  same ('weight distribution of RM(1,5)', py_dict (keys, arrayfun (@(k) sum (ws == k), keys)), ...
        py_dict ([0 16 32], [1 62 1]));
  % Every pair of distinct codewords, as itertools.combinations(words, 2) gives them.
  pairs = nchoosek (1:64, 2);
  pair = unique (weight (bitxor (words(pairs(:, 1)), words(pairs(:, 2)))));
  same ('distances between distinct codewords', int64 (pair(:).'), int64 ([16 32]));
  same ('minimum distance d', int64 (min (pair)), int64 (16));
  same ('errors corrected, (d - 1) // 2', int64 (floor ((16 - 1) / 2)), int64 (7));
  check ('rate 6/32', 6 / 32, '0.1875');
  x1 = rm15_codeword (double (0b000010));
  x2 = rm15_codeword (double (0b000100));
  % format(x, "032b")[::-1] lists the bits from the low end, bit 0 first.
  b1 = fliplr (dec2bin (x1, 32));
  b2 = fliplr (dec2bin (x2, 32));
  same ('codeword for a1 = 1, first eight bits', b1(1:8), '01010101');
  same ('codeword for a2 = 1, first eight bits', b2(1:8), '00110011');
  same ('positions where they differ', int64 (weight (bitxor (x1, x2))), int64 (16));
  % Orthogonality in the plus-and-minus-one picture: the inner product of two
  % codewords of the 32 with a0 = 0 is 32 - 2 * distance, which is zero. The
  % 32 words are built once here, where the Python rebuilds them for each pair.
  even = arrayfun (@rm15_codeword, 2 * (0:31));
  pairs = nchoosek (1:32, 2);
  ortho = all (32 - 2 * weight (bitxor (even(pairs(:, 1)), even(pairs(:, 2)))) == 0);
  same ('the 32 words with a0 = 0 are mutually orthogonal', ortho, true);
  same ('length of the fixed comma-free vector', int64 (numel (MARINER_FIXED)), int64 (32));
  % The fixed vector is a pseudonoise sequence of period 31 from a five-stage
  % shift register, with one bit appended. Here s(n + 1) is the Python's s[n],
  % so the five bits before it, s[n - 1] down to s[n - 5], are s(n:-1:n - 4).
  s = MARINER_FIXED - '0';
  found = [];
  for taps = 1:31
    c = bitand (bitshift (taps, -(0:4)), 1);
    if (c(5) == 0)
      continue;
    end
    ok = true;
    for n = 5:30
      ok = ok && s(n + 1) == mod (sum (c .* s(n:-1:n - 4)), 2);
    end
    if (ok)
      found(end + 1) = taps;
    end
  end
  same ('a five-stage recurrence generates its first 31 bits', ~isempty (found), true);
  seq = s(1:5);
  c = bitand (bitshift (found(1), -(0:4)), 1);
  for k = 1:62
    seq(end + 1) = mod (sum (c .* seq(end:-1:end - 4)), 2);
  end
  % The Python's seq[t:t + 31] is seq(t + 1:t + 31). A slice that runs past the
  % end is shorter in Python and so unequal, and Octave would raise an error.
  period = 0;
  for t = 1:62
    if (t + 31 <= numel (seq) && isequal (seq(t + 1:t + 31), seq(1:31)))
      period = t;
      break;
    end
  end
  same ('period of that recurrence', int64 (period), int64 (31));
  same ('ones among the first 31 bits', int64 (sum (s(1:31))), int64 (16));
  check ('required ST_B/N_0 saved, 5.20 - 3.00', 5.20 - 3.00, '2.20');
  check ('power ratio of 2.20 dB', 10 ^ 0.22, '1.66');
  check ('uncoded BPSK bit error at 5.20 dB', q_func (sqrt (2 * 10 ^ 0.52)), '5e-3');
  check ('bit rate, bits per second', 16200, '16,200');
end


% Voyager and concatenated codes -------------------------------------------------------

function g = G1 ()
  % Octave has no 0o literal, so the octal generators are read with base2dec.
  g = base2dec ('171', 8);
end

function g = G2 ()
  g = base2dec ('133', 8);
end

function [nxt, out] = conv_step (state, bit, g1, g2, memory)
  % The Python's step() is nested inside conv_free_distance. Octave scripts have
  % no nested functions, so it takes the generators and memory as arguments.
  reg = bitor (bitshift (bit, memory), state);
  out = mod (weight (bitand (reg, g1)), 2) + mod (weight (bitand (reg, g2)), 2);
  nxt = bitshift (reg, -1);
end

function heap = heap_push (heap, w, s)
  heap(end + 1, :) = [w, s];
end

function [heap, w, s] = heap_pop (heap)
  % heapq.heappop returns the smallest (weight, state) pair, comparing weights
  % first and states on a tie. The heap here never holds more than a few
  % hundred pairs, so a scan for the smallest row returns the same pair.
  [~, order] = sortrows (heap);
  k = order(1);
  w = heap(k, 1);
  s = heap(k, 2);
  heap(k, :) = [];
end

function d = conv_free_distance (g1, g2, memory)
  % The least output weight of a path that leaves the zero state and returns to it.
  % best(s + 1) holds the best weight found for state s, and 2^30 stands for
  % the missing entry, as best.get(s, 1 << 30) does.
  best = repmat (2 ^ 30, 1, 2 ^ memory);
  [first, w0] = conv_step (0, 1, g1, g2, memory);
  best(first + 1) = w0;
  heap = [w0, first];
  while (~isempty (heap))
    [heap, w, s] = heap_pop (heap);
    if (s == 0)
      d = w;
      return;
    end
    if (w > best(s + 1))
      continue;
    end
    for bit = 0:1
      [nxt, out] = conv_step (s, bit, g1, g2, memory);
      if (w + out < best(nxt + 1))
        best(nxt + 1) = w + out;
        heap = heap_push (heap, w + out, nxt);
      end
    end
  end
  d = -1;
end

function [exp_t, log_t] = gf_tables (poly)
  % exp_t(i + 1) is alpha^i and log_t(x + 1) is the exponent of x, since Octave
  % indexes from 1 and both tables are indexed by values that start at 0.
  exp_t = zeros (1, 510);
  log_t = zeros (1, 256);
  x = 1;
  for i = 0:254
    exp_t(i + 1) = x;
    log_t(x + 1) = i;
    x = bitshift (x, 1);
    if (bitand (x, 256))
      x = bitxor (x, poly);
    end
  end
  for i = 255:509
    exp_t(i + 1) = exp_t(i - 255 + 1);
  end
end

function c = gf_mul (a, b, exp_t, log_t)
  if (a == 0 || b == 0)
    c = 0;
    return;
  end
  c = exp_t(log_t(a + 1) + log_t(b + 1) + 1);
end

function g = ccsds_rs_generator (e)
  % g(x) = prod over j from 128 - E to 127 + E of (x - alpha^(11 j)), lowest degree first.
  [exp_t, log_t] = gf_tables (double (0x187));
  g = 1;
  for j = 128 - e:127 + e
    root = exp_t(mod (11 * j, 255) + 1);
    nxt = zeros (1, numel (g) + 1);
    for i = 0:numel (g) - 1
      c = g(i + 1);
      nxt(i + 2) = bitxor (nxt(i + 2), c);
      nxt(i + 1) = bitxor (nxt(i + 1), gf_mul (c, root, exp_t, log_t));
    end
    g = nxt;
  end
end

function t = VOYAGER_TABLE ()
  % planet, factor in the table, expected rate, achieved rate, improvement.
  % Python's None for Jupiter's improvement is the empty [] here. The Python's
  % VOYAGER_AU is read by nothing, in that script or in ch10_charts.py, so it
  % has no twin.
  t = {'Jupiter', 1, 115200, 115200, []
       'Saturn', 4, 29000, 44800, '1.5'
       'Uranus', 13, 9000, 29900, '3.3'
       'Neptune', 36, 3200, 21600, '6.8'};
end

function check_voyager ()
  section ('Voyager and Concatenated Codes');
  same ('free distance of the (171, 133) code', int64 (conv_free_distance (G1, G2, 6)), int64 (10));
  same ('generators in binary', tup (dec2bin (G1, 7), dec2bin (G2, 7)), tup ('1111001', '1011011'));
  check ('rate of RS(255,223)', 223 / 255, '0.875');
  check ('concatenated rate 223/255 x 1/2', 223 / 255 / 2, '0.437');
  check ('Golay plus convolutional rate 1/2 x 1/2', 0.25, '0.25');
  same ('RS check symbols', int64 (255 - 223), int64 (32));
  same ('RS minimum distance n - k + 1', int64 (255 - 223 + 1), int64 (33));
  same ('RS symbol errors corrected', int64 (floor (32 / 2)), int64 (16));
  check ('RS check overhead 32/223', 32 / 223, '14.3%');
  check ('RS(255,239) rate', 239 / 255, '0.937');
  exp_t = gf_tables (double (0x187));
  % exp_t(i + 1) is alpha^i, so position j of exp_t(2:256) is alpha^j.
  order = find (exp_t(2:256) == 1, 1);
  same ('order of alpha in GF(256) built on 0x187', int64 (order), int64 (255));
  g = ccsds_rs_generator (16);
  same ('degree of the E = 16 generator', int64 (numel (g) - 1), int64 (32));
  same ('generator coefficients read the same both ways', isequal (g, fliplr (g)), true);
  same ('the E = 8 generator is palindromic too', isequal (ccsds_rs_generator (8), fliplr (ccsds_rs_generator (8))), true);
  check ('general bound 1/16!', 1 / factorial (16), '4.78e-14');
  check ('undetected word error, 1e-5 x 3e-14', 1e-5 * 3e-14, '3e-19');
  for row = [5, 80, 640; 8, 128, 1024].'
    depth = row(1);
    symbols = row(2);
    bits = row(3);
    same (sprintf ('burst of symbols cleared at depth %d', depth), int64 (16 * depth), int64 (symbols));
    same (sprintf ('the same burst in bits at depth %d', depth), int64 (16 * depth * 8), int64 (bits));
  end
end

function check_voyager_table ()
  section ('Voyager''s Table 6-1');
  table = VOYAGER_TABLE;
  jupiter = table{1, 3};
  exact = struct ('Saturn', '28,800', 'Uranus', '8,862', 'Neptune', '3,200');
  for k = 2:size (table, 1)
    [planet, factor, expected, achieved, improvement] = table{k, :};
    check (sprintf ('%s: 115,200 / %d', planet, factor), jupiter / factor, exact.(planet));
    near (sprintf ('%s: the table''s rounded expected rate', planet), expected, jupiter / factor, 0.02);
    check (sprintf ('%s: achieved over expected', planet), achieved / expected, improvement);
  end
  check ('(5.2/10)^2 as one over', (10 / 5.2) ^ 2, '3.70');
  check ('(5.2/19)^2 as one over', (19 / 5.2) ^ 2, '13.35');
  check ('(5.2/30)^2 as one over', (30 / 5.2) ^ 2, '33.28');
  check ('table''s 36 against 33.28', 36 / (30 / 5.2) ^ 2 - 1, '8%');
  % Python's / always gives a float, so these two print as 4.0 and 36.0 and
  % still equal the ints 4 and 36.
  same ('the table''s 4 is (10/5)^2 exactly', (10 / 5) ^ 2, int64 (4));
  same ('the table''s 36 is (30/5)^2 exactly', (30 / 5) ^ 2, int64 (36));
  check ('Uranus with Jupiter at 5 AU, (19/5)^2', (19 / 5) ^ 2, '14.4');
  check ('Neptune expected from the distances', 115200 / (30 / 5.2) ^ 2, '3,461');
  check ('Neptune improvement from the distances', 21600 / (115200 / (30 / 5.2) ^ 2), '6.2');
end


% The compact disc ---------------------------------------------------------------------

function r = CHANNEL_RATE ()
  r = 4321800;
end

function d = D_FRAMES ()
  d = 4;
end

function best = cd_max_hits (burst)
  % The most bytes of one C2 codeword a run of erased frames can reach, from the delay model.
  %
  % Byte i of the C2 codeword formed at frame f leaves in frame f + 4 i, so a burst that
  % erases frames s to s + B - 1 erases byte i exactly when f + 4 i falls in that run.
  % Every start frame f is a row and every byte i a column, so one matrix covers both loops.
  f = (-120:39).';
  i = 0:27;
  frame = f + D_FRAMES * i;
  hits = sum (0 <= frame & frame < burst, 2);
  best = max ([0; hits]);
end

function check_cd ()
  section ('The Compact Disc');
  same ('audio bits per second, 44,100 x 16 x 2', int64 (44100 * 16 * 2), int64 (1411200));
  same ('audio bytes per second', int64 (44100 * 2 * 2), int64 (176400));
  same ('bytes in six stereo samples', int64 (6 * 2 * 2), int64 (24));
  same ('channel bits per frame', int64 (24 + 3 + 14 + 3 + 32 * (14 + 3)), int64 (588));
  same ('frames per second, 4,321,800 / 588', int64 (floor (CHANNEL_RATE / 588)), int64 (7350));
  same ('4,321,800 divides evenly by 588', int64 (mod (CHANNEL_RATE, 588)), int64 (0));
  same ('sections per second, 7,350 / 98', int64 (floor (7350 / 98)), int64 (75));
  same ('audio bytes per second from the frame rate', int64 (24 * 7350), int64 (176400));
  check ('CIRC rate 24/32', 24 / 32, '0.75');
  check ('share of channel bits that carry audio, 192/588', 192 / 588, '32.7%');
  same ('longest C2 delay, 27 x 4 frames', int64 (27 * D_FRAMES), int64 (108));
  same ('minimum distance of RS(28,24) and RS(32,28)', int_tuple ([28 - 24 + 1, 32 - 28 + 1]), int_tuple ([5, 5]));
  for b = 0:40
    if (cd_max_hits (b) ~= ceil (b / D_FRAMES))
      printf ('MISMATCH: burst %d reaches %d bytes, not %d\n', b, cd_max_hits (b), ceil (b / D_FRAMES));
      exit (1);
    end
  end
  same ('delay model gives ceil(B/4) for every B up to 40', true, true);
  b = 0:40;
  longest = max (b(arrayfun (@cd_max_hits, b) <= 4));
  same ('longest burst of frames that C2 can fill', int64 (longest), int64 (16));
  check ('16 frames in milliseconds', 16 / 7350 * 1e3, '2.18');
  check ('track length at 1.2 m/s, mm', 16 / 7350 * 1.2 * 1e3, '2.6');
  check ('track length at 1.4 m/s, mm', 16 / 7350 * 1.4 * 1e3, '3.0');
  check ('bound with D = 2: frames', 2 * 4, '8');
  check ('bound with D = 2: milliseconds', 8 / 7350 * 1e3, '1.09');
  same ('bytes in a CD-ROM sector', int64 (98 * 24), int64 (2352));
end


% QR codes ----------------------------------------------------------------------------

function v = QR_LEVELS ()
  % The Python's two dicts keyed by level become rows in the order L, M, Q, H.
  v = 'LMQH';
end

function v = QR_DATA_CODEWORDS ()
  v = [19, 16, 13, 9];
end

function v = QR_CAPACITY ()
  % DENSO WAVE's version 1 table: numeric, alphanumeric, binary, kanji. One row per level.
  v = [41, 25, 17, 10
       34, 20, 14, 8
       27, 16, 11, 7
       17, 10, 7, 4];
end

function s = qr_side (version)
  s = 21 + 4 * (version - 1);
end

function m = qr_version1_data_modules ()
  % used(r + 1, c + 1) marks module (r, c). A logical matrix counts each module
  % once, as the Python's set of (r, c) pairs does.
  n = 21;
  used = false (n, n);
  for corner = [0, 0; 0, n - 8; n - 8, 0].'
    r0 = corner(1);
    c0 = corner(2);
    used(r0 + 1:r0 + 8, c0 + 1:c0 + 8) = true;
  end
  used(6 + 1, (8:n - 9) + 1) = true;
  used((8:n - 9) + 1, 6 + 1) = true;
  c = 0:8;
  c = c(c ~= 6);
  used(8 + 1, c + 1) = true;
  used(c + 1, 8 + 1) = true;
  used(8 + 1, (n - 8:n - 1) + 1) = true;
  used((n - 7:n - 1) + 1, 8 + 1) = true;
  used(4 * 1 + 9 + 1, 8 + 1) = true;
  m = n * n - nnz (used);
end

function caps = qr_capacity (data_codewords)
  % Characters per mode in a version 1 to 9 symbol, from the mode indicator and count lengths.
  bits = 8 * data_codewords;
  numeric_bits = bits - 4 - 10;
  groups = floor (numeric_bits / 10);
  rest = mod (numeric_bits, 10);
  if (rest >= 7)
    extra = 2;
  elseif (rest >= 4)
    extra = 1;
  else
    extra = 0;
  end
  numeric = 3 * groups + extra;
  alnum_bits = bits - 4 - 9;
  pairs = floor (alnum_bits / 11);
  rest = mod (alnum_bits, 11);
  alnum = 2 * pairs + (rest >= 6);
  byte = floor ((bits - 4 - 8) / 8);
  kanji = floor ((bits - 4 - 8) / 13);
  caps = [numeric, alnum, byte, kanji];
end

function check_qr ()
  section ('QR Codes');
  same ('version 1 side', int64 (qr_side (1)), int64 (21));
  same ('version 40 side', int64 (qr_side (40)), int64 (177));
  same ('version 10 side', int64 (qr_side (10)), int64 (57));
  same ('modules in version 1', int64 (21 * 21), int64 (441));
  same ('data modules left in version 1', int64 (qr_version1_data_modules ()), int64 (208));
  same ('codewords in version 1', int64 (floor (qr_version1_data_modules () / 8)), int64 (26));
  levels = QR_LEVELS;
  data_codewords = QR_DATA_CODEWORDS;
  capacity = QR_CAPACITY;
  stated_checks = [7, 10, 13, 17];
  corrected = zeros (1, 4);
  for k = 1:4
    level = levels(k);
    data = data_codewords(k);
    r = 26 - data;
    corrected(k) = floor (r / 2);
    same (sprintf ('level %s: check codewords', level), int64 (r), int64 (stated_checks(k)));
    same (sprintf ('level %s: capacities derived from %d data codewords', level, data), ...
          int_tuple (qr_capacity (data)), int_tuple (capacity(k, :)));
  end
  same ('codeword errors corrected, L M Q H', int_tuple (corrected), int_tuple ([3, 5, 6, 8]));
  shares = {'11.5%', '19.2%', '23.1%', '30.8%'};
  for k = 1:4
    check (sprintf ('level %s: share of 26 codewords corrected', levels(k)), corrected(k) / 26, shares{k});
  end
  check ('level M rate 16/26', 16 / 26, '0.615');
  same ('level M byte capacity, (128 - 12) // 8', int64 (floor ((128 - 12) / 8)), int64 (14));
  check ('DENSO WAVE example: 50 of 200', 50 / 200, '25%');
  same ('level H erasures filled', int64 (26 - 9), int64 (17));
end


% Memory codes after Hsiao ---------------------------------------------------------------

function ok = is_secded (columns)
  % Every column distinct and nonzero, and the sum of any two columns neither
  % zero nor equal to a single column. All pairs are tested at once.
  singles = unique (columns);
  if (numel (singles) ~= numel (columns) || any (singles == 0))
    ok = false;
    return;
  end
  pairs = nchoosek (1:numel (columns), 2);
  s = bitxor (columns(pairs(:, 1)), columns(pairs(:, 2)));
  ok = ~any (s == 0 | ismember (s, singles));
end

function v = subset_masks (n, k)
  % The k-element subsets of range(n), in itertools.combinations order, each as
  % the integer with those bits set. nchoosek lists the subsets in the same order.
  c = nchoosek (0:n - 1, k);
  v = sum (2 .^ c, 2).';
end

function h = hsiao_72_64 ()
  units = 2 .^ (0:7);
  threes = subset_masks (8, 3);
  five = double (0b00011111);
  fives = zeros (1, 8);
  for k = 0:7
    fives(k + 1) = bitand (bitor (bitshift (five, k), bitshift (five, -(8 - k))), 255);
  end
  h = [units, threes, fives];
end

function e = extended_hamming_72_64 ()
  % Seven Hamming rows from the lightest columns, then an all-ones overall-parity row.
  data = [];
  for w = [2, 3, 4]
    for c = subset_masks (7, w)
      if (numel (data) < 64)
        data(end + 1) = c;
      end
    end
  end
  columns = [data, 2 .^ (0:6)];
  e = [bitor(columns, 2 ^ 7), 2 ^ 7];
end

function s = hsiao_39_32 ()
  threes = subset_masks (7, 3);
  threes = threes(1:32);
  s = [2 .^ (0:6), threes];
end

function check_hsiao ()
  section ('Memory Codes after Hsiao');
  h = hsiao_72_64 ();
  same ('Hsiao (72,64): columns', int64 (numel (h)), int64 (72));
  same ('Hsiao (72,64): every column has odd weight', all (mod (weight (h), 2)), true);
  wh = weight (h);
  same ('Hsiao (72,64): columns of weight 1, 3, 5', int_tuple (arrayfun (@(w) sum (wh == w), [1, 3, 5])), ...
        int_tuple ([8, 56, 8]));
  same ('Hsiao (72,64): ones in the check matrix', int64 (sum (weight (h))), int64 (216));
  rows = arrayfun (@(r) sum (bitand (bitshift (h, -r), 1)), 0:7);
  same ('Hsiao (72,64): ones in each row', py_set (rows), py_set (27));
  same ('Hsiao (72,64): corrects one error and detects two', is_secded (h), true);
  e = extended_hamming_72_64 ();
  same ('extended Hamming (72,64): corrects one and detects two', is_secded (e), true);
  same ('extended Hamming: columns', int64 (numel (e)), int64 (72));
  % 0x7F keeps the seven Hamming rows, and e(1:64) is the Python's e[:64].
  low = weight (bitand (e(1:64), 127));
  same ('extended Hamming: data columns of weight 2, 3, 4 in the Hamming rows', ...
        int_tuple (arrayfun (@(w) sum (low == w), [2, 3, 4])), int_tuple ([21, 35, 8]));
  same ('extended Hamming: ones in the seven Hamming rows', int64 (sum (weight (bitand (e, 127)))), int64 (186));
  same ('extended Hamming: ones in all eight rows', int64 (sum (weight (e))), int64 (258));
  same ('ones Hsiao saves', int64 (258 - 216), int64 (42));
  check ('saving as a share', 42 / 258, '16%');
  same ('XOR tree depth for 27 inputs', int64 (ceil (log2 (27))), int64 (5));
  same ('XOR tree depth for 72 inputs', int64 (ceil (log2 (72))), int64 (7));
  s = hsiao_39_32 ();
  same ('Hsiao (39,32): corrects one and detects two', is_secded (s), true);
  same ('Hsiao (39,32): ones', int64 (sum (weight (s))), int64 (103));
  check ('rate of the (72,64) code', 64 / 72, '0.889');
  same ('double-error syndromes in the (72,64) code', int64 (nchoosek (72, 2)), int64 (2556));
end

function check_parity_against_secded ()
  section ('Worked example: byte parity against SEC-DED');
  p = 1e-6;
  [alarm, silent] = byte_parity_word (p);
  check ('byte parity: alarm per read', alarm, '7.2e-5');
  check ('byte parity: silent corruption per read', silent, '2.88e-10');
  near ('silent corruption against 8 x C(9,2) p^2', silent, 8 * nchoosek (9, 2) * p ^ 2, 1e-4);
  unc = secded_uncorrected (p);
  check ('SEC-DED: two or more errors per read', unc, '2.56e-9');
  near ('two or more errors against C(72,2) p^2', unc, nchoosek (72, 2) * p ^ 2, 1e-4);
  check ('ratio of the two stopping rates', alarm / unc, '2.8e4');
  check ('SEC-DED: three or more errors per read', tail_at_least (72, 3, p), '5.96e-14');
  same ('C(9,2) and 8 x C(9,2)', int_tuple ([nchoosek(9, 2), 8 * nchoosek(9, 2)]), int_tuple ([36, 288]));
end


% 5G New Radio --------------------------------------------------------------------------

function sets = LIFTING_SETS ()
  sets = {};
  for a = [2, 3, 5, 7, 9, 11, 13, 15]
    j = 0:9;
    z = a * 2 .^ j;
    sets{end + 1} = z(z <= 384);
  end
end

function z = LIFTING_SIZES ()
  sets = LIFTING_SETS;
  z = sort ([sets{:}]);
end

function v = KCB (bg)
  % The Python's dicts keyed by base graph 1 and 2 become rows indexed by bg.
  t = [8448, 3840];
  v = t(bg);
end

function v = N_OVER_Z (bg)
  t = [66, 50];
  v = t(bg);
end

function v = K_OVER_Z (bg)
  t = [22, 10];
  v = t(bg);
end

function v = RV_NUMERATORS (bg, rv)
  % Row bg, and rv counts from 0, so it reads column rv + 1.
  t = [0, 17, 33, 56
       0, 13, 25, 43];
  v = t(bg, rv + 1);
end

function n = crc_length (a)
  if (a > 3824)
    n = 24;
  else
    n = 16;
  end
end

function g = base_graph (a, r)
  if (a <= 292 || (a <= 3824 && r <= 0.67) || r <= 0.25)
    g = 2;
  else
    g = 1;
  end
end

function c = code_blocks (b, bg)
  if (b <= KCB (bg))
    c = 1;
  else
    c = ceil (b / (KCB (bg) - 24));
  end
end

function v = kb (b, bg)
  if (bg == 1)
    v = 22;
  elseif (b > 640)
    v = 10;
  elseif (b > 560)
    v = 9;
  elseif (b > 192)
    v = 8;
  else
    v = 6;
  end
end

function out = lifting (b, bg)
  % One code block (C = 1): the lifting size, the code block length, fillers and buffer length.
  sizes = LIFTING_SIZES;
  zc = min (sizes(kb (b, bg) * sizes >= b));
  k = K_OVER_Z (bg) * zc;
  n = N_OVER_Z (bg) * zc;
  out = struct ('Kb', kb (b, bg), 'Zc', zc, 'K', k, 'F', k - b, 'N', n);
end

function v = k0 (rv, bg, ncb, zc)
  v = floor (RV_NUMERATORS (bg, rv) * ncb / (N_OVER_Z (bg) * zc)) * zc;
end

function m = SMALL_BASIS ()
  % TS 38.212 Table 5.3.3.3-1 as the Python's text holds it: the row number i,
  % then M(i,0) to M(i,10).
  m = [0 1 1 0 0 0 0 0 0 0 0 1
       1 1 1 1 0 0 0 0 0 0 1 1
       2 1 0 0 1 0 0 1 0 1 1 1
       3 1 0 1 1 0 0 0 0 1 0 1
       4 1 1 1 1 0 0 0 1 0 0 1
       5 1 1 0 0 1 0 1 1 1 0 1
       6 1 0 1 0 1 0 1 0 1 1 1
       7 1 0 0 1 1 0 0 1 1 0 1
       8 1 1 0 1 1 0 0 1 0 1 1
       9 1 0 1 1 1 0 1 0 0 1 1
       10 1 0 1 0 0 1 1 1 0 1 1
       11 1 1 1 0 0 1 1 0 1 0 1
       12 1 0 0 1 0 1 0 1 1 1 1
       13 1 1 0 1 0 1 0 1 0 1 1
       14 1 0 0 0 1 1 0 1 0 0 1
       15 1 1 0 0 1 1 1 1 0 1 1
       16 1 1 1 0 1 1 1 0 0 1 0
       17 1 0 0 1 1 1 0 0 1 0 0
       18 1 1 0 1 1 1 1 1 0 0 0
       19 1 0 0 0 0 1 1 0 0 0 0
       20 1 0 1 0 0 0 1 0 0 0 1
       21 1 1 0 1 0 0 0 0 0 1 1
       22 1 0 0 0 1 0 0 1 1 0 1
       23 1 1 1 0 1 0 0 0 1 1 1
       24 1 1 1 1 1 0 1 1 1 1 0
       25 1 1 0 0 0 1 1 1 0 0 1
       26 1 0 1 1 0 1 0 0 1 1 0
       27 1 1 1 1 0 1 0 1 1 1 0
       28 1 0 1 0 1 1 1 0 1 0 0
       29 1 0 1 1 1 1 1 1 1 0 0
       30 1 1 1 1 1 1 1 1 1 1 1
       31 1 0 0 0 0 0 0 0 0 0 0];
end

function rows = small_basis ()
  % TS 38.212 Table 5.3.3.3-1: row i + 1 holds M(i,0) to M(i,10), so the
  % Python's m[i][n] is rows(i + 1, n + 1). The row numbers are dropped.
  m = SMALL_BASIS;
  rows = m(:, 2:end);
end

function w = small_code_weights (m, k)
  % The weight of the codeword for every message a from 0 to 2^k - 1, one per
  % row. Bit n of a picks column n + 1 of the basis, and the 32 code bits are
  % the mod-2 sums. The Python loops over a and i, this does all at once.
  % bits(a + 1, n + 1) is bit n of a, found by floor and mod because those
  % broadcast a column of messages against a row of bit positions.
  a = (0:2 ^ k - 1).';
  bits = mod (floor (a ./ 2 .^ (0:k - 1)), 2);
  w = sum (mod (bits * m(:, 1:k).', 2), 2);
end

function best = small_code_distance (k)
  m = small_basis ();
  w = small_code_weights (m, k);
  best = min ([32; w(2:end)]);
end

function check_nr ()
  section ('5G New Radio');
  sizes = LIFTING_SIZES;
  same ('number of lifting sizes', int64 (numel (sizes)), int64 (51));
  same ('distinct lifting sizes', int64 (numel (unique (sizes))), int64 (51));
  same ('smallest and largest', int_tuple ([sizes(1), sizes(end)]), int_tuple ([2, 384]));
  check ('BG1 lowest rate 22/66', 22 / 66, '0.333');
  check ('BG2 lowest rate 10/50', 10 / 50, '0.2');
  same ('BG1 largest code block, 22 x 384', int64 (22 * 384), int64 (8448));
  same ('BG2 largest code block', int64 (KCB (2)), int64 (3840));
  a = 984;
  r = 0.5;
  same ('CRC bits for A = 984', int64 (crc_length (a)), int64 (16));
  b = a + crc_length (a);
  same ('B for A = 984', int64 (b), int64 (1000));
  same ('base graph for A = 984, R = 1/2', int64 (base_graph (a, r)), int64 (2));
  same ('code blocks', int64 (code_blocks (b, 2)), int64 (1));
  bg2 = lifting (b, 2);
  same ('BG2: Kb, Zc, K, F, N', int_tuple ([bg2.Kb, bg2.Zc, bg2.K, bg2.F, bg2.N]), ...
        int_tuple ([10, 104, 1040, 40, 5200]));
  same ('BG2: the lifting size below 104 is 96', int64 (max (sizes(sizes < 104))), int64 (96));
  starts = arrayfun (@(rv) k0 (rv, 2, bg2.N, bg2.Zc), 0:3);
  same ('BG2: k0 for rv 0 to 3', int_tuple (starts), int_tuple ([0, 1352, 2600, 4472]));
  same ('BG2: k0 in units of Zc', int_tuple (floor (starts / 104)), int_tuple ([0, 13, 25, 43]));
  bg1 = lifting (b, 1);
  same ('BG1 instead: Zc, K, N', int_tuple ([bg1.Zc, bg1.K, bg1.N]), int_tuple ([48, 1056, 3168]));
  same ('base graph for A = 984, R = 0.8', int64 (base_graph (984, 0.8)), int64 (1));
  same ('base graph for A = 200, R = 0.8', int64 (base_graph (200, 0.8)), int64 (2));
  ex = lifting (500, 2);
  same ('exercise B = 500, BG2: Kb, Zc, K, F', int_tuple ([ex.Kb, ex.Zc, ex.K, ex.F]), int_tuple ([8, 64, 640, 140]));
  same ('exercise: lifting sizes from 63 up', int64 (sizes(60 <= sizes & sizes <= 64)), int64 ([60, 64]));
  same ('segmentation of B = 20,000 with BG1', int64 (code_blocks (20000, 1)), int64 (3));
  same ('polar mother lengths 2^9 and 2^10', int_tuple ([2 ^ 9, 2 ^ 10]), int_tuple ([512, 1024]));
  m = small_basis ();
  same ('(32,K) basis: M(i,0) is all ones', all (m(:, 1) == 1), true);
  same ('(32,K) basis: M(i,1) to M(i,5) take all 32 patterns', int64 (size (unique (m(:, 2:6), 'rows'), 1)), int64 (32));
  distances = arrayfun (@small_code_distance, 3:11);
  same ('(32,K) minimum distance for K = 3 to 11', int_tuple (distances), int_tuple ([16, 16, 16, 16, 12, 12, 12, 12, 10]));
  w6 = small_code_weights (m, 6);
  keys = unique (w6);
  same ('(32,6) weight distribution', py_dict (keys, arrayfun (@(k) sum (w6 == k), keys)), py_dict ([0 16 32], [1 62 1]));
end


% The case summary table -------------------------------------------------------------------

function check_summary ()
  section ('The systems at a glance');
  check ('byte parity', 8 / 9, '0.889');
  check ('Hamming (7,4)', 4 / 7, '0.571');
  check ('Mariner (32,6)', 6 / 32, '0.1875');
  check ('Voyager concatenated', 223 / 510, '0.437');
  check ('compact disc CIRC', 24 / 32, '0.75');
  check ('QR version 1-M', 16 / 26, '0.615');
  check ('Hsiao (72,64)', 64 / 72, '0.889');
end

function out = impulse_response (g1, g2, memory)
  % The output pairs when a single 1 enters an empty register and zeros follow.
  out = '';
  state = 0;
  for bit = [1, zeros(1, memory)]
    reg = bitor (bitshift (bit, memory), state);
    out = [out, sprintf('%d%d', mod (weight (bitand (reg, g1)), 2), mod (weight (bitand (reg, g2)), 2))];
    state = bitshift (reg, -1);
  end
end

function check_prose_extras ()
  section ('Further numbers in the prose, examples and exercises');
  % Hamming's Code at Bell Labs and Memory Codes after Hsiao
  same ('data bits that six check bits protect, 2^6 - 1 - 6', int64 (2 ^ 6 - 1 - 6), int64 (57));
  same ('six check bits are too few for 64 data bits', 2 ^ 6 >= 64 + 6 + 1, false);
  % Voyager
  response = impulse_response (G1, G2, 6);
  same ('impulse response of the (171, 133) code', response, '11101111000111');
  same ('weight of the impulse response', int64 (sum (response == '1')), int64 (10));
  check ('Viterbi output bit error rate quoted by Ludwig and Taylor', 5e-3, '5e-3');
  check ('Reed-Solomon output bit error rate quoted by Ludwig and Taylor', 1e-6, '1e-6');
  % The Compact Disc
  same ('frames spanned by one C2 codeword, 4 x 27 + 1', int64 (4 * 27 + 1), int64 (109));
  % QR Codes
  same ('version 1 function modules, 441 - 208', int64 (21 * 21 - 208), int64 (233));
  same ('finder 147 + separators 45 + timing 10 + format 30 + dark module 1', int64 (147 + 45 + 10 + 30 + 1), int64 (233));
  % ISO/IEC 18004 holds p check codewords back from correction in version 1:
  % p = 3 at level L, 2 at M, 1 at Q and 1 at H, so t = (r - p) // 2.
  reserve = [3, 2, 1, 1];
  specified = floor ((26 - QR_DATA_CODEWORDS - reserve) / 2);
  same ('version 1 specified correction, (r - p) // 2 for L M Q H', int_tuple (specified), int_tuple ([2, 4, 6, 8]));
  same ('level L, (7 - 3) // 2', int64 (floor ((7 - 3) / 2)), int64 (2));
  same ('level M, (10 - 2) // 2', int64 (floor ((10 - 2) / 2)), int64 (4));
  check ('version 1-L specified correction 2/26', specified(1) / 26, '7.7%');
  check ('version 1-M specified correction 4/26', specified(2) / 26, '15.4%');
  check ('version 1-Q specified correction 6/26, below the nominal 25%', specified(3) / 26, '23.1%');
  % Worked example: byte parity against SEC-DED
  p = 1e-6;
  [~, silent] = byte_parity_word (p);
  check ('silent parity error over silent SEC-DED error', silent / tail_at_least (72, 3, p), '4.8e3');
  % Worked example: the Mariner code against five-fold repetition
  p = 0.05;
  check ('six bits sent bare at p = 0.05', 1 - 0.95 ^ 6, '0.2649');
  check ('five-fold repetition, one bit', tail_at_least (5, 3, p), '0.00116');
  near ('five-fold repetition, one bit, by the three-term formula', ...
        10 * p ^ 3 * (1 - p) ^ 2 + 5 * p ^ 4 * (1 - p) + p ^ 5, tail_at_least (5, 3, p), 1e-12);
  check ('five-fold repetition, six bits', repetition5_word_error (p), '0.00693');
  check ('Mariner code, eight or more errors in 32', rm15_word_error (p), '0.000139');
  check ('repetition over Mariner', repetition5_word_error (p) / rm15_word_error (p), '50');
  same ('extra bits the Mariner code spends over repetition', int64 (32 - 30), int64 (2));
  % Worked example: a scratch on a compact disc
  frame_mm = 1200 / 7350;
  check ('frame length on the track at 1.2 m/s, mm', frame_mm, '0.163');
  check ('frame lengths in a 1.5 mm speck', 1.5 / frame_mm, '9.19');
  same ('frames a 9.19-frame stretch can touch, ceil + 1', int64 (ceil (1.5 / frame_mm) + 1), int64 (11));
  same ('bytes of one C2 codeword lost to 11 frames', int64 (cd_max_hits (11)), int64 (3));
  check ('longest fillable dropout at 1.2 m/s, mm', 16 * frame_mm, '2.61');
  % Exercise: decoding a Hamming word
  received = [1, 1, 1, 0, 0, 1, 1];
  s = hamming_syndrome (received);
  same ('syndrome of 1110011', dec2bin (s, 3), '001');
  corrected = received;
  % The syndrome names a position counted from 1, which is the Octave index.
  corrected(s) = bitxor (corrected(s), 1);
  same ('corrected word', sprintf ('%d', corrected), '0110011');
  same ('data bits at positions 3, 5, 6, 7', sprintf ('%d', corrected([3, 5, 6, 7])), '1011');
  % Exercise: silent errors in a nine-bit group
  p = 1e-3;
  exact = even_count_probability (9, p) - (1 - p) ^ 9;
  check ('exact silent probability, p = 1e-3', exact, '3.57e-5');
  check ('approximation 36 p^2', 36 * p ^ 2, '3.6e-5');
  check ('approximation high by', 36 * p ^ 2 / exact - 1, '0.7%');
  check ('(1 - p)^7', (1 - p) ^ 7, '0.993');
  % Exercise: the Mariner code with erasures
  t = 0:7;
  same ('largest t with 2t + 5 < 16', int64 (max (t(2 * t + 5 < 16))), int64 (5));
  same ('errors given up to 5 erasures', int64 (7 - 5), int64 (2));
  % Exercise: a lighter Reed-Solomon code
  same ('RS(255,239) distance', int64 (255 - 239 + 1), int64 (17));
  same ('RS(255,239) errors corrected', int64 (floor ((255 - 239) / 2)), int64 (8));
  check ('1/8!', 1 / factorial (8), '2.48e-5');
  check ('1/8! over 1/16!', factorial (16) / factorial (8), '5e8');
  % Exercise: QR version 1 at levels L and H
  same ('level L: check codewords, corrected, erasures, bytes', ...
        int_tuple ([26 - 19, floor((26 - 19) / 2), 26 - 19, floor((19 * 8 - 12) / 8)]), int_tuple ([7, 3, 7, 17]));
  same ('level H: check codewords, corrected, erasures, bytes', ...
        int_tuple ([26 - 9, floor((26 - 9) / 2), 26 - 9, floor((9 * 8 - 12) / 8)]), int_tuple ([17, 8, 17, 7]));
  same ('bytes given up and codewords gained from L to H', int_tuple ([17 - 7, 8 - 3]), int_tuple ([10, 5]));
  % Exercise: a 5G block of 500 bits
  check ('500 / 8', 500 / 8, '62.5');
  block = lifting (500, 2);
  same ('B = 500: Kb, Zc, K, F, N', int_tuple ([block.Kb, block.Zc, block.K, block.F, block.N]), ...
        int_tuple ([8, 64, 640, 140, 3200]));
  same ('B = 500: start of rv2', int64 (k0 (2, 2, block.N, block.Zc)), int64 (1600));
end


function main ()
  printf ('  %-62s %-16s stated\n', 'quantity', 'computed');
  check_summary ();
  check_parity ();
  check_hamming ();
  check_mariner ();
  check_voyager ();
  check_voyager_table ();
  check_cd ();
  check_qr ();
  check_hsiao ();
  check_parity_against_secded ();
  check_nr ();
  check_prose_extras ();
  printf ('\nAll chapter 10 numbers reproduced.\n');
end

main ();
