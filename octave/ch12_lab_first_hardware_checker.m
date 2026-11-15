% Reproduce every number in chapter 12 of Error Correction Coding, "Lab: A First Hardware Checker".
%
% This is the GNU Octave twin of python/ch12_lab_first_hardware_checker.py. It needs
% core Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch12_lab_first_hardware_checker.m
%
% The file is the software model of the parity bench in hdl/parity/. It holds the
% generator, the channel and the checker as a few lines of integer arithmetic, the
% gate and depth counts of the recursive XOR tree in hdl/parity/xor_tree.sv, the
% closed forms for the alarm, undetected-error and coverage rates of single parity
% on a binary symmetric channel, and the Wilson interval, zero-failure bound and
% trial count that chapter 11 builds.
%
% Run, it replays the seeded streams that tb/parity/test_parity.py drives through
% the Verilator simulation, recounts everything the cocotb tests log, and prints
% each value computed here beside the value the chapter states. It stops with exit
% status 1 at the first pair that disagrees by more than half a unit in the last
% digit the chapter prints, and at the first count that differs at all. The two
% routes print the same lines, and tools/compare-routes.sh checks that against the
% Python. The gate counts and depths Yosys reported are stored below as measured
% by synth/synth_parity.py and checked against the formulas the chapter proves.
%
% Every random number is a draw of Python's random.Random(seed).random(), made bit
% for bit in Octave by random_streams below, the Mersenne Twister of octave/common's
% PyRandom stepped for many seeds at once. The fault-injection runs take hundreds of
% seeds, and stepping their generators together is what keeps this file fast.
%
% Differences from the Python that a side-by-side reader meets:
%   - Octave indexes from 1, so the bit the Python calls d[0] is d(1) here, and a
%     0-based quantity such as a bit position keeps its Python value.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave functions cannot see script variables, so the module constants of the
%     Python (SEED, P_VALUES, MEASURED and the rest) are functions that return them.
%   - Python's bsc_stream is a generator that yields one (data, mask) pair at a
%     time. Here bsc_stream returns the first so many pairs at once, for a list of
%     seeds and a list of crossover probabilities. The bench restarts both
%     generators for every p, so every p sees the same draws, and run_bsc draws
%     them once for all the p it is given.
%   - tree_gates and tree_depth fill in the Python's recursion as a table, from one
%     input up, so a width is worked out once rather than once per call.
%   - exhaustive builds every (data, mask) vector at once, in the Python's loop order.
%   - 2^63 is one past the largest int64, so that count is carried as decimal text.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's.
%   - Python's x ** 2 calls the C library's pow, while Octave's x .^ 2 on an array
%     multiplies x by itself, and the two can differ in the last bit. py_pow below
%     makes Octave call pow as Python does.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

% Constants ----------------------------------------------------------------------

function s = SEED ()
  s = 12;   % the channel generator, which the bench restarts for every p
end

function s = SOURCE_OFFSET ()
  s = 10^9;   % the data generator is seeded with SEED + SOURCE_OFFSET
end

function z = Z95 ()
  z = 1.96;
end

function p = P_VALUES ()
  p = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];
end

function n = WORDS_PER_P ()
  n = 100000;
end

function p = WAVE_P ()
  p = 0.2;
end

function n = WAVE_WORDS ()
  n = 12;
end

function n = STEP_NS ()
  n = 10;
end

function n = BATCH_DRAWS ()
  % How many uniform draws run_bsc holds in memory at once. This is an Octave
  % detail with no Python counterpart, and any value gives the same counts.
  n = 6e6;
end

function m = MEASURED ()
  % What the cocotb run logged for K = 7, seed 12, 100,000 words at each p:
  % p, damaged, flagged, slipped, coverage, Wilson low, Wilson high, formula, verdict.
  m = {
    0.01, 7640, 7379, 261, '0.96584', '0.96152', '0.96968', '0.96587', 'inside';
    0.02, 14829, 13806, 1023, '0.93101', '0.92682', '0.93498', '0.93345', 'inside';
    0.05, 33558, 28362, 5196, '0.84516', '0.84125', '0.84899', '0.84606', 'inside';
    0.1, 56906, 41303, 15603, '0.72581', '0.72213', '0.72946', '0.73062', 'OUTSIDE';
    0.2, 83373, 49271, 34102, '0.59097', '0.58763', '0.59430', '0.59071', 'inside';
    0.3, 94208, 49611, 44597, '0.52661', '0.52342', '0.52980', '0.53024', 'OUTSIDE';
    0.4, 98315, 49962, 48353, '0.50818', '0.50506', '0.51131', '0.50854', 'inside';
    0.5, 99641, 50215, 49426, '0.50396', '0.50085', '0.50706', '0.50196', 'inside';
  };
end

function w = WAVEFORM ()
  % The twelve words of the waveform window, as test_waveform logs them:
  % time in ns, data, error mask, received word, error flag.
  w = {
    0, '0000101', '00011000', '00011101', 0;
    10, '0101010', '01010000', '11111010', 0;
    20, '1101100', '10000101', '11101001', 1;
    30, '1011110', '01001001', '10010111', 1;
    40, '1010000', '01101001', '00111001', 0;
    50, '0110001', '00100000', '10010001', 1;
    60, '0011010', '01000000', '11011010', 1;
    70, '0110101', '00010000', '00100101', 1;
    80, '0001110', '01010000', '11011110', 0;
    90, '0001110', '00010001', '10011111', 0;
    100, '1011110', '10010000', '01001110', 0;
    110, '0100001', '10000100', '10100101', 0;
  };
end

function y = YOSYS ()
  % What synth/synth_parity.py measured with Yosys 0.71:
  % block, n, width W of its tree, gates, XOR cells, XNOR cells, longest path in cells.
  y = {
    'parity_gen', 8, 7, 6, 0, 6, 3;
    'parity_gen', 16, 15, 14, 4, 10, 4;
    'parity_gen', 32, 31, 30, 8, 22, 5;
    'parity_gen', 64, 63, 62, 20, 42, 6;
    'parity_check', 8, 8, 7, 3, 4, 3;
    'parity_check', 16, 16, 15, 5, 10, 4;
    'parity_check', 32, 32, 31, 11, 20, 5;
    'parity_check', 64, 64, 63, 21, 42, 6;
  };
end


% Helpers ------------------------------------------------------------------------
% half_unit comes from octave/common and matches the Python's own half_unit.

function check (label, computed, stated)
  value = str2double (strrep (stated, ',', ''));
  if (stated(end) == '%')
    value = str2double (strrep (stated(1:end - 1), ',', ''));
    computed = 100 * double (computed);
  end
  computed = double (computed);
  printf ('  %-64s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % Python compares the two values with !=. Here the two are compared through
  % the text Python's str() prints for them, which is equal exactly when the
  % values are, for every kind of value this file passes.
  printf ('  %-64s %-16s %s\n', label, py_str (computed), py_str (stated));
  if (~strcmp (py_str (computed), py_str (stated)))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, ...
            py_str (computed, true), py_str (stated, true));
    exit (1);
  end
end

function holds (label, condition)
  if (condition)
    printf ('  %-64s %s\n', label, 'yes');
  else
    printf ('  %-64s %s\n', label, 'NO');
    printf ('MISMATCH: ''%s'' does not hold\n', label);
    exit (1);
  end
end

function t = tuple (varargin)
  % A Python tuple, in the form py_str prints as one.
  t = struct ('tuple', {varargin});
end

function s = py_sum (x)
  % Python's built-in sum() over a row of floats. Since Python 3.12 it carries
  % a Neumaier compensation term c beside the running total f and adds c once
  % at the end, which Octave's sum does not.
  f = 0 + x(1);
  c = 0;
  for j = 2:numel (x)
    xj = x(j);
    t = f + xj;
    if (abs (f) >= abs (xj))
      c = c + ((f - t) + xj);
    else
      c = c + ((xj - t) + f);
    end
    f = t;
  end
  s = f;
  if (c ~= 0 && isfinite (c))
    s = f + c;
  end
end

function y = py_pow (x, e)
  % Python's x ** e for each entry of x. Octave's x .^ 2 multiplies x by
  % itself, but an exponent passed as an array makes Octave call the C
  % library's pow for every entry, which is what Python's ** calls.
  y = x .^ (e * ones (size (x)));
end

function s = pow2_text (e)
  % The decimal digits of 2^e, as Python's exact int prints them. A double
  % holds every integer exactly only below 2^53, and int64 stops at 2^63 - 1,
  % so the digits are doubled e times as a list, least significant first.
  % Doubling a digit gives at most 18, so each carry is 0 or 1 and a digit
  % plus its carry never passes 9.
  digits = 1;
  for i = 1:e
    digits = 2 * digits;
    carry = floor (digits / 10);
    digits = [mod(digits, 10), 0] + [0, carry];
    if (digits(end) == 0)
      digits(end) = [];
    end
  end
  s = char (fliplr (digits) + '0');
end


% Many generators at once --------------------------------------------------------
% random_streams(seeds, n) gives what PyRandom gives for each seed, the Mersenne
% Twister MT19937 that Python's random.Random runs, but it holds the 624-word
% state of every seed as one column of a matrix and steps all the columns
% together, so four hundred short runs cost about as much as one long one. The
% state is held as uint32, whose bit operations Octave runs natively. Octave's
% integer types saturate rather than wrap, so the one product the twist needs is
% a 0 or 1 times a constant, which can never overflow, and the seeding, whose
% products do pass 2^32, is done in doubles split into 16-bit halves, as
% PyRandom does. The same functions serve the chapter 11 twin.

function r = random_streams (seeds, n)
  % R(:, j) holds the first n values random.Random(seeds(j)).random() returns,
  % for whole-number seeds below 2^32.
  if (any (seeds(:) >= 2^32))
    error ('random_streams: every seed must be below 2^32');
  end
  mt = mt_seed (seeds);
  words = 2 * n;
  blocks = ceil (words / 624);
  w = zeros (624 * blocks, numel (seeds), 'uint32');
  for b = 1:blocks
    mt = mt_twist (mt);
    w((b - 1) * 624 + (1:624), :) = mt_temper (mt);
  end
  % Python's random(): 27 bits from one word and 26 from the next, as a
  % fraction of 2^53. Both parts and the sum are exact in a double.
  a = double (bitshift (w(1:2:words, :), -5));
  b = double (bitshift (w(2:2:words, :), -6));
  r = (a * 67108864 + b) / 9007199254740992;
end

function mt = mt_seed (seeds)
  % The state random.Random(seed) starts from, for each seed below 2^32. CPython
  % splits the seed into 32-bit words for init_by_array, so each of these seeds
  % is a key of one word (seed 0 too), and every seed runs the same steps on its
  % own row.
  base = zeros (1, 624);
  base(1) = 19650218;
  for i = 2:624
    prev = base(i - 1);
    base(i) = mod (mulmod (bitxor (prev, floor (prev / 2^30)), 1812433253) + (i - 1), 2^32);
  end
  key = double (seeds(:));
  mt = repmat (base, numel (key), 1);
  i = 2;
  for count = 1:624
    prev = mt(:, i - 1);
    mixed = mulmod (bitxor (prev, floor (prev / 2^30)), 1664525);
    mt(:, i) = mod (bitxor (mt(:, i), mixed) + key, 2^32);
    i = i + 1;
    if (i > 624)
      mt(:, 1) = mt(:, 624);
      i = 2;
    end
  end
  for count = 1:623
    prev = mt(:, i - 1);
    mixed = mulmod (bitxor (prev, floor (prev / 2^30)), 1566083941);
    % mod, unlike rem, brings a negative difference back into 0 to 2^32 - 1.
    mt(:, i) = mod (bitxor (mt(:, i), mixed) - (i - 1), 2^32);
    i = i + 1;
    if (i > 624)
      mt(:, 1) = mt(:, 624);
      i = 2;
    end
  end
  mt(:, 1) = 2^31;
  mt = uint32 (mt.');
end

function c = mulmod (a, b)
  % (a * b) mod 2^32 for words a and b, exactly, in doubles.
  ah = floor (a / 65536);
  al = a - ah * 65536;
  c = mod (mod (ah * b, 65536) * 65536 + al * b, 2^32);
end

function mt = mt_twist (mt)
  % One pass of MT19937's recurrence over all 624 rows of every column. Row k
  % is rebuilt from rows k, k+1 and k+397, counted around the ring, in four
  % slices so that each slice reads only rows a slice before it has rebuilt.
  upper = uint32 (2147483648);
  lower = uint32 (2147483647);
  matrix_a = uint32 (2567483615);
  a = 1:227;
  b = 228:454;
  c = 455:623;
  y = bitor (bitand (mt(a, :), upper), bitand (mt(a + 1, :), lower));
  mt(a, :) = bitxor (bitxor (mt(a + 397, :), bitshift (y, -1)), bitand (y, 1) * matrix_a);
  y = bitor (bitand (mt(b, :), upper), bitand (mt(b + 1, :), lower));
  mt(b, :) = bitxor (bitxor (mt(b - 227, :), bitshift (y, -1)), bitand (y, 1) * matrix_a);
  y = bitor (bitand (mt(c, :), upper), bitand (mt(c + 1, :), lower));
  mt(c, :) = bitxor (bitxor (mt(c - 227, :), bitshift (y, -1)), bitand (y, 1) * matrix_a);
  y = bitor (bitand (mt(624, :), upper), bitand (mt(1, :), lower));
  mt(624, :) = bitxor (bitxor (mt(397, :), bitshift (y, -1)), bitand (y, 1) * matrix_a);
end

function y = mt_temper (y)
  % MT19937's tempering. A left shift of a uint32 drops the bits pushed past
  % bit 31, as the C code's 32-bit arithmetic does.
  y = bitxor (y, bitshift (y, -11));
  y = bitxor (y, bitand (bitshift (y, 7), uint32 (2636928640)));    % 0x9d2c5680
  y = bitxor (y, bitand (bitshift (y, 15), uint32 (4022730752)));   % 0xefc60000
  y = bitxor (y, bitshift (y, -18));
end


% The bench in software ----------------------------------------------------------
% Words are whole numbers held in doubles, and every function below works on each
% entry of an array of them. Octave's bitand, bitor, bitxor and bitshift accept
% doubles that hold whole numbers below 2^53.

function w = weight (x)
  % The number of 1 bits in each entry of x.
  w = zeros (size (x));
  while (any (x(:) > 0))
    w = w + mod (x, 2);
    x = floor (x / 2);
  end
end

function c = codeword (data, k)
  % What parity_gen.sv sends: the k data bits with the even-parity bit on top, at bit k.
  c = bitor (data, bitshift (bitand (weight (data), 1), k));
end

function f = flag (word)
  % What parity_check.sv raises: 1 when the word has odd weight.
  f = bitand (weight (word), 1);
end

function f = or_root_flag (word, n)
  % A broken checker whose top gate is OR where it should be XOR.
  low = floor (n / 2);
  f = bitor (bitand (weight (bitand (word, bitshift (1, low) - 1)), 1), ...
             bitand (weight (bitshift (word, -low)), 1));
end

function g = tree_gates (w)
  % Gates in xor_tree.sv over w inputs, counted the way the module builds them.
  % A tree over one input has no gate, and a wider one joins a tree over half
  % the inputs, rounded down, to a tree over the rest with one gate. The table
  % runs that recursion upward, so table(v) holds the count for v inputs.
  table = zeros (1, max (w(:)));
  for v = 2:numel (table)
    table(v) = table(floor (v / 2)) + table(v - floor (v / 2)) + 1;
  end
  g = table(w);
end

function d = tree_depth (w)
  % Gates on the deepest input-to-output path of the same tree.
  table = zeros (1, max (w(:)));
  for v = 2:numel (table)
    table(v) = 1 + max (table(floor (v / 2)), table(v - floor (v / 2)));
  end
  d = table(w);
end

function [data, mask] = bsc_stream (p, k, seed, words)
  % The first words (data, mask) pairs test_parity.py drives, from the same two
  % generators. DATA(i, j) is word i from seed(j). MASK(i, j, m) is the error
  % mask of that word at crossover probability p(m), and every p reads the same
  % channel draws, as a bench restarted for each p does.
  if (nargin < 3)
    seed = SEED;
  end
  seed = seed(:).';
  source = random_streams (seed + SOURCE_OFFSET, k * words);
  channel = random_streams (seed, (k + 1) * words);
  % Draw i of a word sets bit i, so the draws of one word fill one column of a
  % k-by-words block, and the weights 2^0, 2^1, ... add the bits into the word.
  bits = reshape (source < 0.5, k, words, numel (seed));
  data = reshape (sum (bits .* (2 .^ (0:k - 1)).', 1), words, numel (seed));
  draws = reshape (channel, k + 1, words, numel (seed));
  mask = zeros (words, numel (seed), numel (p));
  for m = 1:numel (p)
    mask(:, :, m) = reshape (sum ((draws < p(m)) .* (2 .^ (0:k)).', 1), words, numel (seed));
  end
end

function [damaged, flagged, slipped] = run_bsc (p, k, words, seed)
  % Damaged, flagged and slipped counts for each crossover probability p(m),
  % in row m, and each seed(j), in column j. The seeds are run a batch at a
  % time so that the draws held at once stay near BATCH_DRAWS.
  if (nargin < 2)
    k = 7;
  end
  if (nargin < 3)
    words = WORDS_PER_P;
  end
  if (nargin < 4)
    seed = SEED;
  end
  damaged = zeros (numel (p), numel (seed));
  flagged = damaged;
  slipped = damaged;
  group = max (1, floor (BATCH_DRAWS / ((2 * k + 1) * words)));
  for first = 1:group:numel (seed)
    cols = first:min (numel (seed), first + group - 1);
    [data, mask] = bsc_stream (p, k, seed(cols), words);
    sent = codeword (data, k);
    for m = 1:numel (p)
      f = flag (bitxor (sent, mask(:, :, m)));
      hit = mask(:, :, m) ~= 0;
      flagged(m, cols) = sum (f, 1);
      damaged(m, cols) = sum (hit, 1);
      slipped(m, cols) = sum (hit & f == 0, 1);
    end
  end
end

function [vectors, flagged, slipped, wrong, first_wrong] = exhaustive (k, checker)
  % Every data word against every mask: vectors, flagged, slipped, checker
  % disagreements, and the first (data, mask) the checker gets wrong, or []
  % where the Python returns None. The vectors are listed data word by data
  % word, each with every mask in turn, the order of the Python's two loops.
  n = k + 1;
  data = kron ((0:2^k - 1).', ones (2^n, 1));
  mask = repmat ((0:2^n - 1).', 2^k, 1);
  received = bitxor (codeword (data, k), mask);
  f = flag (received);
  flagged = sum (f);
  slipped = sum (mask ~= 0 & f == 0);
  wrong = 0;
  first_wrong = [];
  if (nargin > 1)
    bad = checker (received, n) ~= f;
    wrong = sum (bad);
    i = find (bad, 1);
    if (~isempty (i))
      first_wrong = tuple (int64 (data(i)), int64 (mask(i)));
    end
  end
  vectors = 2^k * 2^n;
end


% Statistics and closed forms ----------------------------------------------------
% These work on each entry of their arguments, so one call can score many runs.

function [lo, hi] = wilson (k, n, z)
  if (nargin < 3)
    z = Z95;
  end
  ph = k ./ n;
  center = ph + z .* z ./ (2 .* n);
  half = z .* sqrt (ph .* (1 - ph) ./ n + z .* z ./ (4 .* n .* n));
  scale = 1 + z .* z ./ n;
  lo = (center - half) ./ scale;
  hi = (center + half) ./ scale;
end

function b = zero_failure_bound (n, confidence)
  if (nargin < 2)
    confidence = 0.95;
  end
  b = 1 - (1 - confidence)^(1 / n);
end

function t = trials_needed (p, r, z)
  if (nargin < 3)
    z = Z95;
  end
  t = ceil ((z / r)^2 * (1 - p) / p);
end

function a = p_alarm (p, n)
  % Probability that an odd number of the n bits flip, so the checker raises its flag.
  a = (1 - (1 - 2 * p)^n) / 2;
end

function d = p_damaged (p, n)
  d = 1 - (1 - p)^n;
end

function u = p_undetected (p, n)
  % Probability of an even, nonzero number of flips.
  u = (1 + (1 - 2 * p)^n) / 2 - (1 - p)^n;
end

function c = coverage (p, n)
  % Fraction of damaged words the checker flags.
  c = p_alarm (p, n) / p_damaged (p, n);
end

function z = z_score (caught, damaged, c)
  z = (caught ./ damaged - c) ./ sqrt (c .* (1 - c) ./ damaged);
end

function r = pearson (a, b)
  ma = py_sum (a) / numel (a);
  mb = py_sum (b) / numel (b);
  cov = py_sum ((a - ma) .* (b - mb));
  r = cov / sqrt (py_sum (py_pow (a - ma, 2)) * py_sum (py_pow (b - mb, 2)));
end


% Checks -------------------------------------------------------------------------

function check_tree ()
  disp ('The XOR tree');
  w = 1:1024;
  holds ('gates = W - 1 for every W from 1 to 1,024', all (tree_gates (w) == w - 1));
  holds ('depth = ceil(log2 W) for every W from 1 to 1,024', all (tree_depth (w) == ceil (log2 (w))));
  same ('generator tree for K = 7: gates', int64 (tree_gates (7)), int64 (6));
  same ('generator tree for K = 7: depth', int64 (tree_depth (7)), int64 (3));
  same ('checker tree for n = 8: gates', int64 (tree_gates (8)), int64 (7));
  same ('checker tree for n = 8: depth', int64 (tree_depth (8)), int64 (3));
  % Python's 13 // 2 floors, and so does floor (13 / 2) for these whole numbers.
  same ('W = 13 splits into', tuple (int64 (floor (13 / 2)), int64 (13 - floor (13 / 2))), ...
        tuple (int64 (6), int64 (7)));
  same ('W = 13: gates', int64 (tree_gates (13)), int64 (12));
  same ('W = 13: depth', int64 (tree_depth (13)), int64 (4));
  same ('W = 6: depth', int64 (tree_depth (6)), int64 (3));
  same ('W = 7: depth', int64 (tree_depth (7)), int64 (3));
  same ('n = 64 checker: gates', int64 (tree_gates (64)), int64 (63));
  same ('n = 64 checker: depth', int64 (tree_depth (64)), int64 (6));
  same ('W = 2, 3, 4: depth', tuple (int64 (tree_depth (2)), int64 (tree_depth (3)), int64 (tree_depth (4))), ...
        tuple (int64 (1), int64 (2), int64 (2)));
  % d(i + 1) holds bit i of 1101100, so d(1) is the Python's d[0].
  d = bitand (bitshift (bin2dec ('1101100'), -(0:6)), 1);
  same ('1101100: d6 down to d0', int64 (fliplr (d)), int64 ([1, 1, 0, 1, 1, 0, 0]));
  same ('low group: d1 xor d2, then d0 xor that', ...
        tuple (int64 (bitxor (d(2), d(3))), int64 (bitxor (bitxor (d(1), d(2)), d(3)))), ...
        tuple (int64 (1), int64 (1)));
  same ('high group: d3 xor d4, d5 xor d6, then both', ...
        tuple (int64 (bitxor (d(4), d(5))), int64 (bitxor (d(6), d(7))), ...
               int64 (bitxor (bitxor (bitxor (d(4), d(5)), d(6)), d(7)))), ...
        tuple (int64 (1), int64 (0), int64 (1)));
  root = 0;
  for i = 1:7
    root = bitxor (root, d(i));
  end
  same ('root: parity bit of 1101100', int64 (root), int64 (0));
  disp ('Yosys netlists against the formulas');
  yosys = YOSYS ();
  for i = 1:size (yosys, 1)
    [block, n, w, gates, xor_cells, xnor_cells, depth] = yosys{i, :};
    % The Python chains its comparisons, a == b == c meaning a == b and b == c.
    holds (sprintf ('%s n = %d: %d gates = W - 1, depth %d = ceil(log2 %d)', block, n, gates, depth, w), ...
           gates == w - 1 && w - 1 == tree_gates (w) && xor_cells + xnor_cells == gates ...
           && depth == ceil (log2 (w)) && ceil (log2 (w)) == tree_depth (w));
    holds (sprintf ('%s n = %d: %d XNOR cells, an even number', block, n, xnor_cells), mod (xnor_cells, 2) == 0);
  end
end

function check_exhaustive ()
  disp ('Every error pattern');
  cases = {7, 32768, 16384, 16256; 3, 128, 64, 56};
  for i = 1:size (cases, 1)
    [k, vectors, flagged, slipped] = cases{i, :};
    [v, f, s] = exhaustive (k);
    same (sprintf ('k = %d: vectors', k), int64 (v), int64 (vectors));
    same (sprintf ('k = %d: flagged', k), int64 (f), int64 (flagged));
    same (sprintf ('k = %d: slipped', k), int64 (s), int64 (slipped));
    same (sprintf ('k = %d: flagged per data word', k), int64 (floor (f / 2^k)), int64 (2^k));
    same (sprintf ('k = %d: slipped per data word', k), int64 (floor (s / 2^k)), int64 (2^k - 1));
  end
  [m, d] = meshgrid (0:255, 0:127);
  holds ('every k = 7 vector: flag = weight(mask) mod 2', ...
         all (flag (bitxor (codeword (d(:), 7), m(:))) == mod (weight (m(:)), 2)));
  holds ('every clean k = 7 and k = 3 word is unflagged', ...
         all (flag (codeword (0:127, 7)) == 0) && all (flag (codeword (0:7, 3)) == 0));
  check ('share of nonzero masks that slip, 127/255', 127 / 255, '0.498');
  same ('clean words sent by test_clean_words for k = 7', int64 (2^7), int64 (128));
  same ('vectors for 32-bit codewords, 2^(2k + 1) with k = 31', pow2_text (2 * 31 + 1), pow2_text (63));
  [v, ~, ~, wrong, first] = exhaustive (7, @or_root_flag);
  same ('OR at the root: vectors the mutant gets wrong', int64 (wrong), int64 (8192));
  same ('OR at the root: first failing (data, mask)', first, tuple (int64 (0), int64 (17)));
  same ('mask 17 in binary', dec2bin (17, 8), '00010001');
  w = 0:255;
  both_odd = mod (weight (bitand (w, 15)), 2) & mod (weight (bitshift (w, -4)), 2);
  same ('words of 8 bits with both halves odd', int64 (sum (both_odd)), int64 (64));
  check ('share of vectors the mutant gets wrong', wrong / v, '0.25');
end

function check_waveform ()
  disp ('The waveform window, p = 0.2, seed 12');
  [data, mask] = bsc_stream (WAVE_P, 7, SEED, WAVE_WORDS);
  caught = 0;
  slipped = 0;
  weights = [];
  wave = WAVEFORM ();
  for i = 1:size (wave, 1)
    [t, data_s, mask_s, recv_s, f_s] = wave{i, :};
    received = bitxor (codeword (data(i), 7), mask(i));
    same (sprintf ('t = %3d ns', t), ...
          tuple (dec2bin (data(i), 7), dec2bin (mask(i), 8), dec2bin (received, 8), int64 (flag (received))), ...
          tuple (data_s, mask_s, recv_s, int64 (f_s)));
    weights(end + 1) = weight (mask(i));
    caught = caught + flag (received);
    slipped = slipped + (mask(i) ~= 0 && flag (received) == 0);
  end
  same ('mask weights in order', int64 (weights), int64 ([2, 2, 3, 3, 4, 1, 1, 1, 2, 2, 2, 2]));
  same ('caught in the window', int64 (caught), int64 (5));
  same ('slipped in the window', int64 (slipped), int64 (7));
  same ('window length in ns', int64 (WAVE_WORDS * STEP_NS), int64 (120));
  first = codeword (bin2dec ('1101100'), 7);
  same ('t = 20 ns: codeword of 1101100', dec2bin (first, 8), '01101100');
  same ('t = 20 ns: weight of the received word', int64 (weight (bitxor (first, bin2dec ('10000101')))), int64 (5));
  i = 7:-1:0;
  same ('t = 20 ns: bits the mask flips', int64 (i(bitand (bitshift (bin2dec ('10000101'), -i), 1) == 1)), ...
        int64 ([7, 2, 0]));
  same ('times at which error is high', int64 ([wave{[wave{:, 5}] == 1, 1}]), int64 ([20, 30, 50, 60, 70]));
  slip = p_undetected (WAVE_P, 8);
  check ('probability a word slips at p = 0.2', slip, '0.341');
  check ('slips expected among the 12 words', WAVE_WORDS * slip, '4.1');
  terms = zeros (1, 6);
  for k = 7:12
    terms(k - 6) = nchoosek (12, k) * slip^k * (1 - slip)^(12 - k);
  end
  check ('seven or more slips among 12', py_sum (terms), '7%');
end

function check_fault_injection ()
  disp ('Fault injection, 100,000 words at each p, seed 12 (replaying the bench)');
  % The Python calls run_bsc once for each p. One call here draws the shared
  % stream once and counts it at all eight p.
  p_values = P_VALUES;
  [damaged, flagged, slipped] = run_bsc (p_values);
  measured = MEASURED ();
  for i = 1:numel (p_values)
    p = p_values(i);
    [d_s, f_s, s_s, cov_s, lo_s, hi_s, form_s, verdict_s] = measured{[measured{:, 1}] == p, 2:end};
    caught = damaged(i) - slipped(i);
    [lo, hi] = wilson (caught, damaged(i));
    c = coverage (p, 8);
    label = sprintf ('p = %s', py_str (p));
    same ([label ': damaged, flagged, slipped'], ...
          tuple (int64 (damaged(i)), int64 (flagged(i)), int64 (slipped(i))), ...
          tuple (int64 (d_s), int64 (f_s), int64 (s_s)));
    holds ([label ': every flagged word was damaged'], flagged(i) == caught);
    check ([label ': measured coverage'], caught / damaged(i), cov_s);
    check ([label ': Wilson low'], lo, lo_s);
    check ([label ': Wilson high'], hi, hi_s);
    check ([label ': formula'], c, form_s);
    if (lo <= c && c <= hi)
      verdict = 'inside';
    else
      verdict = 'OUTSIDE';
    end
    same ([label ': formula inside the interval'], verdict, verdict_s);
  end
  damaged_counts = damaged;
  holds ('damaged counts rise with p', all (damaged_counts(1:end - 1) < damaged_counts(2:end)));
  [~, low_mask] = bsc_stream (0.1, 7, SEED, WORDS_PER_P);
  [~, high_mask] = bsc_stream (0.3, 7, SEED, WORDS_PER_P);
  holds ('every p = 0.1 mask sits inside the p = 0.3 mask, word by word', ...
         all (bitand (low_mask, high_mask) == low_mask));
  same ('words sent through the channel in all', int64 (numel (P_VALUES) * WORDS_PER_P), int64 (800000));
  same ('simulated time at the end of the eight runs, ns', int64 (numel (P_VALUES) * WORDS_PER_P * STEP_NS), ...
        int64 (8000000));
  same ('seed of the data generator', int64 (SEED + SOURCE_OFFSET), int64 (1000000012));
  check ('misses expected among eight 95% intervals', 8 * 0.05, '0.4');
  check ('z at p = 0.1', z_score (41303, 56906, coverage (0.1, 8)), '-2.59');
  check ('z at p = 0.3', z_score (49611, 94208, coverage (0.3, 8)), '-2.23');
  check ('two or more misses of 8 if they were independent', 1 - 0.95^8 - 8 * 0.05 * 0.95^7, '0.057');
  check ('that is one chance in', 1 / (1 - 0.95^8 - 8 * 0.05 * 0.95^7), '17');
  disp ('How strongly the shared draws couple the runs, seeds 0 to 399, 5,000 words each');
  % Row m of zs_by_p holds the z of every seed at coupled(m), the Python's zs_by_p[p].
  coupled = [0.01, 0.02, 0.1, 0.3];
  [damaged, ~, slipped] = run_bsc (coupled, 7, 5000, 0:399);
  c = arrayfun (@(p) coverage (p, 8), coupled(:));
  zs_by_p = z_score (damaged - slipped, damaged, c);
  check ('correlation of z at p = 0.01 and p = 0.02', pearson (zs_by_p(1, :), zs_by_p(2, :)), '0.43');
  check ('correlation of z at p = 0.1 and p = 0.3', pearson (zs_by_p(3, :), zs_by_p(4, :)), '0.02');
  disp ('The same experiment at p = 0.1 over seeds 0 to 299, 20,000 words each');
  c = coverage (0.1, 8);
  [damaged, ~, slipped] = run_bsc (0.1, 7, 20000, 0:299);
  [lo, hi] = wilson (damaged - slipped, damaged);
  misses = sum (~(lo <= c & c <= hi));
  zs = z_score (damaged - slipped, damaged, c);
  mean = py_sum (zs) / numel (zs);
  sd = sqrt (py_sum (py_pow (zs - mean, 2)) / (numel (zs) - 1));
  same ('intervals that miss the formula', int64 (misses), int64 (12));
  check ('share that miss', misses / 300, '4%');
  [lo, hi] = wilson (misses, 300);
  check ('Wilson interval on 12 of 300: low', lo, '0.023');
  check ('Wilson interval on 12 of 300: high', hi, '0.069');
  check ('mean z over 300 seeds', mean, '-0.09');
  check ('standard deviation of z over 300 seeds', sd, '0.97');
end

function check_closed_forms ()
  disp ('Closed forms at n = 8, p = 0.01');
  n = 8;
  p = 0.01;
  check ('(1 - 2p)^n', (1 - 2 * p)^n, '0.850763');
  check ('probability of an alarm', p_alarm (p, n), '0.074618');
  check ('probability the word is damaged', p_damaged (p, n), '0.077255');
  check ('probability of an undetected error', p_undetected (p, n), '0.0026368');
  check ('alarm + undetected = damaged', p_alarm (p, n) + p_undetected (p, n), '0.077255');
  check ('coverage', coverage (p, n), '0.96587');
  check ('approximation 1 - (n - 1)p/2', 1 - (n - 1) * p / 2, '0.965');
  check ('coverage at p = 0.5, 128/255', coverage (0.5, n), '0.50196');
  holds ('coverage at p = 0.5 is exactly 128/255', abs (coverage (0.5, n) - 128 / 255) < 1e-15);
  check ('expected damaged words in 100,000 at p = 0.01', 100000 * p_damaged (p, n), '7,725.5');
  check ('(1 - p)^n', (1 - p)^n, '0.922745');
  check ('second flip against first, (n - 1)p/2', (n - 1) * p / 2, '3.5%');
  two_or_more = 1 - (1 - p)^n - n * p * (1 - p)^(n - 1);
  check ('share of damaged words that carry a second flip', two_or_more / p_damaged (p, n), '3.5%');
  sd = sqrt (100000 * p_damaged (p, n) * (1 - p_damaged (p, n)));
  check ('standard deviation of the damaged count', sd, '84.4');
  check ('7,640 damaged in standard deviations from 7,725.5', (7640 - 100000 * p_damaged (p, n)) / sd, '-1.01');
  disp ('The Wilson interval for 7,379 of 7,640, step by step');
  c_hat = 7379 / 7640;
  D = 7640;
  z2 = Z95^2;
  check ('estimate', c_hat, '0.965838');
  check ('z^2', z2, '3.8416');
  check ('center before scaling', c_hat + z2 / (2 * D), '0.966089');
  check ('spread', Z95 * sqrt (c_hat * (1 - c_hat) / D + z2 / (4 * D * D)), '0.0040809');
  check ('scale', 1 + z2 / D, '1.000503');
  disp ('Precision of the estimates');
  check ('Wilson half-width scale at p = 0.01, 1.96 sqrt(c(1-c)/7,640)', ...
         Z95 * sqrt (coverage (p, n) * (1 - coverage (p, n)) / 7640), '0.0041');
  check ('measured undetected rate at p = 0.01, 261/100,000', 261 / 100000, '0.00261');
  [lo, hi] = wilson (261, 100000);
  check ('its Wilson low', lo, '0.00231');
  check ('its Wilson high', hi, '0.00295');
  holds ('the closed form 0.0026368 lies inside', lo <= p_undetected (p, n) && p_undetected (p, n) <= hi);
  P = p_undetected (p, n);
  check ('(1 - P)/P for the undetected rate', (1 - P) / P, '378.245');
  check ('(1.96/0.1)^2', (Z95 / 0.1)^2, '384.16');
  same ('words for 10% precision on the undetected rate', int64 (trials_needed (P, 0.1)), int64 (145307));
  check ('precision 100,000 words give, 1.96 sqrt((1 - P)/(N P))', Z95 * sqrt ((1 - P) / (100000 * P)), '12%');
  check ('extra words 10% precision needs beyond 100,000', trials_needed (P, 0.1) / 100000 - 1, '45%');
end

function check_exercise_answers ()
  disp ('Exercise answers');
  n = 8;
  p = 0.001;
  P = p_undetected (p, n);
  check ('undetected rate at p = 0.001', P, '2.78e-5');
  check ('expected slips in 20,000 words', 20000 * P, '0.56');
  check ('probability of no slip in 20,000 words', (1 - P)^20000, '0.57');
  check ('zero-failure bound at N = 20,000', zero_failure_bound (20000), '1.4978e-4');
  check ('bound over the true rate', zero_failure_bound (20000) / P, '5.4');
  m = 0:15;
  same ('odd masks among the 16 of length 4', int64 (sum (bitand (weight (m), 1))), int64 (8));
  m = 1:15;
  same ('even nonzero masks among the 16 of length 4', int64 (sum (mod (weight (m), 2) == 0)), int64 (7));
  check ('coverage at n = 4, p = 0.5, 8/15', coverage (0.5, 4), '0.5333');
  same ('three XNOR gates invert the parity, 3 mod 2', int64 (mod (3, 2)), int64 (1));
end

function main ()
  check_tree ();
  check_exhaustive ();
  check_waveform ();
  check_closed_forms ();
  check_fault_injection ();
  check_exercise_answers ();
  printf ('\nAll chapter 12 numbers reproduced.\n');
end

main ();
