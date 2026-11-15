% Reproduce every number in chapter 11 of Error Correction Coding, "Project: A Channel Simulator".
%
% This is the GNU Octave twin of python/ch11_project_channel_simulator.py. It needs
% core Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch11_project_channel_simulator.m
%
% The file is the chapter's channel simulator and its reference model at once. It
% holds the three channels (binary symmetric, binary erasure and BPSK over Gaussian
% noise, with the Gaussian samples made by the Box-Muller transform), the
% repetition and single parity check codes with their decoders, the Monte Carlo
% estimator with its Wald and Wilson intervals, the zero-failure bound, the number
% of trials a target precision needs, and the closed forms each simulation is
% checked against.
%
% Run, it repeats every seeded experiment the chapter reports and prints each
% value computed here beside the value the chapter states. It stops with exit
% status 1 at the first pair that disagrees by more than half a unit in the last
% digit the chapter prints, and at the first count that differs at all. The two
% routes print the same lines, and tools/compare-routes.sh checks that against
% the Python.
%
% Every random number is a draw of Python's random.Random(seed).random(), made
% bit for bit in Octave. Where the Python builds one generator and draws from it,
% this file uses PyRandom from octave/common. The simulation loop count_failures
% instead runs whole experiments at once, often a hundred or a thousand seeds in
% one call, so it steps the Mersenne Twisters of all those seeds together in
% random_streams below, which gives the same draws as PyRandom many times faster.
%
% Differences from the Python that a side-by-side reader meets:
%   - Octave indexes from 1, so a Python s[1] is s(2), and a 0-based quantity such
%     as a histogram bin keeps its Python value and gains 1 only to index.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave functions cannot see script variables, so the module constants of the
%     Python (Z95, SOURCE_OFFSET and the rest) are functions that return a value.
%   - Python passes a trial to count_failures as a function that runs one trial
%     at a time. Here a trial is a struct saying how many source and channel draws
%     one trial takes, with a function that runs every trial of an experiment at
%     once on a matrix of draws, one trial per row. A word is a row vector, so a
%     matrix of words holds one word per row.
%   - Octave has no None. An erased bit is NaN, and parity_fill reports a word it
%     cannot fill through a second output.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's.
%   - Python's x ** 2 calls the C library's pow, while Octave's x .^ 2 on an array
%     multiplies x by itself, and the two differ in the last bit for about one
%     value in a thousand. py_pow below makes Octave call pow as Python does.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

% Constants ----------------------------------------------------------------------

function z = Z95 ()
  z = 1.96;
end

function s = SOURCE_OFFSET ()
  s = 10^9;
end

function n = BATCH_DRAWS ()
  % How many uniform draws count_failures holds in memory at once. This is an
  % Octave detail with no Python counterpart, and any value gives the same counts.
  n = 6e6;
end


% Helpers ------------------------------------------------------------------------

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
  % Python's built-in sum() over each row of x, so a single list is a row.
  % Since Python 3.12 it carries a Neumaier compensation term c beside the
  % running total f and adds c once at the end, which Octave's sum does not.
  f = 0 + x(:, 1);
  c = zeros (size (f));
  if (rows (x) == 1)
    for j = 2:columns (x)
      xj = x(j);
      t = f + xj;
      if (abs (f) >= abs (xj))
        c = c + ((f - t) + xj);
      else
        c = c + ((xj - t) + f);
      end
      f = t;
    end
  else
    for j = 2:columns (x)
      xj = x(:, j);
      t = f + xj;
      big = abs (f) >= abs (xj);
      comp = (xj - t) + f;
      comp(big) = (f(big) - t(big)) + xj(big);
      c = c + comp;
      f = t;
    end
  end
  fix = c ~= 0 & isfinite (c);
  s = f;
  s(fix) = f(fix) + c(fix);
end

function y = py_pow (x, e)
  % Python's x ** e for each entry of x. Octave's x .^ 2 multiplies x by
  % itself, but an exponent passed as an array makes Octave call the C
  % library's pow for every entry, which is what Python's ** calls.
  y = x .^ (e * ones (size (x)));
end

function s = commas (n)
  % Python's f"{n:,}": a whole number with a comma between groups of three digits.
  s = sprintf ('%d', n);
  for i = numel (s) - 3:-3:1
    s = [s(1:i) ',' s(i + 1:end)];
  end
end

function s = word_text (word)
  % The text Python prints for a list of bits in which an erasure is None.
  % Octave holds an erased bit as NaN.
  items = cell (1, numel (word));
  for i = 1:numel (word)
    if (isnan (word(i)))
      items{i} = 'None';
    else
      items{i} = sprintf ('%d', word(i));
    end
  end
  s = ['[' strjoin(items, ', ') ']'];
end

function s = fill_text (word, ok)
  % The text Python prints for what parity_fill returns: the word, or None.
  if (ok)
    s = word_text (word);
  else
    s = 'None';
  end
end


% Many generators at once --------------------------------------------------------
% random_streams(seeds, n) gives what PyRandom gives for each seed, the Mersenne
% Twister MT19937 that Python's random.Random runs, but it holds the 624-word
% state of every seed as one column of a matrix and steps all the columns
% together, so a thousand short runs cost about as much as one long one. The
% state is held as uint32, whose bit operations Octave runs natively. Octave's
% integer types saturate rather than wrap, so the one product the twist needs is
% a 0 or 1 times a constant, which can never overflow, and the seeding, whose
% products do pass 2^32, is done in doubles split into 16-bit halves, as
% PyRandom does.

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
  % is a key of one word, and every seed runs the same steps on its own row.
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


% Channels -----------------------------------------------------------------------

function u = uniforms (rng, sz)
  % The uniform draws a channel compares with its threshold. rng is either a
  % PyRandom, drawn from row by row in the order Python's loop draws, or a
  % matrix of draws made ahead of time, one per bit.
  if (isa (rng, 'PyRandom'))
    u = reshape (rng.random (prod (sz)), sz(2), sz(1)).';
  else
    u = rng;
  end
end

function out = bsc (word, p, rng)
  % Binary symmetric channel: each bit flips when a uniform draw falls below p.
  out = double (xor (word, uniforms (rng, size (word)) < p));
end

function out = bec (word, eps, rng)
  % Binary erasure channel: each bit becomes NaN, an erasure, when a draw falls below eps.
  out = double (word);
  out(uniforms (rng, size (word)) < eps) = NaN;
end

function z = gaussian_source (u)
  % The standard normal samples the Python class GaussianSource hands out,
  % made from its generator's uniform draws u by the Box-Muller transform.
  % Each pair of draws gives two independent samples, r cos(theta) first and
  % r sin(theta), the spare the class keeps for its next call, second. Each
  % column of u is one generator's draws in order.
  u1 = 1.0 - u(1:2:end, :);   % in (0, 1], so log(u1) is finite
  u2 = u(2:2:end, :);         % in [0, 1)
  r = sqrt (-2.0 * log (u1));
  theta = 2.0 * pi * u2;
  z = zeros (2 * rows (u1), columns (u1));
  z(1:2:end, :) = r .* cos (theta);
  z(2:2:end, :) = r .* sin (theta);
end

function sigma = awgn_sigma (ebn0_db, rate)
  % Noise standard deviation for pulses of energy Es = 1 carrying `rate` data bits each.
  ebn0 = 10 ^ (ebn0_db / 10);
  sigma = sqrt (1.0 / (2.0 * rate * ebn0));
end

function y = bpsk_awgn (word, sigma, noise)
  % Send 0 as +1 and 1 as -1, then add Gaussian noise of standard deviation
  % sigma. noise holds one standard normal sample per bit.
  y = (1.0 - 2.0 * word) + sigma * noise;
end


% Codes --------------------------------------------------------------------------

function word = repetition_encode (bit, n)
  word = repmat (bit, 1, n);
end

function b = majority (word)
  b = double (sum (word, 2) > floor (columns (word) / 2));
end

function b = hard_repetition (received)
  % Decide each pulse by its sign, then take the majority.
  b = majority (double (received < 0));
end

function b = soft_repetition (received)
  % Add the received values and decide by the sign of the sum.
  b = double (py_sum (received) < 0);
end

function b = repetition_erasures (received)
  % On the erasure channel any copy that arrives is correct. NaN when every copy is erased.
  b = NaN (rows (received), 1);
  for j = columns (received):-1:1
    arrived = ~isnan (received(:, j));
    b(arrived) = received(arrived, j);
  end
end

function word = parity_encode (data)
  word = [data, mod(sum (data, 2), 2)];
end

function ok = parity_ok (word)
  ok = mod (sum (word, 2), 2) == 0;
end

function [word, ok] = parity_fill (received)
  % Fill a single erasure from the parity rule. ok is false where the Python
  % returns None, for a word with two or more positions erased.
  missing = isnan (received);
  ok = sum (missing, 2) <= 1;
  known = received;
  known(missing) = 0;
  word = received;
  [r, c] = find (missing & ok);
  word(sub2ind (size (word), r, c)) = mod (sum (known(r, :), 2), 2);
end


% The simulation loop ------------------------------------------------------------

function failures = count_failures (trials, seed, trial, gaussian)
  % Run `trial` `trials` times and count the runs that return true.
  %
  % The message bits come from their own generator, seeded with seed + 10^9, so
  % the channel sees the same uniform draws whatever the messages are. Every
  % draw an experiment needs is made first, one trial per row, and the trial
  % then runs on all the rows at once. seed may be a vector, which runs one
  % experiment per seed and returns one count each.
  if (nargin < 4)
    gaussian = false;
  end
  seeds = seed(:).';
  failures = zeros (1, numel (seeds));
  ks = trial.source_draws;
  kc = trial.channel_draws;
  batch = max (1, floor (BATCH_DRAWS / (trials * (ks + kc))));
  for first = 1:batch:numel (seeds)
    idx = first:min (first + batch - 1, numel (seeds));
    s = seeds(idx);
    source = reshape (random_streams (s + SOURCE_OFFSET, trials * ks), ks, []).';
    if (gaussian)
      % The noise generator hands out kc samples per trial, two per pair of draws.
      noise = gaussian_source (random_streams (s, 2 * ceil (trials * kc / 2)));
      channel = noise(1:trials * kc, :);
    else
      channel = random_streams (s, trials * kc);
    end
    channel = reshape (channel, kc, []).';
    failures(idx) = sum (reshape (trial.run (source, channel), trials, numel (s)), 1);
  end
end

function bits = random_bits (source)
  % One bit per source draw, 1 when the draw falls below 0.5. Python draws k
  % of them per call, and here each row of source holds one trial's k draws.
  bits = double (source < 0.5);
end

function trial = repetition_bsc_trial (n, p)
  trial = struct ('source_draws', 1, 'channel_draws', n, ...
                  'run', @(source, channel) repetition_bsc_run (source, channel, n, p));
end

function failed = repetition_bsc_run (source, channel, n, p)
  b = random_bits (source);
  failed = majority (bsc (repetition_encode (b, n), p, channel)) ~= b;
end

function trial = repetition_bec_trial (n, eps)
  trial = struct ('source_draws', 1, 'channel_draws', n, ...
                  'run', @(source, channel) repetition_bec_run (source, channel, n, eps));
end

function failed = repetition_bec_run (source, channel, n, eps)
  b = random_bits (source);
  % NaN, the all-erased outcome, differs from every bit, as None does in Python.
  failed = repetition_erasures (bec (repetition_encode (b, n), eps, channel)) ~= b;
end

function trial = parity_undetected_trial (n, p)
  % True when the received word passes the parity check but differs from the word sent.
  trial = struct ('source_draws', n - 1, 'channel_draws', n, ...
                  'run', @(source, channel) parity_undetected_run (source, channel, p));
end

function failed = parity_undetected_run (source, channel, p)
  word = parity_encode (random_bits (source));
  received = bsc (word, p, channel);
  failed = parity_ok (received) & any (received ~= word, 2);
end

function trial = parity_detected_trial (n, p)
  trial = struct ('source_draws', n - 1, 'channel_draws', n, ...
                  'run', @(source, channel) parity_detected_run (source, channel, p));
end

function failed = parity_detected_run (source, channel, p)
  word = parity_encode (random_bits (source));
  failed = ~parity_ok (bsc (word, p, channel));
end

function trial = parity_bec_trial (n, eps)
  trial = struct ('source_draws', n - 1, 'channel_draws', n, ...
                  'run', @(source, channel) parity_bec_run (source, channel, eps));
end

function failed = parity_bec_run (source, channel, eps)
  word = parity_encode (random_bits (source));
  [filled, ok] = parity_fill (bec (word, eps, channel));
  if (any (ok & any (filled ~= word, 2)))
    error ('parity filled an erasure wrongly');
  end
  failed = ~ok;
end

function trial = uncoded_bpsk_trial (ebn0_db)
  sigma = awgn_sigma (ebn0_db, 1.0);
  trial = struct ('source_draws', 1, 'channel_draws', 1, ...
                  'run', @(source, noise) uncoded_bpsk_run (source, noise, sigma));
end

function failed = uncoded_bpsk_run (source, noise, sigma)
  b = random_bits (source);
  failed = double (bpsk_awgn (b, sigma, noise) < 0) ~= b;
end

function trial = repetition_awgn_trial (n, ebn0_db, decoder)
  sigma = awgn_sigma (ebn0_db, 1.0 / n);
  trial = struct ('source_draws', 1, 'channel_draws', n, ...
                  'run', @(source, noise) repetition_awgn_run (source, noise, n, sigma, decoder));
end

function failed = repetition_awgn_run (source, noise, n, sigma, decoder)
  b = random_bits (source);
  failed = decoder (bpsk_awgn (repetition_encode (b, n), sigma, noise)) ~= b;
end


% Estimation ---------------------------------------------------------------------
% Each of these works entry by entry, so k may be a vector of failure counts.

function [lo, hi] = wald (k, n, z)
  if (nargin < 3)
    z = Z95;
  end
  ph = k ./ n;
  half = z .* sqrt (ph .* (1 - ph) ./ n);
  lo = ph - half;
  hi = ph + half;
end

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

function total = exact_coverage (p, n, interval)
  % Probability that the interval after n trials contains the true rate p.
  % nchoosek is exact only below 2^53, which every term this file adds stays
  % under, since the intervals that cover these rates have small K.
  if (nargin < 3)
    interval = @wilson;
  end
  total = 0.0;
  for k = 0:n
    [lo, hi] = interval (k, n);
    if (lo <= p && p <= hi)
      total = total + nchoosek (n, k) * p^k * (1 - p)^(n - k);
    end
  end
end

function b = zero_failure_bound (n, confidence)
  % Largest rate that still shows no failure in n trials with probability 1 - confidence.
  if (nargin < 2)
    confidence = 0.95;
  end
  b = 1 - (1 - confidence)^(1 / n);
end

function t = trials_needed (p, r, z)
  % Trials for the interval's half-width to be at most r times p.
  if (nargin < 3)
    z = Z95;
  end
  t = ceil ((z / r)^2 * (1 - p) / p);
end


% Closed forms -------------------------------------------------------------------

function v = q (x)
  v = 0.5 * erfc (x / sqrt (2));
end

function v = repetition_bsc_exact (n, p)
  t = floor ((n - 1) / 2);
  terms = zeros (1, n - t);
  for i = t + 1:n
    terms(i - t) = nchoosek (n, i) * p^i * (1 - p)^(n - i);
  end
  v = py_sum (terms);
end

function v = repetition_bec_exact (n, eps)
  v = eps^n;
end

function v = uncoded_bpsk_exact (ebn0_db)
  v = q (sqrt (2 * 10^(ebn0_db / 10)));
end

function v = hard_repetition_awgn_exact (n, ebn0_db)
  v = repetition_bsc_exact (n, q (sqrt (2 * 10^(ebn0_db / 10) / n)));
end

function v = soft_repetition_awgn_exact (n, ebn0_db)
  v = uncoded_bpsk_exact (ebn0_db);
end

function v = parity_undetected_exact (n, p)
  v = (1 + (1 - 2 * p)^n) / 2 - (1 - p)^n;
end

function v = parity_detected_exact (n, p)
  v = (1 - (1 - 2 * p)^n) / 2;
end

function v = parity_bec_exact (n, eps)
  v = 1 - (1 - eps)^n - n * eps * (1 - eps)^(n - 1);
end


% The experiments ----------------------------------------------------------------
% Each row: name, seed, trials, the trial, Gaussian or not, exact value.

function e = experiments ()
  e = {
    'R3 on BSC(0.1)', 2026, 10000, repetition_bsc_trial(3, 0.1), false, repetition_bsc_exact(3, 0.1);
    'R3 on BEC(0.2)', 12, 100000, repetition_bec_trial(3, 0.2), false, repetition_bec_exact(3, 0.2);
    'parity (8,7) on BSC(0.01), undetected', 13, 100000, parity_undetected_trial(8, 0.01), false, ...
      parity_undetected_exact(8, 0.01);
    'parity (8,7) on BEC(0.1), unfilled', 14, 20000, parity_bec_trial(8, 0.1), false, ...
      parity_bec_exact(8, 0.1);
    'uncoded BPSK at 4 dB', 15, 100000, uncoded_bpsk_trial(4.0), true, uncoded_bpsk_exact(4.0);
    'hard R3 at 4 dB', 16, 100000, repetition_awgn_trial(3, 4.0, @hard_repetition), true, ...
      hard_repetition_awgn_exact(3, 4.0);
    'soft R3 at 4 dB', 16, 100000, repetition_awgn_trial(3, 4.0, @soft_repetition), true, ...
      soft_repetition_awgn_exact(3, 4.0);
  };
end

function out = run_experiments ()
  % One row per experiment: name, seed, trials, failures, Wilson low and high, exact.
  e = experiments ();
  out = cell (size (e, 1), 7);
  for i = 1:size (e, 1)
    [name, seed, trials, trial, gaussian, exact] = e{i, :};
    k = count_failures (trials, seed, trial, gaussian);
    [lo, hi] = wilson (k, trials);
    out(i, :) = {name, seed, trials, k, lo, hi, exact};
  end
end

function n = COVERAGE_RUNS ()
  n = 1000;
end

function n = COVERAGE_LONG ()
  n = 1000;
end

function n = COVERAGE_SHORT ()
  n = 100;
end

function p = COVERAGE_EXACT ()
  p = repetition_bsc_exact (3, 0.1);
end

function counts = coverage_counts (trials, runs)
  % Failures of R3 on BSC(0.1) in each of `runs` runs, seeds 1 to runs, all in one call.
  if (nargin < 2)
    runs = COVERAGE_RUNS;
  end
  counts = count_failures (trials, 1:runs, repetition_bsc_trial (3, 0.1));
end

function zs = box_muller_samples (count, seed)
  rng = PyRandom (seed);
  zs = gaussian_source (rng.random (2 * ceil (count / 2)));
  zs = zs(1:count).';
end

function s = BOX_MULLER_SEED ()
  s = 17;
end

function n = BOX_MULLER_COUNT ()
  n = 100000;
end

function d = SWEEP_DB ()
  d = 0:7;
end

function n = SWEEP_TRIALS ()
  n = 100000;
end

function s = sweep_seed (db)
  s = 100 + db;
end

function out = awgn_sweep ()
  % Uncoded, hard R3 and soft R3 at 0 to 7 dB. Each point uses seed 100 + dB and
  % 100,000 trials. One row per point: dB, then the three failure counts.
  out = zeros (0, 4);
  for db = SWEEP_DB
    seed = sweep_seed (db);
    ku = count_failures (SWEEP_TRIALS, seed, uncoded_bpsk_trial (db), true);
    kh = count_failures (SWEEP_TRIALS, seed, repetition_awgn_trial (3, db, @hard_repetition), true);
    ks = count_failures (SWEEP_TRIALS, seed, repetition_awgn_trial (3, db, @soft_repetition), true);
    out(end + 1, :) = [db, ku, kh, ks];
  end
end

function s = ZERO_SEED ()
  s = 18;
end

function n = ZERO_TRIALS ()
  n = 20000;
end

function [both, hard_only, soft_only] = paired_hard_soft (trials, seed, ebn0_db, n)
  % Decode the same received words both ways and count the words both get wrong, only hard, only soft.
  if (nargin < 4)
    n = 3;
  end
  source = PyRandom (seed + SOURCE_OFFSET);
  channel = PyRandom (seed);
  noise = gaussian_source (channel.random (2 * ceil (trials * n / 2)));
  noise = reshape (noise(1:trials * n), n, trials).';
  sigma = awgn_sigma (ebn0_db, 1.0 / n);
  b = random_bits (source.random (trials));
  received = bpsk_awgn (repetition_encode (b, n), sigma, noise);
  hard_wrong = hard_repetition (received) ~= b;
  soft_wrong = soft_repetition (received) ~= b;
  both = sum (hard_wrong & soft_wrong);
  hard_only = sum (hard_wrong & ~soft_wrong);
  soft_only = sum (soft_wrong & ~hard_wrong);
end

function out = simulate_repetition (n, p, trials, rng)
  % Copied from python/ch09_error_control_strategies_and_performance.py, which
  % the Python imports. One entry per trial, 1 when majority voting over n
  % copies on BSC(p) fails. Each column of the draws is one trial's n copies.
  t = floor ((n - 1) / 2);
  flips = sum (reshape (rng.random (n * trials), n, trials) < p, 1);
  out = double (flips > t);
end


% Checks -------------------------------------------------------------------------

function check_channels_and_codes ()
  disp ('Channels and decoders');
  rng = PyRandom (1);
  same ('BSC with p = 0 changes nothing', int64 (bsc ([0, 1, 1, 0], 0.0, rng)), int64 ([0, 1, 1, 0]));
  same ('BSC with p = 1 flips every bit', int64 (bsc ([0, 1, 1, 0], 1.0, rng)), int64 ([1, 0, 0, 1]));
  same ('BEC with eps = 1 erases every bit', word_text (bec ([0, 1], 1.0, rng)), '[None, None]');
  same ('parity of 1011001', int64 (parity_encode ([1, 0, 1, 1, 0, 0, 1])), int64 ([1, 0, 1, 1, 0, 0, 1, 0]));
  same ('parity check passes 10110010', parity_ok ([1, 0, 1, 1, 0, 0, 1, 0]), true);
  same ('parity check fails 10100010', parity_ok ([1, 0, 1, 0, 0, 0, 1, 0]), false);
  same ('parity check passes 10100011, two errors', parity_ok ([1, 0, 1, 0, 0, 0, 1, 1]), true);
  [word, ok] = parity_fill ([1, NaN, 1, 1, 0, 0, 1, 0]);
  same ('parity fills 1 ? 1 1 0 0 1 0', fill_text (word, ok), '[1, 0, 1, 1, 0, 0, 1, 0]');
  [word, ok] = parity_fill ([1, NaN, 1, NaN, 0, 0, 1, 0]);
  same ('parity cannot fill two erasures', fill_text (word, ok), 'None');
  same ('soft R3 decision on 0.9, -1.4, 0.2', int64 (soft_repetition ([0.9, -1.4, 0.2])), int64 (1));
  same ('hard R3 decision on 0.9, -1.4, 0.2', int64 (hard_repetition ([0.9, -1.4, 0.2])), int64 (0));
  check ('sum of 0.9, -1.4, 0.2', 0.9 - 1.4 + 0.2, '-0.3');
  check ('Eb/N0 = 4 dB as a ratio', 10^0.4, '2.5119');
  check ('noise sigma for uncoded BPSK at 4 dB', awgn_sigma (4.0, 1.0), '0.4462');
  check ('noise sigma for R3 pulses at 4 dB', awgn_sigma (4.0, 1 / 3), '0.7728');
end

function check_box_muller ()
  disp ('Box-Muller');
  u1 = 0.25;
  u2 = 0.1;
  r = sqrt (-2 * log (u1));
  check ('worked pair: radius for u1 = 0.25', r, '1.6651');
  check ('worked pair: angle for u2 = 0.1, radians', 2 * pi * u2, '0.6283');
  check ('worked pair: z0', r * cos (2 * pi * u2), '1.3471');
  check ('worked pair: z1', r * sin (2 * pi * u2), '0.9787');
  check ('worked pair: cos(2 pi 0.1)', cos (2 * pi * u2), '0.8090');
  check ('worked pair: sin(2 pi 0.1)', sin (2 * pi * u2), '0.5878');
  check ('P(R^2 > 2) = e^(-1)', exp (-1), '0.3679');
  rng = PyRandom (BOX_MULLER_SEED);
  a = rng.random ();
  b = rng.random ();
  same ('first uniform from seed 17, as Python prints it', py_str (a), '0.5219839097124932');
  same ('second uniform from seed 17, as Python prints it', py_str (b), '0.8066907771186791');
  check ('first uniform from seed 17', a, '0.52198');
  check ('second uniform from seed 17', b, '0.80669');
  check ('u1 = 1 - first uniform', 1 - a, '0.47802');
  check ('R^2 from seed 17', -2 * log (1 - a), '1.4762');
  r = sqrt (-2 * log (1 - a));
  check ('radius from seed 17', r, '1.2150');
  check ('angle from seed 17, radians', 2 * pi * b, '5.0686');
  check ('cos of the seed-17 angle', cos (2 * pi * b), '0.3487');
  check ('sin of the seed-17 angle', sin (2 * pi * b), '-0.9372');
  zs = box_muller_samples (BOX_MULLER_COUNT, BOX_MULLER_SEED);
  check ('first sample from seed 17', zs(1), '0.4237');
  check ('second sample from seed 17', zs(2), '-1.1387');
  mean = py_sum (zs) / numel (zs);
  var = py_sum (py_pow (zs - mean, 2)) / numel (zs);
  tail = sum (abs (zs) > Z95);
  check ('mean of 100,000 samples, seed 17', mean, '0.0044');
  check ('that mean in standard errors', mean * sqrt (numel (zs)), '1.4');
  check ('variance of 100,000 samples, seed 17', var, '0.9981');
  same ('samples beyond 1.96 in size, of 100,000', int64 (tail), int64 (5014));
  check ('2Q(1.96)', 2 * q (Z95), '0.0500');
  ph = tail / numel (zs);
  check ('tail fraction', ph, '0.05014');
  check ('one minus the tail fraction', 1 - ph, '0.94986');
  se = sqrt (ph * (1 - ph) / numel (zs));
  check ('estimated standard deviation of the tail fraction', se, '0.00069012');
  check ('Wald half-width for the tail fraction', Z95 * se, '0.0013526');
  [lo, hi] = wald (tail, numel (zs));
  check ('Wald low for the tail count', lo, '0.04879');
  check ('Wald high for the tail count', hi, '0.05149');
  [lo, hi] = wilson (tail, numel (zs));
  check ('Wilson low for the tail count', lo, '0.04880');
  check ('Wilson high for the tail count', hi, '0.05151');
  % Bin i, counted from 0 as in the Python, holds -4 + 0.25 i <= z < -4 + 0.25 (i + 1).
  % Dividing by 0.25, a power of two, is exact, so floor matches Python's //.
  shown = zs(-4 <= zs & zs < 4);
  bins = floor ((shown + 4) / 0.25);
  counts = accumarray (bins(:) + 1, 1, [32, 1]).';
  [~, top] = max (counts);   % the first tallest bar, as Python's max() picks
  peak = top - 1;
  same ('lower edge of the tallest histogram bar', -4 + 0.25 * peak, -0.25);
  same ('samples in the tallest histogram bar', int64 (counts(peak + 1)), int64 (9907));
  same ('samples outside -4 to 4', int64 (numel (zs) - sum (counts)), int64 (3));
  check ('peak of the standard normal density', 1 / sqrt (2 * pi), '0.3989');
end

function yes = inside (interval, k, n)
  % Whether each interval after k failures in n trials covers the true rate.
  % The Python defines this inside check_estimator.
  [a, b] = interval (k, n);
  yes = a <= COVERAGE_EXACT & COVERAGE_EXACT <= b;
end

function check_estimator ()
  disp ('Estimator and intervals');
  rows = run_experiments ();
  [name, seed, trials, k, lo, hi, exact] = rows{1, :};
  same ('R3 on BSC(0.1), seed 2026, 10,000 words: failures', int64 (k), int64 (266));
  same ('chapter 9''s simulator on the same seed', ...
        int64 (sum (simulate_repetition (3, 0.1, 10000, PyRandom (2026)))), int64 (266));
  check ('exact R3 word error at p = 0.1', exact, '0.028');
  check ('Wilson low', lo, '0.02362');
  check ('Wilson high', hi, '0.02994');
  [wl, wh] = wald (k, trials);
  check ('Wald low', wl, '0.02345');
  check ('Wald high', wh, '0.02975');

  disp ('Coverage of the intervals, 1,000 seeded runs');
  longs = coverage_counts (COVERAGE_LONG);
  shorts = coverage_counts (COVERAGE_SHORT);
  estimates = longs / COVERAGE_LONG;
  mean = py_sum (estimates) / numel (estimates);
  sd = sqrt (py_sum (py_pow (estimates - mean, 2)) / (numel (estimates) - 1));
  check ('mean of 1,000 estimates at N = 1,000', mean, '0.0280');
  check ('their standard deviation', sd, '0.00520');
  check ('predicted standard deviation sqrt(p(1-p)/N)', ...
         sqrt (COVERAGE_EXACT * (1 - COVERAGE_EXACT) / COVERAGE_LONG), '0.00522');
  same ('fewest failures in a run of 1,000 words', int64 (min (longs)), int64 (13));
  same ('most failures in a run of 1,000 words', int64 (max (longs)), int64 (47));

  same ('Wilson intervals covering 0.028, N = 1,000', int64 (sum (inside (@wilson, longs, COVERAGE_LONG))), int64 (957));
  same ('Wald intervals covering 0.028, N = 1,000', int64 (sum (inside (@wald, longs, COVERAGE_LONG))), int64 (947));
  same ('Wilson intervals covering 0.028, N = 100', int64 (sum (inside (@wilson, shorts, COVERAGE_SHORT))), int64 (981));
  same ('Wald intervals covering 0.028, N = 100', int64 (sum (inside (@wald, shorts, COVERAGE_SHORT))), int64 (950));
  same ('runs at N = 100 with no failure', int64 (sum (shorts == 0)), int64 (49));
  check ('probability of no failure in 100 words', (1 - COVERAGE_EXACT)^100, '0.0584');
  same ('Wald misses at N = 100 that come from K = 0', ...
        int64 (sum (shorts == 0 & ~inside (@wald, shorts, COVERAGE_SHORT))), int64 (49));
  same ('Wilson misses at N = 100 with K below 7', ...
        int64 (sum (shorts < 7 & ~inside (@wilson, shorts, COVERAGE_SHORT))), int64 (0));
  [lo, ~] = wald (1, 100);
  check ('Wald low at K = 1, N = 100', lo, '-0.0095');
  first50 = longs(1:50);
  same ('misses among the first 50 Wilson intervals at N = 1,000', ...
        int64 (sum (~inside (@wilson, first50, COVERAGE_LONG))), int64 (2));
  % The run with seed s is entry s of longs, so its 1-based position is its seed.
  missed = find (~inside (@wilson, first50, COVERAGE_LONG));
  pairs = arrayfun (@(s) tuple (int64 (s), int64 (first50(s))), missed, 'UniformOutput', false);
  same ('those two misses as (seed, failures)', pairs, ...
        {tuple(int64 (20), int64 (45)), tuple(int64 (50), int64 (16))});
  same ('failure counts of the Wilson misses at N = 100', ...
        int64 (unique (shorts(~inside (@wilson, shorts, COVERAGE_SHORT)))), int64 ([7, 8, 9]));
  same ('Wald misses at N = 100', int64 (sum (~inside (@wald, shorts, COVERAGE_SHORT))), int64 (50));
  check ('exact Wilson coverage at P = 0.028, N = 100', exact_coverage (COVERAGE_EXACT, COVERAGE_SHORT), '0.9774');
  [lo, ~] = wilson (1, 100);
  check ('Wilson low at K = 1, N = 100', lo, '0.0018');
  [lo0, hi0] = wilson (0, 100);
  [lo1, ~] = wilson (1, 100);
  holds ('at P = 0.0017, N = 100 only K = 0 covers', lo0 <= 0.0017 && 0.0017 <= hi0 && lo1 > 0.0017);
  check ('1 - 0.0017', 1 - 0.0017, '0.9983');
  check ('exact Wilson coverage at P = 0.0017, N = 100', exact_coverage (0.0017, COVERAGE_SHORT), '0.8435');

  disp ('Zero failures and the rule of three');
  k0 = count_failures (ZERO_TRIALS, ZERO_SEED, repetition_bsc_trial (5, 0.01));
  same ('R5 on BSC(0.01), seed 18, 20,000 words: failures', int64 (k0), int64 (0));
  p5 = repetition_bsc_exact (5, 0.01);
  check ('exact R5 word error at p = 0.01', p5, '9.8506e-6');
  check ('expected failures in 20,000 words', ZERO_TRIALS * p5, '0.197');
  check ('probability of no failure in 20,000 words', (1 - p5)^ZERO_TRIALS, '0.821');
  check ('-ln 0.05', -log (0.05), '2.9957');
  check ('exact zero-failure bound for N = 20,000', zero_failure_bound (ZERO_TRIALS), '1.4978e-4');
  check ('rule of three for N = 20,000', 3 / ZERO_TRIALS, '1.5e-4');
  check ('rule-of-three bound over the true R5 rate', 3 / ZERO_TRIALS / p5, '15');
  check ('e^-3', exp (-3), '0.0498');
  [~, hi] = wilson (0, ZERO_TRIALS);
  check ('Wilson high at K = 0, N = 20,000', hi, '1.9204e-4');
  check ('-ln 0.025', -log (0.025), '3.6889');
  [lo, hi] = wald (0, ZERO_TRIALS);
  same ('Wald interval at K = 0 has zero width', tuple (lo, hi), tuple (0.0, 0.0));

  disp ('Trials for a target precision');
  check ('(1.96/0.2)^2', (Z95 / 0.2)^2, '96.04');
  check ('(1.96/0.1)^2', (Z95 / 0.1)^2, '384.16');
  pu = parity_undetected_exact (8, 0.01);
  same ('trials for 10% precision on the parity escape rate', int64 (trials_needed (pu, 0.1)), int64 (145307));
  same ('trials for 20% precision on the parity escape rate', int64 (trials_needed (pu, 0.2)), int64 (36327));
  check ('trials for 20% at P = 1e-5, millions', trials_needed (1e-5, 0.2) / 1e6, '9.60');
end

function s = STATED ()
  % What the chapter's summary table states for each experiment: failures, Wilson
  % limits, exact value, and whether the exact value lies inside the interval.
  s = {
    'R3 on BSC(0.1)', 266, '0.02362', '0.02994', '0.028', true;
    'R3 on BEC(0.2)', 741, '0.006897', '0.007961', '0.008', false;
    'parity (8,7) on BSC(0.01), undetected', 265, '0.002350', '0.002988', '0.0026368', true;
    'parity (8,7) on BEC(0.1), unfilled', 3720, '0.18067', '0.19145', '0.18690', true;
    'uncoded BPSK at 4 dB', 1246, '0.011791', '0.013167', '0.012501', true;
    'hard R3 at 4 dB', 2691, '0.025925', '0.027931', '0.026835', true;
    'soft R3 at 4 dB', 1250, '0.011830', '0.013208', '0.012501', true;
  };
end

function check_closed_forms ()
  disp ('Closed-form checks');
  rows = run_experiments ();
  stated = STATED ();
  for i = 1:size (rows, 1)
    [name, seed, trials, k, lo, hi, exact] = rows{i, :};
    [sk, slo, shi, sexact, sinside] = stated{strcmp (stated(:, 1), name), 2:end};
    same (sprintf ('%s: failures in %s, seed %d', name, commas (trials), seed), int64 (k), int64 (sk));
    check ([name ': Wilson low'], lo, slo);
    check ([name ': Wilson high'], hi, shi);
    check ([name ': exact'], exact, sexact);
    same ([name ': exact value inside the interval'], lo <= exact && exact <= hi, sinside);
  end
  lows = [rows{:, 5}];
  highs = [rows{:, 6}];
  exacts = [rows{:, 7}];
  same ('intervals in the table that miss, of 7', int64 (sum (~(lows <= exacts & exacts <= highs))), int64 (1));

  disp ('The erasure miss');
  sd = sqrt (100000 * 0.008 * 0.992);
  check ('standard deviation of K for R3 on BEC(0.2), N = 100,000', sd, '28.17');
  check ('distance of 741 from 800 in standard deviations', (741 - 800) / sd, '-2.09');
  k = count_failures (1000000, 12, repetition_bec_trial (3, 0.2));
  same ('same run continued to 1,000,000 words: failures', int64 (k), int64 (7974));
  [lo, hi] = wilson (k, 1000000);
  check ('its Wilson low', lo, '0.007802');
  check ('its Wilson high', hi, '0.008150');
  % The Python loops over the seeds. count_failures takes all 100 in one call.
  trial = repetition_bec_trial (3, 0.2);
  [lo, hi] = wilson (count_failures (100000, 1:100, trial), 100000);
  covered = sum (lo <= 0.008 & 0.008 <= hi);
  same ('runs with seeds 1 to 100 whose interval covers 0.008', int64 (covered), int64 (98));

  disp ('Parity and Gaussian details');
  check ('0.98^8', 0.98^8, '0.850763');
  check ('0.99^8', 0.99^8, '0.922745');
  check ('0.9^8', 0.9^8, '0.430467');
  check ('8 x 0.1 x 0.9^7', 8 * 0.1 * 0.9^7, '0.382638');
  check ('parity detected rate, n = 8, p = 0.01', parity_detected_exact (8, 0.01), '0.07462');
  same ('parity words flagged, seed 13, 100,000 words', ...
        int64 (count_failures (100000, 13, parity_detected_trial (8, 0.01))), int64 (7565));
  [lo, hi] = wilson (265, 100000);
  check ('relative half-width of the parity escape interval', (hi - lo) / 2 / (265 / 100000), '0.120');
  check ('1.96 / sqrt(265)', Z95 / sqrt (265), '0.120');
  check ('argument of Q for uncoded BPSK at 4 dB', sqrt (2 * 10^0.4), '2.2414');
  check ('argument of Q for each R3 pulse at 4 dB', sqrt (2 * 10^0.4 / 3), '1.2941');
  check ('symbol error of each R3 pulse at 4 dB', q (sqrt (2 * 10^0.4 / 3)), '0.09782');
  [both, hard_only, soft_only] = paired_hard_soft (100000, 16, 4.0);
  same ('words both decoders get wrong, seed 16', int64 (both), int64 (894));
  same ('words only the hard decoder gets wrong', int64 (hard_only), int64 (1797));
  same ('words only the soft decoder gets wrong', int64 (soft_only), int64 (356));
  same ('extra failures of hard over soft', int64 (hard_only - soft_only), int64 (1441));

  disp ('Gaussian sweep, seeds 100 + dB, 100,000 trials per point');
  misses = {};
  sweep = awgn_sweep ();
  for i = 1:size (sweep, 1)
    db = sweep(i, 1);
    points = {'uncoded', sweep(i, 2), uncoded_bpsk_exact(db);
              'hard R3', sweep(i, 3), hard_repetition_awgn_exact(3, db);
              'soft R3', sweep(i, 4), soft_repetition_awgn_exact(3, db)};
    for j = 1:3
      [label, k, exact] = points{j, :};
      [lo, hi] = wilson (k, SWEEP_TRIALS);
      if (~(lo <= exact && exact <= hi))
        misses{end + 1} = sprintf ('%s at %d dB', label, db);
      end
    end
  end
  same ('sweep points whose interval misses the closed form, of 24', int64 (numel (misses)), int64 (1));
  same ('the one miss', misses, {'hard R3 at 1 dB'});
  % Python's sweep[1][2] is row 2, column 3 here.
  same ('hard R3 failures at 1 dB, seed 101', int64 (sweep(2, 3)), int64 (8727));
  check ('expected hard R3 failures at 1 dB', SWEEP_TRIALS * hard_repetition_awgn_exact (3, 1), '8,536.0');
  check ('misses expected among 24 intervals', 24 * 0.05, '1.2');
  check ('chance that 24 independent 95% intervals all cover', 0.95^24, '0.2920');
  same ('uncoded, soft and hard failures at 7 dB, seed 107', ...
        tuple (int64 (sweep(8, 2)), int64 (sweep(8, 4)), int64 (sweep(8, 3))), ...
        tuple (int64 (92), int64 (63), int64 (308)));
end

function check_worked_numbers ()
  disp ('Other numbers the chapter works by hand');
  check ('-2 ln 0.25', -2 * log (0.25), '2.7726');
  check ('standard error of a mean of 100,000 standard normal samples', 1 / sqrt (100000), '0.0032');
  noise = box_muller_samples (2, BOX_MULLER_SEED);
  sigma = awgn_sigma (4.0, 1.0);
  check ('received value for +1 plus the first seed-17 sample at 4 dB', 1 + sigma * noise(1), '1.1890');
  check ('received value for +1 plus the second seed-17 sample at 4 dB', 1 + sigma * noise(2), '0.4919');
  check ('(1 + 0.98^8) / 2', (1 + 0.98^8) / 2, '0.92538151');
  check ('0.99^8 to eight places', 0.99^8, '0.92274469');
  check ('their difference', (1 + 0.98^8) / 2 - 0.99^8, '0.00263682');
  check ('0.95^7', 0.95^7, '0.6983');
  check ('chance that at least one of 7 intervals misses', 1 - 0.95^7, '0.3017');
  check ('misses expected among 7 intervals', 7 * 0.05, '0.35');
  check ('1.96^2', Z95^2, '3.8416');
  check ('-ln(0.025) / 20,000', -log (0.025) / ZERO_TRIALS, '1.8444e-4');
  check ('trials for 100 failures at the R5 rate, millions', 100 / repetition_bsc_exact (5, 0.01) / 1e6, '10.2');
  % 2^19937 - 1 overflows a double, so the digits are counted by logarithm.
  % As whole numbers, 2^19937 - 1 >= 10^6000 exactly when 2^19937 > 10^6000,
  % that is when 19937 log10(2) > 6000. The left side is 6001.6, a margin
  % no rounding in the logarithm can close.
  holds ('2^19937 - 1 has more than 6,000 digits', 19937 * log10 (2) > 6000);
  % Python's math.degrees(x) multiplies x by 180 / pi.
  check ('Box-Muller angle for u2 = 0.1, degrees', (180 / pi) * (2 * pi * 0.1), '36');
  check ('misses expected among 50 intervals', 50 * 0.05, '2.5');
  [lo, hi] = wilson (45, 1000);
  check ('Wilson low for run 20, 45 failures in 1,000', lo, '0.0338');
  check ('Wilson high for run 20', hi, '0.0597');
  [lo, hi] = wilson (16, 1000);
  check ('Wilson low for run 50, 16 failures in 1,000', lo, '0.0099');
  check ('Wilson high for run 50', hi, '0.0258');
  check ('failures expected in 100 words', 100 * repetition_bsc_exact (3, 0.1), '2.8');
  check ('runs of 1,000 expected to see no failure in 100 words', 1000 * (1 - repetition_bsc_exact (3, 0.1))^100, '58');
  [~, hi] = wilson (0, 100);
  check ('Wilson high at K = 0, N = 100', hi, '0.0370');
  [lo, ~] = wilson (7, 100);
  check ('Wilson low at K = 7, N = 100', lo, '0.0343');
  [lo, ~] = wilson (6, 100);
  holds ('Wilson at K = 6, N = 100 still reaches below 0.028', lo < 0.028);
  [lo8, ~] = wald (8, 100);
  [lo9, ~] = wald (9, 100);
  holds ('Wald at N = 100 misses 0.028 from above exactly when K >= 9', lo8 < 0.028 && 0.028 < lo9);
  p = 0.01;
  check ('R5 term 10 p^3 (1-p)^2 at p = 0.01', 10 * p^3 * (1 - p)^2, '9.801e-6');
  check ('R5 term 5 p^4 (1-p) at p = 0.01', 5 * p^4 * (1 - p), '4.95e-8');
  check ('expected parity flags in 100,000 words', 100000 * parity_detected_exact (8, 0.01), '7,461.8');
  pu = parity_undetected_exact (8, 0.01);
  check ('(1 - P)/P for the parity escape rate', (1 - pu) / pu, '378.245');
  check ('expected uncoded failures in 100,000 words at 4 dB', 100000 * uncoded_bpsk_exact (4), '1,250.1');
  check ('expected hard R3 failures in 100,000 words at 4 dB', 100000 * hard_repetition_awgn_exact (3, 4), '2,683.5');
  check ('expected soft R3 failures in 100,000 words at 4 dB', 100000 * soft_repetition_awgn_exact (3, 4), '1,250.1');
  check ('hard over soft at 4 dB', hard_repetition_awgn_exact (3, 4) / soft_repetition_awgn_exact (3, 4), '2.147');
  check ('N + z^2 at N = 20,000', ZERO_TRIALS + Z95^2, '20,003.8416');
  disp ('Exercise answers');
  check ('-2 ln 0.5', -2 * log (0.5), '1.3863');
  check ('radius for u1 = 0.5', sqrt (-2 * log (0.5)), '1.1774');
  check ('(1 - P)/P at P = 1e-4', (1 - 1e-4) / 1e-4, '9,999');
  check ('384.16 x 9,999', (Z95 / 0.1)^2 * (1 - 1e-4) / 1e-4, '3,841,215.8');
  same ('trials for 10% precision at P = 1e-4', int64 (trials_needed (1e-4, 0.1)), int64 (3841216));
  check ('zero-failure bound for N = 1,000,000', zero_failure_bound (1000000), '2.9957e-6');
  check ('(1 + 0.8^4) / 2', (1 + 0.8^4) / 2, '0.7048');
  check ('0.9^4', 0.9^4, '0.6561');
  check ('two flips in four at p = 0.1', nchoosek (4, 2) * 0.1^2 * 0.9^2, '0.0486');
  check ('four flips in four at p = 0.1', 0.1^4, '0.0001');
  check ('parity undetected rate, n = 4, p = 0.1', parity_undetected_exact (4, 0.1), '0.0487');
  same ('hard R3 decision on 0.3, 0.4, -1.2', int64 (hard_repetition ([0.3, 0.4, -1.2])), int64 (0));
  same ('soft R3 decision on 0.3, 0.4, -1.2', int64 (soft_repetition ([0.3, 0.4, -1.2])), int64 (1));
  same ('BSC(0.1) on 1011 with draws 0.31, 0.04, 0.77, 0.09', ...
        int64 (xor ([1, 0, 1, 1], [0.31, 0.04, 0.77, 0.09] < 0.1)), int64 ([1, 1, 1, 0]));
end

function main ()
  check_channels_and_codes ();
  check_box_muller ();
  check_estimator ();
  check_closed_forms ();
  check_worked_numbers ();
  printf ('\nAll chapter 11 numbers reproduced.\n');
end

main ();
