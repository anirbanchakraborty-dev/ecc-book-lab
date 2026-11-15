1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

% CH09_ERROR_CONTROL_STRATEGIES_AND_PERFORMANCE  Reproduce every number in chapter 9 of
% Error Correction Coding, "Error-Control Strategies and Performance".
%
% This is the Octave twin of python/ch09_error_control_strategies_and_performance.py.
% It needs only core GNU Octave, with no pkg load and no Forge packages. From the
% repository root:
%
%     octave --no-window-system --norc --quiet octave/ch09_error_control_strategies_and_performance.m
%
% Each check prints the value computed here beside the value the chapter states, and
% the program stops with exit status 1 at the first pair that disagrees by more than
% half a unit in the last digit the chapter prints. The ARQ throughput formulas are
% checked against a frame-by-frame simulation of each protocol, the binomial word-error
% formula against exhaustive enumeration of every error pattern, and every Monte Carlo
% estimate against its exact formula. The simulations draw from PyRandom, a bit-for-bit
% copy of Python's random.Random, in the same order as the Python script, so every
% seeded count comes out the same. tools/compare-routes.sh runs both scripts and
% requires identical output. A run that reaches the final line has reproduced every
% number in the chapter.
%
% The charts stay in Python (python/ch09_charts.py), so this file holds only the checks.
%
% Octave loops are slow, so the long simulations draw their random numbers in batches
% and count with whole-array operations. Each batch is sized so that it consumes exactly
% the draws the Python loop consumes, in the same order, and the comments at each one
% say why.


% Helpers ------------------------------------------------------------------------
%
% half_unit comes from octave/common and does what the Python's own half_unit does.

function check (label, computed, stated)
  text = strrep (stated, ',', '');
  while (~isempty (text) && text(end) == '%')
    text(end) = [];
  end
  value = str2double (text);
  computed = double (computed);
  if (~isempty (stated) && stated(end) == '%')
    computed = 100 * computed;
  end
  printf ('  %-62s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ("MISMATCH in '%s': computed %s, the chapter states %s\n", label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % A Python int is passed as int64 and a Python float as a double, so that py_str
  % prints each the way Python's str() does.
  printf ('  %-62s %-16s %s\n', label, py_str (computed), py_str (stated));
  if (computed ~= stated)
    printf ("MISMATCH in '%s': computed %s, the chapter states %s\n", label, py_str (computed), py_str (stated));
    exit (1);
  end
end

function near (label, a, b, rel)
  % Two computations of one quantity that must agree to a relative tolerance.
  printf ('  %-62s %-16s %s\n', label, py_g (a), py_g (b));
  if (abs (a - b) > rel * abs (b))
    printf ("MISMATCH in '%s': %s against %s\n", label, py_str (a), py_str (b));
    exit (1);
  end
end

function section (title)
  printf ('\n%s\n', title);
end


% Python's arithmetic, rebuilt ------------------------------------------------------
%
% Four things the Python takes from its language and library come out differently in
% Octave's own functions: the compensated sum() of Python 3.12 and later, CPython's own
% math.lgamma, the exact integer math.comb, and round() to a number of decimal places.
% Each is rebuilt here so that every printed digit matches.

function total = py_sum (x)
  % Python's built-in sum() of floats. Since Python 3.12 it adds with Neumaier's
  % compensated summation, carrying the rounding error of each addition in c and adding
  % c once at the end, while Octave's sum() adds plainly. The two can differ in the last
  % digit, so every Python sum() of floats goes through this function, in the Python's
  % order. A matrix is summed down each column, so many sums run side by side.
  if (isvector (x))
    x = x(:);
  end
  if (isempty (x))
    total = 0;
    return;
  end
  total = x(1, :);
  c = zeros (size (total));
  for i = 2:rows (x)
    xi = x(i, :);
    t = total + xi;
    big = abs (total) >= abs (xi);
    d = (xi - t) + total;
    d(big) = (total(big) - t(big)) + xi(big);
    c = c + d;
    total = t;
  end
  % CPython adds the compensation only when it is nonzero and finite.
  use = c ~= 0 & isfinite (c);
  total(use) = total(use) + c(use);
end

function [s, e] = two_sum (a, b)
  % s = a + b rounded, and e the exact error, so that a + b = s + e (Knuth).
  s = a + b;
  bb = s - a;
  e = (a - (s - bb)) + (b - bb);
end

function [hi, lo] = split_double (a)
  % Veltkamp's split of a into two halves of 26 bits each, a = hi + lo exactly.
  c = 134217729 * a;   % 2^27 + 1
  hi = c - (c - a);
  lo = a - hi;
end

function [p, e] = two_prod (a, b)
  % p = a * b rounded, and e the exact error, so that a * b = p + e (Dekker).
  p = a .* b;
  [ah, al] = split_double (a);
  [bh, bl] = split_double (b);
  e = ((ah .* bh - p) + ah .* bl + al .* bh) + al .* bl;
end

function r = add_round_odd (x, y)
  % x + y rounded to odd: when the sum is inexact, the one of its two neighboring
  % doubles whose last mantissa bit is 1. Ordinary rounding gives one neighbor, s, and
  % the exact error e says on which side the other one lies.
  [s, e] = two_sum (x, y);
  r = s;
  [f, ex] = log2 (s);                 % s = f * 2^ex with 0.5 <= |f| < 1
  odd = mod (abs (f) * 2^53, 2) == 1; % the last of the 53 mantissa bits
  fix = (e ~= 0) & ~odd;
  if (any (fix(:)))
    sf = s(fix);
    ef = e(fix);
    step = 2 .^ (ex(fix) - 53);       % the spacing of doubles at s
    % Below a power of two the spacing halves.
    down = sign (ef) ~= sign (sf) & abs (f(fix)) == 0.5;
    step(down) = step(down) / 2;
    r(fix) = sf + sign (ef) .* step;
  end
end

function z = fused_multiply_add (a, b, c)
  % a * b + c with a single rounding, as a fused multiply-add instruction gives it.
  % Octave has no fma(), so this follows Boldo and Melquiond, "Emulation of FMA and
  % correctly rounded sums: proved algorithms using rounding to odd", IEEE Transactions
  % on Computers 57(4), 2008: the exact product, an exact sum with c, then the two
  % small parts added with rounding to odd before the one final rounding.
  [ph, pl] = two_prod (a, b);
  [th, tl] = two_sum (c, ph);
  z = th + add_round_odd (tl, pl);
end

function s = lanczos_sum (x)
  % CPython's lanczos_sum from Modules/mathmodule.c, for positive x.
  num_c = [23531376880.410759688572007674451636754734846804940, ...
           42919803642.649098768957899047001988850926355848959, ...
           35711959237.355668049440185451547166705960488635843, ...
           17921034426.037209699919755754458931112671403265390, ...
           6039542586.3520280050642916443072979210699388420708, ...
           1439720407.3117216736632230727949123939715485786772, ...
           248874557.86205415651146038641322942321632125127801, ...
           31426415.585400194380614231628318205362874684987640, ...
           2876370.6289353724412254090516208496135991145378768, ...
           186056.26539522349504029498971604569928220784236328, ...
           8071.6720023658162106380029022722506138218516325024, ...
           210.82427775157934587250973392071336271166969580291, ...
           2.5066282746310002701649081771338373386264310793408];
  den_c = [0.0, 39916800.0, 120543840.0, 150917976.0, 105258076.0, 45995730.0, ...
           13339535.0, 2637558.0, 357423.0, 32670.0, 1925.0, 66.0, 1.0];
  num = zeros (size (x));
  den = zeros (size (x));
  small = x < 5.0;
  % Below 5, Horner's rule in x, where each step is one fused multiply-add.
  xs = x(small);
  ns = num(small);
  ds = den(small);
  for i = 13:-1:1
    ns = fused_multiply_add (ns, xs, num_c(i));
    ds = fused_multiply_add (ds, xs, den_c(i));
  end
  % From 5 up, the same rational function in 1/x, with no multiply to fuse.
  xl = x(~small);
  nl = num(~small);
  dl = den(~small);
  for i = 1:13
    nl = nl ./ xl + num_c(i);
    dl = dl ./ xl + den_c(i);
  end
  num(small) = ns;
  den(small) = ds;
  num(~small) = nl;
  den(~small) = dl;
  s = num ./ den;
end

function r = py_lgamma (x)
  % Python's math.lgamma, for positive x. It is not the C library's lgamma, which
  % Octave's gammaln calls: CPython computes it with its own Lanczos formula, and the
  % two disagree in the last bit for about half of all arguments. binom_log_pmf
  % subtracts values near 1.7e9 to leave a few thousand, so a one-bit difference there
  % moves the seventh printed digit of window_prob. This follows the CPython source
  % line by line. The CPython build that tools/compare-routes.sh runs was compiled with
  % multiply-add contraction, so its line "r += a * b" rounds once, and
  % fused_multiply_add reproduces that.
  g = 6.024680040776729583740234375;
  r = log (lanczos_sum (x)) - g;
  r = fused_multiply_add (x - 0.5, log (x + g - 0.5) - 1, r);
  r(x == floor (x) & x <= 2) = 0;     % lgamma(1) = lgamma(2) = 0 exactly
end

function x = limbs_to_double (limbs)
  % The integer sum(limbs(j) * 2^(24 (j - 1))) rounded to the nearest double, ties to
  % even, as Python converts an int to a float.
  L = numel (limbs);
  bits = zeros (1, 24 * L);           % least significant bit first
  for j = 1:L
    bits(24 * (j - 1) + (1:24)) = bitget (limbs(j), 1:24);
  end
  B = find (bits, 1, 'last');         % the bit length
  if (B <= 53)
    x = sum (bits(1:B) .* 2 .^ (0:B - 1));
  else
    m = sum (bits(B - 52:B) .* 2 .^ (0:52));   % the top 53 bits
    guard = bits(B - 53);
    sticky = any (bits(1:B - 54));
    if (guard && (sticky || mod (m, 2) == 1))
      m = m + 1;
    end
    x = m * 2 ^ (B - 53);
  end
end

function row = comb_row (n)
  % float(comb(n, k)) for k = 0 to n, stored at row(k + 1). Python's comb is an exact
  % integer, and comb(n, i) * p ** i rounds it once to a double before multiplying.
  % Octave's nchoosek works in doubles and loses exactness beyond 2^53, which
  % comb(1000, k) passes from k = 7. So each value is built exactly in base-2^24
  % limbs, least significant first, by C(n, k) = C(n, k - 1) (n - k + 1) / k, and
  % rounded only at the end. A limb times n - k + 1 stays far below 2^53.
  base = 2 ^ 24;
  row = zeros (1, n + 1);
  row(1) = 1;
  limbs = 1;
  for k = 1:n
    limbs = limbs * (n - k + 1);
    carry = 0;
    for j = 1:numel (limbs)
      v = limbs(j) + carry;
      limbs(j) = mod (v, base);
      carry = floor (v / base);
    end
    while (carry > 0)
      limbs(end + 1) = mod (carry, base);
      carry = floor (carry / base);
    end
    % The division by k is exact, so long division from the top limb leaves no remainder.
    r = 0;
    for j = numel (limbs):-1:1
      v = r * base + limbs(j);
      limbs(j) = floor (v / k);
      r = v - limbs(j) * k;
    end
    while (numel (limbs) > 1 && limbs(end) == 0)
      limbs(end) = [];
    end
    row(k + 1) = limbs_to_double (limbs);
  end
end

function c = comb_float (n, k)
  % float(comb(n, k)) for a scalar n and a vector k. Each row of Pascal's triangle is
  % built once and kept, since word_error asks for the n = 1000 row hundreds of times.
  persistent rows_by_n;
  if (isempty (rows_by_n))
    rows_by_n = {};
  end
  if (numel (rows_by_n) < n + 1 || isempty (rows_by_n{n + 1}))
    rows_by_n{n + 1} = comb_row (n);
  end
  row = rows_by_n{n + 1};
  c = row(k + 1);                     % k is 0-based, the array 1-based
end

function y = py_round (x, ndigits)
  % Python's round(x, ndigits) for a float: x rounded in decimal to ndigits places,
  % read back as the nearest double. sprintf does that decimal rounding exactly.
  y = str2double (sprintf ('%.*f', ndigits, x));
end

function s = commas (n)
  % Python's format(n, ','): an integer with commas between groups of three digits.
  s = sprintf ('%d', n);
  k = numel (s) - 3;
  while (k > 0)
    s = [s(1:k) ',' s(k + 1:end)];
    k = k - 3;
  end
end


% The Gaussian tail ----------------------------------------------------------------

function y = q_func (x)
  % The Q-function, the probability that a standard Gaussian exceeds x.
  y = 0.5 * erfc (x / sqrt (2.0));
end

function x = q_inv (p)
  % The x with Q(x) = p, found by bisection.
  lo = -40.0;
  hi = 40.0;
  for i = 1:300
    mid = (lo + hi) / 2;
    if (q_func (mid) > p)
      lo = mid;
    else
      hi = mid;
    end
  end
  x = (lo + hi) / 2;
end

function y = db (x)
  y = 10 * log10 (x);
end


% Binomial counting -----------------------------------------------------------------

function y = binom_pmf (n, i, p)
  % i may be a vector. A scalar base raised to a vector of powers calls the C library's
  % pow() for each element, as Python's ** does. A vector raised to the power 2 or 3 is
  % multiplied out instead, which can differ in the last bit, so that form is avoided.
  y = comb_float (n, i) .* p .^ i .* (1 - p) .^ (n - i);
end

function y = word_error (n, t, p)
  % Probability that more than t of n bits flip on BSC(p), the binomial word-error formula.
  y = py_sum (binom_pmf (n, (t + 1):n, p));
end

function y = binom_log_pmf (n, k, p)
  y = py_lgamma (n + 1) - py_lgamma (k + 1) - py_lgamma (n - k + 1) + k * log (p) + (n - k) * log (1 - p);
end

function y = window_prob (n, p, rel)
  % Exact probability that K/n lands within rel*p of p when K is binomial(n, p).
  lo = ceil (n * p * (1 - rel) - 1e-9);
  hi = floor (n * p * (1 + rel) + 1e-9);
  k = max (lo, 0):min (hi, n);
  y = py_sum (exp (binom_log_pmf (n, k, p)));
end


% ARQ throughput ------------------------------------------------------------------

function y = frame_error (n, p)
  y = 1.0 - (1.0 - p) ^ n;
end

function y = eta_sw (P, a, rate)
  y = (1 - P) / (1 + a) * rate;
end

function y = eta_gbn (P, window, rate)
  y = (1 - P) / (1 + (window - 1) * P) * rate;
end

function y = eta_sr (P, rate)
  y = (1 - P) * rate;
end

function s = simulate_protocol (kind, P, a, frames, rng)
  % Slots, in frame times, that a protocol spends delivering a number of frames.
  %
  % Each transmission of a frame fails independently with probability P. A
  % stop-and-wait attempt holds the channel for 1 + a slots. Go-back-N with a window
  % of 1 + a slots loses the failed frame and the a frames sent behind it. Selective
  % repeat loses only the failed frame. Because each failure is charged the slots the
  % chapter's proof charges, this checks the averaging in that proof, not the slot
  % model of the protocols.
  %
  % The Python draws one random() per attempt until frames attempts have succeeded.
  % With "left" frames still to deliver, at least "left" more attempts are certain, so
  % a batch of exactly that many draws is always used up whole, in the Python's order.
  % Repeating until no frame is left consumes the same draws as the Python loop.
  slots = 0;
  left = frames;
  while (left > 0)
    failed = rng.random (left) < P;
    nfail = sum (failed);
    nok = left - nfail;
    switch (kind)
      case 'sw'
        slots = slots + left * (1 + a);
      case 'gbn'
        slots = slots + nfail * (1 + a) + nok;
      otherwise
        slots = slots + left;
    end
    left = left - nok;
  end
  s = slots / frames;
end

function delays = sr_release_delays (a, P, slots, seed)
  % Extra delay, in slots, of each frame under ideal selective repeat with in-order release.
  %
  % The sender fills every slot, resending a failed frame as soon as its NAK arrives a
  % slots after the frame ends and otherwise sending the next new frame. The receiver
  % buffers frames that arrive out of order and releases a frame only once every
  % earlier frame is in, so a frame's delay counts the wait behind earlier damaged
  % frames as well as its own retries. Delay runs from the frame's first transmission
  % to its release.
  %
  % Every slot draws exactly one random(), so all of them are drawn up front. The
  % Python keeps its resends in a heap keyed by (t + a, f). Those keys are pushed in
  % strictly increasing order of t, so the heap always pops the oldest entry, and a
  % first-in first-out queue (two columns, a head and a tail) does the same work.
  % Frame numbers stay 0-based as in the Python, and f + 1 indexes the arrays.
  rng = PyRandom (seed);
  r = rng.random (slots);
  due = zeros (slots, 2);
  head = 1;
  tail = 0;
  first = zeros (slots, 1);
  recv = -ones (slots, 1);             % -1 marks a frame not yet received
  nxt = 0;
  for t = 1:slots
    if (head <= tail && due(head, 1) <= t - 1)
      f = due(head, 2);
      head = head + 1;
    else
      f = nxt;
      nxt = nxt + 1;
      first(f + 1) = t;
    end
    if (r(t) < P)
      tail = tail + 1;
      due(tail, :) = [t + a, f];
    else
      recv(f + 1) = t;
    end
  end
  % Release stops at the first frame never received. The running maximum of the
  % receive times is the slot at which each frame is released.
  m = find (recv(1:nxt) < 0, 1) - 1;
  if (isempty (m))
    m = nxt;
  end
  released = cummax (recv(1:m));
  delays = released - first(1:m);
end

function delays = gbn_delivery_delays (a, P, slots, seed)
  % Extra delay, in slots, of each frame under go-back-N with a window of 1 + a.
  %
  % The receiver accepts only the next frame in order. When that frame fails, the
  % sender learns of it a slots after the frame ends and goes back to it, resending it
  % and every frame sent behind it. Delay runs from the frame's first transmission to
  % its acceptance.
  %
  % As in sr_release_delays, every slot draws one random(), drawn here up front. The
  % Python's goback tuple, or None, becomes goback_t and goback_f, with goback_t < 0
  % standing for None. A first time of 0 marks a frame not yet sent.
  rng = PyRandom (seed);
  r = rng.random (slots);
  first = zeros (slots, 1);
  recv = zeros (slots, 1);
  nxt = 0;
  expect = 0;
  goback_t = -1;
  goback_f = 0;
  for t = 1:slots
    if (goback_t >= 0 && goback_t <= t - 1)
      nxt = goback_f;
      goback_t = -1;
    end
    f = nxt;
    nxt = nxt + 1;
    if (first(f + 1) == 0)
      first(f + 1) = t;
    end
    ok = r(t) >= P;
    if (f == expect && ok)
      recv(f + 1) = t;
      expect = expect + 1;
    elseif (f == expect && goback_t < 0)
      goback_t = t + a;
      goback_f = f;
    end
  end
  % Frames are accepted in order, so the received frames are 0 to expect - 1.
  delays = recv(1:expect) - first(1:expect);
end

function y = hybrid_failure (n, t, p)
  y = word_error (n, t, p);
end

function x = hybrid_crossover (n, rate_crc, rate_hybrid, t)
  % The bit error probability where type-I hybrid ARQ and plain selective repeat tie.
  lo = 1e-7;
  hi = 1e-3;
  for i = 1:200
    mid = sqrt (lo * hi);
    if (eta_sr (frame_error (n, mid), rate_crc) > (1 - hybrid_failure (n, t, mid)) * rate_hybrid)
      lo = mid;
    else
      hi = mid;
    end
  end
  x = sqrt (lo * hi);
end


% Single parity under ARQ -----------------------------------------------------------

function out = parity_outcomes (n, p)
  % Probabilities that one parity-checked transmission is clean, flagged, or wrong but
  % accepted, returned as the row [clean, alarm, slip] in place of the Python's tuple.
  clean = (1 - p) ^ n;
  alarm = (1 - (1 - 2 * p) ^ n) / 2;
  slip = 1 - clean - alarm;
  out = [clean, alarm, slip];
end

function out = parity_outcomes_enumerated (n, p)
  % product((0, 1), repeat=n) runs through the patterns in the order of the binary
  % numbers 0 to 2^n - 1, and the weight of a pattern is the count of its 1 bits.
  % The Python adds with plain +=, so this loop does too, in the same order.
  clean = 0.0;
  alarm = 0.0;
  slip = 0.0;
  for m = 0:2 ^ n - 1
    w = sum (bitget (m, 1:n));
    pr = p ^ w * (1 - p) ^ (n - w);
    if (w == 0)
      clean = clean + pr;
    elseif (mod (w, 2))
      alarm = alarm + pr;
    else
      slip = slip + pr;
    end
  end
  out = [clean, alarm, slip];
end

function [wrong_rate, sends_rate] = simulate_arq_parity (n, p, frames, rng)
  % Deliver frames under ARQ with single parity, counting wrong deliveries and transmissions.
  %
  % Each transmission draws n random() values, and a frame ends at the first
  % transmission of even weight. With "left" frames still to deliver, at least "left"
  % more transmissions are certain, so a batch of that many is used up whole. Column j
  % of the reshaped batch holds the n draws of transmission j, in the Python's order.
  wrong = 0;
  sends = 0;
  left = frames;
  while (left > 0)
    w = sum (reshape (rng.random (n * left), n, left) < p, 1);
    even = mod (w, 2) == 0;
    sends = sends + left;
    wrong = wrong + sum (even & w > 0);
    left = left - sum (even);
  end
  wrong_rate = wrong / frames;
  sends_rate = sends / frames;
end


% Coding gain ---------------------------------------------------------------------

function [x, eb_db] = required_uncoded (target)
  % The Gaussian-channel argument x and the Eb/N0 in dB at which uncoded BPSK meets a target.
  x = q_inv (target);
  eb_db = db (x * x / 2);
end

function [p, x, es_db, eb_db] = required_repetition (n, target)
  % Channel flip probability, Q argument, Es/N0 and Eb/N0 in dB for hard-decision R_n to meet a target.
  t = floor ((n - 1) / 2);
  lo = 1e-15;
  hi = 0.5;
  for i = 1:300
    mid = (lo + hi) / 2;
    if (word_error (n, t, mid) > target)
      hi = mid;
    else
      lo = mid;
    end
  end
  p = (lo + hi) / 2;
  x = q_inv (p);
  es = x * x / 2;
  es_db = db (es);
  eb_db = db (es * n);
end

function y = hard_rep_ber (n, ebn0_db)
  % Decoded bit error rate of hard-decision R_n on the Gaussian channel.
  ebn0 = 10 ^ (ebn0_db / 10);
  p = q_func (sqrt (2 * ebn0 / n));
  y = word_error (n, floor ((n - 1) / 2), p);
end

function y = uncoded_ber (ebn0_db)
  y = q_func (sqrt (2 * 10 ^ (ebn0_db / 10)));
end

function y = asymptotic_hard_gain (n)
  t = floor ((n - 1) / 2);
  y = db ((t + 1) / n);
end

function ber = soft_rep_ber_simulated (n, ebn0_db, trials, rng)
  % Soft-decision R_n: add the n received values and take the sign.
  %
  % All n * trials Gaussian draws come in one batch, in the Python's order, and column
  % j of the reshaped batch holds the n values of trial j. py_sum adds down each
  % column with Python's compensated sum, so a total near zero gets the same sign.
  ebn0 = 10 ^ (ebn0_db / 10);
  es = ebn0 / n;
  sigma = sqrt (0.5);
  received = sqrt (es) + reshape (rng.gauss (0.0, sigma, n * trials), n, trials);
  total = py_sum (received);
  errors = sum (total < 0);
  ber = errors / trials;
end


% Monte Carlo ---------------------------------------------------------------------
%
% Octave functions cannot see a script's variables, so each module-level constant of
% the Python is a function here that returns it.

function v = Z95 ()
  v = 1.96;
end

function out = simulate_repetition (n, p, trials, rng)
  % One entry per trial, 1 when majority voting over n copies on BSC(p) fails.
  % Column j of the reshaped batch holds the n draws of trial j, in the Python's order.
  t = floor ((n - 1) / 2);
  flips = sum (reshape (rng.random (n * trials), n, trials) < p, 1);
  out = double (flips > t);
end

function [lo, hi] = wald (k, n, z)
  if (nargin < 3)
    z = Z95;
  end
  ph = k / n;
  hw = z * sqrt (ph * (1 - ph) / n);
  lo = ph - hw;
  hi = ph + hw;
end

function [lo, hi] = wilson (k, n, z)
  if (nargin < 3)
    z = Z95;
  end
  ph = k / n;
  centre = (ph + z * z / (2 * n)) / (1 + z * z / n);
  hw = z / (1 + z * z / n) * sqrt (ph * (1 - ph) / n + z * z / (4 * n * n));
  lo = centre - hw;
  hi = centre + hw;
end

function y = trials_needed (p, rel, z)
  if (nargin < 3)
    z = Z95;
  end
  y = (z / rel) ^ 2 * (1 - p) / p;
end

function v = MAIN_SEED ()
  v = 2026;
end

function v = MAIN_TRIALS ()
  v = 10000;
end

function v = MAIN_P ()
  v = 0.1;
end

function out = main_run ()
  out = simulate_repetition (3, MAIN_P, MAIN_TRIALS, PyRandom (MAIN_SEED));
end

function v = POINT_PS ()
  v = [0.02, 0.05, 0.1, 0.2, 0.3];
end

function v = POINT_SEED ()
  v = 7;
end

function v = POINT_TRIALS ()
  v = 10000;
end

function rows = mc_points ()
  % Estimates of the R_3 word error rate at several p, with exact values and Wilson
  % intervals, one row [p, exact, k, trials, lo, hi] per p in place of the Python's tuples.
  rng = PyRandom (POINT_SEED);
  rows = zeros (0, 6);
  for p = POINT_PS
    k = sum (simulate_repetition (3, p, POINT_TRIALS, rng));
    [lo, hi] = wilson (k, POINT_TRIALS);
    rows(end + 1, :) = [p, word_error(3, 1, p), k, POINT_TRIALS, lo, hi];
  end
end


% Reliability -----------------------------------------------------------------------

function v = HOURS_PER_YEAR ()
  v = 8760;
end

function y = mtbf_hours (fit)
  y = 1e9 / fit;
end


% The checks ------------------------------------------------------------------------

function check_arq ()
  section ('Stop-and-wait, go-back-N and selective repeat');
  n = 1000;
  k = 968;
  R = 1e6;
  t_rt = 10e-3;
  p = 1e-4;
  t_f = n / R;
  a = t_rt / t_f;
  window = 1 + a;
  rate = k / n;
  q = (1 - p) ^ n;
  P = frame_error (n, p);
  check ('frame time in ms', t_f * 1e3, '1');
  check ('a = T_rt / T_f', a, '10');
  check ('window N = 1 + a', window, '11');
  check ('frame efficiency k/n', rate, '0.968');
  check ('clean frame probability (chapter 1)', q, '0.904833');
  check ('frame error probability P', P, '0.095167');
  check ('stop-and-wait slots per delivered frame', (1 + a) / q, '12.1569');
  check ('go-back-N slots per delivered frame', (1 + (window - 1) * P) / q, '2.15694');
  check ('selective repeat slots per delivered frame', 1 / q, '1.10518');
  check ('go-back-N denominator 1 + (N - 1)P', 1 + (window - 1) * P, '1.95167');
  sw = eta_sw (P, a, rate);
  gbn = eta_gbn (P, window, rate);
  sr = eta_sr (P, rate);
  check ('stop-and-wait efficiency', sw, '0.079625');
  check ('go-back-N efficiency', gbn, '0.44878');
  check ('selective repeat efficiency', sr, '0.87588');
  check ('stop-and-wait goodput kbit/s', sw * R / 1e3, '79.6');
  check ('go-back-N goodput kbit/s', gbn * R / 1e3, '448.8');
  check ('selective repeat goodput kbit/s', sr * R / 1e3, '875.9');
  check ('go-back-N over stop-and-wait', gbn / sw, '5.64');
  check ('selective repeat over go-back-N', sr / gbn, '1.95');
  check ('stop-and-wait ceiling k/(n(1 + a)) with P = 0', rate / (1 + a), '0.088');
  % Toy timeline of figure 9-2: a = 3, 16 slots, frame 2 damaged once. Python's //
  % is floor division, and the counts are ints, so they are passed as int64.
  same ('toy timeline, stop-and-wait frames delivered', int64 (floor (16 / (1 + 3)) - 1), int64 (3));
  same ('toy timeline, go-back-N frames delivered (window 4 lost)', int64 (16 - (1 + 3)), int64 (12));
  same ('toy timeline, ideal selective repeat frames delivered', int64 (16 - 1), int64 (15));
  % Held to a window of 4, selective repeat sends frames 3 to 5, resends 2 in slot 6,
  % idles in slots 7 to 9 until that resend is acknowledged, then sends frames 6 to 12
  % in slots 10 to 16.
  same ('toy timeline, selective repeat held to a window of 4', int64 (5 + (16 - 9)), int64 (12));
  % One generator serves all three protocols in turn, as in the Python. fix() is
  % Python's int(), which truncates toward zero.
  rng = PyRandom (11);
  kinds = {'sw', 'gbn', 'sr'};
  formulas = [(1 + a) / q, (1 + (window - 1) * P) / q, 1 / q];
  for i = 1:3
    near (sprintf ('%s slots per frame, simulated over 200,000 frames', kinds{i}), ...
          simulate_protocol (kinds{i}, P, fix (a), 200000, rng), formulas(i), 0.01);
  end
  check ('mean retries P/(1 - P)', P / q, '0.10518');
  check ('one retry costs (1 + a) frame times, ms', (1 + a) * t_f * 1e3, '11');
  check ('mean extra delay in ms', P / q * (1 + a) * t_f * 1e3, '1.16');
  check ('probability of two or more retries', P ^ 2, '0.0090568');
  check ('one frame in this many waits two retries', 1 / P ^ 2, '110');
  check ('extra wait of two retries, ms', 2 * (1 + a) * t_f * 1e3, '22');
  % In-order release: a frame also waits behind earlier damaged frames. The delays
  % are whole slots, so their sum is exact in a double, as Python's int sum is.
  sr_d = sr_release_delays (fix (a), P, 400000, 7);
  check ('selective repeat, simulated mean extra delay, ms', sum (sr_d) / numel (sr_d) * t_f * 1e3, '5.68');
  check ('selective repeat, one frame in this many waits 22 ms or more', ...
         numel (sr_d) / sum (sr_d >= 22), '58');
  gbn_d = gbn_delivery_delays (fix (a), P, 400000, 7);
  near ('go-back-N simulated throughput against 1/S', numel (gbn_d) / 400000, q / (1 + (window - 1) * P), 0.01);
  % Each resend of a frame trails its previous attempt by exactly 1 + a slots, and a
  % frame is sent S_GBN times on average, so the mean extra delay is (1 + a)(S_GBN - 1) T_f.
  s_gbn = (1 + (window - 1) * P) / q;
  gbn_mean_ms = (1 + a) * (s_gbn - 1) * t_f * 1e3;
  check ('go-back-N, mean extra delay (1 + a)(S_GBN - 1) T_f, ms', gbn_mean_ms, '12.7');
  near ('go-back-N simulated mean extra delay against that', sum (gbn_d) / numel (gbn_d) * t_f * 1e3, gbn_mean_ms, 0.02);
  check ('go-back-N, one frame in this many waits 22 ms or more', ...
         numel (gbn_d) / sum (gbn_d >= 22), '3');
end

function check_satellite ()
  section ('Worked example: a satellite hop');
  n = 1000;
  k = 968;
  R = 1e6;
  t_rt = 0.48;
  p = 1e-5;
  a = t_rt / (n / R);
  window = 1 + a;
  P = frame_error (n, p);
  rate = k / n;
  check ('one-way path, km', 2 * 35786, '71,572');
  check ('one-way time, s', 2 * 35786 / 299792.458, '0.239');
  check ('round trip, s', 4 * 35786 / 299792.458, '0.477');
  check ('a', a, '480');
  check ('window', window, '481');
  check ('frame error probability', P, '0.0099502');
  check ('stop-and-wait efficiency', eta_sw (P, a, rate), '0.0019924');
  check ('go-back-N efficiency', eta_gbn (P, window, rate), '0.16592');
  check ('go-back-N denominator', 1 + (window - 1) * P, '5.7761');
  check ('selective repeat efficiency', eta_sr (P, rate), '0.95837');
  check ('frames held after one failure, about 1 + a, in kilobytes', window * n / 8 / 1000, '60.1');
end

function check_residual ()
  section ('The residual error of a checked link');
  cases = {8, 0.01, {'0.850763', '0.922745', '0.074618', '0.0026368', '1.080635', '0.0028494'}; ...
           1000, 1e-4, {'0.818714', '0.904833', '0.090643', '0.0045243', '1.099678', '0.0049753'}};
  for c = 1:rows (cases)
    [n, p, stated] = cases{c, :};
    o = parity_outcomes (n, p);
    clean = o(1);
    alarm = o(2);
    slip = o(3);
    check (sprintf ('n = %d: (1 - 2p)^n', n), (1 - 2 * p) ^ n, stated{1});
    check (sprintf ('n = %d: clean', n), clean, stated{2});
    check (sprintf ('n = %d: alarm', n), alarm, stated{3});
    check (sprintf ('n = %d: slip', n), slip, stated{4});
    check (sprintf ('n = %d: expected transmissions', n), 1 / (1 - alarm), stated{5});
    check (sprintf ('n = %d: residual per accepted frame', n), slip / (1 - alarm), stated{6});
  end
  enumerated = parity_outcomes_enumerated (8, 0.01);
  formula = parity_outcomes (8, 0.01);
  names = {'clean', 'alarm', 'slip'};
  for i = 1:3
    near (sprintf ('n = 8 %s, all 256 patterns enumerated', names{i}), enumerated(i), formula(i), 1e-9);
  end
  [wrong, sends] = simulate_arq_parity (8, 0.01, 400000, PyRandom (5));
  near ('n = 8 residual, simulated over 400,000 frames', wrong, formula(3) / (1 - formula(2)), 0.05);
  near ('n = 8 transmissions, simulated', sends, 1 / (1 - formula(2)), 0.005);
  check ('n = 8: wrong with no check, 1 - c', 1 - 0.99 ^ 8, '0.077255');
  o = parity_outcomes (1000, 1e-4);
  clean = o(1);
  alarm = o(2);
  slip = o(3);
  check ('n = 1000: one accepted frame in this many is wrong', (1 - alarm) / slip, '201');
  check ('n = 1000: frames accepted per second when 1000 are sent', 1000 * (1 - alarm), '909');
  check ('n = 1000: wrong frames accepted per second when 1000 are sent', 1000 * slip, '4.5');
  check ('the same, as 909 accepted frames over 201', 1000 * (1 - alarm) / 201, '4.5');
  check ('n = 1000: share of corrupted frames that slip', slip / (1 - clean), '0.0475');
  check ('2^-32', 2 ^ -32, '2.33e-10');
  P = frame_error (1000, 1e-4);
  undetected = P * 2 ^ -32;
  check ('heuristic undetected probability with 32 check bits', undetected, '2.22e-11');
  check ('mean time to an undetected frame, years', 1 / (1000 * undetected) / 3600 / 24 / 365.25, '1.43');
  % The exact count behind the heuristic: 2^k - 1 of the 2^n - 1 nonzero patterns are
  % codewords. Python divides the exact integers. In doubles, 2^968 - 1 rounds to
  % 2^968 and 2^1000 - 1 to 2^1000, so the quotient is 2^-32 exactly. The exact
  % quotient lies within 2^-968 of 2^-32 in relative terms, so Python's correctly
  % rounded division gives 2^-32 too.
  nn = 1000;
  kk = 968;
  frac = (2 ^ kk - 1) / (2 ^ nn - 1);
  near ('escape fraction (2^k - 1)/(2^n - 1) against 2^-r', frac, 2.0 ^ -(nn - kk), 1e-9);
end

function check_hybrid ()
  section ('Hybrid ARQ');
  n = 1000;
  t = 3;
  rate_crc = 968 / 1000;
  rate_h = 938 / 1000;
  check ('type-I frame efficiency', rate_h, '0.938');
  p = 1e-3;
  q = (1 - p) ^ n;
  check ('p = 1e-3: clean frame probability', q, '0.367695');
  check ('p = 1e-3: frame error probability', 1 - q, '0.632305');
  pf = hybrid_failure (n, t, p);
  check ('p = 1e-3: probability of four or more errors', pf, '0.018927');
  check ('p = 1e-3: hybrid efficiency', (1 - pf) * rate_h, '0.92025');
  check ('p = 1e-3: selective repeat efficiency', eta_sr (1 - q, rate_crc), '0.35593');
  p = 1e-4;
  pf = hybrid_failure (n, t, p);
  check ('p = 1e-4: probability of four or more errors', pf, '3.8e-6');
  check ('p = 1e-4: hybrid efficiency', (1 - pf) * rate_h, '0.93800');
  cross = hybrid_crossover (n, rate_crc, rate_h, t);
  check ('crossover bit error probability', cross, '3.1e-5');
  check ('BCH parity bits m t for m = 10, t = 3', 10 * 3, '30');
  check ('data bits in the hybrid frame', 1000 - 32 - 30, '938');
end

function check_ber_wer ()
  section ('Bit and word error rates');
  p = 0.01;
  dist = zeros (1, 5);
  for j = 0:4
    dist(j + 1) = binom_pmf (4, j, p);  % dist(j + 1) is the Python's dist[j]
  end
  stated = {'0.96059601', '0.03881196', '0.00058806', '3.96e-6', '1e-8'};
  for j = 0:4
    check (sprintf ('uncoded 4-bit word: P(J = %d)', j), dist(j + 1), stated{j + 1});
  end
  pw = 1 - dist(1);
  ej = py_sum ((0:4) .* dist);
  check ('uncoded word error rate', pw, '0.03940399');
  check ('expected wrong bits E[J]', ej, '0.04');
  check ('uncoded bit error rate E[J]/4', ej / 4, '0.01');
  check ('lower bound P_w/k', pw / 4, '0.0098510');
  terms = {2, '2.97e-4'; 3, '1e-6'};
  for r = 1:rows (terms)
    check (sprintf ('R3 term i = %d', terms{r, 1}), binom_pmf (3, terms{r, 1}, p), terms{r, 2});
  end
  r3 = word_error (3, 1, p);
  check ('R3 word error rate', r3, '2.98e-4');
  terms = {3, '9.801e-6'; 4, '4.95e-8'; 5, '1e-10'};
  for r = 1:rows (terms)
    check (sprintf ('R5 term i = %d', terms{r, 1}), binom_pmf (5, terms{r, 1}, p), terms{r, 2});
  end
  r5 = word_error (5, 2, p);
  check ('R5 word error rate', r5, '9.8506e-6');
  check ('R3 over R5', r3 / r5, '30.3');
  % Exhaustive enumeration of every pattern for R3 and R5. The patterns run in the
  % order of the binary numbers 0 to 2^n - 1, as product((0, 1), repeat=n) runs them,
  % and Python's sum() adds the majority-failure terms in that order.
  for nv = [3, 5; r3, r5]
    n = nv(1);
    value = nv(2);
    w = sum (dec2bin (0:2 ^ n - 1, n) == '1', 2).';
    w = w(w > floor (n / 2));
    total = py_sum (p .^ w .* (1 - p) .^ (n - w));
    near (sprintf ('R%d word error rate, all %d patterns enumerated', n, 2 ^ n), total, value, 1e-9);
  end
  block = 1 - (1 - r3) ^ 4;
  check ('4 bits each sent with R3: word error rate', block, '0.0011915');
  check ('4 bits each sent with R3: lower bound P_w/4', block / 4, '2.9787e-4');
  check ('uncoded word errors per R3 word error, 4-bit message', pw / block, '33.1');
end

function check_coding_gain ()
  section ('Coding gain');
  [x, unc] = required_uncoded (1e-6);
  check ('uncoded: Q argument at 1e-6', x, '4.7534');
  check ('uncoded: Eb/N0 ratio at 1e-6', x * x / 2, '11.2975');
  check ('uncoded: Eb/N0 dB at 1e-6', unc, '10.530');
  [p3, x3, es3, eb3] = required_repetition (3, 1e-6);
  check ('R3: channel flip probability', p3, '5.7746e-4');
  check ('R3: Q argument', x3, '3.2498');
  check ('R3: Es/N0 dB', es3, '7.227');
  check ('R3: Eb/N0 dB', eb3, '11.998');
  check ('R3: coding gain dB', unc - eb3, '-1.468');
  [p5, x5, es5, eb5] = required_repetition (5, 1e-6);
  check ('R5: channel flip probability', p5, '4.6524e-3');
  check ('R5: Q argument', x5, '2.6006');
  check ('R5: Es/N0 dB', es5, '5.291');
  check ('R5: Eb/N0 dB', eb5, '12.281');
  check ('R5: coding gain dB', unc - eb5, '-1.751');
  check ('10 log10 3', db (3), '4.771');
  check ('10 log10 5', db (5), '6.990');
  [~, unc5] = required_uncoded (1e-5);
  [~, ~, ~, eb3_5] = required_repetition (3, 1e-5);
  check ('uncoded Eb/N0 dB at 1e-5 (chapter 3)', unc5, '9.588');
  check ('R3 gain at 1e-5 (chapter 3)', unc5 - eb3_5, '-1.441');
  [~, ~, ~, eb5_5] = required_repetition (5, 1e-5);
  check ('R5 gain at 1e-5', unc5 - eb5_5, '-1.718');
  check ('asymptotic hard gain R3', asymptotic_hard_gain (3), '-1.761');
  check ('asymptotic hard gain R5', asymptotic_hard_gain (5), '-2.218');
  check ('asymptotic hard gain R7', asymptotic_hard_gain (7), '-2.430');
  check ('limit 10 log10(1/2)', db (0.5), '-3.010');
  % The finite-target loss creeps toward the asymptote as the target falls.
  prev = 0.0;
  for target = [1e-5, 1e-6, 1e-9, 1e-12, 1e-15]
    [~, u] = required_uncoded (target);
    [~, ~, ~, e3] = required_repetition (3, target);
    loss = u - e3;
    printf ('  R3 gain at %g: %.4f dB\n', target, loss);
    if (~(loss < prev && loss > asymptotic_hard_gain (3)))
      printf ('MISMATCH: the R3 gain does not fall monotonically toward its asymptote\n');
      exit (1);
    end
    prev = loss;
  end
  % Soft-decision R_n equals uncoded BPSK.
  rng = PyRandom (3);
  sim = soft_rep_ber_simulated (3, 2.0, 200000, rng);
  near ('soft R3 at 2 dB, simulated over 200,000 bits', sim, uncoded_ber (2.0), 0.03);
  check ('rate-1/2 soft headroom at 1e-6, dB', unc - 0.19, '10.34');
  check ('rate-1/2 hard headroom at 1e-6, dB', unc - 1.77, '8.76');
  check ('BER of hard R3 at 12 dB', hard_rep_ber (3, 12.0), '9.9e-7');
end

function check_monte_carlo ()
  section ('Measuring error rates by simulation');
  exact = word_error (3, 1, MAIN_P);
  check ('exact R3 word error at p = 0.1', exact, '0.028');
  outcomes = main_run ();
  k = sum (outcomes);
  ph = k / MAIN_TRIALS;
  same ('failures K in 10,000 seeded trials', int64 (k), int64 (266));
  check ('estimate p-hat', ph, '0.0266');
  se = sqrt (ph * (1 - ph) / MAIN_TRIALS);
  check ('estimated standard error', se, '0.0016091');
  check ('true standard error', sqrt (exact * (1 - exact) / MAIN_TRIALS), '0.0016497');
  check ('1.96 times the standard error', Z95 * se, '0.0031539');
  [lo, hi] = wald (k, MAIN_TRIALS);
  check ('Wald lower', lo, '0.023446');
  check ('Wald upper', hi, '0.029754');
  [wlo, whi] = wilson (k, MAIN_TRIALS);
  check ('Wilson lower', wlo, '0.023623');
  check ('Wilson upper', whi, '0.029940');
  check ('Wilson centre', (ph + Z95 ^ 2 / (2 * MAIN_TRIALS)) / (1 + Z95 ^ 2 / MAIN_TRIALS), '0.026782');
  check ('estimate below truth, in standard errors', (exact - ph) / sqrt (exact * (1 - exact) / MAIN_TRIALS), '0.85');
  if (~(lo < exact && exact < hi && wlo < exact && exact < whi))
    printf ('MISMATCH: an interval misses the exact value\n');
    exit (1);
  end
  % Running estimate checkpoints for the chart. outcomes(1:n) is the Python's outcomes[:n].
  checkpoints = {100, '0.01'; 300, '0.02'; 1000, '0.022'; 3000, '0.026333'};
  for r = 1:rows (checkpoints)
    n = checkpoints{r, 1};
    check (sprintf ('running estimate after %d trials', n), sum (outcomes(1:n)) / n, checkpoints{r, 2});
  end
  % Trials for a target precision.
  check ('1 - 2Q(1.96)', 1 - 2 * q_func (Z95), '0.95');
  check ('(1.96/0.2)^2', (Z95 / 0.2) ^ 2, '96.04');
  check ('relative half-width with 100 errors', Z95 / sqrt (100), '0.196');
  check ('trials for 20% at p = 0.028', trials_needed (exact, 0.2), '3,334');
  same ('ceil(100/p) at p = 0.028', int64 (ceil (100 / exact)), int64 (3572));
  check ('trials for 20% at p = 1e-3', trials_needed (1e-3, 0.2), '95,944');
  check ('trials for 20% at p = 1e-6', trials_needed (1e-6, 0.2), '9.604e7');
  windows = {0.028, 3572, '0.958'; 1e-3, 100000, '0.960'; 1e-6, 100000000, '0.960'};
  for r = 1:rows (windows)
    [p, n, stated] = windows{r, :};
    check (sprintf ('exact chance of landing within 20%%, p = %g, N = %s', p, commas (n)), ...
           window_prob (n, p, 0.2), stated);
  end
  check ('the same at p = 0.028, to five places', window_prob (3572, exact, 0.2), '0.95767');
  for pn = [1e-3, 1e-6; 100000, 100000000]
    p = pn(1);
    n = pn(2);
    same (sprintf ('N p is exactly 100 at p = %g', p), py_round (n * p, 9), 100.0);
    check (sprintf ('exact chance without the boundary counts K = 80 and 120, p = %g', p), ...
           window_prob (n, p, 0.2) - exp (binom_log_pmf (n, 80, p)) - exp (binom_log_pmf (n, 120, p)), '0.949');
  end
  % 400 experiments of 3,572 trials each, drawn in turn from one generator.
  rng = PyRandom (99);
  hits = 0;
  for i = 1:400
    kk = sum (simulate_repetition (3, MAIN_P, 3572, rng));
    hits = hits + (abs (kk / 3572 - exact) <= 0.2 * exact);
  end
  same ('seeded experiments of 3,572 trials within 20%, of 400', int64 (hits), int64 (381));
  near ('that share against the exact chance', hits / 400, window_prob (3572, exact, 0.2), 0.03);
  same ('Chebyshev trials for 20% at 95%, p = 0.028, rounded up', ...
        int64 (ceil ((1 - exact) / (0.2 ^ 2 * exact * 0.05))), int64 (17358));
  check ('Chebyshev factor 1/(0.05 r^2) at r = 0.2', 1 / (0.05 * 0.2 ^ 2), '500');
  check ('Chebyshev over the Gaussian count', (1 / (0.04 * 0.05)) / (Z95 / 0.2) ^ 2, '5.2');
  % Zero failures.
  check ('ln 20', log (20), '2.996');
  check ('exact 95% upper bound after 10^6 clean trials', 1 - 0.05 ^ (1 / 1e6), '2.9957e-6');
  % The five points of the second chart. Each row is [p, exact, k, trials, lo, hi].
  points = mc_points ();
  for r = 1:rows (points)
    p = points(r, 1);
    ex = points(r, 2);
    kk = points(r, 3);
    wl = points(r, 5);
    wh = points(r, 6);
    printf ('  p = %s: exact %s, K = %d, Wilson [%s, %s]\n', py_str (p), py_g (ex, 6), kk, py_g (wl, 6), py_g (wh, 6));
    if (~(wl < ex && ex < wh))
      printf ('MISMATCH: the interval at p = %s misses the exact value\n', py_str (p));
      exit (1);
    end
  end
  pts = mc_points ();
  same ('failures at p = 0.02', int64 (pts(1, 3)), int64 (9));
  check ('exact R3 word error at p = 0.02', pts(1, 2), '0.001184');
  same ('failures at p = 0.3', int64 (pts(5, 3)), int64 (2155));
  check ('exact R3 word error at p = 0.3', pts(5, 2), '0.216');
end

function check_wilson_exercise ()
  section ('Exercise: two failures in 1,000 trials');
  [lo, hi] = wald (2, 1000);
  check ('Wald lower', lo, '-0.000769');
  check ('Wald upper', hi, '0.004769');
  [wlo, whi] = wilson (2, 1000);
  check ('Wilson lower', wlo, '0.00055');
  check ('Wilson upper', whi, '0.00726');
end

function check_reliability ()
  section ('Reliability in FIT, SDC and DUE');
  raw = 40000;
  uncorrected = raw * 0.0025;
  due = uncorrected * 0.95;
  sdc = uncorrected * 0.05;
  check ('raw flip rate per hour', raw * 1e-9, '4e-5');
  check ('uncorrected FIT', uncorrected, '100');
  check ('DUE FIT', due, '95');
  check ('SDC FIT', sdc, '5');
  check ('raw MTBF hours', mtbf_hours (raw), '25,000');
  check ('raw MTBF years', mtbf_hours (raw) / HOURS_PER_YEAR, '2.85');
  check ('DUE MTBF hours', mtbf_hours (due), '1.0526e7');
  check ('DUE MTBF years', mtbf_hours (due) / HOURS_PER_YEAR, '1,202');
  check ('SDC MTBF hours', mtbf_hours (sdc), '2e8');
  check ('SDC MTBF years', mtbf_hours (sdc) / HOURS_PER_YEAR, '22,831');
  servers = 100000;
  check ('fleet raw MTBF minutes', mtbf_hours (raw * servers) * 60, '15');
  check ('fleet DUE MTBF hours', mtbf_hours (due * servers), '105.3');
  check ('fleet DUE MTBF days', mtbf_hours (due * servers) / 24, '4.39');
  check ('fleet SDC MTBF hours', mtbf_hours (sdc * servers), '2,000');
  check ('fleet SDC MTBF days', mtbf_hours (sdc * servers) / 24, '83.3');
  lam = sdc * servers * 1e-9;
  check ('fleet SDC rate per hour', lam, '5e-4');
  check ('fleet lambda t over a year', lam * HOURS_PER_YEAR, '4.38');
  check ('fleet chance of at least one SDC in a year', 1 - exp (-lam * HOURS_PER_YEAR), '0.98747');
  check ('one server''s chance of an SDC in a year', 1 - exp (-sdc * 1e-9 * HOURS_PER_YEAR), '4.38e-5');
  % The survival limit (1 - lambda t/m)^m -> exp(-lambda t), checked at a large m.
  lt = 4.38;
  near ('(1 - 4.38/m)^m at m = 10^7 against exp(-4.38)', (1 - lt / 1e7) ^ 1e7, exp (-lt), 1e-5);
  check ('exercise: 500 FIT MTBF hours', mtbf_hours (500), '2e6');
  check ('exercise: 500 FIT MTBF years', mtbf_hours (500) / HOURS_PER_YEAR, '228.3');
  check ('exercise: 20,000 devices at 500 FIT, hours', mtbf_hours (500 * 20000), '100');
end

function check_exercises ()
  section ('Exercises');
  n = 1000;
  R = 10e6;
  t_rt = 0.05e-3;
  p = 1e-5;
  a = t_rt / (n / R);
  P = frame_error (n, p);
  check ('short link: a', a, '0.5');
  check ('short link: stop-and-wait efficiency', eta_sw (P, a, 0.968), '0.63891');
  o = parity_outcomes (4, 0.1);
  clean = o(1);
  alarm = o(2);
  slip = o(3);
  check ('n = 4, p = 0.1: clean', clean, '0.6561');
  check ('n = 4, p = 0.1: alarm', alarm, '0.2952');
  check ('n = 4, p = 0.1: slip', slip, '0.0487');
  check ('n = 4, p = 0.1: residual', slip / (1 - alarm), '0.069098');
  check ('n = 4, p = 0.1: transmissions', 1 / (1 - alarm), '1.41884');
  check ('trials for 10% at 95%, p = 1e-5', trials_needed (1e-5, 0.1), '3.8416e7');
  check ('(1.96/0.1)^2', (Z95 / 0.1) ^ 2, '384.16');
  pw = word_error (5, 2, 0.01);
  check ('trials to size an R5 run at p = 0.01, 100/P_w', 100 / pw, '1.015e7');
  check ('channel bits for that run', 5 * 100 / pw, '5.08e7');
end

function main ()
  printf ('Chapter 9: Error-Control Strategies and Performance\n');
  check_arq ();
  check_satellite ();
  check_residual ();
  check_hybrid ();
  check_ber_wer ();
  check_coding_gain ();
  check_monte_carlo ();
  check_wilson_exercise ();
  check_reliability ();
  check_exercises ();
  printf ('\nAll chapter 9 numbers reproduced.\n');
end

main ();
