% Reproduce every number in chapter 7 of Error Correction Coding, "Capacity and Shannon's Theorem".
%
% This is the GNU Octave twin of python/ch07_capacity_and_shannons_theorem.py. It needs
% core Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch07_capacity_and_shannons_theorem.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. Binomial sums are
% built from every binomial term, the small random code is drawn from the same
% fixed seed as the Python's and decoded exactly over every received word, and
% the BPSK capacity is a Simpson-rule integral with 4000 intervals. A run that
% reaches the final line has reproduced every number in the chapter. The two
% routes print the same lines, and tools/compare-routes.sh checks that against
% the Python.
%
% The Python module also holds the capacity, typical-set and random-coding functions
% that ch07_charts.py plots. Charts stay Python-only, and the functions are here
% because the checks use them too.
%
% Differences from the Python that a side-by-side reader meets:
%   - Python's comb and the sums of comb are exact integers, and the probabilities
%     of the typical sets are exact Fractions. Most of these integers are far past
%     2^53, the largest range in which a double holds every integer exactly, and
%     some are past the largest double, near 2^1024. Here comb_scaled holds each
%     binomial coefficient as a mantissa times a power of 2, so nothing overflows,
%     and each probability is a sum of doubles. The exact integers are never
%     printed, only their logarithms and the probabilities, and those agree with
%     the exact values far below the nine digits a check prints.
%   - Python takes log2 of an int too large for a double by splitting it into a
%     mantissa and a power of 2, and big_log2 does the same.
%   - Python computes a tail such as 1 - sum (...) exactly with Fractions. In
%     doubles that subtraction would wipe out most digits of a tail near 1e-5, so
%     the twin adds the tail terms directly instead.
%   - math.lgamma is CPython's own Lanczos approximation, not the C library's
%     lgamma that Octave's gammaln calls, and py_lgamma reproduces it bit for bit.
%     CPython for Apple silicon is compiled so that a multiply followed by an add
%     in that code is one fused operation with a single rounding. Core Octave has
%     no fma function, so fma_emul computes the fused result exactly from two
%     error-free transformations. On a machine whose CPython does not fuse, the
%     two routes can differ in the last bit, far below the digits printed.
%   - NormalDist ().inv_cdf is Wichura's algorithm AS241 from Python's statistics
%     module, which inv_cdf below reproduces with the same fused multiply-adds.
%   - The Python's float loops add left to right, and so does Octave's sum, so
%     the Simpson sum and the decoder's sum over received words are vectorized
%     without changing a bit of the result.
%   - Octave indexes from 1, so the Python's tails[k] is tails(k + 1) here.
%   - A function passed to another function, as to bisect, is a handle such as
%     @binary_entropy.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function check (label, computed, stated)
  computed = double (computed);
  value = str2double (strrep (stated, ',', ''));
  printf ('  %-62s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % Python compares the two values with !=. Here they are compared through their
  % repr text, which for the ints, strings, lists and tuples of this chapter is
  % equal exactly when the values are.
  printf ('  %-62s %-16s %s\n', label, py_str (computed), py_str (stated));
  if (~strcmp (py_str (computed, true), py_str (stated, true)))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed, true), py_str (stated, true));
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

function t = tup (varargin)
  % A Python tuple, in the form py_str prints.
  t = struct ('tuple', {varargin});
end

function c = ilist (v)
  % A Python list of ints, as a cell row so that one item keeps its brackets.
  c = num2cell (int64 (v(:).'));
end


% Exact arithmetic the Python gets from its ints ---------------------------------

function [s, err] = two_sum (a, b)
  % s = a + b rounded, and err the exact rounding error, so that s + err = a + b.
  s = a + b;
  bb = s - a;
  err = (a - (s - bb)) + (b - bb);
end

function [hi, lo] = split_double (a)
  % Dekker's split of a double into two halves of 26 bits each, hi + lo = a.
  c = 134217729 * a;
  hi = c - (c - a);
  lo = a - hi;
end

function [p, err] = two_prod (a, b)
  % p = a b rounded, and err the exact rounding error, so that p + err = a b.
  p = a .* b;
  [ah, al] = split_double (a);
  [bh, bl] = split_double (b);
  err = ((ah .* bh - p) + ah .* bl + al .* bh) + al .* bl;
end

function f = fma_emul (a, b, c)
  % a b + c with a single rounding, as the C function fma returns it.
  % The exact sum p + e + c is gathered into f and a remainder d, and the last
  % lines settle the one case that rounding f alone gets wrong: d exactly half a
  % unit of f, with more of the exact value beyond it in the same direction.
  [p, e] = two_prod (a, b);
  [s, t] = two_sum (p, c);
  [u, w] = two_sum (t, e);
  [r, z] = two_sum (s, u);
  [zw, ww] = two_sum (z, w);
  [f, d] = two_sum (r, zw);
  [m, ~] = log2 (abs (f));
  sp = eps (f);
  % Just below a power of 2, the unit toward zero is half as large.
  down = (sign (d) ~= sign (f)) & (m == 0.5);
  sp(down) = sp(down) / 2;
  tie = (d ~= 0) & (abs (d) == sp / 2) & (ww ~= 0) & (sign (ww) == sign (d));
  f(tie) = f(tie) + 2 * d(tie);
end

function r = py_lgamma (x)
  % CPython's math.lgamma for x > 0, the Lanczos sum with the coefficients and
  % the order of operations of CPython's mathmodule.c.
  NUM = [23531376880.410759688572007674451636754734846804940, ...
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
  DEN = [0.0, 39916800.0, 120543840.0, 150917976.0, 105258076.0, 45995730.0, ...
         13339535.0, 2637558.0, 357423.0, 32670.0, 1925.0, 66.0, 1.0];
  g = 6.024680040776729583740234375;
  num = zeros (size (x));
  den = zeros (size (x));
  % Below 5 the sum is Horner's rule in x, each step a fused multiply-add.
  lo = x < 5;
  xl = x(lo);
  nl = zeros (size (xl));
  dl = nl;
  for i = 13:-1:1
    nl = fma_emul (nl, xl, NUM(i));
    dl = fma_emul (dl, xl, DEN(i));
  end
  % From 5 on it is Horner's rule in 1/x, a division and an add, which do not fuse.
  xh = x(~lo);
  nh = zeros (size (xh));
  dh = nh;
  for i = 1:13
    nh = nh ./ xh + NUM(i);
    dh = dh ./ xh + DEN(i);
  end
  num(lo) = nl;
  den(lo) = dl;
  num(~lo) = nh;
  den(~lo) = dh;
  r = log (num ./ den) - g;
  r = fma_emul (x - 0.5, log ((x + g) - 0.5) - 1, r);
  % CPython returns exactly 0 at 1 and 2.
  r(x == floor (x) & x <= 2) = 0;
end

function x = inv_cdf (p)
  % NormalDist ().inv_cdf (p): Wichura's algorithm AS241, with the coefficients
  % of Python's statistics module and each Horner step a fused multiply-add.
  q = p - 0.5;
  x = zeros (size (p));
  c = abs (q) <= 0.425;
  qc = q(c);
  r = fma_emul (-qc, qc, 0.180625);
  A = [2.5090809287301226727e+3, 3.3430575583588128105e+4, 6.7265770927008700853e+4, ...
       4.5921953931549871457e+4, 1.3731693765509461125e+4, 1.9715909503065514427e+3, ...
       1.3314166789178437745e+2, 3.3871328727963666080e+0];
  B = [5.2264952788528545610e+3, 2.8729085735721942674e+4, 3.9307895800092710610e+4, ...
       2.1213794301586595867e+4, 5.3941960214247511077e+3, 6.8718700749205790830e+2, ...
       4.2313330701600911252e+1, 1.0];
  num = A(1) * ones (size (r));
  den = B(1) * ones (size (r));
  for i = 2:8
    num = fma_emul (num, r, A(i));
    den = fma_emul (den, r, B(i));
  end
  x(c) = (num .* qc) ./ den;
  % The tails, where p is within 0.075 of 0 or of 1.
  t = ~c;
  rr = p(t);
  up = q(t) > 0;
  rr(up) = 1.0 - rr(up);
  rr = sqrt (-log (rr));
  C1 = [7.7454501427834140764e-4, 2.2723844989269184583e-2, 2.4178072517745061177e-1, ...
        1.2704582524523683826e+0, 3.6478483247632046050e+0, 5.7694972214606914055e+0, ...
        4.6303378461565452959e+0, 1.4234371107496835773e+0];
  D1 = [1.0507500716444168432e-9, 5.4759380849953449460e-4, 1.5198666563616457197e-2, ...
        1.4810397642748007459e-1, 6.8976733498510000455e-1, 1.6763848301838038494e+0, ...
        2.0531916266377588219e+0, 1.0];
  C2 = [2.0103343992922881326e-7, 2.7115555687434875782e-5, 1.2426609473880784386e-3, ...
        2.6532189526576123093e-2, 2.9656057182850489123e-1, 1.7848265399172913358e+0, ...
        5.4637849111641143699e+0, 6.6579046435011037772e+0];
  D2 = [2.0442631033899397856e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5, ...
        7.8686913114561325910e-4, 1.4875361290850614852e-2, 1.3692988092273580531e-1, ...
        5.9983220655588793769e-1, 1.0];
  small = rr <= 5;
  rr(small) = rr(small) - 1.6;
  rr(~small) = rr(~small) - 5.0;
  % One column of coefficients per tail value, C1 and D1 up to 5 and C2 and D2 beyond.
  CN = repmat (C1(:), 1, numel (rr));
  DN = repmat (D1(:), 1, numel (rr));
  CN(:, ~small) = repmat (C2(:), 1, sum (~small));
  DN(:, ~small) = repmat (D2(:), 1, sum (~small));
  rr = rr(:).';
  num = CN(1, :);
  den = DN(1, :);
  for i = 2:8
    num = fma_emul (num, rr, CN(i, :));
    den = fma_emul (den, rr, DN(i, :));
  end
  xt = num ./ den;
  qt = q(t);
  xt(qt(:).' < 0) = -xt(qt(:).' < 0);
  x(t) = xt;
end

function [M, E] = comb_scaled (n)
  % comb (n, k) for k = 0 to n, as M .* 2 .^ E with each M in [0.5, 1).
  % The Python's comb (n, k) is an exact int, and for n = 10000 it reaches 2^9996,
  % beyond the largest double. Each step of comb (n, k + 1) = comb (n, k) (n - k) / (k + 1)
  % rounds by at most a unit in the 53rd bit, and log2 splits the result back into
  % a mantissa and an exponent so that nothing overflows. While the coefficients
  % stay below 2^53 every step is exact. The rows are kept, because the ball-size
  % check asks for the same n five hundred times.
  persistent memo
  if (isempty (memo))
    memo = containers.Map ('KeyType', 'double', 'ValueType', 'any');
  end
  if (isKey (memo, n))
    v = memo(n);
    M = v{1};
    E = v{2};
    return;
  end
  M = zeros (n + 1, 1);
  E = zeros (n + 1, 1);
  % comb (n, 0) = 1 = 0.5 * 2^1.
  m = 0.5;
  e = 1;
  M(1) = m;
  E(1) = e;
  for k = 0:n - 1
    m = m * (n - k) / (k + 1);
    [m, de] = log2 (m);
    e = e + de;
    M(k + 2) = m;
    E(k + 2) = e;
  end
  memo(n) = {M, E};
end

function b = big_sum (M, E)
  % The sum of M .* 2 .^ E, returned as [mantissa, exponent] in the same form.
  top = max (E);
  s = sum (M .* 2 .^ (E - top));
  [m, de] = log2 (s);
  b = [m, top + de];
end

function v = big_log2 (b)
  % Python's log2 of a positive int, here held as b(1) * 2^b(2).
  if (b(2) <= 1024)
    % The int fits a double, and Python takes log2 of that double. It is written
    % as (2 b(1)) 2^(b(2) - 1) so that 2^1024, which overflows, is never formed.
    v = log2 (pow2 (2 * b(1), b(2) - 1));
  else
    % Beyond the largest double CPython splits the int into a mantissa in [0.5, 1)
    % and a power of 2, and adds the exponent to log2 of the mantissa.
    v = b(2) + log2 (b(1));
  end
end

function v = big_int (b)
  % A number held as b(1) * 2^b(2) that is an integer below 2^53, as a Python int.
  v = int64 (pow2 (b(1), b(2)));
end

function pmf = binomial_pmf (n, p)
  % P(K = k) for k = 0 to n and K ~ Binomial(n, p), with p the Fraction [a b].
  % comb (n, k) p^k (1 - p)^(n - k) is formed through its base-2 logarithm, since
  % for n = 10000 the coefficient overflows a double and 0.9^9000 underflows one.
  % The logarithm is within about 1e-12 of the exact one, so each probability is
  % within about one part in 10^12 of the Python's exact Fraction.
  [M, E] = comb_scaled (n);
  k = (0:n).';
  L = E + k * log2 (p(1) / p(2)) + (n - k) * log2 ((p(2) - p(1)) / p(2));
  pmf = M .* 2 .^ L;
end


% Information measures ----------------------------------------------------------------

function h = binary_entropy (p)
  % Written for a vector of p as well, for the grid search in check_examples.
  h = zeros (size (p));
  in = p > 0 & p < 1;
  q = p(in);
  h(in) = -q .* log2 (q) - (1 - q) .* log2 (1 - q);
end

function v = db (x)
  v = 10 * log10 (x);
end

function v = q_function (x)
  v = 0.5 * erfc (x / sqrt (2));
end

function v = q_inverse (prob)
  v = -inv_cdf (prob);
end

function v = log2_comb (n, k)
  [M, E] = comb_scaled (n);
  v = big_log2 ([M(k + 1), E(k + 1)]);
end


% Capacities of the binary channels ------------------------------------------------

function v = bsc_output_one (q, p)
  % P(Y = 1) for a BSC(p) whose input is 1 with probability q.
  v = q * (1 - p) + (1 - q) * p;
end

function v = bsc_mutual_information (q, p)
  v = binary_entropy (bsc_output_one (q, p)) - binary_entropy (p);
end

function v = bec_mutual_information (q, eps)
  v = (1 - eps) * binary_entropy (q);
end

function v = z_mutual_information (q, p)
  % Z-channel: a 1 turns into 0 with probability p, a 0 is never disturbed.
  % Written for a vector of q, for the grid search in check_examples.
  v = binary_entropy (q .* (1 - p)) - q .* binary_entropy (p);
end

function v = z_capacity (p)
  v = log2 (1 + 2 ^ (-binary_entropy (p) / (1 - p)));
end

function v = z_best_output (p)
  v = 1 / (1 + 2 ^ (binary_entropy (p) / (1 - p)));
end

function v = bsc_capacity (p)
  v = 1 - binary_entropy (p);
end

function x = bisect (f, lo, hi, target, steps)
  % The x in [lo, hi] where the increasing function f reaches target.
  if (nargin < 5)
    steps = 80;
  end
  for i = 1:steps
    mid = (lo + hi) / 2;
    if (f (mid) < target)
      lo = mid;
    else
      hi = mid;
    end
  end
  x = (lo + hi) / 2;
end


% The Gaussian channel ------------------------------------------------------------

function v = gaussian_capacity (snr)
  % Bits per real dimension, snr = P/N.
  v = 0.5 * log2 (1 + snr);
end

function c = bpsk_soft_capacity (snr, intervals)
  % 1 - E[log2(1 + exp(-2Y/s2))] for Y ~ N(1, s2) and s2 = 1/snr, by Simpson's rule on +-12 sigma.
  % The Python loops over the 4001 points. Here every point is computed at once,
  % with the same operations in the same order, and sum adds the terms left to
  % right as the loop does.
  if (nargin < 2)
    intervals = 4000;
  end
  s2 = 1 / snr;
  s = sqrt (s2);
  a = 1 - 12 * s;
  b = 1 + 12 * s;
  step = (b - a) / intervals;
  i = (0:intervals).';
  y = a + i * step;
  weight = 2 * ones (intervals + 1, 1);
  weight(mod (i, 2) == 1) = 4;
  weight([1, end]) = 1;
  z = -2 * y / s2;
  f = zeros (size (z));
  pos = z > 0;
  f(pos) = (z(pos) + log1p (exp (-z(pos)))) / log (2);
  f(~pos) = log1p (exp (z(~pos))) / log (2);
  total = sum (weight .* f .* exp (-(y - 1) .^ 2 / (2 * s2)) / sqrt (2 * pi * s2));
  c = 1 - total * step / 3;
end

function v = bpsk_hard_capacity (snr)
  v = bsc_capacity (q_function (sqrt (snr)));
end

function v = min_ebn0_db (capacity, rate)
  % Smallest Eb/N0 in dB at which capacity(snr) reaches rate, with Eb/N0 = snr / (2 rate).
  snr = bisect (capacity, 1e-6, 1e3, rate, 70);
  v = db (snr / (2 * rate));
end

function v = shannon_hartley (bandwidth, snr)
  v = bandwidth * log2 (1 + snr);
end

function v = ebn0_boundary (eta)
  % Smallest Eb/N0 (as a ratio) at spectral efficiency eta bits per second per hertz.
  % Written for a vector of eta, for the grid in check_shannon_limit. A scalar
  % raised to a vector of powers calls pow for each one, as Python's ** does.
  v = (2 .^ eta - 1) ./ eta;
end


% Typical sequences ---------------------------------------------------------------

function [ks, set_size, prob, atypical] = typical_counts (n, p, window)
  % Ones counts k with |k/n - p| <= window, the number of such sequences and their total probability.
  % p and window are Fractions [numerator denominator]. The Python's size is
  % set_size here, since size is an Octave function, and it comes back as
  % [mantissa, exponent]. The fourth output is the probability outside the set,
  % added directly, for the checks that the Python makes on 1 - prob.
  k = (0:n).';
  % |k/n - a/b| <= c/d multiplied through by n b d, which keeps the test in exact
  % integer arithmetic as the Python's Fractions do.
  inside = abs (k * p(2) - p(1) * n) * window(2) <= window(1) * n * p(2);
  ks = k(inside).';
  [M, E] = comb_scaled (n);
  set_size = big_sum (M(inside), E(inside));
  pmf = binomial_pmf (n, p);
  prob = sum (pmf(inside));
  atypical = sum (pmf(~inside));
end

function total = typical_probability_float (n, p, window)
  % The same probability through log-gamma, fast enough for the chart.
  k = (0:n).';
  k = k(abs (k / n - p) <= window + 1e-12);
  terms = exp (py_lgamma (n + 1) - py_lgamma (k + 1) - py_lgamma (n - k + 1) + k * log (p) + (n - k) * log (1 - p));
  total = sum (terms);
end


% Random coding on the BSC --------------------------------------------------------

function b = ball_size (n, r)
  % sum (comb (n, k) for k = 0 to r), as [mantissa, exponent].
  [M, E] = comb_scaled (n);
  b = big_sum (M(1:r + 1), E(1:r + 1));
end

function t = binomial_tail (n, p, r)
  % P(K > r) for K ~ Binomial(n, p). The Python subtracts the exact head from 1.
  % Here the terms for k = r + 1 to n are added directly, which keeps every digit
  % of a small tail.
  pmf = binomial_pmf (n, p);
  t = sum (pmf(r + 2:end));
end

function c = log_add (a, b)
  if (a == -Inf)
    c = b;
    return;
  end
  if (b == -Inf)
    c = a;
    return;
  end
  m = max (a, b);
  c = m + log1p (exp (-abs (a - b)));
end

function best = sphere_decoder_bound (n, rate, p)
  % Average error of a random code of 2^(nR) words with the radius-r decoder, minimized over r.
  %
  % The bound is P(K > r) + (2^(nR) - 1) |ball(r)| / 2^n, each term computed in the natural-log domain.
  k = (0:n).';
  log_pmf = py_lgamma (n + 1) - py_lgamma (k + 1) - py_lgamma (n - k + 1) + k * log (p) + (n - k) * log (1 - p);
  log_binom = py_lgamma (n + 1) - py_lgamma (k + 1) - py_lgamma (n - k + 1);
  % tails(k + 1) is the Python's tails[k].
  tails = -Inf (n + 2, 1);
  for k = n:-1:0
    tails(k + 1) = log_add (tails(k + 2), log_pmf(k + 1));
  end
  best = 1.0;
  log_ball = -Inf;
  log_words = n * rate * log (2) + log1p (-2.0 ^ (-n * rate));
  for r = 0:floor (n / 2)
    log_ball = log_add (log_ball, log_binom(r + 1));
    % Python's int () truncates toward zero, which fix does.
    if (r < fix (n * p))
      continue;
    end
    union = log_words + log_ball - n * log (2);
    tail = tails(r + 2);
    if (union < 0)
      value = exp (tail) + exp (union);
    else
      value = exp (tail) + 1.0;
    end
    best = min (best, value);
  end
end

function code = random_code (seed, words, n)
  % PyRandom draws the same words as Python's random.Random (seed), and a column
  % of draws is the same as one draw after another.
  rng = PyRandom (seed);
  code = rng.getrandbits (n, words);
end

function w = weight (x)
  % The number of 1 bits in each entry of x, which bin (x).count ("1") counts.
  w = zeros (size (x));
  while (any (x(:)))
    w = w + bitand (x, 1);
    x = bitshift (x, -1);
  end
end

function pe = ml_block_error (code, n, p)
  % Exact block error of maximum-likelihood (minimum-distance) decoding on a BSC(p), ties split evenly.
  % Row y + 1 of dists holds the distances from the received word y to every
  % codeword, so the Python's loop over y becomes one matrix.
  w = 0:n;
  pw = p .^ w .* (1 - p) .^ (n - w);
  y = (0:2 ^ n - 1).';
  dists = weight (bsxfun (@bitxor, y, code(:).'));
  % Each of the w codewords tied at the minimum distance d wins y with probability 1/w, and
  % each has likelihood pw[d], so y adds w * pw[d] / w = pw[d] to the summed correct probability.
  correct = sum (pw(min (dists, [], 2) + 1));
  pe = 1 - correct / numel (code);
end


% Fano and the converse ------------------------------------------------------------

function v = fano_bound (pe, m)
  v = binary_entropy (pe) + pe * log2 (m - 1);
end

function v = converse_bound (capacity, rate, n)
  v = 1 - capacity / rate - 1 / (n * rate);
end

function v = dispersion_bsc (p)
  v = p * (1 - p) * log2 ((1 - p) / p) ^ 2;
end

function v = normal_approximation (n, eps, p)
  v = bsc_capacity (p) - sqrt (dispersion_bsc (p) / n) * q_inverse (eps);
end


% The checks ---------------------------------------------------------------------

function check_counting ()
  section ('Counting noise patterns');
  h = binary_entropy (0.1);
  check ('h(0.1)', h, '0.4690');
  check ('1 - h(0.1)', 1 - h, '0.5310');
  tbl = {100, '43.98', '0.4398', '46.90'; 1000, '464.42', '0.4644', '469.00'; ...
          10000, '4683.72', '0.4684', '4689.96'};
  for i = 1:rows (tbl)
    [n, stated, per, nh] = tbl{i, :};
    check (sprintf ('log2 C(%d, %d)', n, n / 10), log2_comb (n, n / 10), stated);
    check (sprintf ('log2 C(%d, %d) / n', n, n / 10), log2_comb (n, n / 10) / n, per);
    check (sprintf ('%d h(0.1)', n), n * h, nh);
  end
  check ('1000 - 1000 h(0.1), exponent of the codeword count', 1000 - 1000 * h, '531.00');
  check ('chapter 6 source over BSC(0.1), uses per letter', h / (1 - h), '0.8832');
end

function check_binary_channels ()
  section ('The binary symmetric and erasure channels');
  check ('h(0.11)', binary_entropy (0.11), '0.4999');
  check ('capacity of BSC(0.11)', bsc_capacity (0.11), '0.5001');
  check ('crossover with capacity exactly 0.5', bisect (@(p) binary_entropy (p), 0.01, 0.5, 0.5), '0.110028');
  tbl = {0.1, '0.188', '0.6973', '0.1974', '0.2345'; 0.2, '0.266', '0.8357', '0.3358', '0.3610'; ...
          0.3, '0.344', '0.9286', '0.4287', '0.4406'; 0.4, '0.422', '0.9824', '0.4825', '0.4855'; ...
          0.5, '0.500', '1.0000', '0.5001', '0.5000'};
  for i = 1:rows (tbl)
    [q, r, hr, i_bsc, i_bec] = tbl{i, :};
    check (sprintf ('q = %s: P(Y = 1) for BSC(0.11)', py_str (q)), bsc_output_one (q, 0.11), r);
    check (sprintf ('q = %s: h(P(Y = 1))', py_str (q)), binary_entropy (bsc_output_one (q, 0.11)), hr);
    check (sprintf ('q = %s: I(X;Y) for BSC(0.11)', py_str (q)), bsc_mutual_information (q, 0.11), i_bsc);
    check (sprintf ('q = %s: I(X;Y) for BEC(0.5)', py_str (q)), bec_mutual_information (q, 0.5), i_bec);
  end
  check ('BSC(0.1) at q = 0.5', bsc_mutual_information (0.5, 0.1), '0.5310');
end

function check_gaussian ()
  section ('The Gaussian channel');
  % e is Octave's name for the double nearest e, the same double as Python's math.e.
  check ('uniform, unit variance, h in bits', log2 (2 * sqrt (3)), '1.7925');
  check ('Laplace, unit variance, h in bits', log2 (sqrt (2) * e), '1.9427');
  check ('Gaussian, unit variance, h in bits', 0.5 * log2 (2 * pi * e), '2.0471');
  check ('Gaussian minus uniform', 0.5 * log2 (2 * pi * e) - log2 (2 * sqrt (3)), '0.2546');
  check ('half log2(pi e / 6)', 0.5 * log2 (pi * e / 6), '0.2546');
  check ('C at P/N = 3', gaussian_capacity (3), '1');
  check ('C at P/N = 15', gaussian_capacity (15), '2');
  check ('C at P/N = 63', gaussian_capacity (63), '3');
  check ('P/N = 3 in dB', db (3), '4.77');
  check ('P/N = 15 in dB', db (15), '11.76');
  check ('P/N = 63 in dB', db (63), '17.99');
  check ('factor 4 in dB', db (4), '6.02');
  tbl = {-10, '0.0688', '0.0687', '0.0449'; -5, '0.1982', '0.1977', '0.1353'; 0, '0.5000', '0.4859', '0.3689'; ...
          5, '1.0287', '0.8592', '0.7685'; 10, '1.7297', '0.9968', '0.9908'};
  for i = 1:rows (tbl)
    [snr_db, cg, cs, ch] = tbl{i, :};
    snr = 10 ^ (snr_db / 10);
    check (sprintf ('%d dB: Gaussian input', snr_db), gaussian_capacity (snr), cg);
    check (sprintf ('%d dB: BPSK, soft', snr_db), bpsk_soft_capacity (snr), cs);
    check (sprintf ('%d dB: BPSK, hard', snr_db), bpsk_hard_capacity (snr), ch);
  end
  % Python's all () stops at the first SNR that fails, and so does this loop.
  ordered = true;
  for d = -10:10
    snr = 10 ^ (d / 10);
    soft = bpsk_soft_capacity (snr);
    if (~(gaussian_capacity (snr) >= soft && soft >= bpsk_hard_capacity (snr)))
      ordered = false;
      break;
    end
  end
  require ('Gaussian >= soft >= hard at every listed SNR', ordered);
end

function check_minimum_ebn0 ()
  section ('Minimum Eb/N0 by rate');
  % One row per rate, with the Python's three dicts gauss, soft and hard as the
  % second, third and fourth columns.
  tbl = {0.5, '0.00', '0.19', '1.77'; 1 / 3, '-0.55', '-0.50', '1.21'; ...
          0.1, '-1.29', '-1.29', '0.59'; 0.01, '-1.56', '-1.56', '0.39'};
  for i = 1:rows (tbl)
    [rate, gauss, soft, hard] = tbl{i, :};
    check (sprintf ('R = %.3g: Gaussian input, (2^(2R) - 1)/(2R) in dB', rate), db (ebn0_boundary (2 * rate)), gauss);
    check (sprintf ('R = %.3g: BPSK soft', rate), min_ebn0_db (@bpsk_soft_capacity, rate), soft);
    check (sprintf ('R = %.3g: BPSK hard', rate), min_ebn0_db (@bpsk_hard_capacity, rate), hard);
  end
  soft_half = min_ebn0_db (@bpsk_soft_capacity, 0.5);
  hard_half = min_ebn0_db (@bpsk_hard_capacity, 0.5);
  check ('rate 1/2, soft, three places', soft_half, '0.187');
  check ('rate 1/2, hard, three places', hard_half, '1.773');
  check ('rate 1/2, hard minus soft', hard_half - soft_half, '1.585');
  snr_hard = bisect (@bpsk_hard_capacity, 1e-6, 1e3, 0.5, 70);
  check ('rate 1/2, crossover of the hard BSC', q_function (sqrt (snr_hard)), '0.1100');
  check ('low-rate hard-decision penalty, 10 log10(pi/2)', db (pi / 2), '1.96');
  check ('rate 0.01, hard minus soft', min_ebn0_db (@bpsk_hard_capacity, 0.01) - min_ebn0_db (@bpsk_soft_capacity, 0.01), ...
         '1.95');
end

function check_bandwidth ()
  section ('Bandwidth and power');
  snr = 10 ^ 2.5;
  check ('25 dB as a ratio', snr, '316.23');
  check ('log2(1 + 316.23)', log2 (1 + snr), '8.3094');
  check ('20 MHz at 25 dB, Mbit/s', shannon_hartley (20, snr), '166.19');
  check ('half the S/N', snr / 2, '158.11');
  check ('log2(1 + 158.11)', log2 (1 + snr / 2), '7.3139');
  check ('40 MHz at the same power, Mbit/s', shannon_hartley (40, snr / 2), '292.56');
  check ('double the S/N', snr * 2, '632.46');
  check ('log2(1 + 632.46)', log2 (1 + 2 * snr), '9.3071');
  check ('20 MHz at double power, Mbit/s', shannon_hartley (20, 2 * snr), '186.14');
  tbl = {1e5, '0.3459'; 1e6, '1.0000'; 1e7, '1.3750'; 1e8, '1.4355'};
  for i = 1:rows (tbl)
    [b, stated] = tbl{i, :};
    % C's %.0e prints 1e+05, as Python's format (b, '.0e') does.
    check (sprintf ('S/N0 = 1e6, B = %.0e, Mbit/s', b), shannon_hartley (b, 1e6 / b) / 1e6, stated);
  end
  check ('log2 e, the wideband limit in Mbit/s', log2 (e), '1.4427');
  % The Python writes the first five values of eta as ints and the rest as
  % floats, and its labels print them that way, so int64 marks the ints here.
  tbl = {int64(8), '31.875', '15.03'; int64(6), '10.500', '10.21'; int64(4), '3.750', '5.74'; ...
          int64(2), '1.500', '1.76'; int64(1), '1.000', '0.00'; 0.5, '0.828', '-0.82'; ...
          0.25, '0.757', '-1.21'; 0.1, '0.718', '-1.44'; 0.01, '0.696', '-1.58'; 0.001, '0.693', '-1.59'};
  for i = 1:rows (tbl)
    [eta, ratio, in_db] = tbl{i, :};
    check (sprintf ('eta = %s: (2^eta - 1)/eta', py_str (eta)), ebn0_boundary (double (eta)), ratio);
    check (sprintf ('eta = %s: in dB', py_str (eta)), db (ebn0_boundary (double (eta))), in_db);
  end
  check ('eta = 3, in dB', db (ebn0_boundary (3)), '3.68');
  check ('uncoded BPSK at 1e-5, dB (chapter 3)', db (q_inverse (1e-5) ^ 2 / 2), '9.59');
  check ('uncoded BPSK minus the bound at eta = 2', db (q_inverse (1e-5) ^ 2 / 2) - db (ebn0_boundary (2)), '7.83');
  check ('uncoded BPSK minus the Shannon limit', db (q_inverse (1e-5) ^ 2 / 2) - db (log (2)), '11.18');
end

function check_shannon_limit ()
  section ('The Shannon limit');
  check ('ln 2', log (2), '0.6931');
  check ('ln 2, six places', log (2), '0.693147');
  check ('log10(ln 2)', log10 (log (2)), '-0.15917');
  check ('10 log10(ln 2)', db (log (2)), '-1.5917');
  check ('10 log10(ln 2), two places', db (log (2)), '-1.59');
  check ('(ln 2)^2 / 2, slope of the series', log (2) ^ 2 / 2, '0.2402');
  check ('eta = 0.1 from the series, ln 2 + 0.1 (ln 2)^2/2', log (2) + 0.1 * log (2) ^ 2 / 2, '0.7172');
  check ('eta = 0.1 exactly', ebn0_boundary (0.1), '0.7177');
  % boundary(k) is the Python's ebn0_boundary (10 ** (-3 + (k - 1) / 100)) for k = 1 to 401.
  boundary = ebn0_boundary (10 .^ (-3 + (0:400) / 100));
  require ('the boundary decreases as eta falls, on a grid from 10 to 0.001', all (boundary(1:400) < boundary(2:401)));
end

function check_typical ()
  section ('Typical sequences');
  p = [1, 10];
  h = binary_entropy (0.1);
  check ('log2 9', log2 (9), '3.1699');
  check ('delta = 0.05 log2 9', 0.05 * log2 (9), '0.1585');
  tbl = {0, '1', '0.12158', '0.12158', '0.1520'; 1, '20', '0.013509', '0.27017', '0.3105'; ...
          2, '190', '0.0015009', '0.28518', '0.4690'; 3, '1140', '0.00016677', '0.19012', '0.6275'; ...
          4, '4845', '0.000018530', '0.08978', '0.7860'; 5, '15504', '0.0000020589', '0.03192', '0.9445'};
  for i = 1:rows (tbl)
    [k, count, each, total, rate] = tbl{i, :};
    % p ** k * (1 - p) ** (20 - k) is the exact Fraction 9^(20 - k) / 10^20. Its
    % numerator passes 2^53 for small k and is rounded once before the division,
    % which moves the value by about one part in 10^16.
    each_exact = p(1) ^ k * (p(2) - p(1)) ^ (20 - k) / p(2) ^ 20;
    same (sprintf ('n = 20, k = %d: C(20, k)', k), int64 (nchoosek (20, k)), int64 (str2double (count)));
    check (sprintf ('n = 20, k = %d: probability of each', k), each_exact, each);
    check (sprintf ('n = 20, k = %d: total probability', k), nchoosek (20, k) * each_exact, total);
    check (sprintf ('n = 20, k = %d: -(1/n) log2 p(x)', k), -log2 (each_exact) / 20, rate);
    check (sprintf ('n = 20, k = %d: h + (k/n - p) log2 9', k), h + (k / 20 - 0.1) * log2 (9), rate);
  end
  check ('2^(20 h)', 2 ^ (20 * h), '666.25');
  check ('2^(20 (h - delta))', 2 ^ (20 * (h - 0.05 * log2 (9))), '74.03');
  check ('2^(20 (h + delta))', 2 ^ (20 * (h + 0.05 * log2 (9))), '5996.2');
  [ks, set_size, prob] = typical_counts (20, p, [1, 20]);
  same ('n = 20: typical ones counts', ilist (ks), ilist ([1, 2, 3]));
  same ('n = 20: |T|', big_int (set_size), int64 (1350));
  check ('n = 20: log2 |T|', big_log2 (set_size), '10.399');
  check ('n = 20: P(T)', prob, '0.74547');
  check ('n = 20: n (h - delta)', 20 * (h - 0.05 * log2 (9)), '6.21');
  check ('n = 20: n (h + delta)', 20 * (h + 0.05 * log2 (9)), '12.55');
  tbl = {100, 5, 15, '0.93640', '58.086', '31.05', '62.75'; ...
          1000, 50, 150, '0.99999972', '605.30', '310.50', '627.49'};
  for i = 1:rows (tbl)
    [n, lo, hi, pt, lt, b_lo, b_hi] = tbl{i, :};
    [ks, set_size, prob] = typical_counts (n, p, [1, 20]);
    same (sprintf ('n = %d: range of typical k', n), tup (int64 (ks(1)), int64 (ks(end))), tup (int64 (lo), int64 (hi)));
    check (sprintf ('n = %d: P(T)', n), prob, pt);
    check (sprintf ('n = %d: log2 |T|', n), big_log2 (set_size), lt);
    check (sprintf ('n = %d: n (h - delta)', n), n * (h - 0.05 * log2 (9)), b_lo);
    check (sprintf ('n = %d: n (h + delta)', n), n * (h + 0.05 * log2 (9)), b_hi);
  end
  [ks, set_size, prob, atypical] = typical_counts (10000, p, [1, 20]);
  check ('n = 10000: log2 |T|', big_log2 (set_size), '6092.2');
  % prob > 1 - 1e-40 is the same as the probability outside the set being below
  % 1e-40. That probability is near 1e-53 here, far from the line.
  require ('n = 10000: P(T) above 1 - 1e-40', atypical < 1e-40);
  [ks, set_size, prob] = typical_counts (10, p, [1, 20]);
  same ('n = 10: typical ones counts', ilist (ks), ilist (1));
  % p * (1 - p) ** 9 is 9^9 / 10^10, and both fit a double exactly.
  check ('n = 10: probability of each typical sequence', 9 ^ 9 / 10 ^ 10, '0.038742');
  check ('n = 10: P(T)', prob, '0.38742');
  check ('n = 10: 2^(10 h)', 2 ^ (10 * h), '25.81');
end

function check_random_coding ()
  section ('The random coding argument');
  n = 1000;
  r = 140;
  p = [1, 10];
  check ('2^400', 2.0 ^ 400, '2.58e120');
  tail = binomial_tail (n, p, r);
  check ('P(K > 140), n = 1000, p = 0.1', tail, '2.424e-5');
  lg = big_log2 (ball_size (n, r)) - n;
  check ('log2 P(a wrong codeword within 140)', lg, '-420.29');
  check ('1 - h(0.14)', 1 - binary_entropy (0.14), '0.4158');
  check ('bound -1000 (1 - h(0.14))', -n * (1 - binary_entropy (0.14)), '-415.76');
  check ('union exponent, exact', 400 + lg, '-20.29');
  check ('union term, exact', 2.0 ^ (400 + lg), '7.8e-7');
  check ('union exponent from the h bound', 400 - n * (1 - binary_entropy (0.14)), '-15.76');
  check ('union term from the h bound', 2.0 ^ (400 - n * (1 - binary_entropy (0.14))), '1.8e-5');
  total = tail + 2.0 ^ (400 + lg);
  check ('average error bound, both terms', total, '2.50e-5');
  check ('P(K > 125)', binomial_tail (n, p, 125), '0.00453');
  check ('log2 P(wrong within 125) + 400', 400 + big_log2 (ball_size (n, 125)) - n, '-60.93');
  require ('ball bound sum C(n,k), k <= rho n, <= 2^(n h(rho)) for n = 1000, rho = 0.14', ...
           big_log2 (ball_size (1000, 140)) <= 1000 * binary_entropy (0.14));
  below = true;
  for rr = 1:500
    if (~(big_log2 (ball_size (1000, rr)) <= 1000 * binary_entropy (rr / 1000) + 1e-9))
      below = false;
      break;
    end
  end
  require ('same bound for n = 1000 and every rho = r/n up to 1/2', below);
  check ('sphere-decoder bound, n = 1000, R = 0.4', sphere_decoder_bound (1000, 0.4, 0.1), '2.1e-5');
end

function check_small_code ()
  section ('A small random code decoded exactly');
  code = random_code (1, 16, 12);
  same ('seed 1, first codeword', dec2bin (code(1), 12), '001000100110');
  same ('seed 1, last codeword', dec2bin (code(end), 12), '001101011011');
  same ('distinct codewords', int64 (numel (unique (code))), int64 (16));
  % nchoosek lists the index pairs in the order itertools.combinations does.
  pairs_at = nchoosek (1:16, 2);
  dists = weight (bitxor (code(pairs_at(:, 1)), code(pairs_at(:, 2))));
  same ('minimum distance', int64 (min (dists)), int64 (2));
  same ('pairs at distance 2', int64 (sum (dists == 2)), int64 (2));
  close = find (dists == 2);
  pairs = cell (1, numel (close));
  keys = cell (1, numel (close));
  for i = 1:numel (close)
    pair = sort ({dec2bin(code(pairs_at(close(i), 1)), 12), dec2bin(code(pairs_at(close(i), 2)), 12)});
    pairs{i} = tup (pair{:});
    % Both strings have 12 characters, so sorting the joined strings sorts the
    % pairs the way Python sorts tuples.
    keys{i} = [pair{:}];
  end
  [~, order] = sort (keys);
  pairs = pairs(order);
  same ('the two close pairs', pairs, {tup('110000101100', '110000111000'), tup('110010011110', '110011010110')});
  check ('block error of ML decoding, n = 12, M = 16', ml_block_error (code, 12, 0.1), '0.10539');
  check ('uncoded 4 bits, block error', 1 - 0.9 ^ 4, '0.3439');
  check ('R3 per bit, 3p^2(1-p) + p^3', 3 * 0.01 * 0.9 + 0.001, '0.028');
  check ('R3 on 4 bits, block error', 1 - 0.972 ^ 4, '0.10738');
  % A fresh generator from seed 1, as the Python's second call to random_code makes.
  code16 = random_code (1, 16, 16);
  same ('n = 16: minimum distance', int64 (min (weight (bitxor (code16(pairs_at(:, 1)), code16(pairs_at(:, 2)))))), int64 (3));
  check ('n = 16, R = 1/4: block error', ml_block_error (code16, 16, 0.1), '0.05220');
end

function check_fano ()
  section ('Fano''s inequality and the converse');
  check ('log2 3', log2 (3), '1.5850');
  check ('Fano, M = 4, Pe = 0.1', fano_bound (0.1, 4), '0.6275');
  check ('h(0.05)', binary_entropy (0.05), '0.2864');
  check ('0.05 log2 7', 0.05 * log2 (7), '0.1404');
  check ('Fano, M = 8, Pe = 0.05', fano_bound (0.05, 8), '0.4268');
  c = bsc_capacity (0.1);
  check ('1 - C/R, BSC(0.1), R = 0.6', 1 - c / 0.6, '0.1150');
  check ('converse, n = 100', converse_bound (c, 0.6, 100), '0.0983');
  check ('converse, n = 1000', converse_bound (c, 0.6, 1000), '0.1133');
  check ('BSC(0.2) capacity', bsc_capacity (0.2), '0.2781');
  check ('converse limit, BSC(0.2), R = 1/3', 1 - bsc_capacity (0.2) / (1 / 3), '0.1658');
  check ('R3 bit error on BSC(0.2)', 3 * 0.04 * 0.8 + 0.008, '0.104');
  check ('converse limit, BSC(0.11), R = 0.75', 1 - bsc_capacity (0.11) / 0.75, '0.3332');
end

function check_sources ()
  section ('Sources over channels');
  check ('h(0.2)', binary_entropy (0.2), '0.7219');
  check ('source p = 0.2 over BSC(0.1), uses per letter', binary_entropy (0.2) / bsc_capacity (0.1), '1.3596');
  check ('D / sigma^2 at P/N = 15, one use', 1 / 16, '0.0625');
  check ('in dB', db (1 / 16), '-12.04');
  check ('D / sigma^2 at P/N = 15, two uses', 1 / 256, '0.00390625');
  check ('in dB', db (1 / 256), '-24.08');
  check ('uncoded scaling with MMSE: N/(P+N) at P/N = 15', 1 / (1 + 15), '0.0625');
end

function check_open_questions ()
  section ('What the theorem leaves open');
  v = dispersion_bsc (0.1);
  check ('dispersion V of BSC(0.1)', v, '0.90436');
  check ('sqrt V', sqrt (v), '0.95098');
  check ('Q^-1(1e-3)', q_inverse (1e-3), '3.0902');
  check ('Q^-1(1e-6)', q_inverse (1e-6), '4.7534');
  tbl = {100, '0.2371', '0.0790'; 1000, '0.4381', '0.3881'; 10000, '0.5016', '0.4858'};
  for i = 1:rows (tbl)
    [n, a, b] = tbl{i, :};
    check (sprintf ('normal approximation, n = %d, eps = 1e-3', n), normal_approximation (n, 1e-3, 0.1), a);
    check (sprintf ('normal approximation, n = %d, eps = 1e-6', n), normal_approximation (n, 1e-6, 0.1), b);
  end
  check ('dropped term log2(1000) / 2000', log2 (1000) / 2000, '0.0050');
  check ('fraction of capacity at n = 1000, eps = 1e-3', normal_approximation (1000, 1e-3, 0.1) / bsc_capacity (0.1), ...
         '0.825');
end

function check_examples ()
  section ('Worked examples');
  p = 0.1;
  check ('Z-channel: h(p)/(1 - p)', binary_entropy (p) / (1 - p), '0.52111');
  check ('Z-channel: 2^(-h/(1-p))', 2 ^ (-binary_entropy (p) / (1 - p)), '0.69684');
  check ('Z-channel: best P(Y = 1)', z_best_output (p), '0.41067');
  check ('Z-channel: best q', z_best_output (p) / (1 - p), '0.45630');
  check ('Z-channel: capacity, closed form', z_capacity (p), '0.76285');
  check ('Z-channel: I at the best q', z_mutual_information (z_best_output (p) / (1 - p), p), '0.76285')
  check ('Z-channel: I at q = 0.5', z_mutual_information (0.5, p), '0.7583');
  % max returns the first index of the largest value, as Python's max with a key
  % returns the first item, and grid point i is index i since the grid starts at 1.
  [~, best_grid] = max (z_mutual_information ((1:999) / 1000, p));
  same ('Z-channel: best q on a grid of 0.001', int64 (best_grid), int64 (456));
  check ('BEC(0.25) capacity', 0.75, '0.75');
  check ('uses for 1000 bits on BEC(0.25)', 1000 / 0.75, '1333.3');
end

function check_prose_details ()
  % Numbers the prose states in passing, outside the tables and the main derivations.
  section ('Numbers stated in passing');
  h = binary_entropy (0.1);
  % Counting noise patterns: the shortfall n h - log2 C(n, n/10) grows far more slowly than n.
  tbl = {100, '2.92'; 1000, '4.57'; 10000, '6.23'};
  for i = 1:rows (tbl)
    [n, stated] = tbl{i, :};
    check (sprintf ('shortfall n h - log2 C(n, n/10), n = %d', n), n * h - log2_comb (n, n / 10), stated);
  end
  check ('BSC crossover for capacity one half, mirrored', 1 - bisect (@binary_entropy, 0.0, 0.5, 0.5), '0.890');
  check ('BSC(0.11): P(Y = 1) term 0.89', 1 - 0.11, '0.89');
  check ('flip cost h(0.11)/0.11', binary_entropy (0.11) / 0.11, '4.54');
  % Gaussian and BPSK.
  check ('Gaussian capacity at 20 dB', gaussian_capacity (10 ^ 2), '3.33');
  check ('pi/2', pi / 2, '1.571');
  check ('soft over hard BPSK capacity at -10 dB', bpsk_soft_capacity (0.1) / bpsk_hard_capacity (0.1), '1.53');
  check ('table ratio 0.0687/0.0449', 0.0687 / 0.0449, '1.53');
  % Bandwidth against power at 25 dB.
  snr = 10 ^ 2.5;
  base = shannon_hartley (20, snr);
  check ('doubled bandwidth: S/N + 1', snr / 2 + 1, '159.11');
  check ('doubled power: S/N + 1', 2 * snr + 1, '633.46');
  check ('bandwidth gain, percent', 100 * (shannon_hartley (40, snr / 2) / base - 1), '76');
  check ('power gain, percent', 100 * (shannon_hartley (20, 2 * snr) / base - 1), '12');
  % Typical sequences.
  check ('delta = 0.05 log2 9', 0.05 * log2 (9), '0.1585');
  check ('variance p(1 - p)', 0.1 * 0.9, '0.09');
  check ('window squared', 0.05 ^ 2, '0.0025');
  check ('Chebyshev constant 0.09/0.0025', 0.09 / 0.0025, '36');
  check ('Chebyshev bound at n = 1000', 36 / 1000, '0.036');
  % The Python's float (1 - prob1000), added directly as the probability outside the set.
  [~, ~, ~, atypical1000] = typical_counts (1000, [1, 10], [1, 20]);
  check ('exact atypical probability at n = 1000', atypical1000, '2.8e-7');
  check ('n = 10: 2^(10 h) over the 10 typical blocks', 2 ^ (10 * h) / 10, '2.6');
  check ('typical probability at n = 10, window 0.05', typical_probability_float (10, 0.1, 0.05), '0.39');
  check ('typical probability at n = 100, window 0.05', typical_probability_float (100, 0.1, 0.05), '0.94');
  require ('typical probability above 0.999 by n = 400', typical_probability_float (400, 0.1, 0.05) > 0.999);
  check ('typical probability at n = 1000, window 0.02', typical_probability_float (1000, 0.1, 0.02), '0.97');
  % Random coding.
  check ('codebook size in bits, 2^400 x 1000', 2.0 ^ 400 * 1000, '2.58e123');
  pe = sphere_decoder_bound (1000, 0.4, 0.1);
  check ('expurgated rate 0.4 - 1/n', 0.4 - 1 / 1000, '0.399');
  check ('expurgated maximum error 2 x 2.50e-5', 2 * 2.50e-5, '5.0e-5');
  require ('rate 0.30 bound below 1e-16 near n = 1400', sphere_decoder_bound (1400, 0.30, 0.1) < 1e-16 ...
           && sphere_decoder_bound (1200, 0.30, 0.1) > 1e-16);
  check ('rate 0.45 bound at n = 4000', sphere_decoder_bound (4000, 0.45, 0.1), '1.1e-7');
  require ('rate 0.40 bound at n = 1000 is the figure''s 2.1e-5', abs (pe - 2.1e-5) < 0.05e-5);
  check ('small code: 2^12 received words', 2 ^ 12, '4096');
  % Fano and the converse.
  c2 = bsc_capacity (0.2);
  check ('converse on BSC(0.2), R = 1/3, n = 3', converse_bound (c2, 1 / 3, 3), '-0.8342');
  check ('1/(nR) at n = 3, R = 1/3', 1 / (3 * (1 / 3)), '1');
  check ('R3 bit success on BSC(0.2)', 1 - (3 * 0.2 ^ 2 * 0.8 + 0.2 ^ 3), '0.896');
  c1 = bsc_capacity (0.1);
  % The Python's chained a < b < c is a < b and b < c.
  middle = converse_bound (c1, 0.6, 2000);
  require ('converse bound rises with n toward its limit', ...
           converse_bound (c1, 0.6, 1000) < middle && middle < 1 - c1 / 0.6);
  % Sources over channels.
  check ('SNR with two repetitions of SNR 15', 2 * 15, '30');
  check ('uncoded error with two uses', 1 / 31, '0.0323');
  % What the theorem leaves open.
  require ('n = 100 at error 1e-3 reaches less than half of capacity', ...
           normal_approximation (100, 1e-3, 0.1) < 0.5 * c1);
  % Worked examples.
  check ('Z-channel: 1 + 2^(-s)', 1 + 2 ^ (-binary_entropy (0.1) / 0.9), '1.69684');
  check ('Z-channel: loss of the uniform input', z_capacity (0.1) - z_mutual_information (0.5, 0.1), '0.0046');
  same ('uses for 1000 bits on BEC(0.25), rounded up', int64 (ceil (1000 / 0.75)), int64 (1334));
end


function main ()
  check_counting ();
  check_binary_channels ();
  check_gaussian ();
  check_minimum_ebn0 ();
  check_bandwidth ();
  check_shannon_limit ();
  check_typical ();
  check_random_coding ();
  check_small_code ();
  check_fano ();
  check_sources ();
  check_open_questions ();
  check_examples ();
  check_prose_details ();
  printf ('All chapter 7 numbers reproduced.\n');
end


main ();
