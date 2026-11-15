% Reproduce every number in chapter 2 of Error Correction Coding, "Probability for Coding".
%
% This is the GNU Octave twin of python/ch02_probability_for_coding.py. It needs core
% Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch02_probability_for_coding.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. A run that reaches
% the final line has reproduced every number in the chapter. The two routes print
% the same lines, and tools/compare-routes.sh checks that against the Python.
%
% Differences from the Python that a side-by-side reader meets:
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave has one default number type, the double, which holds every integer
%     below 2^53 exactly and no more. Python's ints and Fractions have no limit, and
%     this chapter needs the difference: the exact binomial tails reach denominators
%     such as 10^20000. The section below the helpers builds exact integers and
%     Fractions from rows of base-one-million digits, so those checks stay exact and
%     need no package. Octave cannot overload + and * for them, so the Python's
%     p ** k * (1 - p) ** m reads frac_mul (frac_pow (p, k), ...) here.
%   - Python's math.comb returns an integer of any size. comb below returns a double
%     and stops the script when the result passes 2^53, and binom_tail builds the
%     larger binomial coefficients it needs as exact integers instead.
%   - Octave indexes from 1, so the Python's k-th list entry is entry k + 1. The
%     Python's k values are kept as they are, and 1 is added only where one touches
%     an array.
%   - Python's round (x, n) rounds the exact binary value to n decimals. Octave's
%     round takes no second argument, so rnd rounds through sprintf, which also
%     rounds the exact binary value.
%   - Python's random.Random(seed).random() is reproduced bit for bit by PyRandom
%     in octave/common. The Python draws one number per step. The twin takes all n
%     draws at once, in the same order, which gives the same stream.
%   - Python chains comparisons, as in a < b < c. Octave does not, so the twin joins
%     the separate comparisons with &&. Sets of patterns are cell arrays of strings,
%     combined with intersect, union and setdiff.
%   - The Python defines its own half_unit. The twin uses the one in octave/common,
%     which does the same arithmetic.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Exact integers and Fractions ----------------------------------------------------
%
% A big integer is a row of limbs, least significant first, each below 1e6 apart from
% the last, which carries the sign. The value is the sum of limb(i) * 1e6^(i - 1). A
% limb product is below 1e12, and a convolution adds at most a few thousand of them,
% which stays under 2^53, so doubles hold every intermediate value exactly. That bound
% holds for rows up to about 9000 limbs, and the largest number here is 10^20000, which
% is 3334 limbs. A Fraction is a struct of two big integers and is never reduced. Two
% Fractions are compared by cross-multiplying, which is how Python's exact == sees them.

function v = big_norm (v)
  % Carry so that every limb but the last lies in [0, 1e6), then drop leading zero limbs.
  B = 1e6;
  while (true)
    L = numel (v);
    low = v(1:L - 1);
    r = mod (low, B);                  % mod rounds toward minus infinity, so r is never negative
    c = (low - r) / B;                 % an exact division, since low - r is a multiple of B
    if (~any (c))
      break;
    end
    v(1:L - 1) = r;
    v(2:L) = v(2:L) + c;
  end
  while (abs (v(end)) >= B)
    r = mod (v(end), B);
    c = (v(end) - r) / B;
    v(end) = r;
    v(end + 1) = c;
  end
  while (numel (v) > 1 && v(end) == 0)
    v(end) = [];
  end
end

function v = big_from (x)
  % The big integer for an integer-valued double of magnitude below 2^53.
  assert (x == fix (x) && abs (x) < 2 ^ 53);
  B = 1e6;
  s = sign (x);
  x = abs (x);
  v = 0;
  i = 1;
  while (x > 0)
    v(i) = mod (x, B);
    x = (x - v(i)) / B;
    i = i + 1;
  end
  v = big_norm (s * v);
end

function r = big_sign (a)
  % The sign of a big integer is the sign of its last limb, which is zero only for zero.
  r = sign (a(end));
end

function r = big_add (a, b)
  n = max (numel (a), numel (b));
  a(end + 1:n) = 0;
  b(end + 1:n) = 0;
  r = big_norm (a + b);
end

function r = big_sub (a, b)
  n = max (numel (a), numel (b));
  a(end + 1:n) = 0;
  b(end + 1:n) = 0;
  r = big_norm (a - b);
end

function r = big_mul (a, b)
  r = big_norm (conv (a, b));
end

function r = big_mul_small (a, s)
  % a times an integer s below 1e6, so each limb product stays below 1e12.
  r = big_norm (a * s);
end

function q = big_div_small (a, d)
  % Exact division of a non-negative big integer by an integer d below 1e6. The
  % remainder carried down is below d, so every partial value is below 1e12.
  B = 1e6;
  q = zeros (size (a));
  left = 0;
  for i = numel (a):-1:1
    cur = left * B + a(i);
    q(i) = floor (cur / d);
    left = cur - q(i) * d;
  end
  assert (left == 0);
  q = big_norm (q);
end

function r = big_pow (a, n)
  % a to the non-negative integer power n, by repeated squaring.
  r = 1;
  while (n > 0)
    if (mod (n, 2) == 1)
      r = big_mul (r, a);
    end
    n = floor (n / 2);
    if (n > 0)
      a = big_mul (a, a);
    end
  end
end

function f = Fraction (n, d)
  % Python's Fraction (n, d) for two integer-valued doubles, with d positive.
  f = struct ('num', big_from (n), 'den', big_from (d));
end

function f = frac_add (a, b)
  if (big_sign (big_sub (a.den, b.den)) == 0)      % equal denominators need no cross-multiplying
    f = struct ('num', big_add (a.num, b.num), 'den', a.den);
  else
    f = struct ('num', big_add (big_mul (a.num, b.den), big_mul (b.num, a.den)), ...
                'den', big_mul (a.den, b.den));
  end
end

function f = frac_sub (a, b)
  if (big_sign (big_sub (a.den, b.den)) == 0)
    f = struct ('num', big_sub (a.num, b.num), 'den', a.den);
  else
    f = struct ('num', big_sub (big_mul (a.num, b.den), big_mul (b.num, a.den)), ...
                'den', big_mul (a.den, b.den));
  end
end

function f = frac_mul (a, b)
  f = struct ('num', big_mul (a.num, b.num), 'den', big_mul (a.den, b.den));
end

function f = frac_div (a, b)
  % a / b for a positive b.
  assert (big_sign (b.num) > 0);
  f = struct ('num', big_mul (a.num, b.den), 'den', big_mul (a.den, b.num));
end

function f = frac_pow (a, n)
  f = struct ('num', big_pow (a.num, n), 'den', big_pow (a.den, n));
end

function s = frac_cmp (a, b)
  % The sign of a - b, so -1, 0 or 1.
  s = big_sign (big_sub (big_mul (a.num, b.den), big_mul (b.num, a.den)));
end

function x = big_to_double (a, first)
  % The value of the limbs a(first:end), correctly rounded when it is below 2^53 and
  % within a few units in the last place beyond that. The loop is Horner's rule.
  x = 0;
  for i = numel (a):-1:first
    x = x * 1e6 + a(i);
  end
end

function r = frac_float (f)
  % float (Fraction). Python rounds the exact ratio once. Here a ratio of two numbers
  % below 2^53 is one correctly rounded division, which agrees exactly, and a larger
  % ratio comes from the leading four limbs of each part (a relative error near 1e-16,
  % far below the nine digits any line prints).
  n = f.num;
  d = f.den;
  s = big_sign (n);
  if (s == 0)
    r = 0;
    return;
  end
  n = n * s;                          % the magnitude, once big_norm has carried the flipped limbs
  n = big_norm (n);
  if (numel (n) <= 3 && numel (d) <= 3)
    vn = big_to_double (n, 1);
    vd = big_to_double (d, 1);
    if (vn < 2 ^ 53 && vd < 2 ^ 53)
      r = s * vn / vd;
      return;
    end
  end
  kn = max (numel (n) - 4, 0);
  kd = max (numel (d) - 4, 0);
  r = big_to_double (n, kn + 1) / big_to_double (d, kd + 1);
  e = kn - kd;                        % the limbs dropped from n and from d, in powers of 1e6
  if (e < 0)
    r = r / 10 ^ (-6 * e);
  else
    r = r * 10 ^ (6 * e);
  end
  r = s * r;
end

function r = comb (n, k)
  % math.comb (n, k) for a result below 2^53, so it can join float arithmetic. Each
  % step multiplies before it divides, so every intermediate value is an integer, and
  % int64 holds them all. A larger result stops the script instead of rounding.
  n = int64 (n);
  k = int64 (k);
  r = int64 (1);
  for i = int64 (1):k
    r = idivide (r * (n - k + i), i, 'floor');
  end
  assert (r < 2 ^ 53);
  r = double (r);
end

function r = rnd (value, places)
  % Python's round (value, places) and this sprintf both round the exact binary value.
  r = str2double (sprintf ('%.*f', places, value));
end


% Events on three bits ---------------------------------------------------------

function f = pattern_prob (pattern, p)
  % Probability that a BSC(p) turns 000... into pattern, flips being independent.
  w = sum (pattern == '1');
  f = frac_mul (frac_pow (p, w), frac_pow (frac_sub (Fraction (1, 1), p), numel (pattern) - w));
end

function total = prob (event, p)
  % Probability of a set of received patterns, by additivity over its outcomes.
  total = Fraction (0, 1);
  for i = 1:numel (event)
    total = frac_add (total, pattern_prob (event{i}, p));
  end
end

function r = patterns (n)
  % Every bit string of length n in counting order, as a row of strings.
  r = cellstr (dec2bin (0:2 ^ n - 1, n)).';
end


% The binomial distribution ------------------------------------------------------

function f = binom_pmf (n, k, p)
  % P(X = k) for X the number of flips among n bits of a BSC(p).
  f = frac_mul (Fraction (comb (n, k), 1), ...
                frac_mul (frac_pow (p, k), frac_pow (frac_sub (Fraction (1, 1), p), n - k)));
end

function f = binom_tail (n, t, p)
  % Exact P(X >= t), computed as one minus the head with integer arithmetic.
  % Python's comb (n, k) is exact at any size, and the coefficients reach 10^421 for
  % n = 10000, so they are built here as big integers, each from the one before,
  % through C(n, k) = C(n, k - 1) (n - k + 1) / k, an exact division by a small number.
  % The cells hold the 0-based k at index k + 1.
  num = p.num;
  den = p.den;
  rest = big_sub (den, num);
  cmb = cell (1, t);
  pn = cell (1, t);
  cmb{1} = 1;
  pn{1} = 1;
  for k = 1:t - 1
    cmb{k + 1} = big_div_small (big_mul_small (cmb{k}, n - k + 1), k);
    pn{k + 1} = big_mul (pn{k}, num);
  end
  % The head is summed from k = t - 1 down to 0, because (den - num)^(n - k) then grows
  % by one factor per step and no big division is needed.
  power = big_pow (rest, n - t + 1);
  head = 0;
  for k = t - 1:-1:0
    head = big_add (head, big_mul (big_mul (cmb{k + 1}, pn{k + 1}), power));
    power = big_mul (power, rest);
  end
  total = big_pow (den, n);
  f = struct ('num', big_sub (total, head), 'den', total);
end

function f = rep_error (n, p)
  % Decoded error of the repetition code R_n: more than half of n copies flipped.
  f = binom_tail (n, floor ((n + 1) / 2), p);
end


% Bounds on a tail ---------------------------------------------------------------

function r = markov (mean, a)
  r = mean / a;
end

function r = chebyshev (var, d)
  r = var / d ^ 2;
end

function r = kl (a, p)
  % D(a||p) = a ln(a/p) + (1-a) ln((1-a)/(1-p)), in nats.
  r = a * log (a / p) + (1 - a) * log ((1 - a) / (1 - p));
end

function r = chernoff (n, a, p)
  % The Chernoff bound exp(-n D(a||p)) on P(X >= an), valid for p < a < 1.
  r = exp (-n * kl (a, p));
end


% The Gaussian and the Q-function ---------------------------------------------

function r = phi (x)
  % The standard Gaussian density.
  r = exp (-x * x / 2) / sqrt (2 * pi);
end

function r = Q (x)
  % P(Z > x) for a standard Gaussian Z, through the complementary error function.
  r = 0.5 * erfc (x / sqrt (2));
end

function r = q_lower (x)
  r = x * phi (x) / (1 + x * x);
end

function r = q_upper (x)
  r = phi (x) / x;
end

function r = q_half (x)
  r = 0.5 * exp (-x * x / 2);
end


% The law of large numbers, simulated -------------------------------------------

function out = running_average (seed, n, p)
  % X_k / k for k = 1..n, where X_k counts flips among the first k simulated bits.
  %
  % PyRandom (seed) gives the same sequence as Python's random.Random(seed) for the
  % same integer seed. Its random (n) returns the first n draws as a column, and the
  % cumulative sum of the comparisons is the running count of flips.
  rng = PyRandom (seed);
  flips = cumsum (rng.random (n) < p);
  out = flips ./ (1:n).';
end


% Checks ----------------------------------------------------------------------

function check (label, computed, stated)
  % Print a computed value beside the chapter's, and stop at the first disagreement.
  % A Fraction (a struct) is turned into the nearest double first, as float () does.
  value = str2double (strrep (stated, ',', ''));
  if (isstruct (computed))
    computed = frac_float (computed);
  else
    computed = double (computed);
  end
  printf ('  %-60s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function main ()
  p = Fraction (1, 10);

  printf ('Events on three bits\n');
  space = patterns (3);
  check ('sum of all eight pattern probabilities', prob (space, p), '1');
  check ('P(000)', pattern_prob ('000', p), '0.729');
  check ('P(100), one flip', pattern_prob ('100', p), '0.081');
  check ('P(110), two flips', pattern_prob ('110', p), '0.009');
  check ('P(111)', pattern_prob ('111', p), '0.001');
  A = space(cellfun (@(x) x(1) == '1', space));
  B = space(cellfun (@(x) sum (x == '1') >= 2, space));
  C = space(cellfun (@(x) x(2) == '1', space));
  one = Fraction (1, 1);
  check ('P(A), first bit flipped', prob (A, p), '0.1');
  check ('P(B), decoder fails', prob (B, p), '0.028');
  check ('P(A and B)', prob (intersect (A, B), p), '0.019');
  check ('P(A or B)', prob (union (A, B), p), '0.109');
  check ('P(A) + P(B) - P(A and B)', ...
         frac_sub (frac_add (prob (A, p), prob (B, p)), prob (intersect (A, B), p)), '0.109');
  check ('P(not B)', frac_sub (one, prob (B, p)), '0.972');
  check ('P(B | A)', frac_div (prob (intersect (A, B), p), prob (A, p)), '0.19');
  check ('1 - 0.9^2', frac_sub (one, frac_pow (frac_sub (one, p), 2)), '0.19');
  notA = setdiff (space, A);
  check ('P(B | not A)', frac_div (prob (intersect (B, notA), p), prob (notA, p)), '0.01');
  check ('total probability 0.19(0.1) + 0.01(0.9)', ...
         frac_add (frac_mul (Fraction (19, 100), p), frac_mul (Fraction (1, 100), frac_sub (one, p))), '0.028');
  bayes = frac_div (prob (intersect (A, B), p), prob (B, p));
  assert (frac_cmp (bayes, Fraction (19, 28)) == 0);
  check ('P(A | B) = 19/28', bayes, '0.6786');
  check ('P(A)P(B)', frac_mul (prob (A, p), prob (B, p)), '0.0028');
  check ('P(A and B) / (P(A)P(B))', ...
         frac_div (prob (intersect (A, B), p), frac_mul (prob (A, p), prob (B, p))), '6.79');
  check ('P(A and C), first and second bits flipped', prob (intersect (A, C), p), '0.01');
  assert (frac_cmp (prob (intersect (A, C), p), frac_mul (prob (A, p), prob (C, p))) == 0);

  printf ('The binomial distribution\n');
  three = {'0.729', '0.243', '0.027', '0.001'};
  for k = 0:3
    check (sprintf ('n = 3, p = 0.1: P(X = %d)', k), binom_pmf (3, k, p), three{k + 1});
  end
  q = Fraction (1, 100);
  eight = {'0.922745', '0.074565', '0.0026361', '5.3255e-5', '6.7242e-7'};
  for k = 0:4
    check (sprintf ('n = 8, p = 0.01: P(X = %d)', k), binom_pmf (8, k, q), eight{k + 1});
  end
  check ('n = 8, p = 0.01: P(X >= 5)', binom_tail (8, 5, q), '5.4612e-9');
  check ('n = 8, p = 0.01: P(X >= 1)', binom_tail (8, 1, q), '0.077255');
  odd = Fraction (0, 1);
  for k = 1:2:8
    odd = frac_add (odd, binom_pmf (8, k, q));
  end
  check ('n = 8, p = 0.01: P(X odd)', odd, '0.0746185');
  check ('(1 - 0.98^8) / 2', ...
         frac_div (frac_sub (one, frac_pow (frac_sub (one, frac_mul (Fraction (2, 1), q)), 8)), Fraction (2, 1)), ...
         '0.0746185');
  even2 = Fraction (0, 1);
  for k = 2:2:8
    even2 = frac_add (even2, binom_pmf (8, k, q));
  end
  check ('n = 8, p = 0.01: P(X even and at least 2)', even2, '0.0026368');
  twenty = {'0.1216', '0.2702', '0.2852', '0.1901', '0.0898', '0.0319', '0.0089', '0.0020', '0.0004'};
  for k = 0:8
    check (sprintf ('n = 20, p = 0.1: P(X = %d)', k), binom_pmf (20, k, p), twenty{k + 1});
  end
  check ('n = 20, p = 0.1: (n + 1)p', frac_mul (Fraction (21, 1), p), '2.1');
  check ('n = 20, p = 0.1: P(X >= 5)', binom_tail (20, 5, p), '0.0432');
  ratio = cell (1, 4);
  for k = 0:3
    ratio{k + 1} = frac_div (binom_pmf (20, k + 1, p), binom_pmf (20, k, p));
  end
  check ('n = 20: P(1)/P(0) = (20/1)(1/9)', ratio{1}, '2.22');
  check ('n = 20: P(2)/P(1) = (19/2)(1/9)', ratio{2}, '1.06');
  check ('n = 20: P(3)/P(2) = (18/3)(1/9)', ratio{3}, '0.667');

  printf ('Expectation and variance\n');
  mean3 = Fraction (0, 1);
  sq3 = Fraction (0, 1);
  for k = 0:3
    mean3 = frac_add (mean3, frac_mul (Fraction (k, 1), binom_pmf (3, k, p)));
    sq3 = frac_add (sq3, frac_mul (Fraction (k * k, 1), binom_pmf (3, k, p)));
  end
  check ('n = 3: E[X]', mean3, '0.3');
  check ('n = 3: E[X^2]', sq3, '0.36');
  check ('n = 3: Var X = E[X^2] - E[X]^2', frac_sub (sq3, frac_pow (mean3, 2)), '0.27');
  check ('n = 3: np(1 - p)', frac_mul (frac_mul (Fraction (3, 1), p), frac_sub (one, p)), '0.27');
  direct = Fraction (0, 1);
  for k = 0:3
    direct = frac_add (direct, frac_mul (frac_pow (frac_sub (Fraction (k, 1), mean3), 2), binom_pmf (3, k, p)));
  end
  check ('n = 3: E[(X - 0.3)^2] term by term', direct, '0.27');
  squares = {'0.09', '0.49', '2.89', '7.29'};
  for k = 0:3
    check (sprintf ('n = 3: (k - 0.3)^2 at k = %d', k), frac_pow (frac_sub (Fraction (k, 1), mean3), 2), squares{k + 1});
  end
  check ('n = 3: standard deviation', sqrt (0.27), '0.5196');
  ber = 1e-12;
  hour = 8e9 * 3600;
  check ('bits in an hour at 8 Gb/s', hour, '2.88e13');
  check ('errors per second at 8 Gb/s and p = 1e-12', 8e9 * ber, '0.008');
  check ('hour: mean errors', hour * ber, '28.8');
  check ('hour: standard deviation', sqrt (hour * ber * (1 - ber)), '5.367');
  day = hour * 24;
  check ('bits in a day at 8 Gb/s', day, '6.912e14');
  check ('day: mean errors', day * ber, '691.2');
  check ('day: standard deviation', sqrt (day * ber * (1 - ber)), '26.29');
  check ('terabyte: mean errors, 8e12 bits', 8e12 * ber, '8');
  check ('terabyte: standard deviation', sqrt (8e12 * ber * (1 - ber)), '2.828');
  check ('hour: sd / mean', sqrt (28.8) / 28.8, '0.186');
  check ('day: sd / mean', sqrt (691.2) / 691.2, '0.038');
  check ('terabyte: sd / mean', sqrt (8) / 8, '0.354');

  printf ('Markov, Chebyshev and the law of large numbers\n');
  exact20 = binom_tail (1000, 20, q);
  check ('n = 1000, p = 0.01: mean', 1000 * 0.01, '10');
  check ('n = 1000, p = 0.01: variance', 1000 * 0.01 * 0.99, '9.9');
  check ('n = 1000: Markov bound on P(X >= 20)', markov (10, 20), '0.5');
  check ('n = 1000: Chebyshev bound on P(X >= 20)', chebyshev (9.9, 10), '0.099');
  check ('n = 1000: exact P(X >= 20)', exact20, '0.0032884');
  sd_day = sqrt (691.2 * (1 - ber));
  check ('day: four standard deviations', 4 * sd_day, '105.2');
  check ('day: mean - 4 sd', 691.2 - 4 * sd_day, '586.0');
  check ('day: mean + 4 sd', 691.2 + 4 * sd_day, '796.4');
  check ('Chebyshev at four sd: 1/16', Fraction (1, 16), '0.0625');
  eps_ = 1e-4;                        % Octave owns the name eps (machine epsilon)
  delta = 0.05;
  pm = 1e-3;
  check ('BER p = 1e-3 to 1e-4 at 95%: p(1-p)/(eps^2 delta)', pm * (1 - pm) / (eps_ ^ 2 * delta), '1,998,000');
  check ('n = 1e6, p = 1e-3: sd of the measured fraction', sqrt (pm * (1 - pm) / 1e6), '3.16e-5');
  check ('n = 1e6, p = 1e-3: two sd of the fraction', 2 * sqrt (pm * (1 - pm) / 1e6), '6.3e-5');
  seeds = [1 2 3];
  for i = 1:3
    seed = seeds(i);
    stated = '0.1';
    series = running_average (seed, 10000, 0.1);
    final = series(end);                 % the last entry, which Python writes as [-1]
    assert (abs (final - 0.1) <= 2 * sqrt (0.1 * 0.9 / 10000));
    check (sprintf ('simulated run %d: X/n at n = 10000, to one digit', seed), rnd (final, 1), stated);
  end
  check ('band half-width 2 sqrt(p(1-p)/n) at n = 10000', 2 * sqrt (0.1 * 0.9 / 10000), '0.006');

  printf ('The Gaussian and the Q-function\n');
  xs = [1 2 3 4 5 6];
  qs = {'0.1587', '0.02275', '0.00135', '3.167e-5', '2.867e-7', '9.866e-10'};
  for i = 1:6
    check (sprintf ('Q(%d)', xs(i)), Q (xs(i)), qs{i});
  end
  check ('peak density of N(1, 0.25): 1/(0.5 sqrt(2 pi))', 1 / (0.5 * sqrt (2 * pi)), '0.798');
  check ('phi(0) = 1/sqrt(2 pi)', phi (0), '0.3989');
  check ('sigma = 0.4: Q(2.5)', Q (2.5), '0.00621');
  check ('sigma = 0.25: Q(4)', Q (4), '3.167e-5');
  check ('Q(2) / Q(4)', Q (2) / Q (4), '718');
  % One row per x = 1..6: lower bound, Q(x), upper bound, and the bound exp(-x^2/2)/2.
  table = {'1.210e-1', '1.587e-1', '2.420e-1', '3.033e-1'; ...
           '2.160e-2', '2.275e-2', '2.700e-2', '6.767e-2'; ...
           '1.330e-3', '1.350e-3', '1.477e-3', '5.554e-3'; ...
           '3.149e-5', '3.167e-5', '3.346e-5', '1.677e-4'; ...
           '2.859e-7', '2.867e-7', '2.973e-7', '1.863e-6'; ...
           '9.853e-10', '9.866e-10', '1.013e-9', '7.615e-9'};
  for x = 1:6
    lo = table{x, 1};
    mid = table{x, 2};
    up = table{x, 3};
    half_bound = table{x, 4};
    assert (q_lower (x) < Q (x) && Q (x) < q_upper (x) && Q (x) <= q_half (x));
    check (sprintf ('x = %d: lower bound x phi/(1 + x^2)', x), q_lower (x), lo);
    check (sprintf ('x = %d: Q(x)', x), Q (x), mid);
    check (sprintf ('x = %d: upper bound phi/x', x), q_upper (x), up);
    check (sprintf ('x = %d: bound exp(-x^2/2)/2', x), q_half (x), half_bound);
  end
  check ('x = 6: upper / Q', q_upper (6) / Q (6), '1.026');
  check ('x = 6: Q / lower', Q (6) / q_lower (6), '1.0013');
  n_mean = 10.0;
  n_sd = sqrt (9.9);
  check ('Gaussian stand-in, n = 1000: (20 - 10)/sqrt(9.9)', (20 - n_mean) / n_sd, '3.178');
  check ('Gaussian stand-in, n = 1000: Q(3.178)', Q ((20 - n_mean) / n_sd), '7.41e-4');
  z2 = (200 - 100) / sqrt (99);
  check ('Gaussian stand-in, n = 10000: (200 - 100)/sqrt(99)', z2, '10.05');
  check ('Gaussian stand-in, n = 10000: Q(10.05)', Q (z2), '4.6e-24');
  check ('n = 10000, p = 0.01: exact P(X >= 200)', binom_tail (10000, 200, q), '5.630e-19');

  printf ('The union bound\n');
  m = Fraction (1, 10 ^ 6);
  check ('72 bits, q = 1e-6: union bound 72q', frac_mul (Fraction (72, 1), m), '7.2e-5');
  check ('72 bits: exact P(X >= 1)', binom_tail (72, 1, m), '7.19974e-5');
  check ('72 bits: C(72,2)', comb (72, 2), '2556');
  check ('72 bits: union bound C(72,2) q^2', frac_mul (Fraction (comb (72, 2), 1), frac_pow (m, 2)), '2.556e-9');
  check ('72 bits: exact P(X >= 2)', binom_tail (72, 2, m), '2.55588e-9');
  check ('R3 at p = 0.1: 3p^2', frac_mul (Fraction (3, 1), frac_pow (p, 2)), '0.03');
  check ('R3 at p = 0.1: exact', rep_error (3, p), '0.028');
  check ('R3 at p = 0.1: overcount 2p^3', frac_mul (Fraction (2, 1), frac_pow (p, 3)), '0.002');
  check ('R5 at p = 0.01: 10p^3', frac_mul (Fraction (10, 1), frac_pow (q, 3)), '1e-5');
  check ('R5 at p = 0.01: exact', rep_error (5, q), '9.8506e-6');
  check ('DIMM: 16 x 0.08', 16 * 0.08, '1.28');
  check ('DIMM: 1 - 0.92^16', 1 - 0.92 ^ 16, '0.7366');

  printf ('The Chernoff bound\n');
  a = 0.02;
  es = a * (1 - 0.01) / (0.01 * (1 - a));
  check ('n = 1000, a = 0.02: e^s = a(1-p)/(p(1-a))', es, '2.0204');
  check ('n = 1000, a = 0.02: s', log (es), '0.7033');
  check ('n = 1000, a = 0.02: (1 - p)/(1 - a)', 0.99 / 0.98, '1.0102');
  check ('n = 1000, a = 0.02: D(a||p) in nats', kl (a, 0.01), '0.0039136');
  check ('n = 1000, a = 0.02: n D', 1000 * kl (a, 0.01), '3.9136');
  check ('n = 1000: Chernoff bound on P(X >= 20)', chernoff (1000, a, 0.01), '0.0200');
  check ('n = 10000: Markov bound on P(X >= 200)', markov (100, 200), '0.5');
  check ('n = 10000: Chebyshev bound on P(X >= 200)', chebyshev (99, 100), '0.0099');
  check ('n = 10000: Chernoff bound on P(X >= 200)', chernoff (10000, a, 0.01), '1.008e-17');
  check ('D(1/2 || 0.1)', kl (0.5, 0.1), '0.5108');
  check ('-ln(0.36)/2', -log (0.36) / 2, '0.5108');
  check ('R25 at p = 0.1: bound 0.36^12.5', 0.36 ^ 12.5, '2.843e-6');
  check ('R25 at p = 0.1: exact', rep_error (25, p), '1.621e-7');
  need = 2 * log (1e9) / log (1 / 0.36);
  check ('R_n target 1e-9: 2 ln(1e9)/ln(1/0.36)', need, '40.57');
  check ('R41 bound 0.36^20.5', 0.36 ^ 20.5, '8.02e-10');
  rn = [33 35 41];
  rn_stated = {'2.396e-9', '8.392e-10', '3.636e-11'};
  for i = 1:3
    check (sprintf ('R%d at p = 0.1: exact', rn(i)), rep_error (rn(i), p), rn_stated{i});
  end
  check ('R39 bound 0.36^19.5', 0.36 ^ 19.5, '2.23e-9');

  printf ('Worked examples\n');
  r = Fraction (1, 10 ^ 4);
  check ('sector, 4096 bits at 1e-4: mean', frac_mul (Fraction (4096, 1), r), '0.4096');
  check ('sector: variance', frac_mul (frac_mul (Fraction (4096, 1), r), frac_sub (one, r)), '0.40956');
  check ('sector: P(X = 0)', binom_pmf (4096, 0, r), '0.66390');
  check ('sector: P(X = 1)', binom_pmf (4096, 1, r), '0.27196');
  check ('sector: P(X >= 2)', binom_tail (4096, 2, r), '0.064136');
  check ('sector: Markov bound on P(X >= 2)', markov (0.4096, 2), '0.2048');
  check ('sector: C(4096,2)', comb (4096, 2), '8,386,560');
  check ('sector: union bound C(4096,2) p^2', frac_mul (Fraction (comb (4096, 2), 1), frac_pow (r, 2)), '0.08387');
  check ('phi(3.5)', phi (3.5), '8.727e-4');
  check ('x = 3.5: upper bound phi/x', q_upper (3.5), '2.493e-4');
  check ('x = 3.5: lower bound x phi/(1 + x^2)', q_lower (3.5), '2.305e-4');
  check ('x = 3.5: 3.5/13.25', 3.5 / 13.25, '0.2642');
  check ('Q(3.5)', Q (3.5), '2.326e-4');
  s = Fraction (1, 1000);
  check ('31-bit block, corrects 2, p = 1e-3: C(31,3)', comb (31, 3), '4495');
  check ('31-bit block: union bound C(31,3) p^3', frac_mul (Fraction (comb (31, 3), 1), frac_pow (s, 3)), '4.495e-6');
  check ('31-bit block: exact P(X >= 3)', binom_tail (31, 3, s), '4.4016e-6');
  check ('31-bit block: ratio bound / exact', ...
         frac_div (frac_mul (Fraction (comb (31, 3), 1), frac_pow (s, 3)), binom_tail (31, 3, s)), '1.021');

  printf ('Exercises\n');
  even = space(cellfun (@(x) mod (sum (x == '1'), 2) == 0, space));
  check ('P(even number of flips) at p = 0.1', prob (even, p), '0.756');
  check ('P(A and even)', prob (intersect (A, even), p), '0.018');
  check ('P(A)P(even)', frac_mul (prob (A, p), prob (even, p)), '0.0756');
  half = Fraction (1, 2);
  % Python chains the comparison x == y == z, which means x == y and y == z.
  both = frac_mul (prob (A, half), prob (even, half));
  assert (frac_cmp (prob (intersect (A, even), half), both) == 0 && frac_cmp (both, Fraction (1, 4)) == 0);
  q1000 = frac_pow (frac_sub (one, Fraction (1, 10 ^ 4)), 1000);
  check ('ARQ: q = (1 - 1e-4)^1000', q1000, '0.904833');
  check ('ARQ: E = 1/q', frac_div (one, q1000), '1.105176');
  pe = 1e-6;
  eps6 = 1e-8;
  nbits = pe * (1 - pe) / (eps6 ^ 2 * 0.05);
  check ('BER p = 1e-6 to 1% at 95%: bits', nbits, '199,999,800,000');
  check ('  seconds at 8 Gb/s', nbits / 8e9, '25.0');
  check ('R7 at p = 0.1: union bound C(7,4) p^4', frac_mul (Fraction (comb (7, 4), 1), frac_pow (p, 4)), '0.0035');
  check ('R7 at p = 0.1: exact', rep_error (7, p), '0.002728');
  check ('R7 at p = 0.1: Chernoff 0.36^3.5', 0.36 ^ 3.5, '0.02799');
  mu_tight = 0.4;
  a_tight = 2;
  check ('Markov tight: P(X = a) = mu/a with mu = 0.4, a = 2', mu_tight / a_tight, '0.2');
  % Ratios and small counts the prose states in words.
  patterns5 = 0;
  for k = 5:8
    patterns5 = patterns5 + comb (8, k);
  end
  check ('patterns with 5 or more flips in 8 bits', patterns5, '93');
  check ('one decision error in about this many bits at Q(3.5)', 1 / Q (3.5), '4.30e3');
  exact_1000 = frac_float (binom_tail (1000, 20, q));
  exact_10000 = frac_float (binom_tail (10000, 200, q));
  check ('Markov over exact, n = 1000 (more than a hundred)', 0.5 / exact_1000, '152');
  check ('Chebyshev over Markov improvement, n = 1000', 0.5 / 0.099, '5');
  check ('Chebyshev over exact, n = 1000 (about thirty)', 0.099 / exact_1000, '30');
  check ('exact over Gaussian stand-in, n = 1000 (about four)', exact_1000 / Q (10 / sqrt (9.9)), '4');
  check ('exact over Gaussian stand-in, n = 10000', exact_10000 / Q (100 / sqrt (99)), '1.2e5');
  check ('exact tail ratio between the two columns', exact_1000 / exact_10000, '5.8e15');
  check ('Chernoff over exact, n = 10000 (about 18)', chernoff (10000, 0.02, 0.01) / exact_10000, '18');
  check ('Chebyshev over Chernoff, n = 1000 (five times)', 0.099 / chernoff (1000, 0.02, 0.01), '5');
  check ('n = 10000: exponent nD', 10000 * kl (0.02, 0.01), '39.136');
  check ('R25 bound over exact (about eighteen)', 0.36 ^ 12.5 / frac_float (rep_error (25, p)), '18');
  printf ('All chapter 2 numbers reproduced.\n');
end


main ();
