% Reproduce every number in chapter 1 of Error Correction Coding, "Why Codes Exist".
%
% This is the GNU Octave twin of python/ch01_why_codes_exist.py. It needs core
% Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch01_why_codes_exist.m
%
% Each check prints the value the chapter quotes and stops with an error if the
% computation disagrees with it, so a run that reaches the last line means every
% number in the chapter was reproduced. The two routes print the same lines, and
% tools/compare-routes.sh checks that against the Python.
%
% Differences from the Python that a side-by-side reader meets:
%   - Octave indexes from 1, so a Python level[-1] is level(end). The Python's bit
%     positions and error weights are kept as they are, and 1 is added only where
%     one touches an array.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave already owns the names flip and xor, so the Python's flip is flip_bit
%     here. Octave functions return several values with [a, b] = f (x), where the
%     Python returns a tuple, and as_tuple builds a tuple for show to print.
%   - Octave has one default number type, the double, which holds every integer
%     below 2^53 exactly and no more. Python's ints and Fractions have no limit, and
%     this chapter needs the difference: the exact checks reach denominators such as
%     10^390. The section below the helpers builds exact integers and Fractions from
%     rows of base-one-million digits, so those checks stay exact and need no
%     package.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's.
%   - Python's round (x, n) rounds the exact binary value to n decimals. Octave's
%     round takes no second argument, so rnd rounds through sprintf, which also
%     rounds the exact binary value.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function show (label, value)
  printf ('  %-58s %s\n', label, py_str (value));
end

function t = as_tuple (varargin)
  % A Python tuple, which py_str prints with parentheses.
  t = struct ('tuple', {varargin});
end

function r = isclose (a, b, varargin)
  % math.isclose (a, b, rel_tol=1e-09, abs_tol=0.0), with the tolerances as name-value pairs.
  rel_tol = 1e-9;
  abs_tol = 0.0;
  for i = 1:2:numel (varargin)
    if (strcmp (varargin{i}, 'rel_tol'))
      rel_tol = varargin{i + 1};
    else
      abs_tol = varargin{i + 1};
    end
  end
  r = abs (a - b) <= max (rel_tol * max (abs (a), abs (b)), abs_tol);
end

function total = py_sum (x)
  % Python's sum of floats: Neumaier's compensated summation, as CPython 3.12 and later add.
  total = 0;
  c = 0;
  for i = 1:numel (x)
    v = x(i);
    t = total + v;
    if (abs (total) >= abs (v))
      c = c + ((total - t) + v);
    else
      c = c + ((v - t) + total);
    end
    total = t;
  end
  if (c ~= 0 && isfinite (c))
    total = total + c;
  end
end

function r = rnd (value, places)
  % Python's round (value, places) and this sprintf both round the exact binary value.
  r = str2double (sprintf ('%.*f', places, value));
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


% Exact integers and Fractions ----------------------------------------------------
%
% A big integer is a row of limbs, least significant first, each below 1e6 apart from
% the last, which carries the sign. The value is the sum of limb(i) * 1e6^(i - 1). A
% limb product is below 1e12, and a convolution adds at most a few thousand of them,
% which stays under 2^53, so doubles hold every intermediate value exactly. A Fraction
% is a struct of two big integers and is never reduced. Two Fractions are compared by
% cross-multiplying, which is how Python's exact == and < see them too.

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
  % far below the digits any line prints).
  n = f.num;
  d = f.den;
  s = big_sign (n);
  if (s == 0)
    r = 0;
    return;
  end
  n = n * s;                          % the magnitude, still normalized, since the limbs of -n flip sign
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


% One flipped bit ------------------------------------------------------------

function r = flip_bit (word, bit)
  % Return word with bit number bit inverted (bit 0 is the least significant).
  r = bitxor (word, bitshift (1, bit));
end

function r = as_signed8 (u)
  % Read an 8-bit pattern as a two's-complement integer.
  if (u >= 128)
    r = u - 256;
  else
    r = u;
  end
end


% How often errors happen -----------------------------------------------------

function r = p_at_least_one (p, n)
  % Probability that n independent bits, each wrong with probability p, hold one error or more.
  %
  % Computed as -expm1(n * log1p(-p)) so that p = 1e-12 does not round to zero.
  r = -expm1 (n * log1p (-p));
end


% Repetition codes ------------------------------------------------------------

function r = rep_encode (bits, n)
  if (nargin < 2)
    n = 3;
  end
  r = strjoin (arrayfun (@(b) repmat (b, 1, n), bits, 'UniformOutput', false), ' ');
end

function r = rep_decode_block (block)
  % Majority vote over one block of a repetition code of odd length.
  if (sum (block == '1') > floor (numel (block) / 2))
    r = '1';
  else
    r = '0';
  end
end

function r = decode_blocks (recv)
  % The Python's "".join (rep_decode_block (b) for b in recv.split ()): the space-separated
  % blocks decoded one by one and joined. (cellfun cannot gather char results directly.)
  votes = cellfun (@rep_decode_block, strsplit (recv, ' '), 'UniformOutput', false);
  r = [votes{:}];
end

function r = rep_error (p, n)
  % Probability that majority decoding of R_n returns the wrong bit on a BSC(p).
  %
  % The vote is wrong when more than half of the n copies are flipped. Works on
  % floats and on Fractions, so identities can be checked exactly. A Fraction is
  % a struct here, which is how the function tells the two apart.
  assert (mod (n, 2) == 1);
  ks = (n + 1) / 2:n;
  if (isstruct (p))
    one = Fraction (1, 1);
    r = Fraction (0, 1);
    for k = ks
      term = frac_mul (frac_mul (Fraction (comb (n, k), 1), frac_pow (p, k)), ...
                       frac_pow (frac_sub (one, p), n - k));
      r = frac_add (r, term);
    end
  else
    terms = zeros (1, numel (ks));
    for i = 1:numel (ks)
      k = ks(i);
      terms(i) = comb (n, k) * p ^ k * (1 - p) ^ (n - k);
    end
    r = py_sum (terms);
  end
end

function r = rep_bound (p, n)
  % The bound (4p(1-p))^(n/2), proved in the chapter for p < 1/2.
  r = (4 * p * (1 - p)) ^ (n / 2);
end

function r = majority (a, b, c)
  r = bitor (bitor (bitand (a, b), bitand (b, c)), bitand (c, a));
end


% Information ---------------------------------------------------------------

function r = binary_entropy (p)
  if (p == 0 || p == 1)
    r = 0.0;
    return;
  end
  r = -p * log2 (p) - (1 - p) * log2 (1 - p);
end


% Single parity ---------------------------------------------------------------

function r = parity_undetected_direct (p, n)
  % Sum over every nonzero even error weight: those are the patterns parity cannot see.
  ws = 2:2:n;
  if (isstruct (p))
    one = Fraction (1, 1);
    r = Fraction (0, 1);
    for w = ws
      term = frac_mul (frac_mul (Fraction (comb (n, w), 1), frac_pow (p, w)), ...
                       frac_pow (frac_sub (one, p), n - w));
      r = frac_add (r, term);
    end
  else
    terms = zeros (1, numel (ws));
    for i = 1:numel (ws)
      w = ws(i);
      terms(i) = comb (n, w) * p ^ w * (1 - p) ^ (n - w);
    end
    r = py_sum (terms);
  end
end

function r = parity_undetected_closed (p, n)
  % The closed form [1 + (1-2p)^n] / 2 - (1-p)^n derived in the chapter.
  if (isstruct (p))
    one = Fraction (1, 1);
    two = Fraction (2, 1);
    a = frac_pow (frac_sub (one, frac_mul (two, p)), n);
    r = frac_sub (frac_div (frac_add (one, a), two), frac_pow (frac_sub (one, p), n));
  else
    r = (1 + (1 - 2 * p) ^ n) / 2 - (1 - p) ^ n;
  end
end

function [gates, depth] = xor_tree (n)
  % Build a balanced tree of two-input XOR gates over n inputs and return gates and depth.
  % The Python returns the pair as a tuple.
  level = zeros (1, n);              % depth at which each signal becomes available
  gates = 0;
  while (numel (level) > 1)
    pairs = floor (numel (level) / 2);
    nxt = max (level(1:2:2 * pairs - 1), level(2:2:2 * pairs)) + 1;
    gates = gates + pairs;
    if (mod (numel (level), 2))
      nxt(end + 1) = level(end);
    end
    level = nxt;
  end
  depth = level(1);
end

function r = xor_bits (varargin)
  width = numel (varargin{1});
  acc = 0;
  for i = 1:numel (varargin)
    acc = bitxor (acc, bin2dec (varargin{i}));
  end
  r = dec2bin (acc, width);
end


% Checks ----------------------------------------------------------------------

function main ()
  printf ('One flipped bit\n');
  x = bin2dec ('01100100');          % Python's 0b01100100
  assert (x == 100);
  hi = flip_bit (x, 7);
  lo = flip_bit (x, 0);
  assert (strcmp (dec2bin (hi, 8), '11100100') && hi == 228 && as_signed8 (hi) == -28);
  assert (strcmp (dec2bin (lo, 8), '01100101') && lo == 101);
  show ('01100100 is', int64 (x));
  show ('flip bit 7: 11100100 unsigned, signed', as_tuple (int64 (hi), int64 (as_signed8 (hi))));
  show ('flip bit 0: 01100101', int64 (lo));

  printf ('Where flips come from\n');
  p_dimm = 0.08;
  p16 = 1 - (1 - p_dimm) ^ 16;
  assert (isclose (p16, 0.7366, 'abs_tol', 5e-5));
  show ('P(one of 16 DIMMs errs in a year | 8% each, independent)', sprintf ('%.4f', p16));

  printf ('How often\n');
  rate = 8e9;
  ber = 1e-12;
  per_s = rate * ber;
  assert (isclose (per_s, 0.008) && isclose (per_s * 3600, 28.8));
  assert (isclose (per_s * 86400, 691.2));
  show ('errors per second at 8 Gb/s, BER 1e-12', per_s);
  show ('errors per hour', rnd (per_s * 3600, 6));
  show ('errors per day', rnd (per_s * 86400, 6));
  n_tb = 8e12;                                   % one terabyte, 10^12 bytes
  exact = p_at_least_one (ber, n_tb);
  np_ = n_tb * ber;
  assert (isclose (np_, 8.0));
  assert (isclose (exact, 1 - exp (-8), 'rel_tol', 1e-9));
  assert (isclose (exact, 0.99966, 'abs_tol', 5e-6));
  show ('terabyte: Np', np_);
  show ('terabyte: P(at least one error)', sprintf ('%.6f', exact));
  show ('terabyte: 1 - e^-8', sprintf ('%.6f', 1 - exp (-8)));
  n_pg = 4096 * 8;
  small = p_at_least_one (ber, n_pg);
  assert (n_pg == 32768);
  assert (isclose (small, n_pg * ber, 'rel_tol', 1e-7));
  assert (1 - exp (-n_pg * ber) <= small && small <= n_pg * ber);
  show ('4 KiB page: Np', sprintf ('%.4e', n_pg * ber));
  show ('4 KiB page: P(at least one error)', sprintf ('%.6e', small));
  % The lemma 1 - x <= e^-x and Bernoulli's inequality, spot-checked.
  for t = [0.0, 1e-9, 0.01, 0.5, 1.0, 3.0]
    assert (1 - t <= exp (-t));
  end
  one = Fraction (1, 1);
  ps = {Fraction(1, 10), Fraction(1, 1000)};
  for j = 1:numel (ps)
    p = ps{j};
    for n = 0:59
      % (1 - p)^n >= 1 - n p, compared as Fractions so nothing rounds
      assert (frac_cmp (frac_pow (frac_sub (one, p), n), frac_sub (one, frac_mul (Fraction (n, 1), p))) >= 0);
    end
  end

  printf ('The binary symmetric channel\n');
  % Three bits at p = 0.1: the probability of each error weight, which the
  % chapter's table lists, and the check that the four add to one.
  p = Fraction (1, 10);
  one_pattern = cell (1, 4);             % one_pattern{w + 1} is the entry for weight w
  by_weight = cell (1, 4);
  for w = 0:3
    one_pattern{w + 1} = frac_mul (frac_pow (p, w), frac_pow (frac_sub (one, p), 3 - w));
    by_weight{w + 1} = frac_mul (Fraction (comb (3, w), 1), one_pattern{w + 1});
  end
  op = cellfun (@frac_float, one_pattern);
  bw = cellfun (@frac_float, by_weight);
  assert (isequal (op, [0.729, 0.081, 0.009, 0.001]));
  assert (isequal (bw, [0.729, 0.243, 0.027, 0.001]));
  total = Fraction (0, 1);
  for w = 1:4
    total = frac_add (total, by_weight{w});
  end
  assert (frac_cmp (total, one) == 0);
  show ('one pattern of weight 0..3, n = 3, p = 0.1', op);
  show ('all patterns of weight 0..3', bw);

  printf ('The repetition code R3\n');
  msg = '1011';
  sent = rep_encode (msg);
  assert (strcmp (sent, '111 000 111 111'));
  recv = '101 000 100 110';
  decoded = decode_blocks (recv);
  assert (strcmp (decoded, '1001'));
  show ('encode 1011', sent);
  show ('decode 101 000 100 110', decoded);
  words = cellstr (dec2bin (0:7, 3)).';          % the keys of the Python's table, '000' to '111'
  votes = cellfun (@rep_decode_block, words, 'UniformOutput', false);
  to_zero = words(strcmp (votes, '0'));
  assert (isequal (to_zero, {'000', '001', '010', '100'}));
  show ('words decoded to 0', to_zero);
  for a = 0:1
    for b = 0:1
      for c = 0:1
        assert (majority (a, b, c) == (a + b + c >= 2));
      end
    end
  end
  show ('maj(a,b,c) = ab + bc + ca matches the vote on all 8 rows', true);

  printf ('How much R3 buys\n');
  ps = {Fraction(1, 10), Fraction(1, 100), Fraction(1, 1000), Fraction(1, 3), Fraction(2, 5)};
  for j = 1:numel (ps)
    p = ps{j};
    rhs = frac_sub (frac_mul (Fraction (3, 1), frac_pow (p, 2)), frac_mul (Fraction (2, 1), frac_pow (p, 3)));
    assert (frac_cmp (rep_error (p, 3), rhs) == 0);
  end
  quoted_p = [0.1, 0.01, 0.001];
  quoted_want = [0.028, 2.98e-4, 2.998e-6];
  for j = 1:3
    p = quoted_p(j);
    want = quoted_want(j);
    got = rep_error (p, 3);
    assert (isclose (got, want, 'rel_tol', 1e-9));
    show (sprintf ('P_e(R3) at p = %s', py_str (p)), sprintf ('%.4g  (gain %.1fx)', got, p / got));
  end
  half = Fraction (1, 2);
  ps = {Fraction(1, 10), Fraction(49, 100), Fraction(1, 2), Fraction(51, 100), Fraction(9, 10)};
  for j = 1:numel (ps)
    p = ps{j};
    better = frac_cmp (rep_error (p, 3), p) < 0;
    assert (better == (frac_cmp (p, half) < 0));
  end
  show ('R3 beats uncoded exactly when p < 1/2', true);

  printf ('The price: rate\n');
  r5 = rep_error (0.01, 5);
  assert (isclose (r5, 9.8506e-6, 'rel_tol', 1e-4));
  show ('P_e(R5) at p = 0.01', sprintf ('%.5g', r5));
  ps = {Fraction(1, 100), Fraction(1, 10), Fraction(1, 3)};
  for j = 1:numel (ps)
    p = ps{j};
    rhs = frac_add (frac_sub (frac_mul (Fraction (10, 1), frac_pow (p, 3)), frac_mul (Fraction (15, 1), frac_pow (p, 4))), ...
                    frac_mul (Fraction (6, 1), frac_pow (p, 5)));
    assert (frac_cmp (rep_error (p, 5), rhs) == 0);
  end
  show ('P_e(R5) = 10p^3 - 15p^4 + 6p^5, checked exactly', true);
  for p = [0.001, 0.01, 0.1, 0.3, 0.49]
    for n = 1:2:40
      assert (rep_error (p, n) <= rep_bound (p, n) + 1e-15);
    end
  end
  show ('P_e(R_n) <= (4p(1-p))^(n/2) for odd n < 41, five values of p', true);
  show ('(4p(1-p))^(1/2) at p = 0.1', rnd ((4 * 0.1 * 0.9) ^ 0.5, 12));
  h = binary_entropy (0.1);
  assert (isclose (h, 0.4690, 'abs_tol', 5e-5) && isclose (1 - h, 0.5310, 'abs_tol', 5e-5));
  show ('H(0.1)', sprintf ('%.4f', h));
  show ('capacity 1 - H(0.1)', sprintf ('%.4f', 1 - h));
  for n = [1, 3, 5, 7, 9, 11, 15, 25]
    show (sprintf ('R%d: rate %.4f, P_e at p = 0.1', n, 1 / n), sprintf ('%.3e', rep_error (0.1, n)));
  end

  printf ('Single parity\n');
  ps = {Fraction(1, 100), Fraction(1, 10), Fraction(1, 3)};
  for j = 1:numel (ps)
    p = ps{j};
    for n = 2:16
      assert (frac_cmp (parity_undetected_direct (p, n), parity_undetected_closed (p, n)) == 0);
    end
  end
  pud = parity_undetected_closed (0.01, 8);
  raw = 1 - 0.99 ^ 8;
  assert (isclose (pud, 0.0026368, 'abs_tol', 5e-8));
  assert (isclose (raw, 0.077255, 'abs_tol', 5e-7));
  two = comb (8, 2) * 0.01 ^ 2 * 0.99 ^ 6;
  assert (isclose (two, 0.0026361, 'abs_tol', 5e-8));
  show ('P(undetected) for n = 8, p = 0.01', sprintf ('%.7f', pud));
  show ('of which weight-2 patterns', sprintf ('%.7f', two));
  show ('P(any error) in the same 8 bits', sprintf ('%.6f', raw));
  show ('fraction of erroneous words that slip through', sprintf ('%.4f', pud / raw));
  for n = 1:128
    [g, d] = xor_tree (n);
    if (n > 1)
      want_depth = ceil (log2 (n));
    else
      want_depth = 0;
    end
    assert (g == n - 1 && d == want_depth);
  end
  [g, d] = xor_tree (8);
  show ('XOR tree, 8 inputs: (gates, depth)', as_tuple (int64 (g), int64 (d)));
  [g, d] = xor_tree (64);
  show ('XOR tree, 64 inputs: (gates, depth)', as_tuple (int64 (g), int64 (d)));
  d1 = '1011';
  d2 = '0110';
  d3 = '1100';
  par = xor_bits (d1, d2, d3);
  assert (strcmp (par, '0001'));
  assert (strcmp (xor_bits (d1, d3, par), d2));
  show ('parity block P = D1 xor D2 xor D3', par);
  show ('rebuilt D2 = D1 xor D3 xor P', xor_bits (d1, d3, par));
  word = '1011001';
  pbit = sprintf ('%d', mod (sum (word == '1'), 2));
  assert (strcmp ([word pbit], '10110010'));
  show ('even-parity codeword for 1011001', [word pbit]);

  printf ('Detection, retransmission and correction\n');
  ok = (1 - 1e-4) ^ 1000;
  assert (isclose (ok, 0.904833, 'abs_tol', 5e-7));
  assert (isclose (1 / ok, 1.105176, 'abs_tol', 5e-6));
  show ('P(1000-bit frame arrives clean at p = 1e-4)', sprintf ('%.6f', ok));
  show ('expected transmissions per frame', sprintf ('%.6f', 1 / ok));

  printf ('Worked examples\n');
  % The length at which an error becomes as likely as not, at BER 1e-12.
  ber = 1e-12;
  n_half = log (2) / ber;
  assert (isclose (n_half, 6.931e11, 'rel_tol', 1e-4));
  assert (isclose (n_half / 8 / 1e9, 86.6, 'abs_tol', 0.05));
  assert (isclose (p_at_least_one (ber, n_half), 0.5, 'rel_tol', 1e-9));
  show ('bits until an error is as likely as not, BER 1e-12', sprintf ('%.4e', n_half));
  show ('the same length in gigabytes (10^9 bytes)', sprintf ('%.1f', n_half / 8 / 1e9));
  % R3 on a channel with p = 0.05.
  p = Fraction (1, 20);
  pe = rep_error (p, 3);
  assert (frac_cmp (pe, Fraction (29, 4000)) == 0 && frac_float (pe) == 0.00725);
  assert (isclose (frac_float (frac_div (p, pe)), 6.9, 'abs_tol', 0.05));
  assert (isclose (rep_bound (0.05, 3), 0.0828, 'abs_tol', 5e-5));
  show ('P_e(R3) at p = 0.05, and the gain', as_tuple (frac_float (pe), rnd (frac_float (frac_div (p, pe)), 1)));
  show ('the bound (4p(1-p))^(3/2) at p = 0.05', rnd (rep_bound (0.05, 3), 4));
  % A 64-bit word with one parity bit, so n = 65, at p = 1e-6.
  p = Fraction (1, 1e6);
  any65 = frac_sub (one, frac_pow (frac_sub (one, p), 65));
  ud65 = parity_undetected_closed (p, 65);
  c65 = frac_mul (Fraction (comb (65, 2), 1), frac_pow (p, 2));
  assert (isclose (frac_float (any65), 6.4998e-5, 'rel_tol', 1e-4));
  assert (isclose (frac_float (ud65), 2.0799e-9, 'rel_tol', 1e-4));
  assert (isclose (frac_float (c65), 2.08e-9, 'rel_tol', 1e-12));
  assert (isclose (frac_float (frac_div (ud65, any65)), 3.2e-5, 'rel_tol', 1e-3));
  show ('n = 65, p = 1e-6: P(any error)', sprintf ('%.4e', frac_float (any65)));
  show ('n = 65, p = 1e-6: P(undetected)', sprintf ('%.4e', frac_float (ud65)));
  show ('n = 65, p = 1e-6: C(65,2) p^2', sprintf ('%.4e', frac_float (c65)));
  % A second RAID rebuild, this time losing the third disk.
  d1 = '1110';
  d2 = '0101';
  d3 = '1000';
  par = xor_bits (d1, d2, d3);
  assert (strcmp (par, '0011') && strcmp (xor_bits (d1, d2, par), d3));
  show ('P for 1110, 0101, 1000; rebuilt D3', as_tuple (par, xor_bits (d1, d2, par)));
  % Decoding a received R3 sequence.
  recv = '110 011 001';
  assert (strcmp (decode_blocks (recv), '110'));
  show ('decode 110 011 001', '110');
  printf ('All chapter 1 numbers reproduced.\n');
end


main ();
