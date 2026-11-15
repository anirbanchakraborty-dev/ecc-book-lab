% Reproduce every number in chapter 8 of Error Correction Coding,
% "Distance, Detection and Correction", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch08_distance_detection_and_correction.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch08_distance_detection_and_correction.py and
% prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees by more than half a unit in the last digit the chapter prints.
% Small codes are checked exhaustively: every codeword against every error
% pattern the chapter names, every combination of errors and erasures inside
% the bound, all 2^16 rectangular parity arrays, and every reflected Gray
% code up to 8 bits. The check-digit schemes are tested at every position
% against every digit pair, on numbers built to be valid, and the Luhn and
% ISBN results are confirmed again on every valid number of five digits. A
% run that reaches the final line has reproduced every number in the chapter.
%
% Two things work differently from the Python, by necessity. A word is held
% in a double, which stores every integer below 2^53 exactly, far more than
% the 25-bit words used here, and the bit operations are bitxor, bitand,
% bitor and bitshift. The single-parity section works with exact fractions
% whose numerators and denominators grow past 2^1000, so it carries them in
% a small exact-integer arithmetic of its own (the bi_ functions below) and
% rounds each one to a double exactly as Python's float(Fraction) does.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ----------------------------------------------------------------------

function check (label, computed, stated)
  % The computed value arrives as a double. Where the Python passes a
  % Fraction, the call site has already rounded it the way float() does.
  text = strrep (stated, ',', '');
  percent = ~isempty (stated) && stated(end) == '%';
  if (percent)
    text = text(1:end - 1);
    computed = 100 * computed;
  end
  value = str2double (text);
  printf ('  %-62s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % Python compares the two values with !=. Every value compared here prints
  % differently whenever it differs, so comparing the printed text is the
  % same test.
  [c_str, c_repr] = py_text (computed);
  [s_str, s_repr] = py_text (stated);
  printf ('  %-62s %-16s %s\n', label, c_str, s_str);
  if (~strcmp (c_str, s_str))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, c_repr, s_repr);
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


% Printing Python values ------------------------------------------------------
%
% same() prints with Python's str() and reports a mismatch with repr(). The
% shared py_str covers numbers, text, lists and tuples. A dict, a Fraction
% and None are spelled out here as a struct holding both texts.

function x = shown (str_text, repr_text)
  x = struct ('str', str_text, 'repr', repr_text);
end

function [s, r] = py_text (x)
  if (isstruct (x) && isfield (x, 'str'))
    s = x.str;
    r = x.repr;
  else
    x = as_int (x);
    s = py_str (x);
    r = py_str (x, true);
  end
end

function x = as_int (x)
  % Every number that same() shows in this chapter is a Python int, while
  % Octave computes in doubles. The whole doubles become int64 here, so that
  % py_str prints 3 and not 3.0.
  if (isstruct (x) && isfield (x, 'tuple'))
    x.tuple = cellfun (@as_int, x.tuple, 'UniformOutput', false);
  elseif (iscell (x))
    x = cellfun (@as_int, x, 'UniformOutput', false);
  elseif (isnumeric (x) && ~isinteger (x))
    if (any (x(:) ~= fix (x(:))))
      error ('as_int: %s is not a whole number', py_str (x));
    end
    x = int64 (x);
  end
end

function t = tuple (varargin)
  t = struct ('tuple', {varargin});
end

function list = tuple_list (pairs)
  % The rows of a matrix as a list of tuples, [(0, 9), (9, 0)].
  list = cell (1, rows (pairs));
  for r = 1:rows (pairs)
    list{r} = struct ('tuple', {num2cell(pairs(r, :))});
  end
end

function list = strings_of (values, n)
  % [show(v, n) for v in values], as a cell row of text.
  list = arrayfun (@(v) show (v, n), values(:).', 'UniformOutput', false);
end

function x = py_none ()
  x = shown ('None', 'None');
end

function x = none_or (value)
  % The decoders return NaN where the Python returns None.
  if (isnan (value))
    x = py_none ();
  else
    x = value;
  end
end

function x = dict_text (keys, values)
  parts = arrayfun (@(k, v) sprintf ('%d: %d', k, v), keys, values, 'UniformOutput', false);
  text = ['{' strjoin(parts, ', ') '}'];
  x = shown (text, text);
end

function x = fraction_shown (num, den)
  % str() and repr() of Fraction(num, den), reduced to lowest terms as the
  % Fraction constructor reduces it. A whole Fraction prints without /1.
  % One argument is read as the pair [num, den].
  if (nargin == 1)
    den = num(2);
    num = num(1);
  end
  g = gcd (num, den);
  num = num / g;
  den = den / g;
  if (den == 1)
    text = sprintf ('%d', num);
  else
    text = sprintf ('%d/%d', num, den);
  end
  x = shown (text, sprintf ('Fraction(%d, %d)', num, den));
end

function s = py_sum (values)
  % Python's built-in sum() of floats. Since Python 3.12 it carries a
  % running correction (Neumaier's compensated summation), which Octave's
  % sum() does not, and the two can differ in the last bit.
  s = 0;
  c = 0;
  for x = values(:).'
    t = s + x;
    if (abs (s) >= abs (x))
      c = c + ((s - t) + x);
    else
      c = c + ((x - t) + s);
    end
    s = t;
  end
  if (c ~= 0 && isfinite (c))
    s = s + c;
  end
end

function c = combos (v, k)
  % itertools.combinations(v, k), one combination per row, in the same
  % order. nchoosek alone reads a one-element v as a count rather than a
  % list, so the edge cases are settled first. With k = 0 there is one
  % empty combination, and with k larger than v there is none.
  v = v(:).';
  if (k == 0)
    c = zeros (1, 0);
  elseif (k > numel (v))
    c = zeros (0, k);
  elseif (k == numel (v))
    c = v;
  else
    c = nchoosek (v, k);
  end
end


% Words, distance and codes ----------------------------------------------------

function w = weight (x)
  % bin(x).count("1") for every element of x at once.
  w = zeros (size (x));
  while (any (x(:)))
    w = w + mod (x, 2);
    x = floor (x / 2);
  end
end

function d = distance (x, y)
  % bsxfun lets a single word be measured against a whole list of them.
  d = weight (bsxfun (@bitxor, x, y));
end

function x = word (s)
  x = bin2dec (s);
end

function s = show (x, n)
  s = dec2bin (x, n);
end

function code = code_of (varargin)
  code = cellfun (@word, varargin);
end

function d = min_distance (code)
  % Each row of code is one code, so a whole stack of codes can be measured
  % at once. A single code is a single row.
  pairs = combos (1:columns (code), 2);
  d = min (distance (code(:, pairs(:, 1)), code(:, pairs(:, 2))), [], 2);
end

function [keys, counts] = distance_distribution (code)
  % Python returns a dict sorted by key. Here the keys come back sorted
  % beside their counts.
  pairs = combos (code, 2);
  [keys, ~, j] = unique (distance (pairs(:, 1), pairs(:, 2)));
  keys = keys(:).';
  counts = accumarray (j(:), 1).';
end

function ok = closed_under_xor (code)
  sums = bsxfun (@bitxor, code(:), code(:).');
  ok = all (ismember (sums(:), code));
end

function v = ball_size (n, t)
  v = 0;
  for i = 0:t
    v = v + nchoosek (n, i);
  end
end

function members = ball (center, n, t)
  % Python returns a set. Here the members come back as a sorted row.
  y = 0:2^n - 1;
  members = y(distance (center, y) <= t);
end

function edges = cube_edges (n)
  % One row per edge (x, x ^ (1 << i)), listed with x before i as in Python.
  x = repelem ((0:2^n - 1).', n);
  i = repmat ((0:n - 1).', 2^n, 1);
  y = bitxor (x, 2.^i);
  keep = x < y;
  edges = [x(keep), y(keep)];
end

function winners = nearest (code, r)
  d = distance (code, r);
  winners = code(d == min (d));
end

function c = decode_nearest (code, r)
  % Nearest-codeword decoding. Returns NaN, standing in for None, on a tie.
  winners = nearest (code, r);
  if (numel (winners) == 1)
    c = winners;
  else
    c = NaN;
  end
end

function c = decode_bounded (code, r, t)
  % Correct when exactly one codeword lies within t, otherwise flag (NaN).
  inside = code(distance (code, r) <= t);
  if (numel (inside) == 1)
    c = inside;
  else
    c = NaN;
  end
end

function c = decode_erasures (code, n, received)
  % Nearest codeword on the positions that are not erased ('?'). NaN on a tie.
  kept = find (received ~= '?');
  bits = show (code(:), n);
  gaps = sum (bits(:, kept) ~= received(kept), 2);
  winners = code(gaps == min (gaps));
  if (numel (winners) == 1)
    c = winners;
  else
    c = NaN;
  end
end

function bits = apply_errors_erasures (c, n, flips, erased)
  % flips and erased hold 0-based positions, as in Python, so each one is
  % one higher as an index into the text.
  bits = show (c, n);
  for i = flips
    if (bits(i + 1) == '0')
      bits(i + 1) = '1';
    else
      bits(i + 1) = '0';
    end
  end
  for i = erased
    bits(i + 1) = '?';
  end
end

function code = random_code (seed, words, n)
  rng = PyRandom (seed);
  code = rng.getrandbits (n, words).';
end

function errs = codeword_errors (code, n, p)
  % Exact error probability of each codeword under ML decoding on BSC(p),
  % ties split evenly. The received words are taken in Python's order, so
  % each codeword's total is built from the same additions.
  pw = zeros (1, n + 1);
  for k = 0:n
    pw(k + 1) = p^k * (1 - p)^(n - k);
  end
  correct = zeros (1, numel (code));
  for y = 0:2^n - 1
    dists = distance (code, y);
    best = min (dists);
    winners = find (dists == best);
    correct(winners) = correct(winners) + pw(best + 1) / numel (winners);
  end
  errs = 1 - correct;
end

function v = ml_block_error (code, n, p)
  errs = codeword_errors (code, n, p);
  v = py_sum (errs) / numel (errs);
end

function v = undetected_probability (code, n, p)
  % For a code closed under XOR: the chance the error pattern is a nonzero codeword.
  nonzero = code(code ~= 0);
  terms = zeros (size (nonzero));
  for k = 1:numel (nonzero)
    w = weight (nonzero(k));
    terms(k) = p^w * (1 - p)^(n - w);
  end
  v = py_sum (terms);
end


% Exact integers and fractions -------------------------------------------------
%
% Python's int has no size limit, and the single-parity checks lean on it:
% the coverage of a 64-bit word at p = 1/1000 is a fraction whose terms
% pass 2^1000. A double holds integers exactly only below 2^53. These
% functions carry a nonnegative integer as a row of base-2^24 digits
% ("limbs"), least significant first. A product of two limbs stays below
% 2^48, inside the exact range of a double, and every sum of such products
% is carried back into range before it can pass 2^53.

function x = bi (v)
  % A small nonnegative integer, below 2^53, as limbs.
  x = zeros (1, 0);
  while (v > 0)
    x(end + 1) = mod (v, 2^24);
    v = floor (v / 2^24);
  end
  if (isempty (x))
    x = 0;
  end
end

function x = bi_norm (x)
  % Carry every limb into 0 .. 2^24 - 1 and drop the leading zero limbs.
  % floor() carries a negative limb as a borrow from the next one.
  x = [x, 0];
  carry = floor (x / 2^24);
  while (any (carry))
    if (carry(end) < 0)
      error ('bi_norm: an exact integer went negative');
    end
    x = x - carry * 2^24;
    x(2:end) = x(2:end) + carry(1:end - 1);
    if (carry(end) > 0)
      x(end + 1) = carry(end);
    end
    carry = floor (x / 2^24);
  end
  top = find (x, 1, 'last');
  if (isempty (top))
    x = 0;
  else
    x = x(1:top);
  end
end

function z = bi_add (x, y)
  n = max (numel (x), numel (y));
  z = bi_norm ([x, zeros(1, n - numel (x))] + [y, zeros(1, n - numel (y))]);
end

function z = bi_sub (x, y)
  % x - y, which must not be negative.
  n = max (numel (x), numel (y));
  z = bi_norm ([x, zeros(1, n - numel (x))] - [y, zeros(1, n - numel (y))]);
end

function z = bi_mul (x, y)
  % Schoolbook multiplication, one limb of y at a time, carrying after each
  % row so that no limb passes 2^53.
  n = numel (x) + numel (y);
  z = zeros (1, n);
  for j = 1:numel (y)
    z(j:j + numel (x) - 1) = z(j:j + numel (x) - 1) + x * y(j);
    z = bi_norm (z);
    z(end + 1:n) = 0;
  end
  z = bi_norm (z);
end

function c = bi_cmp (x, y)
  % -1, 0 or 1 as x is less than, equal to or greater than y.
  if (numel (x) ~= numel (y))
    c = sign (numel (x) - numel (y));
    return;
  end
  k = find (x ~= y, 1, 'last');
  if (isempty (k))
    c = 0;
  else
    c = sign (x(k) - y(k));
  end
end

function n = bi_bitlen (x)
  % int.bit_length(). The two-output log2 splits a limb exactly into a
  % fraction in [0.5, 1) and a power of two.
  if (x(end) == 0)
    n = 0;
  else
    [~, e] = log2 (x(end));
    n = 24 * (numel (x) - 1) + e;
  end
end

function z = bi_shift (x, s)
  % x * 2^s for s >= 0: whole limbs first, then the remaining bits.
  z = bi_norm ([zeros(1, floor (s / 24)), x * 2^mod(s, 24)]);
end

function v = bi_double (x)
  % The exact value as a double, allowed only below 2^53.
  if (bi_bitlen (x) > 53)
    error ('bi_double: the integer does not fit a double exactly');
  end
  v = sum (x .* 2.^(24 * (0:numel (x) - 1)));
end

function v = bi_ratio (num, den)
  % num / den rounded to the nearest double, ties to even, which is what
  % Python's int / int gives and so what float(Fraction) gives. The quotient
  % is taken to 55 or 56 bits by long division, one bit at a time, and a
  % nonzero remainder is folded into its lowest bit. Converting that integer
  % to a double then drops at least two bits, and the hardware rounding of
  % the conversion is the correct rounding of the exact ratio. The quotient
  % is held in uint64, which carries all 56 bits where a double would not.
  if (bi_cmp (num, 0) == 0)
    v = 0;
    return;
  end
  s = 55 - (bi_bitlen (num) - bi_bitlen (den));
  if (s >= 0)
    rest = bi_shift (num, s);
    divisor = den;
  else
    rest = num;
    divisor = bi_shift (den, -s);
  end
  q = uint64 (0);
  for k = 55:-1:0
    step = bi_shift (divisor, k);
    if (bi_cmp (rest, step) >= 0)
      rest = bi_sub (rest, step);
      q = q + bitshift (uint64 (1), k);
    end
  end
  if (any (rest))
    q = bitor (q, uint64 (1));
  end
  v = pow2 (double (q), -s);
end

% A fraction is a struct of two such integers, num and den. It is kept
% unreduced, since equality is tested by cross-multiplying and rounding to a
% double does not depend on the reduction.

function f = frac (a, b)
  f = struct ('num', bi (a), 'den', bi (b));
end

function f = frac_add (x, y)
  if (bi_cmp (x.den, y.den) == 0)
    f = struct ('num', bi_add (x.num, y.num), 'den', x.den);
  else
    f = struct ('num', bi_add (bi_mul (x.num, y.den), bi_mul (y.num, x.den)), 'den', bi_mul (x.den, y.den));
  end
end

function f = frac_sub (x, y)
  % x - y, which must not be negative.
  if (bi_cmp (x.den, y.den) == 0)
    f = struct ('num', bi_sub (x.num, y.num), 'den', x.den);
  else
    f = struct ('num', bi_sub (bi_mul (x.num, y.den), bi_mul (y.num, x.den)), 'den', bi_mul (x.den, y.den));
  end
end

function f = frac_mul (x, y)
  f = struct ('num', bi_mul (x.num, y.num), 'den', bi_mul (x.den, y.den));
end

function f = frac_div (x, y)
  f = struct ('num', bi_mul (x.num, y.den), 'den', bi_mul (x.den, y.num));
end

function f = frac_pow (x, n)
  f = frac (1, 1);
  for k = 1:n
    f = frac_mul (f, x);
  end
end

function t = frac_eq (x, y)
  t = bi_cmp (bi_mul (x.num, y.den), bi_mul (y.num, x.den)) == 0;
end

function v = frac_float (x)
  v = bi_ratio (x.num, x.den);
end

function x = frac_text (f)
  % str() of a fraction small enough to reduce in doubles.
  x = fraction_shown (bi_double (f.num), bi_double (f.den));
end


% Single parity -------------------------------------------------------------------
%
% Python's functions take p as a float or a Fraction alike. Here p is a
% fraction struct, except in any_error, which also accepts a double because
% the chapter calls it with one.

function f = parity_alarm (n, p)
  one = frac (1, 1);
  two = frac (2, 1);
  f = frac_div (frac_sub (one, frac_pow (frac_sub (one, frac_mul (two, p)), n)), two);
end

function f = any_error (n, p)
  if (isstruct (p))
    f = frac_sub (frac (1, 1), frac_pow (frac_sub (frac (1, 1), p), n));
  else
    f = 1 - (1 - p)^n;
  end
end

function f = parity_coverage (n, p)
  f = frac_div (parity_alarm (n, p), any_error (n, p));
end

function f = parity_alarm_by_enumeration (n, p)
  % Every pattern of odd weight w adds p^w (1-p)^(n-w), so the 2^n patterns
  % are counted by weight and each weight's term is added that many times.
  % The exact sum is the same as adding them one at a time.
  w = weight (0:2^n - 1);
  count = accumarray (w(:) + 1, 1, [n + 1, 1]);
  q = frac_sub (frac (1, 1), p);
  f = frac (0, 1);
  for k = 1:2:n
    term = frac_mul (frac_pow (p, k), frac_pow (q, n - k));
    f = frac_add (f, frac_mul (frac (count(k + 1), 1), term));
  end
end


% Check digits ---------------------------------------------------------------------
%
% A number is a row of digits. The sums and validity tests also take a
% matrix with one number per row, which lets make_valid try all 100 digit
% pairs at once.

function ds = digits (s)
  s = s(s ~= '-' & s ~= ' ');
  ds = double (s) - double ('0');
  ds(s == 'X' | s == 'x') = 10;
end

function t = isbn10_sum (ds)
  t = ds * (10 - (0:columns (ds) - 1)).';
end

function c = isbn10_check (payload)
  % mod, unlike rem, takes the sign of the divisor, as Python's % does.
  c = mod (-isbn10_sum (payload), 11);
end

function t = isbn10_valid (ds)
  t = columns (ds) == 10 & mod (isbn10_sum (ds), 11) == 0;
end

function w = isbn13_weights (len)
  % Weight 1 at the even 0-based positions, 3 at the odd ones.
  if (nargin < 1)
    len = 13;
  end
  w = ones (1, len);
  w(2:2:end) = 3;
end

function t = isbn13_sum (ds)
  t = ds * isbn13_weights (columns (ds)).';
end

function c = isbn13_check (payload)
  c = mod (-isbn13_sum (payload), 10);
end

function t = isbn13_valid (ds)
  t = columns (ds) == 13 & mod (isbn13_sum (ds), 10) == 0;
end

function y = luhn_double (x)
  y = 2 * x;
  y(y > 9) = y(y > 9) - 9;
end

function t = luhn_terms (ds)
  % The value each digit adds, counting positions from the right, where
  % position 1 is the check digit.
  n = columns (ds);
  doubled = mod (n - (0:n - 1), 2) == 0;
  t = ds;
  t(:, doubled) = luhn_double (ds(:, doubled));
end

function s = luhn_sum (ds)
  s = sum (luhn_terms (ds), 2);
end

function c = luhn_check (payload)
  for c = 0:9
    if (mod (luhn_sum ([payload, c]), 10) == 0)
      return;
    end
  end
  error ('no Luhn check digit');
end

function t = luhn_valid (ds)
  t = mod (luhn_sum (ds), 10) == 0;
end

function ds = make_valid (valid, base, fixed_pos, fixed_val)
  % A number equal to base except at the fixed positions and two free
  % positions, that passes valid(). Python passes the fixed digits as a dict
  % {position: digit}. Here they come as two rows, the 0-based positions
  % and their digits. For each pair of free positions the 100 digit pairs of
  % product(range(10), repeat=2) are tried at once, a the slower digit, and
  % the first that passes is the one Python's loop would return.
  n = numel (base);
  free = setdiff (0:n - 1, fixed_pos);
  start = base;
  start(fixed_pos + 1) = fixed_val;
  ab = (0:99).';
  trials = repmat (start, 100, 1);
  for x = 1:numel (free) - 1
    for y = x + 1:numel (free)
      trial = trials;
      trial(:, free(x) + 1) = floor (ab / 10);
      trial(:, free(y) + 1) = mod (ab, 10);
      hit = find (valid (trial), 1);
      if (~isempty (hit))
        ds = trial(hit, :);
        return;
      end
    end
  end
  error ('no valid number');
end

function t = swap_escapes (valid, base, i, j, a, b)
  ds = make_valid (valid, base, [i, j], [a, b]);
  ds([i, j] + 1) = ds([j, i] + 1);
  t = valid (ds);
end

function t = substitution_escapes (valid, base, i, a, b)
  ds = make_valid (valid, base, i, a);
  ds(i + 1) = b;
  t = valid (ds);
end

function rates = scheme_rates (valid, base)
  % Fractions of single substitutions, adjacent swaps and jump swaps (two
  % apart) that are caught. Each fraction is returned as [caught, cases].
  n = numel (base);
  ab = (0:99).';
  pairs = [floor(ab / 10), mod(ab, 10)];
  pairs = pairs(pairs(:, 1) ~= pairs(:, 2), :);
  subs = false (0, 1);
  adj = false (0, 1);
  jump = false (0, 1);
  for i = 0:n - 1
    for r = 1:rows (pairs)
      subs(end + 1) = substitution_escapes (valid, base, i, pairs(r, 1), pairs(r, 2));
    end
  end
  for i = 0:n - 2
    for r = 1:rows (pairs)
      adj(end + 1) = swap_escapes (valid, base, i, i + 1, pairs(r, 1), pairs(r, 2));
    end
  end
  for i = 0:n - 3
    for r = 1:rows (pairs)
      jump(end + 1) = swap_escapes (valid, base, i, i + 2, pairs(r, 1), pairs(r, 2));
    end
  end
  rate = @(cases) [sum(~cases), numel(cases)];
  rates = struct ('substitution', rate (subs), 'adjacent', rate (adj), 'jump', rate (jump));
end

function pairs = escaping_adjacent_pairs (valid, base)
  % Python returns a set of (a, b). Here it is a matrix of distinct rows.
  pairs = zeros (0, 2);
  for i = 0:numel (base) - 2
    for a = 0:9
      for b = 0:9
        if (a ~= b && swap_escapes (valid, base, i, i + 1, a, b))
          pairs(end + 1, :) = [a, b];
        end
      end
    end
  end
  pairs = unique (pairs, 'rows');
end


% Rectangular parity ----------------------------------------------------------------
%
% An array is a matrix of bits. Python's array[i][j] is array(i + 1, j + 1).

function out = rect_encode (rows_in)
  % Append a parity bit to each row, then a parity row, so every row and
  % column of the result is even. A stack of arrays along the third
  % dimension is encoded page by page.
  with_rows = [rows_in, mod(sum (rows_in, 2), 2)];
  out = [with_rows; mod(sum (with_rows, 1), 2)];
end

function [bad_rows, bad_cols] = rect_failures (array)
  % The 0-based rows and columns whose parity fails, as rows of numbers.
  bad_rows = find (mod (sum (array, 2), 2)).' - 1;
  bad_cols = find (mod (sum (array, 1), 2)) - 1;
end

function x = failures_shown (array)
  % rect_failures as Python prints it, a tuple of two lists.
  [bad_rows, bad_cols] = rect_failures (array);
  x = tuple (num2cell (bad_rows), num2cell (bad_cols));
end

function [fixed, verdict] = rect_decode (array)
  % Correct one error, flag anything else that breaks a parity. Returns the
  % array, or [] where Python returns None, and the verdict.
  [bad_rows, bad_cols] = rect_failures (array);
  if (isempty (bad_rows) && isempty (bad_cols))
    fixed = array;
    verdict = 'clean';
  elseif (numel (bad_rows) == 1 && numel (bad_cols) == 1)
    fixed = array;
    fixed(bad_rows + 1, bad_cols + 1) = bitxor (fixed(bad_rows + 1, bad_cols + 1), 1);
    verdict = 'corrected';
  else
    fixed = [];
    verdict = 'flagged';
  end
end

function value = flatten (array)
  % The bits read row by row, the first one most significant. A stack of
  % arrays gives one value per page.
  [m, n, pages] = size (array);
  bits = reshape (permute (array, [2, 1, 3]), m * n, pages);
  value = 2.^(m * n - 1:-1:0) * bits;
end


% Gray codes ------------------------------------------------------------------------

function g = gray (b)
  g = bitxor (b, bitshift (b, -1));
end

function b = gray_inverse (g)
  b = zeros (size (g));
  while (any (g(:)))
    b = bitxor (b, g);
    g = bitshift (g, -1);
  end
end

function list = reflected (n)
  if (n == 1)
    list = {'0', '1'};
    return;
  end
  prev = reflected (n - 1);
  list = [cellfun(@(s) ['0' s], prev, 'UniformOutput', false), ...
          cellfun(@(s) ['1' s], fliplr (prev), 'UniformOutput', false)];
end


% Checks ----------------------------------------------------------------------------
%
% Python keeps these six codes as module constants. Octave's local functions
% cannot see a script's variables, so each one is a function.

function c = C5 ()
  c = code_of ('00000', '01101', '10110', '11011');
end

function c = C6 ()
  c = code_of ('000000', '011011', '101101', '110110');
end

function c = E6 ()
  c = code_of ('000000', '001111', '110011', '111100');
end

function c = R3 ()
  c = code_of ('000', '111');
end

function c = R4 ()
  c = code_of ('0000', '1111');
end

function c = R5 ()
  c = code_of ('00000', '11111');
end


function check_distance_and_weight ()
  section ('Distance and Weight');
  x = word ('10110');
  y = word ('11011');
  same ('10110 XOR 11011', show (bitxor (x, y), 5), '01101');
  same ('weight of 01101', weight (bitxor (x, y)), 3);
  same ('d(10110, 11011)', distance (x, y), 3);
  same ('weight of 10110', weight (x), 3);
  same ('weight of 11011', weight (y), 4);
  same ('10110 AND 11011', show (bitand (x, y), 5), '10010');
  same ('w(x) + w(y) - 2 w(x AND y)', weight (x) + weight (y) - 2 * weight (bitand (x, y)), 3);
  [a, b] = ndgrid (0:31, 0:31);
  a = a(:);
  b = b(:);
  require ('w(x XOR y) = w(x) + w(y) - 2 w(x AND y) for all 5-bit pairs', ...
           all (weight (bitxor (a, b)) == weight (a) + weight (b) - 2 * weight (bitand (a, b))));
  [a, b, c] = ndgrid (0:63, 0:63, 0:63);
  a = a(:);
  b = b(:);
  c = c(:);
  require ('triangle inequality holds for all 6-bit triples', ...
           all (distance (a, c) <= distance (a, b) + distance (b, c)));
  a = word ('00000');
  b = word ('01101');
  c = word ('11011');
  same ('d(00000, 11011)', distance (a, c), 4);
  same ('d(00000, 01101) + d(01101, 11011)', distance (a, b) + distance (b, c), 6);
  % Python keys its dict by a tuple of the two words. A Map keys by text.
  table = containers.Map ();
  pairs = combos (C5, 2);
  for k = 1:rows (pairs)
    table([show(pairs(k, 1), 5) '-' show(pairs(k, 2), 5)]) = distance (pairs(k, 1), pairs(k, 2));
  end
  same ('C5 distance 00000-01101', table('00000-01101'), 3);
  same ('C5 distance 00000-10110', table('00000-10110'), 3);
  same ('C5 distance 00000-11011', table('00000-11011'), 4);
  same ('C5 distance 01101-10110', table('01101-10110'), 4);
  same ('C5 distance 01101-11011', table('01101-11011'), 3);
  same ('C5 distance 10110-11011', table('10110-11011'), 3);
  same ('pairs in C5', table.Count, 6);
  same ('minimum distance of C5', min_distance (C5), 3);
  require ('C5 is closed under XOR', closed_under_xor (C5));
  same ('01101 XOR 10110', show (bitxor (word ('01101'), word ('10110')), 5), '11011');
  code = C5 ();
  same ('smallest nonzero weight in C5', min (weight (code(code ~= 0))), 3);
end

function check_minimum_distance ()
  section ('Minimum Distance');
  same ('d of all 2-bit words', min_distance (0:3), 1);
  v = 0:7;
  same ('d of 3-bit even-weight words', min_distance (v(mod (weight (v), 2) == 0)), 2);
  same ('d of R3', min_distance (R3), 3);
  code = random_code (1, 16, 12);
  same ('random code, distinct codewords', numel (unique (code)), 16);
  same ('random code, first codeword', show (code(1), 12), '001000100110');
  same ('random code, last codeword', show (code(end), 12), '001101011011');
  same ('random code, pairs', nchoosek (16, 2), 120);
  [keys, counts] = distance_distribution (code);
  same ('random code, distance distribution', dict_text (keys, counts), ...
        dict_text ([2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 6, 15, 23, 17, 29, 22, 5, 1]));
  same ('random code, minimum distance', min_distance (code), 2);
  % Sorting each pair and then the list, as Python sorts the text. Binary
  % words of one width sort as text the way they sort as numbers.
  pairs = combos (code, 2);
  close = sortrows (sort (pairs(distance (pairs(:, 1), pairs(:, 2)) == 2, :), 2));
  close_text = cell (1, rows (close));
  for k = 1:rows (close)
    close_text{k} = tuple (show (close(k, 1), 12), show (close(k, 2), 12));
  end
  same ('random code, pairs at distance 2', close_text, ...
        {tuple('110000101100', '110000111000'), tuple('110010011110', '110011010110')});
  % A ratio of two exact integers rounds once, as float(Fraction) does.
  mean_distance = sum (keys .* counts) / 120;
  check ('random code, mean pairwise distance', mean_distance, '6.11');
  check ('random pair, expected distance n/2', 12 / 2, '6');
  check ('P(two random 12-bit words within distance 2)', (1 + 12 + 66) / 4096, '0.0193');
  same ('words within distance 2 of a 12-bit word, 1 + 12 + 66', ball_size (12, 2), 79);
  check ('expected pairs within distance 2 among 120', 120 * 79 / 4096, '2.31');
  expected = 120 * arrayfun (@(k) nchoosek (12, k), 0:12) / 4096;
  check ('expected pairs at distance 6', expected(6 + 1), '27.1');
  check ('expected pairs at distance 7', expected(7 + 1), '23.2');
  check ('ML block error, BSC(0.1)', ml_block_error (code, 12, 0.1), '0.10539');
  errs = codeword_errors (code, 12, 0.01);
  block = py_sum (errs) / 16;
  check ('ML block error, BSC(0.01)', block, '0.0032306');
  in_pairs = unique (close(:));
  share = py_sum (errs(ismember (code, in_pairs))) / 16;
  check ('share of the four close-pair codewords at p = 0.01', share / block, '0.82');
  check ('first-order estimate p/4 at p = 0.01', 0.01 / 4, '0.0025');
  check ('uncoded 4 bits, BSC(0.1)', 1 - 0.9^4, '0.3439');
end

function check_cube ()
  section ('The Geometry of the Binary Cube');
  same ('vertices of the 3-cube', bitshift (1, 3), 8);
  same ('edges of the 3-cube', rows (cube_edges (3)), 12);
  same ('edges of the 4-cube', rows (cube_edges (4)), 32);
  ok = true;
  for n = 1:8
    ok = ok && rows (cube_edges (n)) == n * 2^(n - 1);
  end
  require ('edges of the n-cube = n 2^(n-1), n = 1..8', ok);
  require ('graph distance on the 4-cube equals Hamming distance', cube_graph_distance_matches (4));
  b0 = ball (word ('000'), 3, 1);
  b1 = ball (word ('111'), 3, 1);
  same ('B(000, 1)', sort (strings_of (b0, 3)), {'000', '001', '010', '100'});
  same ('B(111, 1)', sort (strings_of (b1, 3)), {'011', '101', '110', '111'});
  require ('the two radius-1 balls of R3 are disjoint and cover the 3-cube', ...
           isempty (intersect (b0, b1)) && numel (union (b0, b1)) == 8);
  same ('V(3, 1)', ball_size (3, 1), 4);
  same ('V(4, 1)', ball_size (4, 1), 5);
  c0 = ball (word ('0000'), 4, 1);
  c1 = ball (word ('1111'), 4, 1);
  same ('words in the two radius-1 balls of R4', numel (union (c0, c1)), 10);
  v = 0:15;
  outside = sort (strings_of (v(~ismember (v, union (c0, c1))), 4));
  same ('words of the 4-cube in neither ball', outside, {'0011', '0101', '0110', '1001', '1010', '1100'});
  require ('each of them is at distance 2 from both codewords', ...
           all (cellfun (@(s) distance (word (s), 0) == 2 && distance (word (s), 15) == 2, outside)));
  same ('4-cube words by weight', arrayfun (@(k) sum (weight (v) == k), 0:4), [1, 4, 6, 4, 1]);
  same ('V(5, 1)', ball_size (5, 1), 6);
  same ('C5: M V(5, 1)', 4 * ball_size (5, 1), 24);
  v = 0:31;
  stray = v(arrayfun (@(x) min (distance (x, C5)) >= 2, v));
  same ('5-bit words in no radius-1 ball of C5', numel (stray), 8);
  require ('each lies at distance 2 from exactly two codewords of C5', ...
           all (arrayfun (@(x) isequal (sort (distance (x, C5)), [2, 2, 3, 3]), stray)));
  ok = true;
  for n = 1:6
    for t = 0:n
      ok = ok && numel (ball (0, n, t)) == ball_size (n, t);
    end
  end
  require ('V(n, t) counts the ball for n <= 6, t <= n', ok);
  same ('V(6, 1)', ball_size (6, 1), 7);
  same ('16 V(6, 1)', 16 * ball_size (6, 1), 112);
  same ('2^6', 2^6, 64);
  same ('V(7, 1)', ball_size (7, 1), 8);
  same ('16 V(7, 1)', 16 * ball_size (7, 1), 128);
  hamming = hamming_7_4 ();
  same ('Hamming (7,4) codewords', numel (unique (hamming)), 16);
  same ('Hamming (7,4) minimum distance', min_distance (hamming), 3);
  balls = arrayfun (@(c) ball (c, 7, 1), hamming, 'UniformOutput', false);
  require ('Hamming (7,4) radius-1 balls tile the 7-cube', numel (unique ([balls{:}])) == 128);
end

function ok = cube_graph_distance_matches (n)
  % Breadth-first search from every vertex. Row v + 1 of adjacency holds the
  % neighbors of v, and seen(v + 1) is -1 until v is reached.
  adjacency = bsxfun (@bitxor, (0:2^n - 1).', 2.^(0:n - 1));
  for start = 0:2^n - 1
    seen = -ones (1, 2^n);
    seen(start + 1) = 0;
    frontier = start;
    depth = 0;
    while (~isempty (frontier))
      depth = depth + 1;
      frontier = reshape (adjacency(frontier + 1, :).', 1, []);
      frontier = frontier(seen(frontier + 1) < 0);
      seen(frontier + 1) = depth;
    end
    if (any (seen ~= distance (start, 0:2^n - 1)))
      ok = false;
      return;
    end
  end
  ok = true;
end

function out = hamming_7_4 ()
  rows_g = code_of ('1000110', '0100101', '0010011', '0001111');
  out = zeros (1, 16);
  for m = 0:15
    c = 0;
    for i = 0:3
      if (bitand (bitshift (m, -(3 - i)), 1))
        c = bitxor (c, rows_g(i + 1));
      end
    end
    out(m + 1) = c;
  end
end

function check_detection ()
  section ('Detecting Errors');
  sent = word ('01101');
  pattern = word ('00110');
  r = bitxor (sent, pattern);
  same ('01101 with errors 00110', show (r, 5), '01011');
  same ('distances from 01011 to C5', distance (r, C5), [3, 2, 4, 1]);
  require ('01011 is not a codeword', ~ismember (r, C5));
  same ('01101 with errors 10110', show (bitxor (sent, word ('10110')), 5), '11011');
  e = 1:31;
  e = e(weight (e) <= 2);
  received = bsxfun (@bitxor, C5.', e);
  require ('C5 flags every pattern of weight 1 or 2 from every codeword', ~any (ismember (received(:), C5)));
  require ('detection bound: d - 1 detected, some weight-d pattern missed, on C5, C6, R4, R5', ...
           detects_exactly (C5, 5) && detects_exactly (C6, 6) && detects_exactly (R4, 4) && detects_exactly (R5, 5));
  same ('weight distribution of C5', arrayfun (@(k) sum (weight (C5) == k), 0:5), [1, 0, 0, 2, 1, 0]);
  check ('C5 undetected-error probability, p = 0.01', undetected_probability (C5, 5, 0.01), '1.9701e-6');
  check ('2 p^3 (1-p)^2 + p^4 (1-p)', 2 * 0.01^3 * 0.99^2 + 0.01^4 * 0.99, '1.9701e-6');
  check ('C5 chance of any error, p = 0.01', any_error (5, 0.01), '0.049010');
  same ('error patterns of weight 3 on 5 bits', nchoosek (5, 3), 10);
  e = 0:31;
  e = e(weight (e) == 3);
  escapes = arrayfun (@(c) sum (ismember (bitxor (c, e), C5)), C5);
  same ('weight-3 patterns that escape, from each codeword of C5', num2cell (unique (escapes)), {2});
end

function ok = detects_exactly (code, n)
  d = min_distance (code);
  e = 1:2^n - 1;
  hits = ismember (bsxfun (@bitxor, code(:), e), code);
  w = weight (e);
  caught = ~any (any (hits(:, w < d)));
  missed = any (any (hits(:, w == d)));
  ok = caught && missed;
end

function check_correction ()
  section ('Correcting Errors');
  r = word ('10111');
  same ('distances from 10111 to C5', distance (r, C5), [4, 3, 1, 2]);
  same ('10111 decodes to', show (decode_nearest (C5, r), 5), '10110');
  r = word ('10000');
  same ('distances from 10000 to C5', distance (r, C5), [1, 4, 2, 3]);
  same ('10000 decodes to', show (decode_nearest (C5, r), 5), '00000');
  d = 1:7;
  % floor of a quotient matches Python's // for these exact doubles.
  same ('t for d = 1..7', floor ((d - 1) / 2), [0, 0, 1, 1, 2, 2, 3]);
  same ('d - 1 for d = 1..7', d - 1, [0, 1, 2, 3, 4, 5, 6]);
  % A NaN from a decoder compares unequal to every codeword, as None does.
  ok = true;
  for c = C5
    for i = 0:4
      ok = ok && decode_nearest (C5, bitxor (c, 2^i)) == c;
    end
  end
  require ('C5 corrects every single error from every codeword', ok);
  e = 0:31;
  ok = true;
  for c = R5
    for x = e(weight (e) <= 2)
      ok = ok && decode_nearest (R5, bitxor (c, x)) == c;
    end
  end
  require ('R5 corrects every pattern of weight 2 or less', ok);
  fails_c5 = false;
  for c = C5
    for x = e(weight (e) == 2)
      fails_c5 = fails_c5 || decode_nearest (C5, bitxor (c, x)) ~= c;
    end
  end
  fails_r5 = false;
  for c = R5
    for x = e(weight (e) == 3)
      fails_r5 = fails_r5 || decode_nearest (R5, bitxor (c, x)) ~= c;
    end
  end
  require ('correction bound is tight: some weight-(t+1) pattern fails, C5 and R5', fails_c5 && fails_r5);
  same ('R4: 0001 with t = 1', show (decode_bounded (R4, word ('0001'), 1), 4), '0000');
  same ('R4: 0011 with t = 1', none_or (decode_bounded (R4, word ('0011'), 1)), py_none ());
  require ('R4 and C6 correct 1 and flag 2 from every codeword', ...
           correct_and_detect (R4, 4, 1, 2) && correct_and_detect (C6, 6, 1, 2));
  require ('R5 corrects 1 and flags 3 (t + s = 4 = d - 1)', correct_and_detect (R5, 5, 1, 3));
  require ('R5 with t = 2 cannot also flag 3', ~correct_and_detect (R5, 5, 2, 3));
end

function ok = correct_and_detect (code, n, t, s)
  for c = code
    for e = 1:2^n - 1
      out = decode_bounded (code, bitxor (c, e), t);
      w = weight (e);
      if (w <= t && out ~= c)
        ok = false;
        return;
      end
      if (t < w && w <= s && ~isnan (out))
        ok = false;
        return;
      end
    end
  end
  ok = true;
end

function check_erasures ()
  section ('Errors and Erasures Together');
  same ('minimum distance of E6', min_distance (E6), 4);
  received = apply_errors_erasures (word ('110011'), 6, 5, 2);
  same ('110011, position 3 erased, position 6 flipped', received, '11?010');
  kept = [0, 1, 3, 4, 5];
  bits = show (E6.', 6);
  gaps = sum (bits(:, kept + 1) ~= received(kept + 1), 2);
  same ('distances on the five kept positions', gaps, [3, 4, 1, 2]);
  same ('11?010 decodes to', show (decode_erasures (E6, 6, received), 6), '110011');
  pairs_d4 = zeros (0, 2);
  for t = 0:2
    for e = 0:4
      if (2 * t + e < 4)
        pairs_d4(end + 1, :) = [t, e];
      end
    end
  end
  same ('(t, e) allowed by d = 4', tuple_list (pairs_d4), tuple_list ([0 0; 0 1; 0 2; 0 3; 1 0; 1 1]));
  pairs_d5 = zeros (0, 2);
  for t = 0:2
    for e = 0:5
      if (2 * t + e < 5)
        pairs_d5(end + 1, :) = [t, e];
      end
    end
  end
  same ('(t, e) allowed by d = 5', tuple_list (pairs_d5), ...
        tuple_list ([0 0; 0 1; 0 2; 0 3; 0 4; 1 0; 1 1; 1 2; 2 0]));
  require ('E6 decodes every combination with 2t + e < 4', erasure_bound_holds (E6, 6));
  require ('R5 decodes every combination with 2t + e < 5', erasure_bound_holds (R5, 5));
  require ('C5 decodes every combination with 2t + e < 3', erasure_bound_holds (C5, 5));
  require ('some combination with 2t + e = d fails, for E6, R5 and C5', ...
           erasure_bound_tight (E6, 6) && erasure_bound_tight (R5, 5) && erasure_bound_tight (C5, 5));
  v = 0:15;
  require ('one erasure in an even-parity word is always filled (d = 2, t = 0, e = 1)', ...
           erasure_bound_holds (v(mod (weight (v), 2) == 0), 4));
end

function ok = erasure_bound_holds (code, n)
  d = min_distance (code);
  for c = code
    for e_count = 0:n
      for t = 0:n
        if (2 * t + e_count >= d)
          continue;
        end
        erased_sets = combos (0:n - 1, e_count);
        for k = 1:rows (erased_sets)
          erased = erased_sets(k, :);
          flip_sets = combos (setdiff (0:n - 1, erased), t);
          for m = 1:rows (flip_sets)
            if (decode_erasures (code, n, apply_errors_erasures (c, n, flip_sets(m, :), erased)) ~= c)
              ok = false;
              return;
            end
          end
        end
      end
    end
  end
  ok = true;
end

function ok = erasure_bound_tight (code, n)
  d = min_distance (code);
  for t = 0:floor (d / 2)
    if (~some_decoding_fails (code, n, t, d - 2 * t))
      ok = false;
      return;
    end
  end
  ok = true;
end

function found = some_decoding_fails (code, n, t, e_count)
  % The any(...) of Python's erasure_bound_tight: some codeword, some
  % e_count erasures and some t flips elsewhere that decode wrongly.
  erased_sets = combos (0:n - 1, e_count);
  for c = code
    for k = 1:rows (erased_sets)
      erased = erased_sets(k, :);
      flip_sets = combos (setdiff (0:n - 1, erased), t);
      for m = 1:rows (flip_sets)
        if (decode_erasures (code, n, apply_errors_erasures (c, n, flip_sets(m, :), erased)) ~= c)
          found = true;
          return;
        end
      end
    end
  end
  found = false;
end

function check_parity ()
  section ('Single Parity');
  p = frac (1, 100);
  ok = true;
  for n = 1:10
    ok = ok && frac_eq (parity_alarm (n, p), parity_alarm_by_enumeration (n, p));
  end
  require ('closed form equals enumeration for n = 1..10, p = 1/100', ok);
  require ('closed form equals enumeration for n = 8, p = 1/3', ...
           frac_eq (parity_alarm (8, frac (1, 3)), parity_alarm_by_enumeration (8, frac (1, 3))));
  check ('n = 8, p = 0.01: any error', frac_float (any_error (8, p)), '0.077255');
  check ('n = 8, p = 0.01: 0.98^8', frac_float (frac_pow (frac (98, 100), 8)), '0.850763');
  check ('n = 8, p = 0.01: alarm', frac_float (parity_alarm (8, p)), '0.074618');
  check ('n = 8, p = 0.01: undetected', frac_float (frac_sub (any_error (8, p), parity_alarm (8, p))), '0.0026368');
  check ('n = 8, p = 0.01: coverage', frac_float (parity_coverage (8, p)), '0.96587');
  check ('n = 8, p = 0.01: coverage in percent', frac_float (parity_coverage (8, p)), '96.6%');
  check ('small-p estimate 1 - (n-1)p/2, n = 8', 1 - 7 * 0.01 / 2, '0.965');
  check ('n = 8, p = 1/2: coverage', frac_float (parity_coverage (8, frac (1, 2))), '0.50196');
  same ('n = 8, p = 1/2: coverage as a fraction', frac_text (parity_coverage (8, frac (1, 2))), fraction_shown (128, 255));
  check ('n = 32, p = 0.01: coverage', frac_float (parity_coverage (32, p)), '0.86561');
  check ('n = 32, p = 0.01: alarm', frac_float (parity_alarm (32, p)), '0.23806');
  check ('n = 32, p = 0.01: any error', frac_float (any_error (32, p)), '0.27502');
  check ('n = 32, p = 0.01: corrupted words that slip through', ...
         frac_float (frac_sub (frac (1, 1), parity_coverage (32, p))), '13%');
  check ('n = 8, p = 0.01: corrupted words that slip through', ...
         frac_float (frac_sub (frac (1, 1), parity_coverage (8, p))), '3.4%');
  check ('n = 64, p = 0.01: coverage', frac_float (parity_coverage (64, p)), '0.76469');
  check ('n = 64, p = 0.001: coverage', frac_float (parity_coverage (64, frac (1, 1000))), '0.96944');
  x = word ('10110110');
  e = 0:255;
  odd = mod (weight (e), 2) == 1;
  require ('injected flips: odd counts always caught, even never, n = 8', ...
           all (mod (weight (bitxor (x, e(odd))), 2) ~= mod (weight (x), 2)) ...
           && all (mod (weight (bitxor (x, e(~odd))), 2) == mod (weight (x), 2)));
end

function check_isbn ()
  section ('Check Digits: ISBN-10');
  payload = digits ('0-13-042672');
  products = (10 - (0:numel (payload) - 1)) .* payload;
  same ('ISBN-10 products', products, [0, 9, 24, 0, 24, 10, 24, 21, 4]);
  same ('ISBN-10 payload sum', sum (products), 116);
  same ('116 mod 11', mod (116, 11), 6);
  same ('ISBN-10 check digit', isbn10_check (payload), 5);
  require ('0-13-042672-5 is valid', isbn10_valid (digits ('0-13-042672-5')));
  swapped = digits ('0-13-042627-5');
  same ('swapped ISBN-10 weighted sum', isbn10_sum (swapped), 116);
  same ('swapped ISBN-10 sum mod 11', mod (isbn10_sum (swapped), 11), 6);
  require ('the swapped ISBN-10 is caught', ~isbn10_valid (swapped));
  same ('payload 0-12-345672: sum', isbn10_sum (digits ('0-12-345672')), 144);
  same ('payload 0-12-345672: check', isbn10_check (digits ('0-12-345672')), 10);
  require ('0-12-345672-X is valid', isbn10_valid (digits ('0-12-345672-X')));
  rates = scheme_rates (@isbn10_valid, digits ('0-13-042672-5'));
  same ('ISBN-10 substitutions caught', fraction_shown (rates.substitution), fraction_shown (1, 1));
  same ('ISBN-10 adjacent swaps caught', fraction_shown (rates.adjacent), fraction_shown (1, 1));
  same ('ISBN-10 jump swaps caught', fraction_shown (rates.jump), fraction_shown (1, 1));
  base = digits ('0-13-042672-5');
  ok = true;
  for i = 0:8
    for j = i + 1:9
      for a = 0:9
        for b = 0:9
          if (a ~= b)
            ok = ok && ~swap_escapes (@isbn10_valid, base, i, j, a, b);
          end
        end
      end
    end
  end
  require ('ISBN-10 catches swaps at every distance', ok);

  section ('Check Digits: ISBN-13');
  payload = digits ('978-0-13-042672');
  terms = isbn13_weights (12) .* payload;
  same ('ISBN-13 terms', terms, [9, 21, 8, 0, 1, 9, 0, 12, 2, 18, 7, 6]);
  same ('ISBN-13 payload sum', sum (terms), 93);
  same ('ISBN-13 check digit', isbn13_check (payload), 7);
  require ('978-0-13-042672-7 is valid', isbn13_valid (digits ('978-0-13-042672-7')));
  swapped = digits ('978-0-13-042627-7');
  same ('swapped ISBN-13 sum', isbn13_sum (swapped), 110);
  require ('the swapped ISBN-13 slips through', isbn13_valid (swapped));
  same ('3 x 7 mod 10', mod (3 * 7, 10), 1);
  rates = scheme_rates (@isbn13_valid, digits ('978-0-13-042672-7'));
  same ('ISBN-13 substitutions caught', fraction_shown (rates.substitution), fraction_shown (1, 1));
  same ('ISBN-13 adjacent swaps caught', fraction_shown (rates.adjacent), fraction_shown (80, 90));
  % 80/90 and 8/9 are one rational number, so both round to one double.
  check ('ISBN-13 adjacent swaps caught, percent', rates.adjacent(1) / rates.adjacent(2), '88.9%');
  same ('ISBN-13 jump swaps caught', fraction_shown (rates.jump), fraction_shown (0, 1));
  escaping = escaping_adjacent_pairs (@isbn13_valid, digits ('978-0-13-042672-7'));
  [b, a] = meshgrid (0:9, 0:9);
  five_apart = sortrows ([a(abs (a - b) == 5), b(abs (a - b) == 5)]);
  same ('ISBN-13 adjacent swaps that escape', tuple_list (sortrows (escaping)), tuple_list (five_apart));
  same ('count of escaping ordered pairs', rows (escaping), 10);
end

function check_luhn ()
  section ('Check Digits: the Luhn formula');
  same ('doubling map D(0..9)', luhn_double (0:9), [0, 2, 4, 6, 8, 1, 3, 5, 7, 9]);
  require ('D is a permutation of 0..9', isequal (sort (luhn_double (0:9)), 0:9));
  same ('f(x) = D(x) - x', luhn_double (0:9) - (0:9), [0, 1, 2, 3, 4, -4, -3, -2, -1, 0]);
  payload = digits ('314159265');
  same ('Luhn check digit of 314159265', luhn_check (payload), 3);
  number = digits ('3141592653');
  same ('Luhn terms of 3141592653, left to right', luhn_terms (number), [6, 1, 8, 1, 1, 9, 4, 6, 1, 3]);
  same ('Luhn sum of 3141592653', luhn_sum (number), 40);
  sub = digits ('3141582653');
  same ('Luhn sum after 9 -> 8', luhn_sum (sub), 39);
  require ('the substitution is caught', ~luhn_valid (sub));
  swap = digits ('3141952653');
  same ('Luhn terms after swapping 5 and 9', luhn_terms (swap), [6, 1, 8, 1, 9, 5, 4, 6, 1, 3]);
  same ('Luhn sum after swapping 5 and 9', luhn_sum (swap), 44);
  require ('the swap is caught', ~luhn_valid (swap));
  bases = {digits('3141592653'), digits('4000123412341234')};
  for k = 1:numel (bases)
    base = make_valid (@luhn_valid, bases{k}, [], []);
    rates = scheme_rates (@luhn_valid, base);
    n = numel (base);
    same (sprintf ('Luhn substitutions caught, length %d', n), fraction_shown (rates.substitution), fraction_shown (1, 1));
    same (sprintf ('Luhn adjacent swaps caught, length %d', n), fraction_shown (rates.adjacent), fraction_shown (88, 90));
    same (sprintf ('Luhn jump swaps caught, length %d', n), fraction_shown (rates.jump), fraction_shown (0, 1));
    same (sprintf ('Luhn adjacent swaps that escape, length %d', n), ...
          tuple_list (sortrows (escaping_adjacent_pairs (@luhn_valid, base))), tuple_list ([0 9; 9 0]));
  end
  check ('Luhn adjacent swaps caught, percent', 88 / 90, '97.8%');
  exercise = digits ('50961');
  same ('Luhn check digit of 50961', luhn_check (exercise), 2);
  same ('Luhn sum of 509612', luhn_sum (digits ('509612')), 20);
  same ('Luhn sum of 590612', luhn_sum (digits ('590612')), 20);
  require ('509612 and 590612 both pass', luhn_valid (digits ('509612')) && luhn_valid (digits ('590612')));

  section ('Check digits confirmed on every valid number of five digits');
  % Every five-digit string at once, one per row, most significant digit first.
  all5 = mod (floor ((0:99999).' ./ 10.^(4:-1:0)), 10);
  luhn5 = all5(luhn_valid (all5), :);
  same ('valid 5-digit Luhn numbers', rows (luhn5), 10000);
  missed_subs = 0;
  for i = 0:4
    for b = 0:9
      changed = luhn5(luhn5(:, i + 1) ~= b, :);
      changed(:, i + 1) = b;
      missed_subs = missed_subs + sum (luhn_valid (changed));
    end
  end
  same ('single substitutions that pass', missed_subs, 0);
  escaped = zeros (0, 2);
  for i = 0:3
    swapped = luhn5;
    swapped(:, [i, i + 1] + 1) = luhn5(:, [i + 1, i] + 1);
    hit = luhn5(:, i + 1) ~= luhn5(:, i + 2) & luhn_valid (swapped);
    escaped = [escaped; luhn5(hit, [i, i + 1] + 1)];
  end
  same ('adjacent swaps that pass', tuple_list (unique (escaped, 'rows')), tuple_list ([0 9; 9 0]));
  w13 = [1, 3, 1, 3, 1];
  toy13 = all5(mod (all5 * w13.', 10) == 0, :);
  escaped13 = zeros (0, 2);
  for i = 0:3
    swapped = toy13;
    swapped(:, [i, i + 1] + 1) = toy13(:, [i + 1, i] + 1);
    hit = toy13(:, i + 1) ~= toy13(:, i + 2) & mod (swapped * w13.', 10) == 0;
    escaped13 = [escaped13; toy13(hit, [i, i + 1] + 1)];
  end
  same ('weights 1,3,1,3,1: adjacent swaps that pass', rows (unique (escaped13, 'rows')), 10);
end

function a = with_errors (zero, cells)
  % Python's with_errors is a closure over zero. Here zero is passed in.
  % Each row of cells is a 0-based (row, column) to flip.
  a = zero;
  for k = 1:rows (cells)
    i = cells(k, 1) + 1;
    j = cells(k, 2) + 1;
    a(i, j) = bitxor (a(i, j), 1);
  end
end

function check_rectangular ()
  section ('Rectangular Parity');
  data = [1 0 1 1; 0 1 1 0; 1 1 0 1; 0 0 1 0];
  array = rect_encode (data);
  same ('row parities', array(1:4, 5), [1, 0, 1, 1]);
  same ('column parities', array(5, 1:4), [0, 0, 1, 0]);
  same ('corner bit', array(5, 5), 1);
  same ('parity of the row parities', mod (sum (array(1:4, 5)), 2), 1);
  same ('parity of the column parities', mod (sum (array(5, 1:4)), 2), 1);
  hit = array;
  hit(3, 2) = bitxor (hit(3, 2), 1);
  same ('failures after an error at row 3, column 2', failures_shown (hit), tuple ({2}, {1}));
  [fixed, verdict] = rect_decode (hit);
  require ('the decoder flips row 3, column 2 back', strcmp (verdict, 'corrected') && isequal (fixed, array));
  same ('bits in the array', 25, 25);
  check ('rate 16/25', 16 / 25, '0.64');
  % All 2^16 arrays at once, one per page of a 4 x 4 x 65536 stack. Bit
  % 15 - 4i - j of value goes to row i, column j, counted from 0.
  value = reshape (0:2^16 - 1, 1, 1, []);
  stack = zeros (4, 4, 2^16);
  for i = 0:3
    for j = 0:3
      stack(i + 1, j + 1, :) = bitand (bitshift (value, -(15 - 4 * i - j)), 1);
    end
  end
  k = weight (flatten (rect_encode (stack)));
  weights = accumarray (k(:) + 1, 1);   % weights(k + 1) arrays have weight k
  present = find (weights).' - 1;
  same ('smallest nonzero weight over all 2^16 arrays', min (present(present > 0)), 4);
  same ('arrays of weight 4', weights(4 + 1), 100);
  same ('C(5,2)^2 rectangles', nchoosek (5, 2)^2, 100);
  zero = rect_encode (zeros (4, 4));
  [j, i] = meshgrid (0:4, 0:4);
  cells = sortrows ([i(:), j(:)]);   % (0,0), (0,1), ... as Python lists them
  singles = cell (1, 25);
  for c = 1:25
    [~, singles{c}] = rect_decode (with_errors (zero, cells(c, :)));
  end
  pairs = combos (1:25, 2);
  doubles = cell (1, rows (pairs));
  for c = 1:rows (pairs)
    [~, doubles{c}] = rect_decode (with_errors (zero, cells(pairs(c, :), :)));
  end
  triples = combos (1:25, 3);
  broken = 0;
  for c = 1:rows (triples)
    [r, col] = rect_failures (with_errors (zero, cells(triples(c, :), :)));
    broken = broken + (~isempty (r) || ~isempty (col));
  end
  same ('single errors corrected', sum (strcmp (singles, 'corrected')), 25);
  same ('double errors flagged', sum (strcmp (doubles, 'flagged')), 300);
  same ('triple errors that break some parity', broken, nchoosek (25, 3));
  same ('two errors in row 1', failures_shown (with_errors (zero, [0 0; 0 2])), tuple ({}, {0, 2}));
  same ('errors at (1,1) and (3,4)', failures_shown (with_errors (zero, [0 0; 2 3])), tuple ({0, 2}, {0, 3}));
  mimic = failures_shown (with_errors (zero, [0 0; 0 2; 2 0]));
  same ('three errors that mimic one', mimic, tuple ({2}, {2}));
  check ('R4 rate', 1 / 4, '0.25');
  check ('8 x 8 array rate 64/81', 64 / 81, '0.790');
end

function check_gray ()
  section ('Gray Codes');
  g3 = reflected (3);
  same ('3-bit reflected Gray code', g3, {'000', '001', '011', '010', '110', '111', '101', '100'});
  g4 = reflected (4);
  same ('4-bit reflected Gray code', g4, ...
        {'0000', '0001', '0011', '0010', '0110', '0111', '0101', '0100', ...
         '1100', '1101', '1111', '1110', '1010', '1011', '1001', '1000'});
  ok = true;
  for n = 1:8
    ok = ok && isequal (reflected (n), strings_of (gray (0:2^n - 1), n));
  end
  require ('reflected code equals b XOR (b >> 1), n = 1..8', ok);
  ok = true;
  for n = 1:8
    b = 0:2^n - 1;
    ok = ok && all (distance (gray (b), gray (mod (b + 1, 2^n))) == 1);
  end
  require ('neighbors differ in one bit, wrap included, n = 1..8', ok);
  b = 0:254;
  require ('the bit that flips from b to b+1 is the count of trailing ones of b', ...
           all (bitxor (gray (b), gray (b + 1)) == 2.^trailing_ones (b)));
  b = 0:4095;
  require ('gray_inverse undoes gray for all 12-bit values', all (gray_inverse (gray (b)) == b));
  same ('1011 XOR 0101', show (bitxor (word ('1011'), word ('0101')), 4), '1110');
  same ('Gray code of 1011', show (gray (word ('1011')), 4), '1110');
  same ('position 11 of the 4-bit list', g4{11 + 1}, '1110');
  same ('binary of Gray 1110', show (gray_inverse (word ('1110')), 4), '1011');
  same ('Gray code of 1101', show (gray (word ('1101')), 4), '1011');
  same ('binary of Gray 0111', show (gray_inverse (word ('0111')), 4), '0101');
  same ('position 5 of the 4-bit list', g4{5 + 1}, '0111');
  same ('bits that change from 011 to 100 in binary', distance (word ('011'), word ('100')), 3);
  same ('bits that change from 3 to 4 in Gray', distance (gray (3), gray (4)), 1);
  same ('Gray codes of 3 and 4', tuple (show (gray (3), 3), show (gray (4), 3)), tuple ('010', '110'));
  % product(("0", "1"), repeat=3): every 3-bit string, one per row.
  [x, y, z] = ndgrid ('01', '01', '01');
  readings = unique (bin2dec ([x(:), y(:), z(:)]));
  same ('binary readings possible mid-transition 3 -> 4', numel (readings), 8);
  same ('Gray readings possible mid-transition 3 -> 4', sort (gray_inverse ([word('010'), word('110')])), [3, 4]);
  same ('binary positions of 111 and 000', tuple (word ('111'), word ('000')), tuple (7, 0));
  same ('Gray reading 111 decodes to', gray_inverse (word ('111')), 5);
  low = gray_inverse (bitand (gray (6), gray (7)));
  high = gray_inverse (bitor (gray (6), gray (7)));
  same ('Gray mid-transition readings, 6 -> 7', unique ([low, high]), [6, 7]);
end

function count = trailing_ones (b)
  % The loop of Python's trailing_ones, run on every element at once: an
  % element drops out once its lowest bit is 0.
  count = zeros (size (b));
  active = bitand (b, 1) == 1;
  while (any (active))
    count(active) = count(active) + 1;
    b(active) = bitshift (b(active), -1);
    active = active & bitand (b, 1) == 1;
  end
end

function check_examples ()
  section ('Worked example: from distance three to four');
  extended = bitor (bitshift (C5, 1), mod (weight (C5), 2));
  same ('C5 with an overall parity bit', sort (strings_of (extended, 6)), sort (strings_of (C6, 6)));
  pairs = combos (C6, 2);
  same ('C6 distances', sort (distance (pairs(:, 1), pairs(:, 2))), [4, 4, 4, 4, 4, 4]);
  same ('minimum distance of C6', min_distance (C6), 4);
  require ('appending parity to any code of odd d gives d + 1, all 4-word codes from 5-bit words with d = 3', ...
           parity_extension_holds ());
  same ('C6 corrects 000001 to', show (decode_bounded (C6, word ('000001'), 1), 6), '000000');
  same ('C6 with two errors 000011', none_or (decode_bounded (C6, word ('000011'), 1)), py_none ());

  section ('Worked example: R5 under errors and erasures');
  same ('1?01? decodes to', show (decode_erasures (R5, 5, '1?01?'), 5), '11111');
  same ('0?01? decodes to', show (decode_erasures (R5, 5, '0?01?'), 5), '00000');
  same ('????1 decodes to', show (decode_erasures (R5, 5, '????1'), 5), '11111');
  same ('2t + e for one error, two erasures', 2 * 1 + 2, 4);
  same ('2t + e for two errors, two erasures', 2 * 2 + 2, 6);
  same ('10?0? decodes to', show (decode_erasures (R5, 5, '10?0?'), 5), '00000');

  section ('Worked example: rectangular parity with two errors');
  data = [1 0 1 1; 0 1 1 0; 1 1 0 1; 0 0 1 0];
  array = rect_encode (data);
  hit = with_errors (array, [0 1; 0 3]);
  same ('errors at (1,2) and (1,4): failures', failures_shown (hit), tuple ({}, {1, 3}));
  [~, verdict] = rect_decode (hit);
  same ('verdict', verdict, 'flagged');
  hit = with_errors (array, [0 1; 2 3]);
  same ('errors at (1,2) and (3,4): failures', failures_shown (hit), tuple ({0, 2}, {1, 3}));
  [~, verdict] = rect_decode (hit);
  same ('verdict', verdict, 'flagged');
  other = with_errors (array, [0 3; 2 1]);
  [rows_o, cols_o] = rect_failures (other);
  [rows_h, cols_h] = rect_failures (hit);
  require ('errors at (1,4) and (3,2) give the same failures', isequal (rows_o, rows_h) && isequal (cols_o, cols_h));
end

function ok = parity_extension_holds ()
  % All 35,960 four-word codes from 5-bit words at once, one code per row.
  % Python stops at the first code that breaks the rule. Testing them all
  % gives the same verdict.
  codes = combos (0:31, 4);
  d = min_distance (codes);
  odd = mod (d, 2) == 1;
  extended = bitor (bitshift (codes(odd, :), 1), mod (weight (codes(odd, :)), 2));
  ok = all (min_distance (extended) == d(odd) + 1);
end

function check_exercises ()
  section ('Exercises');
  code = code_of ('00000', '00111', '11001', '11110');
  pairs = combos (code, 2);
  same ('exercise code distances', sort (distance (pairs(:, 1), pairs(:, 2))), [3, 3, 3, 3, 4, 4]);
  same ('exercise code minimum distance', min_distance (code), 3);
  same ('V(7, 2)', ball_size (7, 2), 29);
  same ('floor(128 / 29)', floor (128 / 29), 4);
  check ('128 / 29', 128 / 29, '4.41');
  % All 341,376 triples of 7-bit words at once, one per row.
  triples = combos (0:127, 3);
  d_ab = distance (triples(:, 1), triples(:, 2));
  d_bc = distance (triples(:, 2), triples(:, 3));
  d_ac = distance (triples(:, 1), triples(:, 3));
  sums = unique (d_ab + d_bc + d_ac);
  require ('three 7-bit words: the three distances sum to an even number at most 14', ...
           all (mod (sums, 2) == 0) && max (sums) == 14);
  require ('no three 7-bit words lie pairwise at distance 5 or more', ...
           ~any (min ([d_ab, d_bc, d_ac], [], 2) >= 5));
  allowed = zeros (4, 2);
  e = 0:6;
  for t = 0:3
    allowed(t + 1, :) = [t, max(e(2 * t + e < 7))];
  end
  same ('(t, e) allowed by d = 7', tuple_list (allowed), tuple_list ([0 6; 1 4; 2 2; 3 0]));
  check ('8 x 8 rectangular rate', 64 / 81, '0.790');
  same ('V(12, 1)', ball_size (12, 1), 13);
  same ('floor(4096 / 13)', floor (4096 / 13), 315);
  require ('three words of any length have d(x,y) + d(y,z) + d(x,z) <= 2n, n = 7', ...
           max (d_ab + d_bc + d_ac) == 14);
end


% Main -------------------------------------------------------------------------

check_distance_and_weight ();
check_minimum_distance ();
check_cube ();
check_detection ();
check_correction ();
check_erasures ();
check_parity ();
check_isbn ();
check_luhn ();
check_rectangular ();
check_gray ();
check_examples ();
check_exercises ();
printf ('\nAll chapter 8 numbers reproduced.\n');
