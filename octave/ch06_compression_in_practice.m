% Reproduce every number in chapter 6 of Error Correction Coding, "Compression in Practice".
%
% This is the GNU Octave twin of python/ch06_compression_in_practice.py. It needs
% core Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch06_compression_in_practice.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. Codewords, parses,
% tables and decoded strings must match the chapter exactly. Every encoder here
% has a decoder, and each pair is run as a round trip. A run that reaches the
% final line has reproduced every number in the chapter. The two routes print
% the same lines, and tools/compare-routes.sh checks that against the Python.
%
% The Python module also holds the arithmetic-coding, LZ78 and quantizer functions
% that ch06_charts.py plots. Charts stay Python-only, and the functions are here
% because the checks use them too.
%
% Differences from the Python that a side-by-side reader meets:
%   - Python's fractions.Fraction is exact. Here a Fraction is a struct with the
%     fields num and den in lowest terms, built by frac and combined by frac_add,
%     frac_sub, frac_mul and frac_div. Each numerator and denominator must stay
%     below 2^53, the largest range in which a double holds every integer exactly,
%     and frac stops the program with an error rather than round one silently.
%   - Two places in the chapter need integers far beyond 2^53: the plain rANS
%     state after 300 letters, near 2^350, and the denominator of Laplace's
%     probability for 1000 letters, near 2^1166. Those are carried as rows of
%     24-bit limbs, least significant first, by the big_ functions below.
%   - A Python dict keyed by letter becomes an array in the order of the letters,
%     so PROBS["a"] is PROBS(LETTERS == 'a') here.
%   - Octave has one number type where Python has two, and py_str prints a double
%     as a Python float. Python ints are converted to int64 before printing, and
%     int_list turns a vector into the cell row that prints as a Python list, so
%     that a list of one item still prints with its brackets.
%   - A Fraction prints as 3/2 on its own and as Fraction(3, 2) inside a tuple,
%     as Python's str and repr do. py_text below adds that, and None, to py_str.
%   - Octave indexes from 1, so the Python's steps[k] is steps(k + 1) here.
%     Quantities the Python counts from 0, such as offsets, codes and row numbers,
%     keep their values, and 1 is added only where one indexes an array.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's.
%     Python's round sends halves to the even neighbor, and py_round does too.
%   - math.lgamma is gammaln here. CPython computes lgamma with its own code
%     rather than the C library's, so the two can differ in the last bit, which
%     never reaches the nine digits a check line prints.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main runs last. The
%     module constants of the Python are globals, declared in every function
%     that reads them.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function check (label, computed, stated)
  % A Fraction is divided once, which rounds the exact quotient the way
  % Python's float () of a Fraction does.
  computed = to_float (computed);
  value = str2double (strrep (stated, ',', ''));
  printf ('  %-62s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % Python compares the two values with !=. Here two numbers are compared by
  % value, so that the float 9.0 equals the int 9 as it does in Python, and
  % anything else through its repr text, which for the strings, ints, floats,
  % Fractions, tuples and lists of this chapter is equal exactly when the
  % values are.
  printf ('  %-62s %-16s %s\n', label, py_text (computed), py_text (stated));
  if (~py_equal (computed, stated))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_text (computed, true), py_text (stated, true));
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
    printf ('MISMATCH in ''%s''\n', label);
    exit (1);
  end
end

function h = entropy (probs)
  % probs is a row of floats, or a struct row of Fractions, which float () turns
  % into one division each.
  if (isstruct (probs))
    probs = [probs.num] ./ [probs.den];
  end
  probs = probs(probs > 0);
  h = py_sum (probs .* log2 (1 ./ probs));
end

function h = binary_entropy (p)
  if (p == 0 || p == 1)
    h = 0.0;
  else
    h = -p * log2 (p) - (1 - p) * log2 (1 - p);
  end
end


% Octave stand-ins for Python's numbers and printing ------------------------------

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

function r = py_round (x)
  % Python's round (): a half goes to the even neighbor, where Octave's round
  % takes it away from zero.
  r = round (x);
  if (abs (x - fix (x)) == 0.5)
    r = 2 * round (x / 2);
  end
end

function n = bit_length (x)
  % int.bit_length () for a nonnegative integer below 2^53.
  if (x == 0)
    n = 0;
  else
    n = numel (dec2bin (x));
  end
end

function v = to_float (x)
  if (isstruct (x))
    v = x.num / x.den;
  else
    v = double (x);
  end
end

function x = py_none ()
  % Python's None, which prints as None.
  x = struct ('none', true);
end

function t = tup (varargin)
  % A Python tuple. Its integer items arrive as doubles and become int64, so that
  % py_str prints them as Python ints. Every number in this chapter's tuples is an int.
  items = varargin;
  for i = 1:numel (items)
    if (isa (items{i}, 'double') && isscalar (items{i}))
      items{i} = int64 (items{i});
    end
  end
  t = struct ('tuple', {items});
end

function c = int_list (v)
  % A list of Python ints, as the cell row that py_str prints with brackets.
  c = num2cell (int64 (v(:).'));
end

function s = py_text (x, quoted)
  % py_str, extended to Fractions, None, and tuples and lists that hold them.
  if (nargin < 2)
    quoted = false;
  end
  if (isstruct (x) && isfield (x, 'num'))
    if (quoted)
      s = sprintf ('Fraction(%d, %d)', x.num, x.den);
    elseif (x.den == 1)
      s = sprintf ('%d', x.num);
    else
      s = sprintf ('%d/%d', x.num, x.den);
    end
  elseif (isstruct (x) && isfield (x, 'none'))
    s = 'None';
  elseif (isstruct (x) && isfield (x, 'tuple'))
    items = cellfun (@(v) py_text (v, true), x.tuple, 'UniformOutput', false);
    if (numel (items) == 1)
      s = ['(' items{1} ',)'];
    else
      s = ['(' strjoin(items, ', ') ')'];
    end
  elseif (iscell (x))
    items = cellfun (@(v) py_text (v, true), x, 'UniformOutput', false);
    s = ['[' strjoin(items, ', ') ']'];
  else
    s = py_str (x, quoted);
  end
end

function eq = py_equal (a, b)
  if ((isnumeric (a) || islogical (a)) && isscalar (a) && (isnumeric (b) || islogical (b)) && isscalar (b))
    eq = double (a) == double (b);
  else
    eq = strcmp (py_text (a, true), py_text (b, true));
  end
end


% Fractions ----------------------------------------------------------------------

function exact_int (v)
  % A double holds every integer exactly only below 2^53. Past that a product
  % would be rounded without a word, so the program stops instead.
  if (abs (v) >= flintmax)
    error ('ch06: the integer %g lies beyond 2^53, where a double is no longer exact', v);
  end
end

function f = frac (n, d)
  % Fraction(n, d) in lowest terms, with a positive denominator. Adding 0 turns
  % a negative zero, which a product can leave, into the 0 Python prints.
  if (nargin < 2)
    d = 1;
  end
  exact_int (n);
  exact_int (d);
  if (d < 0)
    n = -n;
    d = -d;
  end
  g = gcd (n, d);
  f = struct ('num', n / g + 0, 'den', d / g);
end

function f = frac_dec (text)
  % Fraction("0.033614"): the decimal digits over a power of ten.
  dot = find (text == '.', 1);
  if (isempty (dot))
    f = frac (str2double (text));
  else
    digits = text(dot + 1:end);
    f = frac (str2double ([text(1:dot - 1) digits]), 10 ^ numel (digits));
  end
end

function f = frac_add (a, b)
  % Over the least common denominator, so that the 10^10 denominators of the
  % sparse example stay well below 2^53.
  g = gcd (a.den, b.den);
  left = a.num * (b.den / g);
  right = b.num * (a.den / g);
  exact_int (left);
  exact_int (right);
  f = frac (left + right, a.den * (b.den / g));
end

function f = frac_sub (a, b)
  f = frac_add (a, struct ('num', -b.num, 'den', b.den));
end

function f = frac_mul (a, b)
  % Each numerator is cancelled against the other denominator first.
  g1 = gcd (a.num, b.den);
  g2 = gcd (b.num, a.den);
  f = frac ((a.num / g1) * (b.num / g2), (a.den / g2) * (b.den / g1));
end

function f = frac_div (a, b)
  f = frac_mul (a, frac (b.den, b.num));
end

function f = frac_pow (a, k)
  f = frac (a.num ^ k, a.den ^ k);
end

function c = frac_cmp (a, b)
  % The sign of a - b, from the cross products a.num * b.den and b.num * a.den.
  left = a.num * b.den;
  right = b.num * a.den;
  exact_int (left);
  exact_int (right);
  c = sign (left - right);
end

function eq = frac_eq (a, b)
  % Both are in lowest terms, so equal Fractions have equal parts.
  eq = a.num == b.num && a.den == b.den;
end

function q = frac_ceil (a)
  % math.ceil of a Fraction. The quotient num / den is rounded, so the
  % remainder corrects it to the exact floor before the ceiling is taken.
  q = floor (a.num / a.den);
  r = a.num - q * a.den;
  if (r < 0)
    q = q - 1;
    r = r + a.den;
  elseif (r >= a.den)
    q = q + 1;
    r = r - a.den;
  end
  if (r > 0)
    q = q + 1;
  end
end


% Integers beyond 2^53, as rows of 24-bit limbs -----------------------------------

function x = big_from (v)
  % The limbs of a nonnegative integer v below 2^53, least significant first.
  x = mod (v, 16777216);
  v = floor (v / 16777216);
  while (v > 0)
    x(end + 1) = mod (v, 16777216);
    v = floor (v / 16777216);
  end
end

function x = big_mul_small (x, m)
  % x * m for m below 2^29, so that a limb times m plus the carry stays below 2^53.
  carry = 0;
  for i = 1:numel (x)
    t = x(i) * m + carry;
    carry = floor (t / 16777216);
    x(i) = t - carry * 16777216;
  end
  while (carry > 0)
    x(end + 1) = mod (carry, 16777216);
    carry = floor (carry / 16777216);
  end
end

function x = big_add_small (x, a)
  carry = a;
  i = 1;
  while (carry > 0)
    if (i > numel (x))
      x(i) = 0;
    end
    t = x(i) + carry;
    carry = floor (t / 16777216);
    x(i) = t - carry * 16777216;
    i = i + 1;
  end
end

function [q, r] = big_divmod_small (x, d)
  % divmod (x, d) for a small d, by long division from the top limb down.
  q = zeros (size (x));
  r = 0;
  for i = numel (x):-1:1
    t = r * 16777216 + x(i);
    q(i) = floor (t / d);
    r = t - q(i) * d;
  end
  while (numel (q) > 1 && q(end) == 0)
    q(end) = [];
  end
end

function bits = big_bits (x)
  % The binary digits of x, most significant first, as bin () prints them.
  bits = reshape (dec2bin (x(end:-1:1), 24).', 1, []);
  bits = bits(find (bits == '1', 1):end);
  if (isempty (bits))
    bits = '0';
  end
end

function v = big_to_double (x)
  % Exact only while x is below 2^53, which the caller checks.
  v = 0;
  for i = numel (x):-1:1
    v = v * 16777216 + x(i);
  end
end

function y = big_log2 (x)
  % math.log2 of a Python int. CPython rounds the int to 53 bits, half to even,
  % and takes log2 of that double. An int of 2^1024 or more overflows a double,
  % so CPython splits it as m * 2^e with m in [0.5, 1) and returns log2 (m) + e.
  bits = big_bits (x);
  nb = numel (bits);
  if (nb <= 53)
    y = log2 (bin2dec (bits));
    return;
  end
  hi = bin2dec (bits(1:53));
  half = bits(54) == '1';
  sticky = any (bits(55:end) == '1');
  if (half && (sticky || mod (hi, 2) == 1))
    hi = hi + 1;
  end
  m = hi / 2 ^ 53;
  e = nb;
  if (m == 1)
    m = 0.5;
    e = e + 1;
  end
  if (e <= 1024)
    y = log2 (pow2 (m, e));
  else
    y = log2 (m) + e;
  end
end

function x = big_from_primes (p, e)
  % The product of p(i) ^ e(i). Each power is applied in pieces below 2^29,
  % which big_mul_small can take.
  x = 1;
  for i = 1:numel (p)
    left = e(i);
    while (left > 0)
      piece = p(i);
      j = 1;
      while (j < left && piece * p(i) < 2 ^ 29)
        piece = piece * p(i);
        j = j + 1;
      end
      x = big_mul_small (x, piece);
      left = left - j;
    end
  end
end

function e = factorial_exponents (p, k)
  % The power of each prime p(i) in k!, by Legendre's formula: the sum of
  % floor (k / p(i) ^ j) over j = 1, 2, ...
  e = zeros (size (p));
  q = p;
  while (any (q <= k))
    e = e + floor (k ./ q);
    q = q .* p;
  end
end


% The running source: a 0.7, b 0.2, c 0.1 ---------------------------------------

function out = cumulative (probs, order)
  % Start of each letter's slice of [0, 1), letters taken in the given order.
  % probs holds one value per letter of order, Fractions or integer counts.
  if (isstruct (probs))
    run = frac (0);
    out = repmat (run, 1, numel (order));
    for k = 1:numel (order)
      out(k) = run;
      run = frac_add (run, probs(k));
    end
  else
    run = 0;
    out = zeros (1, numel (order));
    for k = 1:numel (order)
      out(k) = run;
      run = run + probs(k);
    end
  end
end


% Arithmetic coding ---------------------------------------------------------------

function steps = arith_intervals (message, probs, order)
  % The interval (low, width) after each letter, starting from (0, 1), in exact fractions.
  % steps(k + 1) is the Python's steps[k], a struct with the fields low and width.
  cum = cumulative (probs, order);
  low = frac (0);
  width = frac (1);
  steps = struct ('low', {low}, 'width', {width});
  for s = message
    j = find (order == s);
    low = frac_add (low, frac_mul (width, cum(j)));
    width = frac_mul (width, probs(j));
    steps(end + 1) = struct ('low', {low}, 'width', {width});
  end
end

function code = arith_codeword (low, width)
  % Shannon-Fano-Elias tag: l = ceil(log2(1/width)) + 1 bits of ceil(low * 2^l) / 2^l.
  l0 = 0;
  while (frac_cmp (frac (1, 2 ^ l0), width) > 0)   % smallest l0 with 2^-l0 <= width, which is ceil(log2(1/width))
    l0 = l0 + 1;
  end
  l = l0 + 1;
  v = frac_ceil (frac_mul (low, frac (2 ^ l)));
  code = dec2bin (v, l);
end

function code = arith_encode (message, probs, order)
  steps = arith_intervals (message, probs, order);
  code = arith_codeword (steps(end).low, steps(end).width);
end

function [decoded, values] = arith_decode (bits, n, probs, order)
  % Decode n letters from a tag, returning the letters and the scaled value before each step.
  cum = cumulative (probs, order);
  t = frac (bin2dec (bits), 2 ^ numel (bits));
  decoded = '';
  values = {};
  for step = 1:n
    values{end + 1} = t;
    for j = 1:numel (order)
      if (frac_cmp (cum(j), t) <= 0 && frac_cmp (t, frac_add (cum(j), probs(j))) < 0)
        decoded(end + 1) = order(j);
        t = frac_div (frac_sub (t, cum(j)), probs(j));
        break;
      end
    end
  end
end

function rate = arith_expected_rate (p, n)
  % Expected bits per letter of the arithmetic code on n letters of a binary source with P(1) = p.
  q = 1 - p;
  total = 0.0;
  for k = 0:n
    log_prob = gammaln (n + 1) - gammaln (k + 1) - gammaln (n - k + 1) + k * log (p) + (n - k) * log (q);
    info = -(k * log2 (p) + (n - k) * log2 (q));
    total = total + exp (log_prob) * (ceil (info - 1e-12) + 1);
  end
  rate = total / n;
end

function [out, low, high] = settle_bits (low, high)
  % Emit the leading bits shared by every number in [low, high), rescaling as they go.
  out = '';
  half = frac (1, 2);
  while (true)
    if (frac_cmp (high, half) <= 0)
      out(end + 1) = '0';
      low = frac_mul (frac (2), low);
      high = frac_mul (frac (2), high);
    elseif (frac_cmp (low, half) >= 0)
      out(end + 1) = '1';
      low = frac_sub (frac_mul (frac (2), low), frac (1));
      high = frac_sub (frac_mul (frac (2), high), frac (1));
    else
      return;
    end
  end
end

function prob = laplace_probability (message, alphabet)
  % Probability of a message under the adaptive model whose counts start at 1 for every letter.
  counts = ones (1, numel (alphabet));
  prob = frac (1);
  for s = message
    k = find (alphabet == s);
    prob = frac_mul (prob, frac (counts(k), sum (counts)));
    counts(k) = counts(k) + 1;
  end
end

function f = laplace_closed_form (message, alphabet)
  % (m - 1)! times the factorial of each letter's count, over (n + m - 1)!.
  % Python builds the factorials and lets Fraction cancel their common factor.
  % Here each prime's power in every factorial comes from Legendre's formula, the
  % common powers cancel, and only what is left is multiplied out, as limbs. The
  % result is a Fraction when both parts fit below 2^53, as they do for six
  % letters, and otherwise a struct with the limbs big_num and big_den, which
  % only bits_of reads.
  m = numel (alphabet);
  n = numel (message);
  p = primes (max (n + m - 1, 2));
  e = factorial_exponents (p, m - 1);
  for s = alphabet
    e = e + factorial_exponents (p, sum (message == s));
  end
  e = e - factorial_exponents (p, n + m - 1);
  num = big_from_primes (p, max (e, 0));
  den = big_from_primes (p, max (-e, 0));
  if (numel (big_bits (num)) <= 53 && numel (big_bits (den)) <= 53)
    f = frac (big_to_double (num), big_to_double (den));
  else
    f = struct ('big_num', num, 'big_den', den);
  end
end


% Range coding with integers ------------------------------------------------------

function [low, rng] = range_step (low, rng, s)
  global LETTERS FREQ TOTAL
  cum = cumulative (FREQ, LETTERS);
  j = find (LETTERS == s);
  % Python's // floors, and floor of a quotient of two small integers is exact.
  low = low + floor (rng * cum(j) / TOTAL);
  rng = floor (rng * FREQ(j) / TOTAL);
end


% rANS ------------------------------------------------------------------------------

function y = rans_push (x, s, M)
  % x may be a vector, which the checks below push all at once.
  global LETTERS FREQ CUM TOTAL
  if (nargin < 3)
    M = TOTAL;
  end
  j = find (LETTERS == s);
  f = FREQ(j);
  c = CUM(j);
  y = floor (x / f) * M + c + mod (x, f);
end

function [s, y] = rans_pop (x, M)
  % x may be a vector, and then s is the row of letters popped from each state.
  global LETTERS FREQ CUM TOTAL
  if (nargin < 2)
    M = TOTAL;
  end
  slot = mod (x, M);
  k = zeros (size (x));
  for j = 1:numel (LETTERS)
    k(CUM(j) <= slot & slot < CUM(j) + FREQ(j)) = j;
  end
  if (any (k(:) == 0))
    error ('slot outside every letter');
  end
  s = reshape (LETTERS(k), size (x));
  y = reshape (FREQ(k), size (x)) .* floor (x / M) + slot - reshape (CUM(k), size (x));
end

function states = rans_encode_plain (message, x0)
  % States while encoding the message backward with no renormalization: [x0, x1, ...].
  states = x0;
  for s = fliplr (message)
    states(end + 1) = rans_push (states(end), s);
  end
end

function [decoded, states] = rans_decode_plain (x, n)
  decoded = '';
  states = x;
  for i = 1:n
    [s, x] = rans_pop (x);
    decoded(end + 1) = s;
    states(end + 1) = x;
  end
end

function [x, digits, log_] = rans_encode_renorm (message, L, b)
  % rANS keeping x in [L, bL): before encoding s, emit digits while x >= (L // M) * b * f_s.
  % Each entry (s, x, emitted, new) of the Python's log is a struct with those fields.
  global LETTERS FREQ TOTAL
  if (nargin < 2)
    L = 10;
  end
  if (nargin < 3)
    b = 10;
  end
  x = L;
  digits = [];
  log_ = struct ('s', {}, 'x', {}, 'emitted', {}, 'new', {});
  for s = fliplr (message)
    emitted = [];
    while (x >= floor (L / TOTAL) * b * FREQ(LETTERS == s))
      emitted(end + 1) = mod (x, b);
      digits(end + 1) = mod (x, b);
      x = floor (x / b);
    end
    new = rans_push (x, s);
    log_(end + 1) = struct ('s', s, 'x', x, 'emitted', emitted, 'new', new);
    x = new;
  end
end

function [decoded, x] = rans_decode_renorm (x, digits, n, L, b)
  if (nargin < 4)
    L = 10;
  end
  if (nargin < 5)
    b = 10;
  end
  stack = digits;
  decoded = '';
  for i = 1:n
    [s, x] = rans_pop (x);
    decoded(end + 1) = s;
    while (x < L)
      x = x * b + stack(end);
      stack(end) = [];
    end
  end
end

% The round trips on 300-letter messages run plain rANS to states near 2^350,
% far past the 2^53 a double holds exactly. These three run the same push and
% pop on limbs, and the Python needs nothing extra because its ints have no limit.

function y = rans_push_big (x, s)
  global LETTERS FREQ CUM TOTAL
  j = find (LETTERS == s);
  [q, r] = big_divmod_small (x, FREQ(j));
  y = big_add_small (big_mul_small (q, TOTAL), CUM(j) + r);
end

function [s, y] = rans_pop_big (x)
  global LETTERS FREQ CUM TOTAL
  [q, slot] = big_divmod_small (x, TOTAL);
  j = find (CUM <= slot & slot < CUM + FREQ, 1);
  s = LETTERS(j);
  y = big_add_small (big_mul_small (q, FREQ(j)), slot - CUM(j));
end

function x = rans_encode_plain_big (message, x0)
  % The last state of rans_encode_plain, the only one the round trip sends.
  x = big_from (x0);
  for s = fliplr (message)
    x = rans_push_big (x, s);
  end
end

function decoded = rans_decode_plain_big (x, n)
  decoded = '';
  for i = 1:n
    [s, x] = rans_pop_big (x);
    decoded(end + 1) = s;
  end
end

function x = kth_state (s, k, limit)
  % The k-th number, counting from 0, whose slot x mod M belongs to letter s.
  global LETTERS FREQ CUM TOTAL
  if (nargin < 3)
    limit = 1000;
  end
  j = find (LETTERS == s);
  xs = 0:limit - 1;
  hits = xs(CUM(j) <= mod (xs, TOTAL) & mod (xs, TOTAL) < CUM(j) + FREQ(j));
  if (numel (hits) <= k)
    error ('kth_state: fewer than k + 1 states below the limit');
  end
  x = hits(k + 1);   % k counts from 0
end


% LZ77 ------------------------------------------------------------------------------

function out = lz77_encode (text)
  % Greedy triples (offset, length, next letter) with an unbounded window. Ties go to the nearest match.
  % Every offset is tried at once: alive marks the offsets whose match is still
  % growing, and a match stops one letter short of the end so that a next letter remains.
  % max returns the first of equal maxima, which is the smallest offset.
  n = numel (text);
  i = 0;   % position, counted from 0 as in the Python
  out = {};
  while (i < n)
    offs = 1:i;
    len = zeros (1, i);
    alive = true (1, i);
    l = 0;
    while (i + l < n - 1 && any (alive))
      idx = find (alive);
      hit = text(i + l - offs(idx) + 1) == text(i + l + 1);
      alive(idx(~hit)) = false;
      len(idx(hit)) = l + 1;
      l = l + 1;
    end
    best_len = 0;
    best_off = 0;
    if (i > 0)
      [top, where] = max (len);
      if (top > 0)
        best_len = top;
        best_off = offs(where);
      end
    end
    out{end + 1} = tup (best_off, best_len, text(i + best_len + 1));
    i = i + best_len + 1;
  end
end

function text = lz77_decode (triples)
  text = '';
  for k = 1:numel (triples)
    off = double (triples{k}.tuple{1});
    len = double (triples{k}.tuple{2});
    for j = 1:len
      text(end + 1) = text(end - off + 1);   % the Python's out[-off]
    end
    text(end + 1) = triples{k}.tuple{3};
  end
end


% LZ78 and LZW ----------------------------------------------------------------------

function [out, phrases] = lz78_encode (text)
  % The Python dict of phrases is a containers.Map from phrase to number.
  phrases = containers.Map ({''}, {0});
  out = {};
  w = '';
  for ch = text
    if (isKey (phrases, [w ch]))
      w = [w ch];
    else
      out{end + 1} = tup (phrases(w), ch);
      phrases([w ch]) = double (phrases.Count);
      w = '';
    end
  end
  if (~isempty (w))
    out{end + 1} = tup (phrases(w), '');
  end
end

function text = lz78_decode (pairs)
  table = {''};   % phrase number i is table{i + 1}
  out = {};
  for k = 1:numel (pairs)
    phrase = [table{double(pairs{k}.tuple{1}) + 1} pairs{k}.tuple{2}];
    table{end + 1} = phrase;
    out{end + 1} = phrase;
  end
  text = [out{:}];
end

function rate = lz78_bits (bits)
  % Bits per letter of LZ78 on a binary sequence: phrase i costs ceil(log2 i) + 1 bits.
  % The Python's trie dict, keyed by (node, bit), is the matrix child, where
  % child(node + 1, bit + 1) is the node reached and 0 means no such node yet.
  child = zeros (numel (bits) + 1, 2);
  node = 0;
  nxt = 1;
  phrases = 0;
  total = 0;
  for k = 1:numel (bits)
    bit = bits(k);
    to = child(node + 1, bit + 1);
    if (to)
      node = to;
    else
      phrases = phrases + 1;
      if (phrases > 1)
        total = total + ceil (log2 (phrases)) + 1;
      else
        total = total + 1;
      end
      child(node + 1, bit + 1) = nxt;
      nxt = nxt + 1;
      node = 0;
    end
  end
  if (node)
    phrases = phrases + 1;
    total = total + ceil (log2 (phrases));
  end
  rate = total / numel (bits);
end

function bits = binary_source (p, n, seed)
  if (nargin < 3)
    seed = 2026;
  end
  rng = PyRandom (seed);
  bits = double (rng.random (n).' < p);
end

function [out, steps] = lzw_encode (text, alphabet)
  % The Python dict of codes is a containers.Map from string to code.
  table = containers.Map (num2cell (alphabet), num2cell (0:numel (alphabet) - 1));
  out = [];
  steps = {};
  w = '';
  for ch = text
    if (isKey (table, [w ch]))
      w = [w ch];
    else
      out(end + 1) = table(w);
      table([w ch]) = double (table.Count);
      steps{end + 1} = tup (w, ch, table(w), [w ch], table([w ch]));
      w = ch;
    end
  end
  out(end + 1) = table(w);
  steps{end + 1} = tup (w, '', table(w), '', py_none ());
end

function [decoded, special, table] = lzw_decode (codes, alphabet)
  % The Python dict from code to string is the cell row table, where code i is
  % table{i + 1}. Its keys are always 0 to len - 1, so "code in table" is code < numel (table).
  table = num2cell (alphabet);
  prev = table{codes(1) + 1};
  out = {prev};
  special = [];
  for code = codes(2:end)
    if (code < numel (table))
      entry = table{code + 1};
    else
      entry = [prev prev(1)];
      special(end + 1) = code;
    end
    table{end + 1} = [prev entry(1)];
    out{end + 1} = entry;
    prev = entry;
  end
  decoded = [out{:}];
end


% Run-length coding -----------------------------------------------------------------

function runs = run_lengths (bits)
  runs = [];
  cur = '0';
  count = 0;
  for b = bits
    if (b == cur)
      count = count + 1;
    else
      runs(end + 1) = count;
      cur = b;
      count = 1;
    end
  end
  runs(end + 1) = count;
end

function bits = run_decode (runs)
  % Run i, counted from 0, is a run of 0s when i is even and of 1s when i is odd.
  bits = '';
  for i = 0:numel (runs) - 1
    bits = [bits repmat(char ('0' + mod (i, 2)), 1, runs(i + 1))];
  end
end


% Burrows-Wheeler and move-to-front -------------------------------------------------

function [last, rows] = bwt (text)
  n = numel (text);
  rows = cell (1, n);
  for i = 0:n - 1
    rows{i + 1} = [text(i + 1:end) text(1:i)];
  end
  rows = sort (rows);
  last = cellfun (@(r) r(end), rows, 'UniformOutput', false);
  last = [last{:}];
end

function out = lf_map (last)
  % Row numbers count from 0, as in the Python.
  first = sort (last);
  out = zeros (1, numel (last));
  for j = 1:numel (last)
    ch = last(j);
    start = find (first == ch, 1) - 1;
    out(j) = start + sum (last(1:j - 1) == ch);
  end
end

function [text, path] = inverse_bwt (last, end_mark)
  % end_mark is the Python's end, a word Octave reserves.
  if (nargin < 2)
    end_mark = '$';
  end
  lf = lf_map (last);
  row = 0;
  letters = '';
  path = 0;
  while (last(row + 1) ~= end_mark)
    letters(end + 1) = last(row + 1);
    row = lf(row + 1);
    path(end + 1) = row;
  end
  text = [fliplr(letters) end_mark];
end

function [out, lists] = mtf_encode (text, alphabet)
  lst = alphabet;
  out = [];
  lists = {};
  for ch = text
    lists{end + 1} = lst;
    i = find (lst == ch, 1) - 1;   % the Python's lst.index (ch), from 0
    out(end + 1) = i;
    lst = [lst(i + 1) lst([1:i, i + 2:end])];
  end
end

function out = mtf_decode (indices, alphabet)
  lst = alphabet;
  out = '';
  for i = indices
    ch = lst(i + 1);
    lst(i + 1) = [];
    out(end + 1) = ch;
    lst = [ch lst];
  end
end

function h = empirical_entropy (seq)
  % The Python iterates set (seq), whose order for letters changes from run to
  % run with string hashing. Here the distinct values are taken in sorted order,
  % and the compensated sum makes the result the same either way at the digits printed.
  n = numel (seq);
  vals = unique (seq);
  counts = zeros (1, numel (vals));
  for i = 1:numel (vals)
    counts(i) = sum (seq == vals(i));
  end
  h = entropy (counts / n);
end


% Quantization, transforms and rate-distortion ---------------------------------------

function [k, r] = quantize (x, step)
  % k is a Python int, and so is k * step when step is an int.
  k = py_round (x / step);
  r = k * step;
end

function v = phi (x)
  % Elementwise, so that a whole quantizer's cells are computed at once.
  v = zeros (size (x));
  keep = abs (x) < 40;
  v(keep) = exp (-x(keep) .* x(keep) / 2) / sqrt (2 * pi);
end

function v = Phi (x)
  v = 0.5 * (1 + erf (x / sqrt (2)));
end

function m = cell_mse (a, b, y)
  % Integral of (x - y)^2 phi(x) over [a, b] for the standard Gaussian density phi.
  afa = zeros (size (a));
  keep = abs (a) < 40;
  afa(keep) = a(keep) .* phi (a(keep));
  bfb = zeros (size (b));
  keep = abs (b) < 40;
  bfb(keep) = b(keep) .* phi (b(keep));
  m = (1 + y .* y) .* (Phi (b) - Phi (a)) - (bfb - afa) + 2 * y .* (phi (b) - phi (a));
end

function total = fixed_rate_mse (levels, step)
  i = 0:levels - 1;
  a = (i - levels / 2) * step;
  a(1) = -50.0;
  b = (i + 1 - levels / 2) * step;
  b(end) = 50.0;
  cells = cell_mse (a, b, (i - (levels - 1) / 2) * step);
  total = 0.0;
  for k = 1:levels
    total = total + cells(k);   % added one at a time, in the Python's order
  end
end

function x = golden_min (f, lo, hi, iterations)
  if (nargin < 4)
    iterations = 200;
  end
  g = (sqrt (5) - 1) / 2;
  a = lo;
  b = hi;
  for it = 1:iterations
    c = b - g * (b - a);
    d = a + g * (b - a);
    if (f (c) < f (d))
      b = d;
    else
      a = c;
    end
  end
  x = (a + b) / 2;
end

function [step, dist] = best_uniform (levels)
  step = golden_min (@(d) fixed_rate_mse (levels, d), 0.01, 3.0);
  dist = fixed_rate_mse (levels, step);
end

function [rate, mse] = entropy_coded_uniform (step)
  % Rate (index entropy) and mean squared error of a uniform quantizer on a unit Gaussian.
  K = fix (12 / step) + 2;   % int () truncates toward zero
  k = -K:K;
  a = (k - 0.5) * step;
  b = (k + 0.5) * step;
  p = Phi (b) - Phi (a);
  cells = cell_mse (a, b, k * step);
  rate = 0.0;
  mse = 0.0;
  for i = 1:numel (k)
    if (p(i) > 0)
      rate = rate - p(i) * log2 (p(i));
    end
    mse = mse + cells(i);
  end
end

function v = db (x)
  v = 10 * log10 (x);
end


% Base-delta-immediate compression ----------------------------------------------------

function ok = fits_signed (v, nbytes)
  ok = -(2 ^ (8 * nbytes - 1)) <= v & v < 2 ^ (8 * nbytes - 1);
end

function r = bdi_compress (words, word_bytes, delta_bytes, immediate)
  % One explicit base (the first word not small enough for the zero base) and signed deltas.
  % Returns a struct with the fields base, deltas, mask and size (in bytes), or []
  % where the Python returns None because some word does not fit.
  base = [];
  for w = words
    if (~(immediate && fits_signed (w, delta_bytes)))
      base = w;
      break;
    end
  end
  if (isempty (base))
    base = words(1);
  end
  deltas = [];
  mask = [];
  r = [];
  for w = words
    if (immediate && fits_signed (w, delta_bytes))
      deltas(end + 1) = w;
      mask(end + 1) = 0;
    elseif (fits_signed (w - base, delta_bytes))
      deltas(end + 1) = w - base;
      mask(end + 1) = 1;
    else
      return;
    end
  end
  nbytes = word_bytes + delta_bytes * numel (words);
  if (immediate)
    nbytes = nbytes + ceil (numel (words) / 8);
  end
  r = struct ('base', base, 'deltas', deltas, 'mask', mask, 'size', nbytes);
end

function words = bdi_decompress (base, deltas, mask)
  words = deltas + base * (mask == 1);
end


% Huffman's code on blocks, copied from python/ch05_source_coding_and_prefix_codes.py --
% The Python imports block_rate from chapter 5. One Octave script cannot call
% another script's functions, so block_rate and the three functions it calls are
% copied here, for a list of floats, the only kind block_rate passes.

function [lengths, merges] = huffman (probs)
  % Huffman's algorithm by the sorted-list method, with each merged entry placed
  % above every entry of equal probability. Octave's sort keeps equal entries in
  % their first order, as Python's sorted does.
  v = probs(:);
  n = numel (v);
  [p, order] = sort (v, 'descend');
  members = num2cell (order(:).');
  lengths = zeros (1, n);
  merges = zeros (n - 1, 1);
  for step = 1:n - 1
    weight = p(end - 1) + p(end);
    group = [members{end - 1}, members{end}];
    p(end - 1:end) = [];
    members(end - 1:end) = [];
    lengths(group) = lengths(group) + 1;
    merges(step) = weight;
    spot = find (p <= weight, 1);
    if (isempty (spot))
      spot = numel (p) + 1;
    end
    p = [p(1:spot - 1); weight; p(spot:end)];
    members = [members(1:spot - 1), {group}, members(spot:end)];
  end
end

function L = huffman_length (probs)
  % Average Huffman codeword length, as the sum of the merged probabilities.
  [~, merges] = huffman (probs);
  L = py_sum (merges);
end

function probs = block_probs (p_one, n)
  % Probabilities of the 2^n blocks of n independent letters, each 1 with probability p_one.
  % itertools.product lists the blocks as the binary numbers 0 to 2^n - 1, so
  % block y holds as many 1s as y has 1 bits. A scalar raised to a vector of
  % powers calls pow for each one, as Python's float ** int does.
  ones_in = sum (dec2bin ((0:2 ^ n - 1).', n) == '1', 2);
  probs = (p_one .^ ones_in) .* ((1 - p_one) .^ (n - ones_in));
end

function r = block_rate (p_one, n)
  % Huffman bits per letter when blocks of n letters are coded.
  r = huffman_length (block_probs (p_one, n)) / n;
end


% The chapter's numbers, section by section ----------------------------------------------

function b = bits_of (prob)
  % log2(1/prob) for an exact fraction too small for a float.
  if (isfield (prob, 'big_num'))
    b = big_log2 (prob.big_den) - big_log2 (prob.big_num);
  else
    b = log2 (prob.den) - log2 (prob.num);
  end
end

function check_intervals ()
  global LETTERS PROBS
  printf ('Coding with Intervals\n');
  check ('entropy of the source a 0.7, b 0.2, c 0.1', entropy (PROBS), '1.1568');
  huff = {'0', '10', '11'};   % the codewords of a, b and c
  terms = zeros (1, 3);
  for k = 1:3
    terms(k) = to_float (PROBS(k)) * numel (huff{k});
  end
  check ('average length of Huffman''s code', py_sum (terms), '1.3');
  msg = 'aabaaa';
  prob = frac_mul (frac_pow (PROBS(LETTERS == 'a'), 5), PROBS(LETTERS == 'b'));
  same ('probability of aabaaa', prob, frac_dec ('0.033614'));
  % math.log2 of a Fraction converts it to a float first, and 1 / prob as a
  % float is den / num, rounded once.
  check ('information in aabaaa', log2 (prob.den / prob.num), '4.8948');
  spent = 0;
  for s = msg
    spent = spent + numel (huff{LETTERS == s});
  end
  same ('bits Huffman''s code spends on aabaaa', int64 (spent), int64 (7));
  check ('information in one a', log2 (10 / 7), '0.5146');
  steps = arith_intervals (msg, PROBS, LETTERS);
  table = {'0', '0.7', '0.7'; '0', '0.49', '0.49'; '0.343', '0.441', '0.098'; ...
           '0.343', '0.4116', '0.0686'; '0.343', '0.39102', '0.04802'; '0.343', '0.376614', '0.033614'};
  for k = 1:6
    low = steps(k + 1).low;
    width = steps(k + 1).width;
    same (sprintf ('interval and width after step %d', k), tup (low, frac_add (low, width), width), ...
          tup (frac_dec (table{k, 1}), frac_dec (table{k, 2}), frac_dec (table{k, 3})));
  end
  same ('b keeps 70 to 90 percent of [0, 0.49)', tup (frac_mul (frac_dec ('0.49'), frac_dec ('0.7')), frac_mul (frac_dec ('0.49'), frac_dec ('0.9'))), ...
        tup (frac_dec ('0.343'), frac_dec ('0.441')));

  low = steps(end).low;
  width = steps(end).width;
  tag = arith_codeword (low, width);
  same ('tag of aabaaa', tag, '010110');
  same ('tag length l', int64 (numel (tag)), int64 (6));
  check ('u times 2^6', to_float (low) * 64, '21.952');
  same ('tag integer v', int64 (bin2dec (tag)), int64 (22));
  v = bin2dec (tag);
  same ('dyadic interval of the tag', tup (frac (v, 64), frac (v + 1, 64)), ...
        tup (frac_dec ('0.34375'), frac_dec ('0.359375')));
  require ('the dyadic interval lies inside [0.343, 0.376614)', frac_cmp (low, frac (v, 64)) <= 0 && frac_cmp (frac (v + 1, 64), frac_add (low, width)) <= 0);
  check ('tag bits beyond the information', 6 - log2 (prob.den / prob.num), '1.1052');
  [decoded, values] = arith_decode (tag, 6, PROBS, LETTERS);
  same ('decoding the tag', decoded, msg);
  stated = {'0.34375', '0.491071', '0.701531', '0.007653', '0.010933', '0.015618'};
  for k = 0:5
    check (sprintf ('decoder value z at step %d', k + 1), values{k + 1}, stated{k + 1});
  end

  for n = 1:5
    tags = cell (1, 3 ^ n);
    total = frac (0);
    inside = true;
    decodes = true;
    % itertools.product (LETTERS, repeat=n) runs through the base-3 numerals in order.
    combos = dec2base (0:3 ^ n - 1, 3, n);
    for row = 1:3 ^ n
      m = LETTERS(combos(row, :) - '0' + 1);
      steps = arith_intervals (m, PROBS, LETTERS);
      lo = steps(end).low;
      w = steps(end).width;
      t = arith_codeword (lo, w);
      vv = bin2dec (t);
      ll = numel (t);
      inside = inside && frac_cmp (lo, frac (vv, 2 ^ ll)) <= 0 && frac_cmp (frac (vv + 1, 2 ^ ll), frac_add (lo, w)) <= 0;
      decodes = decodes && strcmp (arith_decode (t, n, PROBS, LETTERS), m);
      tags{row} = t;
      total = frac_add (total, w);
    end
    tags = sort (tags);
    prefix_free = true;
    for i = 1:numel (tags) - 1
      % strncmp is false when tags{i + 1} is shorter than tags{i}, as startswith is.
      prefix_free = prefix_free && ~strncmp (tags{i + 1}, tags{i}, numel (tags{i}));
    end
    require (sprintf ('n = %d: intervals tile [0, 1), tags inside, prefix-free, decodable', n), ...
             frac_eq (total, frac (1)) && inside && prefix_free && decodes);
  end
end

function check_rates (block_rate)
  printf ('Why Arithmetic Coding Reaches the Entropy\n');
  hb = binary_entropy (0.1);
  check ('h(0.1)', hb, '0.4690');
  table = {1, '2.3', '2.4690'; 2, '1.3', '1.4690'; 4, '0.8009', '0.9690'; ...
           8, '0.6751', '0.7190'; 100, '0.4841', '0.4890'; 1000, '0.4705', '0.4710'};
  for i = 1:rows (table)
    n = table{i, 1};
    r = arith_expected_rate (0.1, n);
    check (sprintf ('arithmetic coding, n = %d', n), r, table{i, 2});
    check (sprintf ('bound h(0.1) + 2/n, n = %d', n), hb + 2 / n, table{i, 3});
    require (sprintf ('h(0.1) <= rate < bound, n = %d', n), hb <= r && r < hb + 2 / n);
  end
  blocks = {1, '1'; 2, '0.6450'; 4, '0.4926'; 8, '0.4758'};
  for i = 1:rows (blocks)
    check (sprintf ('Huffman''s code on blocks of %d', blocks{i, 1}), block_rate (0.1, blocks{i, 1}), blocks{i, 2});
  end
  same ('tag length of the letter 0 alone', int64 (numel (arith_codeword (frac (0), frac (9, 10)))), int64 (2));
  same ('tag length of the letter 1 alone', int64 (numel (arith_codeword (frac (9, 10), frac (1, 10)))), int64 (5));
  check ('0.9 x 2 + 0.1 x 5', 0.9 * 2 + 0.1 * 5, '2.3');
  r6 = arith_expected_rate (0.1, 6);
  r7 = arith_expected_rate (0.1, 7);
  check ('arithmetic coding, n = 6', r6, '0.7114');
  check ('arithmetic coding, n = 7', r7, '0.7286');
  require ('the rate rises from n = 6 to n = 7', r7 > r6);
  check ('distance from the entropy at n = 1000', arith_expected_rate (0.1, 1000) - hb, '0.0015');
  h01 = binary_entropy (0.01);
  check ('h(0.01)', h01, '0.0808');
  check ('1 bit per letter over h(0.01)', 1 / h01, '12.4');
  check ('arithmetic coding at p = 0.01, n = 1000', arith_expected_rate (0.01, 1000), '0.0823');
  check ('bound h(0.01) + 2/1000', h01 + 2 / 1000, '0.0828');
end

function check_range_coding ()
  global LETTERS FREQ CUM TOTAL PROBS
  printf ('Range Coding and Adaptive Models\n');
  pa = PROBS(LETTERS == 'a');
  pb = PROBS(LETTERS == 'b');
  [bits, lo, hi] = settle_bits (frac (0), frac_dec ('0.49'));
  same ('bits settled after aa', tup (bits, lo, hi), tup ('0', frac (0), frac_dec ('0.98')));
  w = frac_sub (hi, lo);
  % Both new ends are computed from the old lo, as the Python's tuple assignment does.
  [lo, hi] = deal (frac_add (lo, frac_mul (w, pa)), frac_add (lo, frac_mul (w, frac_add (pa, pb))));
  same ('interval after the b', tup (lo, hi), tup (frac_dec ('0.686'), frac_dec ('0.882')));
  [bits, lo, hi] = settle_bits (lo, hi);
  same ('bits settled after the b', tup (bits, lo, hi), tup ('1', frac_dec ('0.372'), frac_dec ('0.764')));
  for k = 0:2
    hi = frac_add (lo, frac_mul (frac_sub (hi, lo), pa));
    [bits, lo, hi] = settle_bits (lo, hi);
    same (sprintf ('bits settled by a number %d', k + 4), bits, '');
  end
  same ('interval after aabaaa', tup (lo, hi), tup (frac_dec ('0.372'), frac_dec ('0.506456')));
  require ('it straddles 1/2 inside [1/4, 3/4)', frac_cmp (frac (1, 4), lo) <= 0 && frac_cmp (lo, frac (1, 2)) < 0 ...
           && frac_cmp (frac (1, 2), hi) < 0 && frac_cmp (hi, frac (3, 4)) <= 0);
  require ('the settled bits 01 begin the tag 010110', strncmp ('010110', '01', 2));

  msg = 'aabaaa';
  low = 0;
  rng = 1000;
  ranges = 1000;
  for s = msg
    [low, rng] = range_step (low, rng, s);
    ranges(end + 1) = rng;
  end
  same ('ranges for aabaaa', int_list (ranges), int_list ([1000, 700, 490, 98, 68, 47, 32]));
  same ('final interval in thousandths', tup (low, low + rng), tup (343, 375));
  check ('cost of truncation', log2 (33.614 / 32), '0.0710');
  low = 0;
  rng = 1000;
  for s = 'aaba'
    [low, rng] = range_step (low, rng, s);
  end
  same ('interval after aaba', tup (low, low + rng), tup (343, 411));
  first = sprintf ('%d', low);
  final = sprintf ('%d', low + rng - 1);
  require ('[343, 411) holds numbers beginning 3 and 4', first(1) == '3' && final(1) == '4');
  [low, rng] = range_step (low, rng, 'a');
  same ('interval after aabaa', tup (low, low + rng), tup (343, 390));
  first = sprintf ('%d', low);
  final = sprintf ('%d', low + rng - 1);
  require ('every number in [343, 390) begins with 3', first(1) == '3' && final(1) == '3');
  r = (1:4999).';
  no_overlap = true;
  for j = 1:numel (LETTERS)
    no_overlap = no_overlap && all (floor (r * CUM(j) / TOTAL) + floor (r * FREQ(j) / TOTAL) <= floor (r * (CUM(j) + FREQ(j)) / TOTAL));
  end
  require ('truncated slices never overlap, R = 1 to 4999', no_overlap);

  same ('Laplace probability of aab', laplace_probability ('aab', LETTERS), frac (1, 30));
  pl = laplace_probability (msg, LETTERS);
  same ('Laplace probability of aabaaa', pl, frac (1, 168));
  same ('the closed form agrees', laplace_closed_form (msg, LETTERS), pl);
  check ('bits under Laplace''s rule', log2 (168), '7.3923');
  same ('tag length under Laplace''s rule', int64 (numel (arith_codeword (frac (0), pl))), int64 (9));
  holds = true;
  for n = 1:6
    combos = dec2base (0:3 ^ n - 1, 3, n);   % itertools.product order
    for row = 1:3 ^ n
      m = LETTERS(combos(row, :) - '0' + 1);
      p_true = frac (1);
      for s = m
        p_true = frac_mul (p_true, PROBS(LETTERS == s));
      end
      holds = holds && frac_cmp (frac_mul (laplace_probability (m, LETTERS), frac (nchoosek (n + 2, 2))), p_true) >= 0;
      holds = holds && frac_eq (laplace_closed_form (m, LETTERS), laplace_probability (m, LETTERS));
    end
  end
  require ('price of learning holds for every message of 1 to 6 letters', holds);
  check ('price bound at n = 100', log2 (nchoosek (102, 2)), '12.331');
  check ('price bound at n = 1000', log2 (nchoosek (1002, 2)), '18.936');
end

function check_rans ()
  global LETTERS FREQ TOTAL
  printf ('Asymmetric Numeral Systems\n');
  msg = 'aabaaa';
  check ('log2 10', log2 (10), '3.32');
  same ('a-integer number 10', int64 (kth_state ('a', 10)), int64 (13));
  same ('pushing a onto 10', int64 (rans_push (10, 'a')), int64 (13));
  states = rans_encode_plain (msg, 10);
  same ('states while encoding aabaaa backward', int_list (states), int_list ([10, 13, 16, 22, 117, 165, 234]));
  rev = fliplr (msg);
  qr = cell (1, numel (rev));
  for i = 1:numel (rev)
    f = FREQ(LETTERS == rev(i));
    qr{i} = tup (floor (states(i) / f), mod (states(i), f));   % divmod
  end
  same ('(q, r) at each push', qr, {tup(1, 3), tup(1, 6), tup(2, 2), tup(11, 0), tup(16, 5), tup(23, 4)});
  stated = {'1.3', '1.2308', '1.375', '5.3182', '1.4103', '1.4182'};
  for i = 1:6
    check (sprintf ('growth factor of push %d', i), states(i + 1) / states(i), stated{i});
  end
  check ('ideal factor for a', 10 / 7, '1.4286');
  check ('ideal factor for b', 10 / 2, '5');
  same ('b-integer number 22', int64 (kth_state ('b', 22)), int64 (117));
  check ('total growth log2(234/10)', log2 (23.4), '4.5484');
  same ('bits to write the state 234', int64 (bit_length (234)), int64 (8));
  same ('slot of 234', int64 (mod (234, TOTAL)), int64 (4));
  [s, x] = rans_pop (234);
  same ('popping 234', tup (s, x), tup ('a', 165));
  same ('7 x 23 + 4', int64 (7 * 23 + 4), int64 (165));
  [decoded, back] = rans_decode_plain (234, 6);
  same ('plain decoding from 234', tup (decoded, int_list (back)), tup (msg, int_list ([234, 165, 117, 22, 16, 13, 10])));
  % The two checks over x = 0 to 2999 push and pop every x of a letter at once.
  inverse = true;
  for s = LETTERS
    xs = 0:2999;
    [t, y] = rans_pop (rans_push (xs, s));
    inverse = inverse && all (t == s) && all (y == xs);
  end
  require ('D(C(s, x)) = x for x = 0 to 2999', inverse);
  growth = true;
  for s = LETTERS
    f = FREQ(LETTERS == s);
    xs = 1:2999;
    growth = growth && all (abs (rans_push (xs, s) * f ./ (TOTAL * xs) - 1) < f ./ xs);
  end
  require ('|delta| < f_s / x for x = 1 to 2999', growth);

  [x, digits, steps] = rans_encode_renorm (msg);
  same ('renormalized states after each push', int_list ([steps.new]), int_list ([13, 16, 22, 17, 23, 32]));
  emitted = cellfun (@int_list, {steps.emitted}, 'UniformOutput', false);
  same ('digits emitted before each push', emitted, {int_list([]), int_list([]), int_list([]), int_list(2), int_list([]), int_list([])});
  same ('state sent and digits sent', tup (x, int_list (digits)), tup (32, int_list (2)));
  same ('C(b, 2)', int64 (rans_push (2, 'b')), int64 (17));
  [decoded, x_end] = rans_decode_renorm (x, digits, 6);
  same ('renormalized decoding', tup (decoded, x_end), tup (msg, 10));
  in_range = true;
  for s = LETTERS
    f = FREQ(LETTERS == s);
    for x0 = 10:99
      y = x0;
      emitted = 0;
      while (y >= 10 * f)
        y = floor (y / 10);
        emitted = emitted + 1;
      end
      [~, popped] = rans_pop (rans_push (y, s));
      in_range = in_range && emitted <= 1 && 10 <= rans_push (y, s) && rans_push (y, s) < 100 && (popped < 10) == (emitted == 1);
    end
  end
  require ('renormalization keeps every state in [10, 100)', in_range);
  rng = PyRandom (6);
  for trial = 0:19
    m = rng.choices (LETTERS, 300, [7, 2, 1]);
    [x, digits] = rans_encode_renorm (m);
    [decoded, x_end] = rans_decode_renorm (x, digits, numel (m));
    ok = strcmp (decoded, m) && x_end == 10;
    ok = ok && strcmp (rans_decode_plain_big (rans_encode_plain_big (m, 10), numel (m)), m);
    if (~ok)
      require (sprintf ('rANS round trip %d', trial), false);
    end
  end
  require ('rANS round trips, plain and renormalized, 20 messages', true);

  same ('states while encoding ca', int_list (rans_encode_plain ('ca', 10)), int_list ([10, 13, 139]));
  same ('10 x 13 + 9 + 0', int64 (10 * 13 + 9 + 0), int64 (139));
  check ('growth for ca', log2 (13.9), '3.797');
  check ('information in ca', log2 (1 / 0.07), '3.8365');
  [s, x] = rans_pop (139);
  same ('popping 139', tup (s, x), tup ('c', 13));
  [s, x] = rans_pop (13);
  same ('popping 13', tup (s, x), tup ('a', 10));
end

function check_dictionary ()
  printf ('Dictionary Coding with LZ77\n');
  text = 'abcbcbcabcd';
  triples = lz77_encode (text);
  same ('LZ77 parse of abcbcbcabcd', triples, {tup(0, 0, 'a'), tup(0, 0, 'b'), tup(0, 0, 'c'), tup(2, 4, 'a'), tup(3, 2, 'd')});
  same ('LZ77 decoding', lz77_decode (triples), text);
  offsets = [];
  for o = 1:8
    if (strcmp (text(9 - o:10 - o), 'bc'))   % the Python's text[8 - o:10 - o]
      offsets(end + 1) = o;
    end
  end
  same ('offsets back to bc from position 9', int_list (offsets), int_list ([3, 5, 7]));
  require ('the copy at position 4 is longer than its offset', triples{4}.tuple{2} > triples{4}.tuple{1});
  same ('LZ77 parse of aaaaaaab', lz77_encode ('aaaaaaab'), {tup(0, 0, 'a'), tup(1, 6, 'b')});
  same ('its decoding', lz77_decode ({tup(0, 0, 'a'), tup(1, 6, 'b')}), 'aaaaaaab');

  printf ('LZ78 and LZW\n');
  pairs = lz78_encode ('aababcabcd');
  same ('LZ78 pairs for aababcabcd', pairs, {tup(0, 'a'), tup(1, 'b'), tup(2, 'c'), tup(3, 'd')});
  same ('LZ78 decoding', lz78_decode (pairs), 'aababcabcd');
  same ('B(4) on the phrases 1, 0, 11, 00', lz78_bits ([1, 0, 1, 1, 0, 0]) * 6, int64 (1 + 2 + 3 + 3));
  hb = binary_entropy (0.1);
  lengths = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000];
  bits = binary_source (0.1, lengths(end));
  rates = zeros (size (lengths));
  for i = 1:numel (lengths)
    rates(i) = lz78_bits (bits(1:lengths(i)));
  end
  check ('LZ78 after 100 letters', rates(lengths == 100), '0.7600');
  check ('LZ78 after 10,000 letters', rates(lengths == 10000), '0.6237');
  check ('LZ78 after 1,000,000 letters', rates(lengths == 1000000), '0.5586');
  check ('percent above the entropy at 1,000,000', 100 * (rates(lengths == 1000000) / hb - 1), '19');
  scaled = (rates - hb) .* log2 (lengths);
  check ('smallest excess times log2 n', min (scaled), '1.68');
  check ('largest excess times log2 n', max (scaled), '2.34');

  [codes, steps] = lzw_encode ('abababa', 'ab');
  same ('LZW numbers for abababa', int_list (codes), int_list ([0, 1, 2, 4]));
  matches = cellfun (@(st) tup (st.tuple{1}, st.tuple{4}, st.tuple{5}), steps, 'UniformOutput', false);
  same ('LZW matches and new entries', matches, ...
        {tup('a', 'ab', 2), tup('b', 'ba', 3), tup('ab', 'aba', 4), tup('aba', '', py_none ())});
  [decoded, special] = lzw_decode (codes, 'ab');
  same ('LZW decoding of 0 1 2 4', tup (decoded, int_list (special)), tup ('abababa', int_list (4)));
  codes = lzw_encode ('aaaaaaa', 'ab');
  same ('LZW numbers for aaaaaaa', int_list (codes), int_list ([0, 2, 3, 0]));
  [decoded, special] = lzw_decode (codes, 'ab');
  same ('its decoding and special cases', tup (decoded, int_list (special)), tup ('aaaaaaa', int_list ([2, 3])));
  codes = lzw_encode (repmat ('a', 1, 1000), 'ab');
  same ('numbers for a run of 1000 a''s', int64 (numel (codes)), int64 (45));
  same ('44 x 45 / 2', int64 (floor (44 * 45 / 2)), int64 (990));
  same ('the 45th phrase', int64 (codes(45)), int64 (10));
  [decoded, special, table] = lzw_decode (codes, 'ab');
  same ('decoding the run', decoded, repmat ('a', 1, 1000));
  same ('special cases are the 2nd to the 44th numbers', int_list (special), int_list (codes(2:44)));   % codes[1:44]
  same ('entries in the dictionary at the end', int64 (numel (table)), int64 (46));
  require ('every number fits in 6 bits', max (codes) < 2 ^ 6 && numel (table) <= 2 ^ 6);
  same ('bits for the run', int64 (45 * 6), int64 (270));
  same ('bits to send the length 1000', int64 (bit_length (1000)), int64 (10));
  [decoded, special, table] = lzw_decode ([1, 0, 2, 4, 1], 'ab');
  same ('decoding 1 0 2 4 1', tup (decoded, int_list (special)), tup ('babababb', int_list (4)));
  same ('entries 2, 3 and 4', tup (table{3}, table{4}, table{5}), tup ('ba', 'ab', 'bab'));   % table[2], table[3], table[4]
  same ('encoding babababb', int_list (lzw_encode ('babababb', 'ab')), int_list ([1, 0, 2, 4, 1]));

  rng = PyRandom (7);
  for trial = 0:29
    t = rng.choices ('abc', 400, [6, 3, 1]);
    ok = strcmp (lz77_decode (lz77_encode (t)), t) && strcmp (lz78_decode (lz78_encode (t)), t);
    ok = ok && strcmp (lzw_decode (lzw_encode (t, 'abc'), 'abc'), t);
    if (~ok)
      require (sprintf ('dictionary round trip %d', trial), false);
    end
  end
  require ('LZ77, LZ78 and LZW round trips, 30 texts', true);
end

function check_runs_bwt_mtf ()
  printf ('Run-Length Coding\n');
  row = [repmat('0', 1, 12) repmat('1', 1, 3) repmat('0', 1, 9) repmat('1', 1, 2)];
  same ('pixels in the row', int64 (numel (row)), int64 (26));
  same ('runs of the row', int_list (run_lengths (row)), int_list ([12, 3, 9, 2]));
  require ('each run fits a 4-bit field', all (run_lengths (row) < 16));
  same ('bits for four 4-bit runs', int64 (4 * numel (run_lengths (row))), int64 (16));
  same ('run decoding', run_decode ([12, 3, 9, 2]), row);
  rng = PyRandom (8);
  ok = true;
  for trial = 1:30
    t = blanks (200);
    for j = 1:200
      t(j) = rng.choice ('01');   % one draw at a time, as each choice may draw more than once
    end
    ok = ok && strcmp (run_decode (run_lengths (t)), t);
  end
  require ('run-length round trips, 30 rows', ok);
  p = 0.1;
  hb = binary_entropy (0.1);
  r = 1:2999;
  gap_entropy = -py_sum ((1 - p) .^ (r - 1) * p .* ((r - 1) * log2 (1 - p) + log2 (p)));
  check ('entropy of a gap, summed term by term', gap_entropy, '4.6900');
  check ('h(p) / p', hb / p, '4.6900');
  check ('mean gap E[R]', py_sum (r .* (1 - p) .^ (r - 1) * p), '10');
  check ('h(p) + p', hb + p, '0.5690');

  printf ('The Burrows-Wheeler Transform\n');
  [last, rows] = bwt ('couscous$');
  same ('sorted rotations', rows, {'$couscous', 'cous$cous', 'couscous$', 'ous$cousc', 'ouscous$c', ...
                                   's$couscou', 'scous$cou', 'us$cousco', 'uscous$co'});
  same ('transform of couscous', last, 'ss$ccuuoo');
  firsts = cellfun (@(r) r(1), rows, 'UniformOutput', false);
  same ('first column', [firsts{:}], '$ccoossuu');
  same ('last-to-first map', int_list (lf_map (last)), int_list ([5, 6, 0, 1, 2, 7, 8, 3, 4]));
  [text, path] = inverse_bwt (last);
  same ('rows visited from row 0', int_list (path), int_list ([0, 5, 7, 3, 1, 6, 8, 4, 2]));
  same ('letters collected', last(path(1:end - 1) + 1), 'suocsuoc');
  same ('inverse transform', text, 'couscous$');
  first = sort ('ii$bbrree');
  same ('first column of ii$bbrree', first, '$bbeeiirr');
  same ('C for $, b, e, i, r', int_list (arrayfun (@(ch) find (first == ch, 1) - 1, '$beir')), int_list ([0, 1, 3, 5, 7]));
  same ('its last-to-first map', int_list (lf_map ('ii$bbrree')), int_list ([5, 6, 0, 1, 2, 7, 8, 3, 4]));
  [text, path] = inverse_bwt ('ii$bbrree');
  word = 'ii$bbrree';
  same ('letters collected', word(path(1:end - 1) + 1), 'irebireb');
  same ('inverse of ii$bbrree', text, 'beriberi$');
  same ('forward transform of beriberi', bwt ('beriberi$'), 'ii$bbrree');
  rng = PyRandom (9);
  ok = true;
  for trial = 1:30
    t = rng.choices ('abcd', 60, []);
    ok = ok && strcmp (inverse_bwt (bwt ([t '$'])), [t '$']);
  end
  require ('Burrows-Wheeler round trips, 30 blocks', ok);

  printf ('Move-to-Front\n');
  require ('the output has the letter counts of couscous', isequal (sort (last), sort ('couscous$')));
  [idx, lists] = mtf_encode (last, '$cosu');
  same ('move-to-front output', int_list (idx), int_list ([3, 0, 1, 2, 0, 4, 0, 4, 0]));
  same ('lists before each step', lists, {'$cosu', 's$cou', 's$cou', '$scou', 'c$sou', 'c$sou', ...
                                          'uc$so', 'uc$so', 'ouc$s'});
  same ('move-to-front decoding', mtf_decode (idx, '$cosu'), last);
  check ('empirical entropy of the transform output', empirical_entropy (last), '2.2810');
  check ('empirical entropy of the move-to-front output', empirical_entropy (idx), '2.0588');
  same ('decoding 1 0 2 0 0 2 from abc', mtf_decode ([1, 0, 2, 0, 0, 2], 'abc'), 'bbccca');
  same ('encoding bbccca', int_list (mtf_encode ('bbccca', 'abc')), int_list ([1, 0, 2, 0, 0, 2]));
  ok = true;
  for trial = 1:30
    t = rng.choices ('abcd', 60, []);   % the same generator, continued from the blocks above
    ok = ok && strcmp (mtf_decode (mtf_encode (t, 'abcd'), 'abcd'), t);
  end
  require ('move-to-front round trips, 30 texts', ok);
end

function check_lossy ()
  printf ('Quantization\n');
  xs = [0.37, -1.12, 2.05, 0.61, -0.28];
  ks = zeros (size (xs));
  rs = zeros (size (xs));
  for i = 1:numel (xs)
    [ks(i), rs(i)] = quantize (xs(i), 0.5);
  end
  same ('indices', int_list (ks), int_list ([1, -2, 4, 1, -1]));
  same ('reconstructions', num2cell (rs), {0.5, -1.0, 2.0, 0.5, -0.5});
  sq_err = zeros (size (xs));
  for i = 1:numel (xs)
    sq_err(i) = (xs(i) - rs(i)) ^ 2;
  end
  check ('mean squared error', py_sum (sq_err) / numel (xs), '0.01886');
  check ('Delta^2 / 12 at Delta = 0.5', 0.5 ^ 2 / 12, '0.020833');
  check ('decibels per bit', db (4), '6.0206');
  check ('signal-to-noise ratio at 8 bits', db (2 ^ 16), '48.16');
  check ('signal-to-noise ratio at 16 bits', db (2 ^ 32), '96.33');

  printf ('Transform Coding\n');
  x1 = 100;
  x2 = 104;
  s = (x1 + x2) / sqrt (2);
  d = (x1 - x2) / sqrt (2);
  check ('s', s, '144.2498');
  check ('d', d, '-2.8284');
  same ('energy x1^2 + x2^2', int64 (x1 * x1 + x2 * x2), int64 (20816));
  check ('energy s^2 + d^2', s * s + d * d, '20816');
  check ('fraction of the energy in s', s * s / (x1 * x1 + x2 * x2), '0.99962');
  [~, sq] = quantize (s, 4);
  [~, dq] = quantize (d, 8);
  same ('s and d quantized', tup (sq, dq), tup (144, 0));
  y1 = (sq + dq) / sqrt (2);
  y2 = (sq - dq) / sqrt (2);
  check ('reconstructed pixels', y1, '101.8234');
  require ('both pixels come back equal', y1 == y2);
  check ('squared error over the pixels', (x1 - y1) ^ 2 + (x2 - y2) ^ 2, '8.0624');
  check ('squared error over the coefficients', (s - sq) ^ 2 + (d - dq) ^ 2, '8.0624');

  printf ('The Rate-Distortion Trade\n');
  y = sqrt (2 / pi);
  check ('best one-bit level E|X|', y, '0.7979');
  check ('its distortion 1 - 2/pi', 1 - 2 / pi, '0.36338');
  check ('the same, integrated over the density', fixed_rate_mse (2, 2 * y), '0.36338');
  check ('in decibels', db (1 - 2 / pi), '-4.40');
  check ('the bound at one bit', db (2 ^ -2), '-6.02');
  table = {1, '1.5958', '0.36338', '-4.40', '-6.02', '1.62'; 2, '0.9957', '0.11885', '-9.25', '-12.04', '2.79'; ...
           3, '0.5860', '0.03744', '-14.27', '-18.06', '3.80'; 4, '0.3352', '0.01154', '-19.38', '-24.08', '4.71'; ...
           5, '0.1881', '0.00350', '-24.57', '-30.10', '5.54'};
  for i = 1:rows (table)
    rate = table{i, 1};
    [step, dist] = best_uniform (2 ^ rate);
    check (sprintf ('best step, %d levels', 2 ^ rate), step, table{i, 2});
    check (sprintf ('distortion, %d levels', 2 ^ rate), dist, table{i, 3});
    check (sprintf ('distortion in dB, %d levels', 2 ^ rate), db (dist), table{i, 4});
    check (sprintf ('bound in dB at %d bits', rate), db (2 ^ (-2 * rate)), table{i, 5});
    check (sprintf ('gap in dB at %d bits', rate), db (dist) - db (2 ^ (-2 * rate)), table{i, 6});
  end
  check ('differential entropy of a unit Gaussian', 0.5 * log2 (2 * pi * exp (1)), '2.0471');
  check ('limiting gap in bits', 0.5 * log2 (pi * exp (1) / 6), '0.2546');
  check ('limiting distortion ratio pi e / 6', pi * exp (1) / 6, '1.4233');
  check ('limiting gap in dB', db (pi * exp (1) / 6), '1.5329');
  check ('limiting gap in dB, as the caption rounds it', db (pi * exp (1) / 6), '1.53');
  % The first step is the Python int 1, which prints as 1, and the others are floats.
  steps = {int64(1), '1.88'; 0.5, '1.62'; 0.25, '1.56'};
  for i = 1:rows (steps)
    [rate, mse] = entropy_coded_uniform (double (steps{i, 1}));
    check (sprintf ('entropy-coded gap in dB at step %s', py_str (steps{i, 1})), db (mse) - db (2 ^ (-2 * rate)), steps{i, 2});
  end
  check ('index entropy at step 0.25', entropy_coded_uniform (0.25), '4.0508');
  check ('R(0.25) for variance 4', 0.5 * log2 (4 / 0.25), '2');
  check ('D(3) for variance 4', 4 * 2 ^ -6, '0.0625');
  check ('a quarter of D(2)', (4 * 2 ^ -6) / (4 * 2 ^ -4), '0.25');
end

function check_memory ()
  printf ('Compression Inside Memory Systems\n');
  % A hexadecimal literal is an integer type in Octave, which saturates rather
  % than wraps, so each one is made a double, exact for these values below 2^53.
  ptrs = double ([0x00403A10, 0x00403A18, 0x00403A20, 0x00403A28, 0x00403A40, 0x00403A48, 0x00403A60, 0x00403A70]);
  require ('only the last byte varies', all (floor (ptrs / 256) == floor (ptrs(1) / 256)));   % w >> 8
  c = bdi_compress (ptrs, 4, 1, false);
  same ('base', sprintf ('%08X', c.base), '00403A10');
  same ('deltas', int_list (c.deltas), int_list ([0, 8, 16, 24, 48, 56, 80, 96]));
  same ('deltas in hexadecimal', arrayfun (@(x) sprintf ('%02X', x), c.deltas, 'UniformOutput', false), ...
        {'00', '08', '10', '18', '30', '38', '50', '60'});
  same ('compressed bytes', int64 (c.size), int64 (12));
  check ('compression ratio', 32 / c.size, '2.67');
  same ('decompression', int_list (bdi_decompress (c.base, c.deltas, c.mask)), int_list (ptrs));
  flipped = bitxor (c.base, 2 ^ 4);   % base ^ (1 << 4)
  same ('base after bit 4 flips', sprintf ('%08X', flipped), '00403A00');
  same ('error in each word', int_list (bdi_decompress (flipped, c.deltas, c.mask) - ptrs), int_list (repmat (-16, 1, 8)));
  same ('stored bits', int64 (8 * c.size), int64 (96));
  same ('bits in the line', int64 (8 * 32), int64 (256));

  mixed = double ([0, 0x00403A10, 3, 0x00403A28, 1, 0x00403A18, 0xC, 0x00403A40]);
  require ('no 1-byte or 2-byte delta reaches every word from one base', ...
           isempty (bdi_compress (mixed, 4, 1, false)) && isempty (bdi_compress (mixed, 4, 2, false)));
  c = bdi_compress (mixed, 4, 1, true);
  same ('deltas with zero as a second base', int_list (c.deltas), int_list ([0, 0, 3, 24, 1, 8, 12, 48]));
  same ('mask', sprintf ('%d', c.mask), '01010101');
  same ('compressed bytes', int64 (c.size), int64 (4 + 8 + 1));
  check ('compression ratio', 32 / c.size, '2.46');
  same ('decompression', int_list (bdi_decompress (c.base, c.deltas, c.mask)), int_list (mixed));

  printf ('Worked Example: Base-Delta Compression of Timestamps\n');
  start = double (0x17F3A2C1000);
  stamps = start + double ([0, 0x120, 0x250, 0x3A8, 0x400, 0x5F0, 0x710, 0x7F8]);
  require ('1-byte deltas fail', isempty (bdi_compress (stamps, 8, 1, false)));
  c = bdi_compress (stamps, 8, 2, false);
  same ('deltas', int_list (c.deltas), int_list ([0, 288, 592, 936, 1024, 1520, 1808, 2040]));
  require ('the second delta exceeds 127', c.deltas(2) > 127);
  same ('compressed bytes', int64 (c.size), int64 (24));
  check ('compression ratio', 64 / c.size, '2.67');
  c4 = bdi_compress (stamps, 8, 4, false);
  same ('bytes with 4-byte deltas', int64 (c4.size), int64 (40));
  same ('decompression', int_list (bdi_decompress (c.base, c.deltas, c.mask)), int_list (stamps));
  require ('bit 12 is set in the base', mod (floor (c.base / 2 ^ 12), 2) == 1);   % (base >> 12) & 1
  flipped = bitxor (c.base, 2 ^ 12);
  same ('base after bit 12 flips', sprintf ('%X', flipped), '17F3A2C0000');
  same ('error in each timestamp', int_list (bdi_decompress (flipped, c.deltas, c.mask) - stamps), int_list (repmat (-4096, 1, 8)));
  same ('stored bits', int64 (8 * c.size), int64 (192));
  same ('bits in the line', int64 (8 * 64), int64 (512));
end

function check_examples ()
  global LETTERS PROBS
  printf ('Worked Example: Arithmetic Coding of a Sparse Binary Message\n');
  bp = [frac(9, 10), frac(1, 10)];   % the Python's {"0": 9/10, "1": 1/10}, in the order "01"
  msg = '0000100000';
  steps = arith_intervals (msg, bp, '01');
  same ('interval after four 0s', tup (steps(5).low, frac_add (steps(5).low, steps(5).width)), tup (frac (0), frac_dec ('0.6561')));
  same ('interval after the 1', tup (steps(6).low, frac_add (steps(6).low, steps(6).width)), tup (frac_dec ('0.59049'), frac_dec ('0.6561')));
  low = steps(11).low;
  width = steps(11).width;
  check ('final upper end', frac_add (low, width), '0.629232');
  check ('final width', width, '0.038742');
  check ('information in the message', log2 (width.den / width.num), '4.690');
  check ('10 h(0.1)', 10 * binary_entropy (0.1), '4.690');
  tag = arith_codeword (low, width);
  same ('tag', tag, '100110');
  check ('u times 2^6', to_float (low) * 64, '37.79');
  v = bin2dec (tag);
  same ('tag integer v', int64 (v), int64 (38));
  same ('dyadic interval', tup (frac (v, 64), frac (v + 1, 64)), tup (frac_dec ('0.59375'), frac_dec ('0.609375')));
  require ('inside the message interval', frac_cmp (low, frac (v, 64)) <= 0 && frac_cmp (frac (v + 1, 64), frac_add (low, width)) <= 0);
  same ('decoding', arith_decode (tag, 10, bp, '01'), msg);

  printf ('Worked Example: Learning the Model as It Goes\n');
  cases = {100, '115.68', '6.096', '0.061', '12.331'; 1000, '1156.78', '9.396', '0.0094', '18.936'};
  for i = 1:rows (cases)
    n = cases{i, 1};
    ka = floor (7 * n / 10);
    kb = floor (2 * n / 10);
    kc = floor (n / 10);
    true_bits = ka * log2 (10 / 7) + kb * log2 (10 / 2) + kc * log2 (10);
    check (sprintf ('bits from the true model, n = %d', n), true_bits, cases{i, 2});
    check (sprintf ('n H(X), n = %d', n), n * entropy (PROBS), cases{i, 2});
    message = [repmat('a', 1, ka) repmat('b', 1, kb) repmat('c', 1, kc)];
    laplace_bits = bits_of (laplace_closed_form (message, LETTERS));
    price = laplace_bits - true_bits;
    check (sprintf ('price of learning, n = %d', n), price, cases{i, 3});
    check (sprintf ('price per letter, n = %d', n), price / n, cases{i, 4});
    bound = log2 (nchoosek (n + 2, 2));
    check (sprintf ('guarantee, n = %d', n), bound, cases{i, 5});
    require (sprintf ('price is about half the guarantee, n = %d', n), 0.45 < price / bound && price / bound < 0.55);
  end
  p100 = bits_of (laplace_closed_form ([repmat('a', 1, 70) repmat('b', 1, 20) repmat('c', 1, 10)], LETTERS)) - 100 * entropy (PROBS);
  p1000 = bits_of (laplace_closed_form ([repmat('a', 1, 700) repmat('b', 1, 200) repmat('c', 1, 100)], LETTERS)) - 1000 * entropy (PROBS);
  check ('growth of the price from 100 to 1000 letters', p1000 - p100, '3.30');

  printf ('Exercise: Arithmetic coding of ba\n');
  steps = arith_intervals ('ba', PROBS, LETTERS);
  low = steps(end).low;
  width = steps(end).width;
  same ('interval of ba', tup (low, frac_add (low, width), width), tup (frac_dec ('0.7'), frac_dec ('0.84'), frac_dec ('0.14')));
  check ('information in ba', log2 (1 / 0.14), '2.8365');
  tag = arith_codeword (low, width);
  same ('tag', tag, '1100');
  check ('u times 2^4', to_float (low) * 16, '11.2');
  v = bin2dec (tag);
  same ('dyadic interval', tup (frac (v, 16), frac (v + 1, 16)), tup (frac_dec ('0.75'), frac_dec ('0.8125')));
  require ('inside the interval of ba', frac_cmp (low, frac (v, 16)) <= 0 && frac_cmp (frac (v + 1, 16), frac_add (low, width)) <= 0);
  same ('Huffman bits on ba', int64 (numel (['10' '0'])), int64 (3));

  printf ('Looking Ahead\n');
  same ('110110 decodes to', arith_decode ('110110', 6, PROBS, LETTERS), 'bbaaaa');
  same ('011110 decodes to', arith_decode ('011110', 6, PROBS, LETTERS), 'aacaba');
  every = true;
  for t = 0:63
    every = every && numel (arith_decode (dec2bin (t, 6), 6, PROBS, LETTERS)) == 6;
  end
  require ('every 6-bit string decodes to some message', every);
  hb = binary_entropy (0.1);
  check ('capacity 1 - h(0.1)', 1 - hb, '0.5310');
  require ('0.4690 < 0.5310', hb < 1 - hb);
end


% Main --------------------------------------------------------------------------------

global LETTERS PROBS FREQ TOTAL CUM
LETTERS = 'abc';
PROBS = [frac(7, 10), frac(2, 10), frac(1, 10)];   % a, b and c, in the order of LETTERS
FREQ = [7, 2, 1];   % integer counts for range coding and rANS
TOTAL = 10;
CUM = cumulative (FREQ, LETTERS);

check_intervals ();
check_rates (@block_rate);
check_range_coding ();
check_rans ();
check_dictionary ();
check_runs_bwt_mtf ();
check_lossy ();
check_memory ();
check_examples ();
printf ('All chapter 6 numbers reproduced.\n');
