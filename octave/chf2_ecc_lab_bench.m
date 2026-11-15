% Reproduce every number in front-matter chapter f2 of Error Correction Coding, "Setting Up the Lab Bench".
%
% This is the GNU Octave twin of python/chf2_ecc_lab_bench.py. It needs core Octave
% only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/chf2_ecc_lab_bench.m
%
% The Python takes a --bench option that inspects the tools installed on the machine
% (git, Verilator, Yosys, cocotb, NumPy, galois and the rest). That part belongs to the
% Python route, since it asks Python for its own version and for the Python packages
% installed. This twin reproduces the default run only, which is every number the
% chapter prints.
%
% The chapter installs the book's tools and then teaches the Python the listings
% use, so its numbers are of three kinds. The first is integer and bit arithmetic:
% binary forms, shifts, masks, floor division and remainders with negative
% operands. The second is floating point: the exact value stored for 0.1, the
% spacing of doubles, the sum 0.1 + 0.2, rounding half to even, and the half-unit
% test that every companion script's check() helper applies. The third is
% version numbers, compared as tuples, and the compatible-release rule that pip
% reads in "cocotb~=2.1". A short NumPy-style preview closes the chapter with a
% mod-2 matrix product and one product in the field of 256 elements, both done
% here with plain loops and integer bit operations. Each theorem the chapter proves
% is also checked by brute force over a range of inputs.
%
% Every value computed here is printed beside the value the chapter states. The
% script stops with exit status 1 at the first pair that disagrees by more than
% half a unit in the last digit the chapter prints, and at the first count, list
% or string that differs at all. tools/compare-routes.sh checks that this twin and
% the Python print the same lines.
%
% Octave has one number type where Python has two, so this file marks the difference.
% A Python int is an int64 here, so that py_str prints it without a decimal point,
% and a Python float is a double. Three habits follow from that and are worth knowing
% before reading on.
%   - Integer classes saturate instead of growing, and mixing one with a double rounds
%     the result to the integer class. A division of two int64 values therefore rounds
%     to nearest, so this file divides with idivide (a, b, 'floor') where Python
%     writes a // b, and mod gives Python's %, since mod takes the sign of the divisor.
%   - A double holds every integer below 2^53 exactly and no more. Python's 2 ** 100
%     is built by hand as decimal text, and Fraction is a pair of int64 values that
%     all stay below 2^63 in this chapter. Python's Decimal(0.1) is built the same way.
%   - Octave indexes from 1. Python's G[i][j] is G(i + 1, j + 1) where an index
%     touches an array, and a bit position keeps its Python value everywhere else.
%
% Local functions stand in for the pieces of Python that Octave lacks: Fraction, Decimal
% text, nextafter, ulp, the half-to-even round, tuple order and the bit helpers.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function t = agrees (computed, stated)
  % The test check() applies: within half a unit of the last digit stated.
  value = str2double (regexprep (strrep (stated, ',', ''), '%+$', ''));
  computed = to_float (computed);
  if (stated(end) == '%')
    computed = 100 * computed;
  end
  t = abs (computed - value) <= half_unit (stated) * (1 + 1e-9);
end

function check (label, computed, stated)
  value = str2double (regexprep (strrep (stated, ',', ''), '%+$', ''));
  computed = to_float (computed);
  if (stated(end) == '%')
    computed = 100 * computed;
  end
  printf ('  %-64s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  printf ('  %-64s %-16s %s\n', label, py_text (computed), py_text (stated));
  if (~py_equal (computed, stated))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_text (computed, true), py_text (stated, true));
    exit (1);
  end
end

function x = to_float (x)
  % Python's float (x): a Fraction divides its two integers, anything else converts.
  if (isstruct (x) && isfield (x, 'num'))
    x = double (x.num) / double (x.den);
  else
    x = double (x);
  end
end

function s = py_text (x, quoted)
  % str (x), or repr (x) when quoted, for the Python types py_str does not cover.
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
  elseif (isstruct (x) && isfield (x, 'dec'))
    if (quoted)
      s = ['Decimal(''' x.dec ''')'];
    else
      s = x.dec;
    end
  elseif (isstruct (x) && isfield (x, 'set'))
    items = arrayfun (@(v) sprintf ('%d', v), x.set, 'UniformOutput', false);
    s = ['{' strjoin(items, ', ') '}'];
  else
    s = py_str (x, quoted);
  end
end

function t = py_equal (a, b)
  % Python's a == b for the values this script compares. Decimals compare by their text,
  % which is canonical here, and the rest compare by value, element by element.
  if (isstruct (a) && isstruct (b) && isfield (a, 'dec') && isfield (b, 'dec'))
    t = strcmp (a.dec, b.dec);
  else
    t = isequal (a, b);
  end
end

function t = tup (v)
  % A row of integers marked as a Python tuple, so py_str prints (5, 54, 0).
  t = struct ('tuple', {num2cell(int64 (v))});
end

function c = tuple_cmp (a, b)
  % Python's order on tuples: -1, 0 or 1. The first difference decides, and when one
  % tuple is a prefix of the other, the shorter one is smaller.
  c = 0;
  for i = 1:min (numel (a), numel (b))
    if (a(i) < b(i))
      c = -1;
      return;
    elseif (a(i) > b(i))
      c = 1;
      return;
    end
  end
  c = sign (numel (a) - numel (b));
end

function t = str_lt (a, b)
  % Python's a < b for strings: the first differing character decides, and a prefix is smaller.
  t = tuple_cmp (double (a), double (b)) < 0;
end

function s = bin (n)
  % Python's bin (n) for a nonnegative integer.
  s = ['0b' dec2bin(double (n))];
end

function n = bit_length (x)
  % Python's int.bit_length ().
  if (x == 0)
    n = 0;
  else
    n = numel (dec2bin (double (x)));
  end
end

function s = dotted (v)
  % '.'.join (map (str, v))
  s = strjoin (arrayfun (@(x) sprintf ('%d', x), v, 'UniformOutput', false), '.');
end

function f = Fraction (a, b)
  % fractions.Fraction, as a reduced pair of int64 values with a positive denominator.
  % Fraction (x) takes the exact value of a double and Fraction (p, q) takes two integers.
  if (nargin == 1)
    [mant, ex] = log2 (a);               % a = mant * 2^ex, with 0.5 <= abs (mant) < 1
    num = int64 (mant * 2 ^ 53);         % an integer below 2^53, so the product is exact
    ex = ex - 53;
    if (ex >= 0)
      num = num * int64 (2) ^ ex;
      den = int64 (1);
    else
      den = int64 (2) ^ (-ex);
    end
  else
    num = int64 (a);
    den = int64 (b);
  end
  if (den < 0)
    num = -num;
    den = -den;
  end
  g = gcd (num, den);
  f = struct ('num', idivide (num, g, 'fix'), 'den', idivide (den, g, 'fix'));
end

function f = fr_sub (a, b)
  % a - b for two Fractions. Every product stays far below 2^63 in this chapter.
  f = Fraction (a.num * b.den - b.num * a.den, a.den * b.den);
end

function c = fr_cmp (a, b)
  % -1, 0 or 1 as a is below, equal to or above b, by cross multiplication.
  c = double (sign (a.num * b.den - b.num * a.den));
end

function r = fraction_round (x)
  % Python's round (Fraction): the nearest integer, a tie going to the even one.
  fl = idivide (x.num, x.den, 'floor');
  twice = 2 * (x.num - fl * x.den);      % twice the fractional part, in units of 1/den
  if (twice < x.den)
    r = fl;
  elseif (twice > x.den)
    r = fl + 1;
  else
    r = fl + mod (fl, 2);
  end
end

function r = py_round (x, ndigits)
  % Python 3's round. With no digits it returns an int, rounding a tie to the even integer,
  % where Octave's round sends a tie away from zero. With digits it returns a float, and
  % Python rounds the exact stored value to the nearest decimal, which is what printf does.
  if (nargin < 2)
    f = floor (x);
    d = x - f;
    if (d < 0.5)
      r = f;
    elseif (d > 0.5)
      r = f + 1;
    else
      r = f + mod (f, 2);                % a tie: keep f if it is even, else step up
    end
    r = int64 (r);
  else
    r = str2double (sprintf ('%.*f', ndigits, x));
  end
end

function u = ulp (x)
  % math.ulp (x): the gap from x to the next double above it, which is Octave's eps (x).
  u = eps (x);
end

function y = nextafter (x, toward)
  % math.nextafter for a positive normal double, the only kind used here. Core Octave has
  % no nextafter, so the neighbor is built from eps. Above x the gap is eps (x). Below x
  % it is the same, except when x is a power of 2, where the gap below is half as wide.
  if (toward > x)
    y = x + eps (x);
  else
    [f, e] = log2 (x);                   % x = f * 2^e with 0.5 <= f < 1
    if (f == 0.5)
      y = x - eps (x) / 2;
    else
      y = x - eps (x);
    end
  end
end

function d = digits_times (d, k)
  % A decimal number held as a row of digits, least significant first, times a small integer.
  carry = 0;
  for j = 1:numel (d)
    v = d(j) * k + carry;
    d(j) = mod (v, 10);
    carry = floor (v / 10);
  end
  while (carry > 0)
    d(end + 1) = mod (carry, 10);
    carry = floor (carry / 10);
  end
end

function t = pow2_text (k)
  % The exact decimal text of 2^k, beyond what a double or an int64 holds, by doubling.
  d = 1;
  for i = 1:k
    d = digits_times (d, 2);
  end
  t = char (fliplr (d) + '0');
end

function t = str_decimal (x)
  % str (Decimal (x)) for a positive double: its exact decimal expansion. A double is
  % m * 2^e with m an integer, and for e < 0 that is m * 5^(-e) shifted -e places right.
  [f, e] = log2 (x);
  m = f * 2 ^ 53;
  e = e - 53;
  while (mod (m, 2) == 0 && e < 0)
    m = m / 2;
    e = e + 1;
  end
  d = fliplr (sprintf ('%d', m) - '0');
  if (e < 0)
    for i = 1:-e
      d = digits_times (d, 5);
    end
  else
    for i = 1:e
      d = digits_times (d, 2);
    end
  end
  text = char (fliplr (d) + '0');
  if (e < 0)
    text = [repmat('0', 1, max (0, -e + 1 - numel (text))) text];
    t = [text(1:end + e) '.' text(end + e + 1:end)];
  else
    t = text;
  end
end

function d = Decimal (text)
  % decimal.Decimal of a string: the text itself, in canonical form.
  d = struct ('dec', text);
end

function d = exact (x)
  % A terminating fraction as an exact decimal, so it prints the way the chapter writes it.
  % The denominator holds only 2s and 5s, so some power of ten 10^k is a multiple of it,
  % and x is then an integer shifted k places.
  k = 0;
  pow10 = int64 (1);
  while (mod (pow10, x.den) ~= 0)
    pow10 = pow10 * 10;
    k = k + 1;
  end
  digits = sprintf ('%d', x.num * idivide (pow10, x.den, 'fix'));
  if (k > 0)
    digits = [repmat('0', 1, max (0, k + 1 - numel (digits))) digits];
    digits = [digits(1:end - k) '.' digits(end - k + 1:end)];
  end
  d = struct ('dec', digits);
end

function s = set_of (v)
  % A Python set of small integers, listed in increasing order as Python prints it.
  s = struct ('set', int64 (sort (v)));
end


% Integers and bits ------------------------------------------------------------

function check_integers_and_bits ()
  printf ('Integers and bits\n');
  same ('bin(100)', bin (100), '0b1100100');
  same ("format(100, '08b')", dec2bin (100, 8), '01100100');
  same ("int('01100100', 2)", int64 (bin2dec ('01100100')), int64 (100));
  same ('100 ^ 0b10000000', bitxor (int64 (100), int64 (0b10000000)), int64 (228));
  same ('100 >> 2', bitshift (int64 (100), -2), int64 (25));
  same ('100 << 1', bitshift (int64 (100), 1), int64 (200));
  same ('100 & 1, the lowest bit', bitand (int64 (100), int64 (1)), int64 (0));
  same ('(100 >> 2) & 1, bit 2', bitand (bitshift (int64 (100), -2), int64 (1)), int64 (1));
  same ('100 | 1', bitor (int64 (100), int64 (1)), int64 (101));
  same ("ones in 100, bin(100).count('1')", int64 (sum (bin (100) == '1')), int64 (3));
  same ('0x64, the same number in hexadecimal', int64 (0x64), int64 (100));
  same ('2 ** 100, exactly', pow2_text (100), '1267650600228229401496703205376');
  same ('digits in 2 ** 100', int64 (numel (pow2_text (100))), int64 (31));
  same ('100 ^ 29', bitxor (int64 (100), int64 (29)), int64 (121));
  same ('121 ^ 29, undoing the XOR', bitxor (int64 (121), int64 (29)), int64 (100));
end

function check_shift_theorem ()
  % n >> k is floor(n / 2^k) and n << k is n * 2^k, for every integer n and k >= 0.
  printf ('Shifts are floor division and multiplication by powers of 2\n');
  cases = int64 (0);
  % The Python loops over n and k. Here one pass over k handles every n at once. bitshift
  % on a signed int64 shifts the two's complement bits, so it is Python's arithmetic shift.
  n = int64 (-1024:1024);
  for k = 0:11
    assert (all (bitshift (n, -k) == idivide (n, int64 (2) ^ k, 'floor')), 'n >> k differs from n // 2 ** k at k = %d', k);
    assert (all (bitshift (n, k) == n * int64 (2) ^ k), 'n << k differs from n * 2 ** k at k = %d', k);
    cases = cases + numel (n);
  end
  same ('cases checked, n from -1024 to 1024, k from 0 to 11', cases, int64 (24588));
  same ('-7 >> 1', bitshift (int64 (-7), -1), int64 (-4));
  same ('13 >> 2', bitshift (int64 (13), -2), int64 (3));
end

function check_floor_division ()
  printf ('Floor division and remainder\n');
  % idivide with 'floor' is Python's //, and mod is Python's %. Both follow the divisor's sign.
  same ('7 // 2', idivide (int64 (7), int64 (2), 'floor'), int64 (3));
  same ('7 % 2', mod (int64 (7), int64 (2)), int64 (1));
  same ('-7 // 2', idivide (int64 (-7), int64 (2), 'floor'), int64 (-4));
  same ('-7 % 2', mod (int64 (-7), int64 (2)), int64 (1));
  same ('7 // -2', idivide (int64 (7), int64 (-2), 'floor'), int64 (-4));
  same ('7 % -2', mod (int64 (7), int64 (-2)), int64 (-1));
  same ('int(-7 / 2), truncation toward zero', int64 (fix (-7 / 2)), int64 (-3));
  same ('divmod(17, 5)', tup ([idivide(int64 (17), int64 (5), 'floor'), mod(int64 (17), int64 (5))]), tup ([3 2]));
  same ('-17 % 5', mod (int64 (-17), int64 (5)), int64 (3));
  same ('1000 % 7', mod (int64 (1000), int64 (7)), int64 (6));
  pairs = int64 (0);
  % The Python loops over a and b. Here one pass over b handles every a at once.
  a = -200:200;
  for b = 1:20
    q = floor (a / b);
    r = mod (a, b);
    assert (all (a == b * q + r & r >= 0 & r < b), 'a = b q + r with 0 <= r < b fails at b = %d', b);
    % The quotient and remainder are the only pair with that property.
    count = zeros (size (a));
    for qq_offset = -3:3
      qq = q + qq_offset;
      count = count + (a - b * qq >= 0 & a - b * qq < b);
    end
    assert (all (count == 1), 'the pair is not unique at b = %d', b);
    pairs = pairs + numel (a);
  end
  same ('pairs checked, a from -200 to 200, b from 1 to 20', pairs, int64 (8020));
end


% Floating point ---------------------------------------------------------------

function out = binary_digits (x, n)
  % The first n bits after the binary point of a fraction in [0, 1).
  num = x.num;
  den = x.den;
  out = repmat ('0', 1, n);
  for i = 1:n
    num = num * 2;
    bit = num >= den;
    out(i) = char ('0' + bit);
    num = num - int64 (bit) * den;
  end
end

function t = terminates (x, limit)
  % Does repeated doubling clear the fraction within the limit? Here x.num may also be a
  % row of numerators over the one denominator x.den, which tests many fractions at once.
  if (nargin < 2)
    limit = 64;
  end
  num = mod (x.num, x.den);              % x -= int(x) for x >= 0
  for i = 1:limit
    num = mod (2 * num, x.den);          % x *= 2, then drop the whole part
  end
  % Once the numerator reaches 0 it stays 0, so the Python's early return and its final
  % test x == 0 both amount to this single test after limit doublings.
  t = (num == 0);
end

function check_finite_expansion_theorem ()
  % p/q in lowest terms has a finite binary expansion exactly when q is a power of 2.
  printf ('Finite binary expansions\n');
  tried = int64 (0);
  for q = 1:256
    power = bitand (q, q - 1) == 0;
    p = 0:q - 1;
    % Fraction(p, q) has denominator q exactly when p/q is already in lowest terms.
    p = p(gcd (p, q) == 1);
    assert (all (terminates (struct ('num', p, 'den', q)) == power), 'p/q fails at q = %d', q);
    tried = tried + numel (p);
  end
  same ('fractions p/q in lowest terms checked, q up to 256', tried, int64 (19948));
  same ('first 24 bits of 1/10', binary_digits (Fraction (1, 10), 24), '000110011001100110011001');
  same ('first 8 bits of 1/10, as a fraction', Fraction (int64 (bin2dec (binary_digits (Fraction (1, 10), 8))), 256), Fraction (25, 256));
  check ('25/256', 25 / 256, '0.09765625');
  check ('error of the 8-bit truncation', 0.1 - 25 / 256, '0.00234375');
  same ('does 3/8 terminate', terminates (Fraction (3, 8)), true);
  same ('3/8 in binary', binary_digits (Fraction (3, 8), 3), '011');
  same ('does 1/3 terminate', terminates (Fraction (1, 3)), false);
  same ('first 8 bits of 1/3', binary_digits (Fraction (1, 3), 8), '01010101');
end

function check_point_one ()
  printf ('The double nearest 0.1\n');
  stored = Fraction (0.1);
  same ('Fraction(0.1)', stored, Fraction (int64 (3602879701896397), int64 (36028797018963968)));
  same ('its denominator is 2 ** 55', stored.den == int64 (2) ^ 55, true);
  same ('2 ** 55', int64 (2) ^ 55, int64 (36028797018963968));
  same ('2 ** 56 / 10, as a fraction', Fraction (int64 (2) ^ 56, 10), Fraction (int64 (36028797018963968), 5));
  check ('2 ** 56 / 10', Fraction (int64 (2) ^ 56, 10), '7205759403792793.6');
  same ('nearest integer to 2 ** 56 / 10', fraction_round (Fraction (int64 (2) ^ 56, 10)), int64 (7205759403792794));
  same ('7205759403792794 / 2 ** 56 in lowest terms', Fraction (int64 (7205759403792794), int64 (2) ^ 56), stored);
  same ('2 ** 52 <= 7205759403792794 < 2 ** 53', int64 (2) ^ 52 <= int64 (7205759403792794) && int64 (7205759403792794) < int64 (2) ^ 53, true);
  same ('stored minus 1/10, exactly', fr_sub (stored, Fraction (1, 10)), Fraction (1, 5 * int64 (2) ^ 55));
  check ('stored minus 1/10', fr_sub (stored, Fraction (1, 10)), '5.551e-18');
  same ('Decimal(0.1)', str_decimal (0.1), '0.1000000000000000055511151231257827021181583404541015625');
  same ('repr(0.1)', py_str (0.1), '0.1');
  same ('0.1 + 0.2', py_str (0.1 + 0.2), '0.30000000000000004');
  same ('0.1 + 0.2 == 0.3', 0.1 + 0.2 == 0.3, false);
  same ('Decimal(0.1 + 0.2)', str_decimal (0.1 + 0.2), '0.3000000000000000444089209850062616169452667236328125');
  same ('Decimal(0.3)', str_decimal (0.3), '0.299999999999999988897769753748434595763683319091796875');
  same ('(0.1 + 0.2) - 0.3 is 2 ** -54', (0.1 + 0.2) - 0.3 == 2.0 ^ -54, true);
  check ('2 ** -54', 2.0 ^ -54, '5.551e-17');
  same ('0.1 + 0.2 is the double just above 0.3', nextafter (0.3, inf) == 0.1 + 0.2, true);
  same ('abs((0.1 + 0.2) - 0.3) < 1e-9', abs ((0.1 + 0.2) - 0.3) < 1e-9, true);
end

function check_spacing_theorem ()
  % Consecutive doubles in [2^e, 2^(e+1)) are 2^(e-52) apart.
  printf ('The spacing of doubles\n');
  for e = -60:60
    x = 2.0 ^ e;
    gap = 2.0 ^ (e - 52);
    assert (ulp (x) == gap, 'ulp (x) is not the gap at e = %d', e);
    assert (nextafter (x, inf) - x == gap, 'the gap above x is wrong at e = %d', e);
    y = nextafter (2.0 ^ (e + 1), 0.0);
    assert (y - nextafter (y, 0.0) == gap, 'the gap below y is wrong at e = %d', e);
  end
  same ('intervals [2**e, 2**(e+1)) checked, e from -60 to 60', int64 (numel (-60:60)), int64 (121));
  same ('ulp(1.0) is 2 ** -52', ulp (1.0) == 2.0 ^ -52, true);
  check ('ulp(1.0)', ulp (1.0), '2.220e-16');
  same ('ulp(0.1) is 2 ** -56', ulp (0.1) == 2.0 ^ -56, true);
  check ('ulp(0.1)', ulp (0.1), '1.388e-17');
  same ('ulp(1000.0) is 2 ** -43', ulp (1000.0) == 2.0 ^ -43, true);
  check ('ulp(1000.0)', ulp (1000.0), '1.137e-13');
  same ('ulp(2.0 ** 53) is 2', ulp (2.0 ^ 53), 2.0);
  same ('2 ** 53', int64 (2) ^ 53, int64 (9007199254740992));
  % A Python int is exact, so int64 stands in for it, and double () rounds it to even as float () does.
  same ('float(2 ** 53) + 1.0 == float(2 ** 53)', double (int64 (2) ^ 53) + 1.0 == double (int64 (2) ^ 53), true);
  same ('float(2 ** 53 + 1) == float(2 ** 53)', double (int64 (2) ^ 53 + 1) == double (int64 (2) ^ 53), true);
  same ('float(2 ** 53 + 2) - float(2 ** 53)', double (int64 (2) ^ 53 + 2) - double (int64 (2) ^ 53), 2.0);
  same ('doubles in [1, 2)', int64 (2) ^ 52, int64 (4503599627370496));
end

function check_rounding ()
  printf ('Rounding half to even\n');
  same ('round(2.5)', py_round (2.5), int64 (2));
  same ('round(3.5)', py_round (3.5), int64 (4));
  same ('round(-2.5)', py_round (-2.5), int64 (-2));
  same ('round(0.5)', py_round (0.5), int64 (0));
  same ('round(1.5)', py_round (1.5), int64 (2));
  same ('round(0.125, 2)', py_round (0.125, 2), 0.12);
  same ('round(0.375, 2)', py_round (0.375, 2), 0.38);
  same ('0.125 is stored exactly', Fraction (0.125), Fraction (1, 8));
  same ('round(2.675, 2)', py_round (2.675, 2), 2.67);
  same ('Decimal(2.675)', str_decimal (2.675), '2.67499999999999982236431605997495353221893310546875');
  same ("f'{2.675:.2f}'", sprintf ('%.2f', 2.675), '2.67');
  same ("f'{0.96583769633507:.5f}'", sprintf ('%.5f', 0.96583769633507), '0.96584');
  % Ties to even leave no drift on average: the four ties 0.5, 1.5, 2.5, 3.5
  % round to 0, 2, 2, 4, which sum to 8, the same as the ties themselves.
  ties = [0.5, 1.5, 2.5, 3.5];
  total = int64 (0);
  for t = ties
    total = total + py_round (t);
  end
  same ('sum of round(t) for t in 0.5, 1.5, 2.5, 3.5', total, int64 (8));
  same ('sum of 0.5, 1.5, 2.5, 3.5', 0.5 + 1.5 + 2.5 + 3.5, 8.0);
  same ('sum when every tie rounds up', int64 (sum (fix (ties + 0.5))), int64 (10));
end


% The half-unit test -----------------------------------------------------------

function check_half_unit_theorem ()
  % A d-place decimal s is a correct rounding of x exactly when |x - s| <= 10^-d / 2.
  printf ('The half-unit test\n');
  same ("half_unit('0.96584')", half_unit ('0.96584'), 0.5 * 10.0 ^ -5);
  check ("half_unit('0.96584')", half_unit ('0.96584'), '0.000005');
  check ("half_unit('7,379')", half_unit ('7,379'), '0.5');
  check ("half_unit('2.220e-16')", half_unit ('2.220e-16'), '5e-20');
  check ("half_unit('96.58%')", half_unit ('96.58%'), '0.005');
  check ('7379 / 7640', 7379 / 7640, '0.96584');
  same ("7379 / 7640 agrees with '0.96584'", agrees (7379 / 7640, '0.96584'), true);
  same ("7379 / 7640 agrees with '0.9658'", agrees (7379 / 7640, '0.9658'), true);
  same ("7379 / 7640 agrees with '0.9659'", agrees (7379 / 7640, '0.9659'), false);
  same ("7379 / 7640 agrees with '96.58%'", agrees (7379 / 7640, '96.58%'), true);
  check ('distance from 0.9659', abs (7379 / 7640 - 0.9659), '0.0000623');
  % Brute force: every x on a fine grid, every d from 0 to 3. A d-place decimal
  % passes the test exactly when it is a nearest d-place decimal to x.
  % The Python uses Fraction. Here every distance is counted in whole units of
  % 1/(160 * 10^d), the unit in which x = i/160 and s = j/10^d are both integers, so a
  % double holds each one exactly. In that unit the half unit, 1/(2 * 10^d), is 80.
  tried = int64 (0);
  for d = 0:3
    scale = 10 ^ d;
    for i = -400:400
      xs = i * scale;                                 % x = i/160, in those units
      js = -30 * scale:30 * scale;
      nearest = min (abs (xs - 160 * js));            % min over j of |x - unit * j|
      low = floor (i * scale / 160);                  % floor (x / unit)
      for j = [low - 1, low, low + 1, low + 2]
        dist = abs (xs - 160 * j);                    % |x - s| for s = unit * j
        passes = dist <= 80;
        assert (passes == (dist == nearest), 'half-unit test fails at i = %d, d = %d, j = %d', i, d, j);
      end
      tried = tried + 1;
    end
  end
  same ('grid points checked, 801 points times 4 precisions', tried, int64 (3204));
end


% Version numbers --------------------------------------------------------------

function v = version_tuple (text)
  % The first dotted number in a tool's output, padded to three parts. This is the Python's
  % version (), renamed because Octave already has a built-in function called version.
  tokens = regexp (text, '(\d+)\.(\d+)(?:\.(\d+))?', 'tokens', 'once');
  if (isempty (tokens))
    error ('no version number in ''%s''', text);
  end
  v = zeros (1, 3);
  for i = 1:numel (tokens)         % an unmatched third part is left out, so it stays 0
    v(i) = str2double (tokens{i});
  end
end

function t = compatible (v, spec)
  % PEP 440's ~= operator: at least spec, and the same release series.
  parts = str2double (strsplit (spec, '.'));
  prefix = parts(1:end - 1);
  padded = [v, zeros(1, max (0, numel (parts) - numel (v)))];
  t = tuple_cmp (padded, parts) >= 0 && isequal (padded(1:min (numel (prefix), numel (padded))), prefix);
end

function check_versions ()
  printf ('Version numbers\n');
  same ("version('Verilator 5.054 2026-11-01')", tup (version_tuple ('Verilator 5.054 2026-11-01 rev vUNKNOWN')), tup ([5 54 0]));
  same ("version('Yosys 0.71+post (git sha1 ...)')", tup (version_tuple ('Yosys 0.71+post (git sha1 5628c13)')), tup ([0 71 0]));
  same ("version('SBY v0.71')", tup (version_tuple ('SBY v0.71')), tup ([0 71 0]));
  same ("version('surfer 0.7.0 (git: v0.7.0)')", tup (version_tuple ('surfer 0.7.0 (git: v0.7.0)')), tup ([0 7 0]));
  same ("version('Icarus Verilog version 13.0 (stable)')", tup (version_tuple ('Icarus Verilog version 13.0 (stable) (v13_0)')), tup ([13 0 0]));
  same ("version('git version 2.56.1')", tup (version_tuple ('git version 2.56.1')), tup ([2 56 1]));
  same ('Verilator 5.054 meets the floor 5.036', tuple_cmp (version_tuple ('5.054'), version_tuple ('5.036')) >= 0, true);
  same ('Verilator 5.020 meets the floor 5.036', tuple_cmp (version_tuple ('5.020'), version_tuple ('5.036')) >= 0, false);
  same ("'3.10' < '3.9' as strings", str_lt ('3.10', '3.9'), true);
  same ('3.10 < 3.9 as decimals', 3.10 < 3.9, true);
  same ('(3, 10) < (3, 9) as tuples', tuple_cmp ([3 10], [3 9]) < 0, false);
  same ("'0.71' < '0.9' as strings", str_lt ('0.71', '0.9'), true);
  same ('(0, 71) > (0, 9) as tuples', tuple_cmp ([0 71], [0 9]) > 0, true);
  same ("(3, 10) < (3, 10, 0), Python's prefix rule", tuple_cmp ([3 10], [3 10 0]) < 0, true);
  same ('padded, (3, 10, 0) == (3, 10, 0)', isequal (version_tuple ('3.10'), version_tuple ('3.10.0')), true);
  % Octave has no sys.version_info, so its own version, read as a tuple, stands in for it.
  version_info = sscanf (OCTAVE_VERSION, '%d.%d.%d').';
  same ('sys.version_info compares as a tuple', tuple_cmp (version_info, [3 10]) >= 0, true);
  cases = {[2 0 0], '2.1',   false; [2 1 0], '2.1',   true;  [2 1 5], '2.1',   true;
           [2 9 3], '2.1',   true;  [3 0 0], '2.1',   false;
           [2 1 0], '2.1.0', true;  [2 1 7], '2.1.0', true;  [2 2 0], '2.1.0', false};
  for row = 1:rows (cases)
    v = cases{row, 1};
    spec = cases{row, 2};
    want = cases{row, 3};
    same ([dotted(v) ' satisfies ~=' spec], compatible (v, spec), want);
  end
  % Lexicographic order on padded tuples is total and transitive.
  triples = zeros (27, 3);
  for i = 0:26
    triples(i + 1, :) = [floor(i / 9), mod(floor (i / 3), 3), mod(i, 3)];
  end
  n = rows (triples);
  less = false (n, n);
  for i = 1:n
    for j = 1:n
      c = tuple_cmp (triples(i, :), triples(j, :));
      assert ((c < 0) + (c == 0) + (c > 0) == 1, 'the order is not total at (%d, %d)', i, j);
      less(i, j) = c < 0;
    end
  end
  % Transitivity: whenever a < b and b < c, a < c. The matrix product counts the b for each (a, c).
  chained = (double (less) * double (less)) > 0;
  assert (all (less(chained)), 'the order is not transitive');
  same ('triples checked for totality and transitivity', int64 (n), int64 (27));
end


% A preview of NumPy and galois, in plain Python --------------------------------

function G = generator_matrix ()
  % The 4 by 7 generator matrix of the chapter's small code. The Python keeps it as the
  % module constant G. Octave functions cannot see script variables, so this returns it.
  G = [1 0 0 0 1 1 0;
       0 1 0 0 1 0 1;
       0 0 1 0 0 1 1;
       0 0 0 1 1 1 1];
end

function c = encode (m)
  % The row vector m times G, reduced mod 2: what (m @ G) % 2 computes in NumPy.
  c = int64 (mod (m * generator_matrix (), 2));
end

function out = clmul (a, b)
  % Carry-less product: polynomial multiplication with coefficients mod 2.
  out = int64 (0);
  while (b ~= 0)
    if (bitand (b, int64 (1)) ~= 0)
      out = bitxor (out, a);
    end
    a = bitshift (a, 1);
    b = bitshift (b, -1);
  end
end

function x = reduce (x, poly)
  % The remainder of x after division by the degree-8 polynomial poly, mod 2.
  if (nargin < 2)
    poly = int64 (0x11D);
  end
  while (bit_length (x) > 8)
    x = bitxor (x, bitshift (poly, bit_length (x) - 9));
  end
end

function r = poly_add (p, q)
  % Coefficient-wise sum mod 2 of two coefficient lists, lowest power first. When q has
  % several rows, each row is added to p, which lets one call handle many polynomials.
  r = mod (p + q, 2);
end

function b = bits_of (v)
  % One row of the eight low bits per entry of v, lowest first: (v >> i) & 1 for i in range(8).
  b = mod (floor (v(:) ./ 2 .^ (0:7)), 2);
end

function check_preview ()
  printf ('A preview of NumPy and galois\n');
  G = generator_matrix ();
  same ('m = [1, 0, 1, 1] times G, mod 2', encode ([1 0 1 1]), int64 ([1 0 1 1 0 1 0]));
  % G(logical (m), :) keeps the rows of G where m is 1, and sum (..., 1) adds them down the columns.
  same ('the sum before reduction', int64 (sum (G(logical ([1 0 1 1]), :), 1)), int64 ([1 0 1 1 2 3 2]));
  same ('[1, 1, 0, 0] times G, mod 2', encode ([1 1 0 0]), int64 ([1 1 0 0 0 1 1]));
  messages = dec2bin (0:15, 4) - '0';             % product((0, 1), repeat=4)
  words = zeros (16, 7, 'int64');
  for i = 1:16
    words(i, :) = encode (messages(i, :));
  end
  codewords = unique (words, 'rows');             % the Python collects them in a set
  same ('distinct codewords from 16 messages', int64 (rows (codewords)), int64 (16));
  nonzero = codewords(any (codewords, 2), :);
  same ('smallest weight of a nonzero codeword', int64 (min (sum (nonzero, 2))), int64 (3));
  same ('0x11D in binary', bin (0x11D), '0b100011101');
  same ('100 in binary', bin (100), '0b1100100');
  same ('29 in binary', bin (29), '0b11101');
  same ('carry-less 100 times 29', clmul (int64 (100), int64 (29)), int64 (1172));
  same ('1172 in binary', bin (1172), '0b10010010100');
  same ('ordinary 100 times 29', int64 (100) * 29, int64 (2900));
  same ('0x11D shifted left by 2', bitshift (int64 (0x11D), 2), int64 (1140));
  same ('1172 ^ 1140', bitxor (int64 (1172), int64 (1140)), int64 (224));
  same ('1172 reduced by 0x11D', reduce (int64 (1172)), int64 (224));
  same ('100 plus 29 in GF(2^8), an XOR', bitxor (int64 (100), int64 (29)), int64 (121));
  % Polynomial addition mod 2 is XOR of the integers that hold the coefficients.
  % The Python loops over a and b. Here one pass over a handles all 256 values of b at once.
  pairs = int64 (0);
  b = 0:255;
  pb = bits_of (b);
  for a = 0:255
    pa = bits_of (a);
    assert (isequal (poly_add (pa, pb), bits_of (bitxor (a, b))), 'polynomial sum differs from XOR at a = %d', a);
    pairs = pairs + numel (b);
  end
  same ('pairs of degree-7 polynomials checked', pairs, int64 (65536));
end


% Worked examples and exercises --------------------------------------------------

function check_examples_and_exercises ()
  G = generator_matrix ();
  printf ('The truncation error of 1/10, worked examples and exercises\n');
  tenth = Fraction (1, 10);
  for k = 1:24
    two_k = int64 (2) ^ k;
    e_k = Fraction (mod (two_k, 10), 10 * two_k);
    assert (isequal (e_k, fr_sub (tenth, Fraction (idivide (two_k, int64 (10), 'floor'), two_k))), 'e_k formula fails at k = %d', k);
    assert (fr_cmp (Fraction (1, 5 * two_k), e_k) <= 0 && fr_cmp (e_k, Fraction (8, 10 * two_k)) <= 0 ...
            && fr_cmp (Fraction (8, 10 * two_k), Fraction (1, two_k)) < 0, 'e_k bounds fail at k = %d', k);
  end
  same ('e_k formula checked, k from 1 to 24', int64 (numel (1:24)), int64 (24));
  same ('last digit of 2**k, k = 1 to 8', int64 (mod (2 .^ (1:8), 10)), int64 ([2 4 8 6 2 4 8 6]));
  same ('e_8', exact (Fraction (mod (int64 (2) ^ 8, 10), 10 * int64 (2) ^ 8)), Decimal ('0.00234375'));
  same ('e_4', exact (Fraction (mod (int64 (2) ^ 4, 10), 10 * int64 (2) ^ 4)), Decimal ('0.0375'));
  same ('e_4 as 6/160', exact (Fraction (6, 160)), Decimal ('0.0375'));
  same ('1/10 cut to four places', exact (Fraction (idivide (int64 (2) ^ 4, int64 (10), 'floor'), int64 (2) ^ 4)), Decimal ('0.0625'));
  same ('2**56 / 10', exact (Fraction (int64 (2) ^ 56, 10)), Decimal ('7205759403792793.6'));
  same ('2**52 <= 2**56/10 < 2**53', fr_cmp (Fraction (int64 (2) ^ 52, 1), Fraction (int64 (2) ^ 56, 10)) <= 0 ...
        && fr_cmp (Fraction (int64 (2) ^ 56, 10), Fraction (int64 (2) ^ 53, 1)) < 0, true);
  same ('2**-4 <= 0.1 < 2**-3', 2.0 ^ -4 <= 0.1 && 0.1 < 2.0 ^ -3, true);
  same ('ulp(0.3) is 2 ** -54', ulp (0.3) == 2.0 ^ -54, true);

  ties = [0.5, 1.5, 2.5, 3.5];
  up = int64 (fix (ties + 0.5));
  same ('the ties rounded half up', up, int64 ([1 2 3 4]));
  same ('their sum', int64 (sum (up)), int64 (10));

  printf ('Example: the byte that holds 29\n');
  same ('29 == 16 + 8 + 4 + 1', 29 == 16 + 8 + 4 + 1, true);
  same ("format(29, '08b')", dec2bin (29, 8), '00011101');
  same ('29 >> 3', bitshift (int64 (29), -3), int64 (3));
  same ('29 // 8', idivide (int64 (29), int64 (8), 'floor'), int64 (3));
  same ('29 & 1', bitand (int64 (29), int64 (1)), int64 (1));

  printf ('Example: dividing -17 by 5\n');
  same ('-17 // 5', idivide (int64 (-17), int64 (5), 'floor'), int64 (-4));
  same ('-17 % 5', mod (int64 (-17), int64 (5)), int64 (3));
  same ('(-4) * 5 + 3', int64 (-4) * 5 + 3, int64 (-17));
  same ('int(-17 / 5)', int64 (fix (-17 / 5)), int64 (-3));
  same ('-17 - 5 * int(-17 / 5)', int64 (-17) - 5 * int64 (fix (-17 / 5)), int64 (-2));

  printf ('Example: which fractions terminate in binary\n');
  same ('3/8 terminates', terminates (Fraction (3, 8)), true);
  same ('5/12 terminates', terminates (Fraction (5, 12)), false);
  same ('7/64 terminates', terminates (Fraction (7, 64)), true);
  same ('1/5 terminates', terminates (Fraction (1, 5)), false);
  same ('6/12 terminates', terminates (Fraction (6, 12)), true);
  same ('6/12 in lowest terms', Fraction (6, 12), Fraction (1, 2));

  printf ('Example: the gap at 3.0 and at 1024.0\n');
  same ('ulp(3.0) is 2 ** -51', ulp (3.0) == 2.0 ^ -51, true);
  same ('ulp(1024.0) is 2 ** -42', ulp (1024.0) == 2.0 ^ -42, true);
  below = 1024.0 - nextafter (1024.0, 0.0);
  same ('gap just below 1024.0 is 2 ** -43', below == 2.0 ^ -43, true);
  same ('gap just below 1024.0 equals ulp(1000.0)', below == ulp (1000.0), true);

  printf ('Example: reading the pin cocotb~=2.1\n');
  pins = {[2 0 5], false; [2 1 3], true; [2 4 0], true; [3 1 0], false};
  for row = 1:rows (pins)
    v = pins{row, 1};
    verdict = pins{row, 2};
    same ([dotted(v) ' accepted by ~=2.1'], compatible (v, '2.1'), verdict);
  end

  printf ('Exercises\n');
  same ('[(100 >> i) & 1 for i in range(8)]', int64 (bitand (bitshift (100, -(0:7)), 1)), int64 ([0 0 1 0 0 1 1 0]));
  same ('round(4.5)', py_round (4.5), int64 (4));
  same ('round(5.5)', py_round (5.5), int64 (6));
  same ('round(-0.5) == 0', py_round (-0.5) == 0, true);
  same ("7379/7640 agrees with '0.966'", agrees (7379 / 7640, '0.966'), true);
  check ('distance from 0.966', abs (7379 / 7640 - 0.966), '0.000162');
  same ("half unit of '0.966'", half_unit ('0.966'), 0.5 * 10.0 ^ -3);
  % Python's rows 0 and 3 of G are rows 1 and 4 here.
  sums = int64 (sum (G([1 4], :), 1));
  same ('column sums of rows 0 and 3', sums, int64 ([1 0 0 1 2 2 1]));
  same ('encode([1, 0, 0, 1])', encode ([1 0 0 1]), int64 ([1 0 0 1 0 0 1]));
  same ('0x11D', int64 (0x11D), int64 (285));
  same ('256 + 16 + 13', int64 (256) + 16 + 13, int64 (285));
  same ("f'{100:08b}'", dec2bin (100, 8), '01100100');
  same ('(0, 71) > (0, 9)', tuple_cmp (version_tuple ('Yosys 0.71'), version_tuple ('Yosys 0.9')) > 0, true);
  same ("'0.71' < '0.9' as strings", str_lt ('0.71', '0.9'), true);
  same ("'3.10' < '3.9' as strings", str_lt ('3.10', '3.9'), true);

  printf ('Version table cells\n');
  table = {[2 2 0], '2.1',   true;  [2 0 0], '2.1.0', false; [2 1 5], '2.1.0', true;
           [2 9 3], '2.1.0', false; [3 0 0], '2.1.0', false; [2 0 0], '2.1',   false;
           [3 0 0], '2.1',   false};
  for row = 1:rows (table)
    v = table{row, 1};
    spec = table{row, 2};
    verdict = table{row, 3};
    same ([dotted(v) ' accepted by ~=' spec], compatible (v, spec), verdict);
  end

  printf ('Bits named in the preview\n');
  % Python's bit position i is the bit worth 2^i, and find returns it plus 1.
  same ('bits set in 100', set_of (find (bits_of (100)) - 1), set_of ([6 5 2]));
  same ('bits set in 29', set_of (find (bits_of (29)) - 1), set_of ([4 3 2 0]));
  same ('(1172).bit_length()', int64 (bit_length (1172)), int64 (11));
end


% The checks -------------------------------------------------------------------

function main ()
  check_integers_and_bits ();
  check_shift_theorem ();
  check_floor_division ();
  check_finite_expansion_theorem ();
  check_point_one ();
  check_spacing_theorem ();
  check_rounding ();
  check_half_unit_theorem ();
  check_versions ();
  check_preview ();
  check_examples_and_exercises ();
  printf ('\nAll chapter f2 numbers reproduced.\n');
end

main ();
