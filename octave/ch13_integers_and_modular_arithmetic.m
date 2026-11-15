% Reproduce every number in chapter 13 of Error Correction Coding,
% "Integers and Modular Arithmetic", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch13_integers_and_modular_arithmetic.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch13_integers_and_modular_arithmetic.py and
% prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. All the arithmetic is on whole numbers small enough for a double
% to hold exactly, so every comparison is exact. Beside the worked numbers,
% each theorem the chapter proves is also tested by brute force on small
% cases, over the same ranges the Python uses. Where the Python leans on its
% own built-ins, // and % and math.fmod and pow(a, -1, n), this file leans on
% Octave's floor, mod, rem and the three-output gcd, so the two routes check
% each other as well as the chapter.

1;


% Helpers ----------------------------------------------------------------------

function s = num (v)
  s = sprintf ('%d', v);
end

function s = lst (v)
  % A row of whole numbers written the way Python prints a list of ints.
  if (isempty (v))
    s = '[]';
  else
    parts = arrayfun (@(x) sprintf ('%d', x), v, 'UniformOutput', false);
    s = ['[' strjoin(parts, ', ') ']'];
  end
end

function same (label, c, s)
  printf ('  %-62s %-16s %s\n', label, c, s);
  if (~strcmp (c, s))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, c, s);
    exit (1);
  end
end

function sn (label, c, s)
  same (label, num (c), num (s));
end

function sl (label, c, s)
  same (label, lst (c), lst (s));
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


% The arithmetic ---------------------------------------------------------------

function r = floor_divmod (a, b)
  q = floor (a / b);
  r = [q, a - q * b];
end

function r = trunc_divmod (a, b)
  q = fix (a / b);
  r = [q, a - q * b];
end

function r = euclid_divmod (a, b)
  rem_ = mod (a, abs (b));
  r = [(a - rem_) / b, rem_];
end

function d = divisors (n)
  d = find (mod (n, 1:n) == 0);
end

function S = euclid_steps (a, b)
  % Each division of Euclid's algorithm as a row [dividend, quotient, divisor, remainder].
  S = zeros (0, 4);
  while (b ~= 0)
    q = floor (a / b);
    r = a - q * b;
    S(end + 1, :) = [a, q, b, r];
    a = b;
    b = r;
  end
end

function n = step_bound (b)
  % 2(floor(log2 b) + 1). log2 with two outputs writes b = f 2^e with f in
  % [0.5, 1), so e is floor(log2 b) + 1 exactly.
  [~, e] = log2 (b);
  n = 2 * e;
end

function steps = euclid_count (A, B)
  % The number of divisions Euclid's algorithm takes on every pair at once.
  steps = zeros (size (A));
  while (any (B(:) > 0))
    k = B > 0;
    steps(k) += 1;
    R = A;
    R(k) = mod (A(k), B(k));
    A(k) = B(k);
    B(k) = R(k);
  end
end

function R = extended_rows (a, b)
  % Rows [r, q, s, t] of the extended algorithm, with r = s a + t b on every row.
  R = [a, 0, 1, 0; b, 0, 0, 1];
  while (R(end, 1) ~= 0)
    q = floor (R(end - 1, 1) / R(end, 1));
    R(end, 2) = q;
    R(end + 1, :) = [R(end - 1, 1) - q * R(end, 1), 0, R(end - 1, 3) - q * R(end, 3), ...
                     R(end - 1, 4) - q * R(end, 4)];
  end
end

function g = bezout (a, b)
  % [gcd, s, t] with gcd = s a + t b, read from the last nonzero row.
  R = extended_rows (a, b);
  g = R(end - 1, [1, 3, 4]);
end

function x = inverse (a, n)
  % The inverse of a modulo n in 0..n-1, or -1 when gcd(a, n) > 1.
  if (mod (a, n) == 0)
    g = [n, 0, 0];
  else
    g = bezout (mod (a, n), n);
  end
  if (g(1) ~= 1)
    x = -1;
  else
    x = mod (g(2), n);
  end
end

function s = inverse_text (a, n)
  x = inverse (a, n);
  if (x < 0)
    s = 'none';
  else
    s = num (x);
  end
end

function r = power_mod (a, e, n)
  % a^e mod n by repeated squaring, one squaring per bit of e.
  r = mod (1, n);
  sq = mod (a, n);
  while (e > 0)
    if (mod (e, 2) == 1)
      r = mod (r * sq, n);
    end
    sq = mod (sq * sq, n);
    e = floor (e / 2);
  end
end

function out = squares (a, count, n)
  out = zeros (1, count);
  x = mod (a, n);
  for k = 1:count
    out(k) = x;
    x = mod (x * x, n);
  end
end

function u = units (n)
  if (n > 1)
    u = find (gcd (1:n, n) == 1);
  else
    u = 1;
  end
end

function p = phi (n)
  p = numel (units (n));
end

function t = is_prime (n)
  t = n >= 2 && all (mod (n, 2:floor (sqrt (n))) ~= 0);
end

function x = crt (a, m, b, n)
  g = bezout (m, n);
  x = mod (a * g(3) * n + b * g(2) * m, m * n);
end

function x = solutions (c, d, n)
  x = 0:n - 1;
  x = x(mod (c * x - d, n) == 0);
end

function k = order (a, n)
  k = 1;
  x = mod (a, n);
  while (x ~= 1)
    x = mod (x * a, n);
    k += 1;
  end
end

function t = is_field (n)
  T = mod ((1:n - 1)' * (1:n - 1), n);
  t = n >= 2 && all (any (T == 1, 2));
end

function T = mult_table (n)
  T = mod ((0:n - 1)' * (0:n - 1), n);
end

function t = isbn10_terms (digits)
  t = (10:-1:1) .* (digits - '0');
end

function y = luhn_double (x)
  y = 2 * x - 9 * (x >= 5);
end

function t = luhn_terms (number)
  % The Luhn terms from the right: digit 1 as it is, digit 2 doubled, and so on.
  d = fliplr (number - '0');
  t = d;
  t(2:2:end) = luhn_double (d(2:2:end));
end

function f = fibonacci (k)
  a = 0;
  b = 1;
  for i = 1:k
    [a, b] = deal (b, a + b);
  end
  f = a;
end


% The checks -------------------------------------------------------------------

function check_division ()
  section ('Division with remainder');
  sl ('17 divided by 5, quotient and remainder', floor_divmod (17, 5), [3, 2]);
  cases = {[17, 5], [3, 2], [3, 2], [3, 2];
           [-17, 5], [-4, 3], [-3, -2], [-4, 3];
           [17, -5], [-4, -3], [-3, 2], [-3, 2];
           [-17, -5], [3, -2], [3, -2], [4, 3]};
  for i = 1:rows (cases)
    a = cases{i, 1}(1);
    b = cases{i, 1}(2);
    sl (sprintf ('%d by %d, floored q and r', a, b), floor_divmod (a, b), cases{i, 2});
    sl (sprintf ('%d by %d, truncated q and r', a, b), trunc_divmod (a, b), cases{i, 3});
    sl (sprintf ('%d by %d, Euclidean q and r', a, b), euclid_divmod (a, b), cases{i, 4});
  end
  sl ('-17 by 5 with the built-in floored operators', [floor(-17 / 5), mod(-17, 5)], [-4, 3]);
  sn ('C-style remainder of -17 by 5 (built-in)', rem (-17, 5), -2);
  unique_ = true;
  q = -60:60;
  for a = -50:50
    for b = 1:12
      r = a - q * b;
      unique_ = unique_ && sum (r >= 0 & r < b) == 1;
    end
  end
  require ('one q and r with 0 <= r < b, a in -50..50, b in 1..12', unique_);
  lemma = true;
  for a = -30:30
    for b = [-30:-1, 1:30]
      rf = floor_divmod (a, b)(2);
      rt = trunc_divmod (a, b)(2);
      if (rt == 0 || (rt > 0) == (b > 0))
        expect = rt;
      else
        expect = rt + b;
      end
      lemma = lemma && rf == expect;
    end
  end
  require ('floored r = truncated r, plus b when the signs differ', lemma);
end

function check_gcd ()
  section ('The greatest common divisor');
  sl ('divisors of 1147', divisors (1147), [1, 31, 37, 1147]);
  sl ('divisors of 899', divisors (899), [1, 29, 31, 899]);
  d = divisors (1147);
  sl ('common divisors of 1147 and 899', d(mod (899, d) == 0), [1, 31]);
  S = euclid_steps (1147, 899);
  sl ('Euclid on 1147 and 899, quotients', S(:, 2)', [1, 3, 1, 1, 1, 2]);
  sl ('Euclid on 1147 and 899, remainders', S(:, 4)', [248, 155, 93, 62, 31, 0]);
  sn ('gcd(1147, 899)', gcd (1147, 899), 31);
  sn ('divisions taken', rows (S), 6);
  sl ('1147 / 31 and 899 / 31', [1147 / 31, 899 / 31], [37, 29]);
  sn ('step bound 2(floor(log2 899) + 1)', step_bound (899), 20);
  halving = true;
  for ab = [1147, 899; 4181, 2584]'
    S = euclid_steps (ab(1), ab(2));
    rs = [ab(1), ab(2), S(:, 4)'];
    halving = halving && all (2 * rs(3:end) < rs(1:end - 2));
  end
  require ('each remainder under half the one two places before', halving);
  [A, B] = meshgrid (1:299, 1:299);
  [~, e] = log2 (B);
  require ('divisions <= 2(floor(log2 b) + 1) for a, b below 300', all (euclid_count (A, B)(:) <= 2 * e(:)));
  [A, B] = meshgrid (2:233, 1:232);
  keep = B < A;
  A = A(keep);
  B = B(keep);
  N = euclid_count (A, B);
  lame = true;
  for n = unique (N)'
    a = A(N == n);
    b = B(N == n);
    f2 = fibonacci (n + 2);
    lame = lame && min (a) == f2 && isequal (b(a == f2), fibonacci (n + 1));
  end
  require ('smallest pair needing n divisions is F(n+2), F(n+1)', lame);
  sl ('Fibonacci numbers F(1) to F(12)', arrayfun (@fibonacci, 1:12), [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]);
end

function check_bezout ()
  section ('Bezout''s identity and the extended algorithm');
  sl ('11 * 1147 and 14 * 899', [11 * 1147, 14 * 899], [12617, 12586]);
  sn ('11 * 1147 - 14 * 899', 11 * 1147 - 14 * 899, 31);
  R = extended_rows (1147, 899);
  sl ('extended table, r', R(:, 1)', [1147, 899, 248, 155, 93, 62, 31, 0]);
  sl ('extended table, q (0 where none)', R(:, 2)', [0, 1, 3, 1, 1, 1, 2, 0]);
  sl ('extended table, s', R(:, 3)', [1, 0, 1, -3, 4, -7, 11, -29]);
  sl ('extended table, t', R(:, 4)', [0, 1, -1, 4, -5, 9, -14, 37]);
  require ('r = 1147 s + 899 t on every row', all (R(:, 1) == 1147 * R(:, 3) + 899 * R(:, 4)));
  sn ('-29 * 1147 + 37 * 899', -29 * 1147 + 37 * 899, 0);
  [Sg, Tg] = meshgrid (-60:60, -60:60);
  V = 1147 * Sg + 899 * Tg;
  sn ('smallest positive 1147 s + 899 t, |s|, |t| <= 60', min (V(V > 0)), 31);
  % The extended algorithm run on every pair at once.
  [A, B] = meshgrid (1:79, 1:79);
  r0 = A; r1 = B; s0 = ones (size (A)); s1 = zeros (size (A)); t0 = s1; t1 = s0;
  while (any (r1(:) ~= 0))
    k = r1 ~= 0;
    q = zeros (size (A));
    q(k) = floor (r0(k) ./ r1(k));
    [r0(k), r1(k)] = deal (r1(k), r0(k) - q(k) .* r1(k));
    [s0(k), s1(k)] = deal (s1(k), s0(k) - q(k) .* s1(k));
    [t0(k), t1(k)] = deal (t1(k), t0(k) - q(k) .* t1(k));
  end
  require ('bezout gives gcd = s a + t b for a, b below 80', ...
           all (r0(:) == gcd (A(:), B(:))) && all (r0(:) == s0(:) .* A(:) + t0(:) .* B(:)));
  sl ('36 mod 6, 4 mod 6, 9 mod 6', mod ([36, 4, 9], 6), [0, 4, 3]);
  [A, B] = meshgrid (0:59, 0:59);
  ok = true;
  for p = 2:39
    if (is_prime (p))
      ok = ok && all (mod (A(:) .* B(:), p) ~= 0 | mod (A(:), p) == 0 | mod (B(:), p) == 0);
    end
  end
  require ('p | ab forces p | a or p | b, primes below 40', ok);
end

function check_congruences ()
  section ('Congruences');
  sl ('9 + 5 and its remainder modulo 12', [9 + 5, mod(9 + 5, 12)], [14, 2]);
  sn ('-17 - 3', -17 - 3, -20);
  sl ('remainders of -17 and 3 modulo 5', mod ([-17, 3], 5), [3, 3]);
  sl ('37 mod 7 and 52 mod 7', mod ([37, 52], 7), [2, 3]);
  sl ('37 * 52 and its remainder modulo 7', [37 * 52, mod(37 * 52, 7)], [1924, 6]);
  sn ('1924 - 6', 1924 - 6, 274 * 7);
  respects = true;
  for n = 1:12
    [Av, Cv] = meshgrid (-12:12, -12:12);
    k = mod (Av - Cv, n) == 0;
    a = Av(k); c = Cv(k);
    [Bv, Dv] = meshgrid (-6:6, -6:6);
    k = mod (Bv - Dv, n) == 0;
    b = Bv(k)'; d = Dv(k)';
    sums = (a + b) - (c + d);
    prods = a * b - c * d;
    respects = respects && all (mod (sums(:), n) == 0) && all (mod (prods(:), n) == 0);
  end
  require ('congruence respects + and x for n up to 12', respects);
  same ('13 in binary', dec2bin (13), '1101');
  sl ('3, 3^2, 3^4, 3^8 modulo 11', squares (3, 4, 11), [3, 9, 4, 5]);
  sn ('3^8 * 3^4 * 3 modulo 11', mod (5 * 4 * 3, 11), 5);
  sn ('5 * 4 * 3', 5 * 4 * 3, 60);
  sn ('3^13', 3 ^ 13, 1594323);
  sn ('3^13 mod 11', power_mod (3, 13, 11), 5);
  terms = isbn10_terms ('0130426725');
  sl ('ISBN-10 terms of 0-13-042672-5', terms, [0, 9, 24, 0, 24, 10, 24, 21, 4, 5]);
  sl ('ISBN-10 weighted sum and its remainder modulo 11', [sum(terms), mod(sum (terms), 11)], [121, 0]);
  lt = luhn_terms ('3141592653');
  sl ('Luhn sum of 3141592653 and its remainder modulo 10', [sum(lt), mod(sum (lt), 10)], [40, 0]);
  sl ('2 * 3 mod 6 and 2 * 0 mod 6', [mod(2 * 3, 6), mod(2 * 0, 6)], [0, 0]);
  cancel = true;
  for n = 2:30
    [X, Y] = meshgrid (0:n - 1, 0:n - 1);
    for c = 1:n - 1
      if (gcd (c, n) == 1)
        cancel = cancel && all (mod (c * X(:) - c * Y(:), n) ~= 0 | mod (X(:) - Y(:), n) == 0);
      end
    end
  end
  require ('c x = c y forces x = y when gcd(c, n) = 1, n up to 30', cancel);
end

function check_inverses ()
  section ('Inverses modulo n');
  sn ('3 * 5 mod 7', mod (3 * 5, 7), 1);
  sl ('inverses of 1 to 6 modulo 7', arrayfun (@(a) inverse (a, 7), 1:6), [1, 4, 5, 2, 3, 6]);
  S = euclid_steps (37, 29);
  sl ('Euclid on 37 and 29, quotients', S(:, 2)', [1, 3, 1, 1, 1, 2]);
  sl ('Euclid on 37 and 29, remainders', S(:, 4)', [8, 5, 3, 2, 1, 0]);
  sl ('gcd, s, t for 37 and 29', bezout (37, 29), [1, 11, -14]);
  sn ('11 * 37 - 14 * 29', 11 * 37 - 14 * 29, 1);
  sn ('-14 mod 37', mod (-14, 37), 23);
  sl ('29 * 23 and 18 * 37', [29 * 23, 18 * 37], [667, 666]);
  sn ('inverse of 29 modulo 37', inverse (29, 37), 23);
  [~, s] = gcd (29, 37);
  sn ('inverse of 29 modulo 37 (built-in)', mod (s, 37), 23);
  terms = isbn10_terms ('0130426725');
  sn ('ISBN-10 sum without position 5', sum (terms) - terms(5), 97);
  sn ('-97 mod 11', mod (-97, 11), 2);
  sn ('inverse of 6 modulo 11', inverse (6, 11), 2);
  sn ('recovered digit 2 * 2 mod 11', mod (2 * 2, 11), 4);
  inv6 = arrayfun (@(a) inverse (a, 6), 1:5);
  sl ('Z_6 elements with inverses', find (inv6 >= 0), [1, 5]);
  sn ('5 * 5 mod 6', mod (5 * 5, 6), 1);
  sl ('2 * 3 mod 6 and 3 * 4 mod 6', [mod(2 * 3, 6), mod(3 * 4, 6)], [0, 0]);
  sn ('gcd(899, 1147)', gcd (899, 1147), 31);
  same ('inverse of 899 modulo 1147', inverse_text (899, 1147), 'none');
  exists = true;
  for n = 2:40
    counts = sum (mod ((0:n - 1)' * (0:n - 1), n) == 1, 2)';
    exists = exists && all ((counts == 1) == (gcd (0:n - 1, n) == 1)) && all (counts <= 1);
  end
  require ('one inverse when gcd(a, n) = 1, none otherwise, n up to 40', exists);
end

function check_fermat_euler ()
  section ('Fermat''s and Euler''s theorems');
  P = zeros (6, 6);
  for a = 1:6
    for k = 1:6
      P(a, k) = power_mod (a, k, 7);
    end
  end
  sl ('powers of 1 modulo 7', P(1, :), [1, 1, 1, 1, 1, 1]);
  sl ('powers of 2 modulo 7', P(2, :), [2, 4, 1, 2, 4, 1]);
  sl ('powers of 3 modulo 7', P(3, :), [3, 2, 6, 4, 5, 1]);
  sl ('powers of 4 modulo 7', P(4, :), [4, 2, 1, 4, 2, 1]);
  sl ('powers of 5 modulo 7', P(5, :), [5, 4, 6, 2, 3, 1]);
  sl ('powers of 6 modulo 7', P(6, :), [6, 1, 6, 1, 6, 1]);
  sl ('3k mod 7 for k = 1 to 6', mod (3 * (1:6), 7), [3, 6, 2, 5, 1, 4]);
  sl ('6! and its remainder modulo 7', [factorial(6), mod(factorial (6), 7)], [720, 6]);
  fermat = true;
  little = true;
  for p = 2:59
    if (is_prime (p))
      for a = 1:p - 1
        fermat = fermat && power_mod (a, p - 1, p) == 1;
      end
      for a = -20:39
        little = little && power_mod (a, p, p) == mod (a, p);
      end
    end
  end
  require ('a^(p-1) = 1 mod p for primes below 60', fermat);
  require ('a^p = a mod p for every a, primes below 60', little);
  sl ('units modulo 12', units (12), [1, 5, 7, 11]);
  sl ('phi(12), phi(7), phi(9), phi(8)', [phi(12), phi(7), phi(9), phi(8)], [4, 6, 6, 4]);
  sl ('units modulo 9', units (9), [1, 2, 4, 5, 7, 8]);
  prime_power = true;
  for p = 2:1000
    if (is_prime (p))
      for k = 1:10
        if (p ^ k <= 1000)
          prime_power = prime_power && phi (p ^ k) == p ^ k - p ^ (k - 1);
        end
      end
    end
  end
  require ('phi(p^k) = p^k - p^(k-1) for p^k up to 1000', prime_power);
  sl ('2^6 and its remainder modulo 9', [2 ^ 6, mod(2 ^ 6, 9)], [64, 1]);
  sl ('squares of 1, 5, 7, 11 modulo 12', mod (units (12) .^ 2, 12), [1, 1, 1, 1]);
  euler = true;
  for n = 2:60
    for a = units (n)
      euler = euler && power_mod (a, phi (n), n) == 1;
    end
  end
  require ('a^phi(n) = 1 mod n when gcd(a, n) = 1, n up to 60', euler);
  sl ('100 = 6q + r', [floor(100 / 6), mod(100, 6)], [16, 4]);
  sl ('3^4 and its remainder modulo 7', [3 ^ 4, mod(3 ^ 4, 7)], [81, 4]);
  sn ('3^100 mod 7', power_mod (3, 100, 7), 4);
  sl ('6, 6^2, 6^4, 6^8 modulo 11', squares (6, 4, 11), [6, 3, 9, 4]);
  sl ('6^9 and its remainder modulo 11', [6 ^ 9, mod(6 ^ 9, 11)], [10077696, 2]);
  sn ('11 * 31', 11 * 31, 341);
  sl ('2^10 and its remainder modulo 341', [2 ^ 10, mod(2 ^ 10, 341)], [1024, 1]);
  sn ('3 * 341', 3 * 341, 1023);
  sn ('2^340 mod 341', power_mod (2, 340, 341), 1);
  sn ('3^340 mod 341', power_mod (3, 340, 341), 56);
end

function check_crt ()
  section ('The Chinese remainder theorem');
  sl ('x mod 3 for x = 3, 8, 13', mod ([3, 8, 13], 3), [0, 2, 1]);
  sn ('x = 2 mod 3 and x = 3 mod 5', crt (2, 3, 3, 5), 8);
  x = 0:14;
  require ('each pair (x mod 3, x mod 5) once for x = 0 to 14', isequal (sort (mod (x, 3) * 5 + mod (x, 5)), 0:14));
  sl ('gcd, s, t for 3 and 5', bezout (3, 5), [1, 2, -1]);
  sn ('2 * (-1) * 5 + 3 * 2 * 3', 2 * (-1) * 5 + 3 * 2 * 3, 8);
  sl ('gcd, s, t for 4 and 9', bezout (4, 9), [1, -2, 1]);
  sn ('3 * 1 * 9 + 5 * (-2) * 4', 3 * 1 * 9 + 5 * (-2) * 4, -13);
  sn ('-13 mod 36', mod (-13, 36), 23);
  sl ('23 mod 4 and 23 mod 9', mod (23, [4, 9]), [3, 5]);
  sn ('x = 3 mod 4 and x = 5 mod 9', crt (3, 4, 5, 9), 23);
  bijective = true;
  for m = 1:15
    for n = 1:15
      if (gcd (m, n) == 1)
        x = 0:m * n - 1;
        bijective = bijective && isequal (sort (mod (x, m) * n + mod (x, n)), x);
        g = bezout (m, n);
        [a, b] = meshgrid (0:m - 1, 0:n - 1);
        X = mod (a * g(3) * n + b * g(2) * m, m * n);
        bijective = bijective && all (mod (X(:), m) == a(:)) && all (mod (X(:), n) == b(:));
      end
    end
  end
  require ('x -> (x mod m, x mod n) one to one when gcd(m, n) = 1', bijective);
  lemma = true;
  x = 0:399;
  for m = 1:12
    for n = 1:12
      if (gcd (m, n) == 1)
        hit = x(mod (x, m) == 0 & mod (x, n) == 0);
        lemma = lemma && all (mod (hit, m * n) == 0);
      end
    end
  end
  require ('m | x and n | x give mn | x when gcd(m, n) = 1', lemma);
  sl ('phi(4), phi(3), phi(12)', [phi(4), phi(3), phi(12)], [2, 2, 4]);
  sl ('phi(4), phi(9), phi(36)', [phi(4), phi(9), phi(36)], [2, 6, 12]);
  ph = arrayfun (@phi, 1:900);
  mult = true;
  for m = 1:30
    for n = 1:30
      if (gcd (m, n) == 1)
        mult = mult && ph(m * n) == ph(m) * ph(n);
      end
    end
  end
  require ('phi(mn) = phi(m) phi(n) when gcd(m, n) = 1, up to 30', mult);
  x = 0:11;
  sl ('x = 1 mod 4 and x = 2 mod 6, x in 0..11', x(mod (x, 4) == 1 & mod (x, 6) == 2), []);
  sl ('x = 1 mod 4 and x = 3 mod 6, x in 0..11', x(mod (x, 4) == 1 & mod (x, 6) == 3), [9]);
end

function check_fields ()
  section ('The integers modulo a prime');
  sl ('elements of Z_7 and Z_6 with inverses', [phi(7), phi(6)], [6, 2]);
  T7 = mult_table (7);
  require ('each nonzero row of Z_7 holds 1 exactly once', all (sum (T7(2:end, :) == 1, 2) == 1));
  T6 = mult_table (6);
  pairs = {};
  for a = 1:5
    for b = 1:5
      if (T6(a + 1, b + 1) == 0)
        pairs{end + 1} = lst ([a, b]);
      end
    end
  end
  same ('pairs of nonzero elements of Z_6 with product 0', ['[' strjoin(pairs, ', ') ']'], ...
        '[[2, 3], [3, 2], [3, 4], [4, 3]]');
  no_one = [];
  for a = 1:5
    if (~any (T6(a + 1, :) == 1))
      no_one(end + 1) = a;
    end
  end
  sl ('rows of Z_6 holding no 1', no_one, [2, 3, 4]);
  field = true;
  for n = 2:60
    field = field && is_field (n) == is_prime (n);
  end
  require ('Z_n is a field exactly when n is prime, n = 2 to 60', field);
  [A, B] = meshgrid (0:1, 0:1);
  require ('a + b mod 2 is a XOR b for all four pairs', isequal (mod (A + B, 2), double (xor (A, B))));
  require ('a b mod 2 is a AND b for all four pairs', isequal (mod (A .* B, 2), double (A & B)));
  sl ('1 + 0 + 1 + 1 and its remainder modulo 2', [1 + 0 + 1 + 1, mod(1 + 0 + 1 + 1, 2)], [3, 1]);
  sn ('2 * 128 mod 256', mod (2 * 128, 256), 0);
  sn ('phi(256)', phi (256), 128);
  require ('257 is prime', is_prime (257));
end

function check_examples ()
  section ('Worked examples');
  S = euclid_steps (4181, 2584);
  sl ('F(19) and F(18)', [fibonacci(19), fibonacci(18)], [4181, 2584]);
  sn ('Euclid on 4181 and 2584, divisions', rows (S), 17);
  sl ('Euclid on 4181 and 2584, quotients', S(:, 2)', [ones(1, 16), 2]);
  sl ('Euclid on 4181 and 2584, remainders', S(:, 4)', [arrayfun(@fibonacci, 17:-1:2), 0]);
  sn ('gcd(4181, 2584)', gcd (4181, 2584), 1);
  sn ('step bound 2(floor(log2 2584) + 1)', step_bound (2584), 24);
  sl ('Luhn doubling of 0 to 9', luhn_double (0:9), [0, 2, 4, 6, 8, 1, 3, 5, 7, 9]);
  lt = luhn_terms ('3141592653');
  sl ('Luhn terms of 3141592653 from the right', lt, [3, 1, 6, 4, 9, 1, 1, 8, 1, 6]);
  sn ('Luhn sum without digit 2 from the right', sum (lt) - lt(2), 39);
  sn ('-39 mod 10', mod (-39, 10), 1);
  x = 0:9;
  sl ('digits x with D(x) = 1', x(luhn_double (x) == 1), [5]);
  sl ('solutions of 2x = 1 mod 10', solutions (2, 1, 10), []);
  x = 0:8;
  require ('D(x) = 2x mod 9 for x = 0 to 8', all (mod (luhn_double (x), 9) == mod (2 * x, 9)));
  sn ('inverse of 2 modulo 9', inverse (2, 9), 5);
  sn ('gcd(6, 10)', gcd (6, 10), 2);
  sl ('solutions of 6x = 4 mod 10', solutions (6, 4, 10), [4, 9]);
  sn ('inverse of 3 modulo 5', inverse (3, 5), 2);
  sn ('2 * 2 mod 5', mod (2 * 2, 5), 4);
  sl ('6 * 4 and 6 * 9', [6 * 4, 6 * 9], [24, 54]);
  sl ('solutions of 6x = 5 mod 10', solutions (6, 5, 10), []);
  sn ('x = 1 mod 4 and x = 2 mod 5', crt (1, 4, 2, 5), 17);
  sl ('17 mod 7 and 20 mod 7', mod ([17, 20], 7), [3, 6]);
  sn ('inverse of 6 modulo 7', inverse (6, 7), 6);
  sl ('k from 6k = 2 mod 7', solutions (6, 2, 7), [5]);
  sn ('17 + 20 * 5', 17 + 20 * 5, 117);
  sl ('117 mod 4, 5, 7', mod (117, [4, 5, 7]), [1, 2, 5]);
  sn ('x = 1 mod 4, 2 mod 5, 5 mod 7', crt (17, 20, 5, 7), 117);
  same ('117 in binary', dec2bin (117), '1110101');
  sl ('5^(2^j) mod 19 for j = 0 to 6', squares (5, 7, 19), [5, 6, 17, 4, 16, 9, 5]);
  sl ('5 * 9, 7 * 16, 17 * 17, 4 * 5 modulo 19', mod ([5 * 9, 7 * 16, 17 * 17, 4 * 5], 19), [7, 17, 4, 1]);
  sn ('5^117 mod 19', power_mod (5, 117, 19), 1);
  sl ('117 = 18q + r', [floor(117 / 18), mod(117, 18)], [6, 9]);
  sn ('5^9 mod 19', power_mod (5, 9, 19), 1);
  sn ('smallest k with 5^k = 1 mod 19', order (5, 19), 9);
end

function check_exercises ()
  section ('Exercises');
  sl ('-23 by 6, floored q and r', floor_divmod (-23, 6), [-4, 1]);
  sl ('-23 by 6, truncated q and r', trunc_divmod (-23, 6), [-3, -5]);
  S = euclid_steps (252, 198);
  sl ('Euclid on 252 and 198, quotients', S(:, 2)', [1, 3, 1, 2]);
  sl ('Euclid on 252 and 198, remainders', S(:, 4)', [54, 36, 18, 0]);
  sl ('gcd, s, t for 252 and 198', bezout (252, 198), [18, 4, -5]);
  sn ('4 * 252 - 5 * 198', 4 * 252 - 5 * 198, 18);
  sn ('inverse of 7 modulo 26', inverse (7, 26), 15);
  sl ('7 * 15 and 4 * 26', [7 * 15, 4 * 26], [105, 104]);
  same ('gcd(13, 26) and the inverse of 13 modulo 26', ...
        sprintf ('[%d, ''%s'']', gcd (13, 26), inverse_text (13, 26)), '[13, ''none'']');
  terms = isbn10_terms ('0130426725');
  sn ('ISBN-10 sum without position 8', sum (terms) - terms(8), 100);
  sn ('-100 mod 11', mod (-100, 11), 10);
  sn ('inverse of 3 modulo 11', inverse (3, 11), 4);
  sn ('4 * 10 mod 11', mod (4 * 10, 11), 7);
  sl ('100 = 12q + r', [floor(100 / 12), mod(100, 12)], [8, 4]);
  sn ('2^100 mod 13', power_mod (2, 100, 13), 3);
  sl ('phi(4), phi(25), phi(100)', [phi(4), phi(25), phi(100)], [2, 20, 40]);
  sn ('2026 mod 40', mod (2026, 40), 26);
  sl ('3, 3^2, 3^4, 3^8, 3^16 modulo 100', squares (3, 5, 100), [3, 9, 81, 61, 21]);
  sl ('21 * 61 mod 100 and 81 * 9 mod 100', [mod(21 * 61, 100), mod(81 * 9, 100)], [81, 29]);
  sl ('3^26 mod 100 and 3^2026 mod 100', [power_mod(3, 26, 100), power_mod(3, 2026, 100)], [29, 29]);
  sn ('inverse of 7 modulo 11', inverse (7, 11), 8);
  sl ('k from 7k = 2 mod 11', solutions (7, 2, 11), [5]);
  sn ('4 + 7 * 5', 4 + 7 * 5, 39);
  sn ('x = 4 mod 7 and x = 6 mod 11', crt (4, 7, 6, 11), 39);
  u = units (9);
  same ('units modulo 9 and their inverses', ...
        ['[' lst(u) ', ' lst(arrayfun (@(a) inverse (a, 9), u)) ']'], ...
        '[[1, 2, 4, 5, 7, 8], [1, 5, 7, 2, 4, 8]]');
  sn ('3 * 3 mod 9', mod (3 * 3, 9), 0);
  sn ('3 * 11 * 17', 3 * 11 * 17, 561);
  sl ('560 mod 2, 10, 16', mod (560, [2, 10, 16]), [0, 0, 0]);
  sn ('phi(561)', phi (561), 320);
  u = units (561);
  require ('a^560 = 1 mod 561 for all 320 a coprime to 561', all (arrayfun (@(a) power_mod (a, 560, 561), u) == 1));
end


check_division ();
check_gcd ();
check_bezout ();
check_congruences ();
check_inverses ();
check_fermat_euler ();
check_crt ();
check_fields ();
check_examples ();
check_exercises ();
printf ('\nAll chapter 13 numbers reproduced.\n');
