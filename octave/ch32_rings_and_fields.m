% Reproduce every number in chapter 32 of Error Correction Coding,
% "Rings and Fields", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch32_rings_and_fields.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch32_rings_and_fields.py and prints the same
% lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. All the arithmetic is on whole numbers small enough for a double
% to hold exactly, so every comparison is exact. Every finite ring the chapter
% uses is held as an addition table and a multiplication table, and the ring
% axioms are checked on every pair and every triple of elements, one row of
% triples at a time: where the Python looks rows up with bytes.translate, this
% file indexes one table by another, so the two routes check each other as
% well as the chapter.

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

function s = yes (flag)
  if (flag)
    s = 'yes';
  else
    s = 'no';
  end
end


% Whole numbers ----------------------------------------------------------------

function t = is_prime_power (n)
  t = n >= 2 && numel (unique (factor (n))) == 1;
end

function p = phi (n)
  p = sum (gcd (1:n, n) == 1);
end

function row = pascal_row (n)
  row = 1;
  for k = 1:n
    row = [1, row(1:end - 1) + row(2:end), 1];
  end
end


% The finite rings of the chapter ----------------------------------------------
% Each ring is a struct holding its addition and multiplication tables on the
% elements 0 .. n-1, with 0 the zero element. Entry (a+1, b+1) is a + b or ab.

function F = f4_mul ()
  F = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
end

function r = ring_zn (n)
  v = 0:n - 1;
  r = struct ('name', sprintf ('Z_%d', n), 'kind', 'zn', 'n', n, 'k', 0, 'one', 1, ...
              'add', mod (v' + v, n), 'mul', mod (v' * v, n));
end

function r = ring_words (k)
  n = 2 ^ k;
  [X, Y] = meshgrid (0:n - 1);
  r = struct ('name', sprintf ('%d-bit words', k), 'kind', 'words', 'n', n, 'k', k, 'one', n - 1, ...
              'add', bitxor (Y, X), 'mul', bitand (Y, X));
end

function z = mat_mul (x, y)
  % Multiply 2 x 2 matrices over GF(2), each packed as the bits a b c d of [[a, b], [c, d]].
  a = bitand (bitshift (x, -3), 1); b = bitand (bitshift (x, -2), 1);
  c = bitand (bitshift (x, -1), 1); d = bitand (x, 1);
  e = bitand (bitshift (y, -3), 1); f = bitand (bitshift (y, -2), 1);
  g = bitand (bitshift (y, -1), 1); h = bitand (y, 1);
  z = mod (a .* e + b .* g, 2) * 8 + mod (a .* f + b .* h, 2) * 4 ...
      + mod (c .* e + d .* g, 2) * 2 + mod (c .* f + d .* h, 2);
end

function r = ring_matrices ()
  [X, Y] = meshgrid (0:15);
  r = struct ('name', 'M2(GF(2))', 'kind', 'mat', 'n', 16, 'k', 0, 'one', 9, ...
              'add', bitxor (Y, X), 'mul', mat_mul (Y, X));
end

function r = ring_f4 ()
  [X, Y] = meshgrid (0:3);
  r = struct ('name', 'F_4', 'kind', 'f4', 'n', 4, 'k', 0, 'one', 1, ...
              'add', bitxor (Y, X), 'mul', f4_mul ());
end

function s = show (r, a)
  switch (r.kind)
    case 'words'
      s = dec2bin (a, r.k);
    case 'mat'
      s = sprintf ('%d%d|%d%d', bitand (bitshift (a, -3), 1), bitand (bitshift (a, -2), 1), ...
                   bitand (bitshift (a, -1), 1), bitand (a, 1));
    case 'f4'
      names = '01ab';
      s = names(a + 1);
    otherwise
      s = num (a);
  end
end

function s = show_all (r, v)
  if (isempty (v))
    s = 'none';
  else
    s = strjoin (arrayfun (@(a) show (r, a), v, 'UniformOutput', false), ' ');
  end
end

function R = chapter_rings ()
  R = {ring_zn(4), ring_zn(12), ring_zn(13), ring_zn(256), ring_words(2), ring_words(8), ...
       ring_matrices(), ring_f4()};
end

function ok = associative (T)
  % (ab)c against a(bc) for every b and c, one a at a time.
  n = rows (T);
  ok = true;
  for a = 1:n
    ta = T(a, :);
    if (~isequal (T(ta + 1, :), ta(T + 1)))
      ok = false;
      return;
    end
  end
end

function ok = distributive (A, M)
  n = rows (A);
  ok = true;
  for a = 1:n
    ma = M(a, :);
    col = M(:, a)';
    if (~isequal (ma(A + 1), A(ma + 1, ma + 1)) || ~isequal (col(A + 1), A(col + 1, col + 1)))
      ok = false;
      return;
    end
  end
end

function t = commutative (T)
  t = isequal (T, T');
end

function u = units (r)
  U = (r.mul == r.one) & (r.mul' == r.one);
  u = find (any (U, 2))' - 1;
end

function v = inverse_in (r, u)
  U = (r.mul(u + 1, :) == r.one) & (r.mul(:, u + 1)' == r.one);
  v = find (U, 1) - 1;
end

function z = zero_divisors (r)
  Z = (r.mul == 0) | (r.mul' == 0);
  Z(1, :) = false;
  Z(:, 1) = false;
  z = find (any (Z, 2))' - 1;
end

function k = characteristic (r)
  s = r.one;
  k = 1;
  while (s ~= 0)
    s = r.add(s + 1, r.one + 1);
    k += 1;
  end
end

function m = multiples_of_one (r)
  m = 0;
  s = r.one;
  while (s ~= 0)
    m(end + 1) = s;
    s = r.add(s + 1, r.one + 1);
  end
  m = sort (m);
end

function u = zn_units (n)
  M = mod ((1:n - 1)' * (0:n - 1), n);
  u = find (any (M == 1, 2))';
end

function z = zn_zero_divisors (n)
  M = mod ((1:n - 1)' * (1:n - 1), n);
  z = find (any (M == 0, 2))';
end

function c = unit_and_zero_divisor_counts (n)
  c = [numel(zn_units (n)), numel(zn_zero_divisors (n))];
end

function x = solutions (a, c, n)
  v = 0:n - 1;
  x = v(mod (a * v, n) == mod (c, n));
end

function m = times_map (a, n)
  m = mod (a * (0:n - 1), n);
end


% Polynomials ------------------------------------------------------------------
% Integer coefficient lists are written constant term first. A polynomial over
% GF(2) is a whole number whose bit i is the coefficient of x^i.

function r = clmul (a, b)
  % Carry-less product, elementwise over arrays of the same size.
  r = zeros (size (a));
  for i = 0:15
    bit = bitand (bitshift (b, -i), 1);
    r = bitxor (r, bit .* bitshift (a, i));
  end
end

function r = spread (a)
  r = zeros (size (a));
  for i = 0:15
    r = bitor (r, bitshift (bitand (bitshift (a, -i), 1), 2 * i));
  end
end

function d = degree (a)
  [~, e] = log2 (a);
  d = e - 1;
end

function s = bits (a)
  s = dec2bin (a);
end

function r = powmod (A, e, m)
  r = ones (size (A));
  for i = 1:e
    r = mod (r .* A, m);
  end
end


% The checks, in chapter order -------------------------------------------------

function check_opening ()
  section ('The opening: two rings inside an 8-bit processor');
  sn ('3 * 171', 3 * 171, 513);
  sn ('2 * 256 + 1', 2 * 256 + 1, 513);
  sn ('3 * 171 mod 256', mod (3 * 171, 256), 1);
  x = 0:255;
  require ('multiplying by 171 undoes multiplying by 3, all 256 bytes', all (mod (3 * x * 171, 256) == x));
  sn ('2 * 128 mod 256', mod (2 * 128, 256), 0);
  sn ('2 * 0 mod 256', mod (2 * 0, 256), 0);
  sn ('different values of 3x mod 256', numel (unique (mod (3 * x, 256))), 256);
  sn ('different values of 2x mod 256', numel (unique (mod (2 * x, 256))), 128);
  sn ('0x0F AND 0xF0', bitand (15, 240), 0);
  sn ('0x0F XOR 0xF0', bitxor (15, 240), 255);
end

function check_axioms (R, stated)
  section ('The ring axioms, checked on every pair and every triple');
  for i = 1:numel (R)
    r = R{i};
    n = r.n;
    name = r.name;
    v = 0:n - 1;
    require (sprintf ('%s: + associative on all %d triples', name, n ^ 3), associative (r.add));
    require (sprintf ('%s: + commutative, 0 an identity, negatives exist', name), ...
             commutative (r.add) && all (r.add(1, :) == v) && all (any (r.add == 0, 2)));
    require (sprintf ('%s: x associative on all %d triples', name, n ^ 3), associative (r.mul));
    require (sprintf ('%s: both distributive laws on all %d triples', name, n ^ 3), distributive (r.add, r.mul));
    require (sprintf ('%s: %s is a two-sided identity for x', name, show (r, r.one)), ...
             all (r.mul(r.one + 1, :) == v) && all (r.mul(:, r.one + 1)' == v));
    same (sprintf ('%s: number of elements, commutative', name), ...
          sprintf ('%d %s', n, yes (commutative (r.mul))), stated{i});
  end
end

function check_gallery ()
  section ('A gallery of rings');
  sl ('(x^2 + x + 1)(x + 1) in Z[x], constant first', conv ([1 1 1], [1 1]), [1, 2, 2, 1]);
  same ('(x^2 + x + 1)(x + 1) in GF(2)[x], bits', bits (clmul (7, 3)), '1001');
  w8 = ring_words (8);
  require ('a AND a = a for every 8-bit word a', all (diag (w8.mul)' == 0:255));
  require ('a XOR a = 0 for every 8-bit word a', all (diag (w8.add) == 0));
  A = 13;
  B = 11;
  m = ring_matrices ();
  same ('A', show (m, A), '11|01');
  same ('B', show (m, B), '10|11');
  same ('AB over GF(2)', show (m, mat_mul (A, B)), '01|11');
  same ('BA over GF(2)', show (m, mat_mul (B, A)), '11|10');
  sn ('number of 2 x 2 matrices over GF(2)', 2 ^ 4, 16);
end

function check_consequences (R)
  section ('First consequences of the axioms');
  ok_zero = true; ok_neg = true; ok_ids = true; ok_unique_neg = true;
  for i = 1:numel (R)
    r = R{i};
    n = r.n;
    v = 0:n - 1;
    [~, j] = max (r.add == 0, [], 2);
    neg = j' - 1;
    ok_zero = ok_zero && all (r.mul(1, :) == 0) && all (r.mul(:, 1) == 0);
    ok_neg = ok_neg && all (all (r.mul(neg + 1, :) == neg(r.mul + 1))) ...
             && all (all (r.mul(:, neg + 1) == neg(r.mul + 1)));
    zeros_ = find (all (r.add == v, 2))' - 1;
    ones_ = find (all (r.mul == v, 2) & all (r.mul' == v, 2))' - 1;
    ok_ids = ok_ids && isequal (zeros_, 0) && isequal (ones_, r.one);
    ok_unique_neg = ok_unique_neg && all (sum (r.add == 0, 2) == 1);
  end
  require ('0a = a0 = 0 for every a, all eight rings', ok_zero);
  require ('(-a)b = a(-b) = -(ab) for every a and b, all eight rings', ok_neg);
  require ('exactly one zero and one multiplicative identity, all eight', ok_ids);
  require ('every element has exactly one negative, all eight rings', ok_unique_neg);
end

function check_homomorphism ()
  section ('Reduction modulo n is a ring homomorphism');
  [B, A] = meshgrid (-30:30);
  ok_add = true; ok_mul = true;
  for n = 2:12
    ok_add = ok_add && all (all (mod (A + B, n) == mod (mod (A, n) + mod (B, n), n)));
    ok_mul = ok_mul && all (all (mod (A .* B, n) == mod (mod (A, n) .* mod (B, n), n)));
  end
  require ('(a + b) mod n = (a mod n + b mod n) mod n, n = 2..12, |a|, |b| <= 30', ok_add);
  require ('(ab) mod n = (a mod n)(b mod n) mod n, n = 2..12, |a|, |b| <= 30', ok_mul);
  sn ('47 * 86', 47 * 86, 4042);
  sl ('47 mod 9, 86 mod 9', [mod(47, 9), mod(86, 9)], [2, 5]);
  sn ('(2 * 5) mod 9', mod (2 * 5, 9), 1);
  sn ('4042 mod 9', mod (4042, 9), 1);
  sn ('4052 mod 9, a wrong product', mod (4052, 9), 2);
end

function check_units (R)
  section ('Units and zero divisors');
  v = -20:20;
  sl ('units of Z among -20..20', v(v ~= 0 & mod (1, v) == 0), [-1, 1]);
  group_ok = true; disjoint_ok = true; cover_ok = true;
  for i = 1:numel (R)
    r = R{i};
    u = units (r);
    z = zero_divisors (r);
    inv = arrayfun (@(a) inverse_in (r, a), u);
    group_ok = group_ok && all (all (ismember (r.mul(u + 1, u + 1), u))) && all (ismember (inv, u));
    disjoint_ok = disjoint_ok && isempty (intersect (u, z));
    cover_ok = cover_ok && isequal (sort ([u, z]), 1:r.n - 1);
  end
  require ('units closed under x and under inverses, all eight rings', group_ok);
  require ('no unit is a zero divisor, all eight rings', disjoint_ok);
  require ('every nonzero element is a unit or a zero divisor, all eight', cover_ok);

  z12 = ring_zn (12);
  u12 = units (z12);
  sl ('units of Z_12', u12, [1, 5, 7, 11]);
  sl ('zero divisors of Z_12', zero_divisors (z12), [2, 3, 4, 6, 8, 9, 10]);
  sl ('squares of the units of Z_12', mod (u12 .* u12, 12), [1, 1, 1, 1]);
  sl ('4 * 3 mod 12, 6 * 2 mod 12', [mod(4 * 3, 12), mod(6 * 2, 12)], [0, 0]);

  z256 = ring_zn (256);
  u256 = units (z256);
  zd256 = zero_divisors (z256);
  same ('units of Z_256, all odd', sprintf ('%d %s', numel (u256), yes (all (mod (u256, 2) == 1))), '128 yes');
  same ('zero divisors of Z_256, all even', sprintf ('%d %s', numel (zd256), yes (all (mod (zd256, 2) == 0))), '127 yes');
  sn ('inverse of 3 in Z_256', inverse_in (z256, 3), 171);

  w2 = ring_words (2);
  w8 = ring_words (8);
  same ('units of the 2-bit words', show_all (w2, units (w2)), '11');
  same ('zero divisors of the 2-bit words', show_all (w2, zero_divisors (w2)), '01 10');
  same ('units of the 8-bit words', show_all (w8, units (w8)), '11111111');
  sn ('zero divisors of the 8-bit words', numel (zero_divisors (w8)), 254);
  x = 0:255;
  require ('a AND (NOT a) = 0 for every 8-bit word a', all (bitand (x, bitxor (255, x)) == 0));

  m = ring_matrices ();
  um = units (m);
  zm = zero_divisors (m);
  sn ('invertible 2 x 2 matrices over GF(2)', numel (um), 6);
  sn ('(2^2 - 1)(2^2 - 2), the count from the matrix chapter', (4 - 1) * (4 - 2), 6);
  sn ('nonzero matrices that are not invertible', numel (zm), 9);
  both = true;
  for a = zm
    both = both && any (mat_mul (a, 1:15) == 0) && any (mat_mul (1:15, a) == 0);
  end
  require ('each of the nine has a nonzero partner on each side', both);

  gcd_ok = true; phi_ok = true;
  for n = 2:40
    zu = zn_units (n);
    zz = zn_zero_divisors (n);
    a = 1:n - 1;
    gcd_ok = gcd_ok && isequal (zz, a(gcd (a, n) > 1));
    phi_ok = phi_ok && isequal (zu, a(gcd (a, n) == 1)) && numel (zu) == phi (n);
    gcd_ok = gcd_ok && numel (zu) + numel (zz) == n - 1;
  end
  require ('units of Z_n are the residues coprime to n, phi(n), n = 2..40', phi_ok);
  require ('zero divisors of Z_n: nonzero a with gcd(a, n) > 1, n = 2..40', gcd_ok);
  sl ('x -> 3x on Z_7', times_map (3, 7), [0, 3, 6, 2, 5, 1, 4]);
  sl ('x -> 2x on Z_6', times_map (2, 6), [0, 2, 4, 0, 2, 4]);
end

function check_solving ()
  section ('Solving ax = c');
  sl ('4x = 8 in Z_12', solutions (4, 8, 12), [2, 5, 8, 11]);
  sl ('4x = 6 in Z_12', solutions (4, 6, 12), []);
  c = 0:11;
  solvable = arrayfun (@(t) ~isempty (solutions (4, t, 12)), c);
  sl ('c in Z_12 for which 4x = c has a solution', c(solvable), [0, 4, 8]);
  sn ('gcd(4, 12)', gcd (4, 12), 4);
  counts = [];
  for a = 1:12
    for t = 0:12
      counts(end + 1) = numel (solutions (a, t, 13));
    end
  end
  sn ('equations ax = c in Z_13 with a nonzero', numel (counts), 156);
  require ('each of them has exactly one solution', all (counts == 1));
  sl ('4x = 6 in Z_13', solutions (4, 6, 13), [8]);
  b = 0:12;
  sn ('inverse of 4 in Z_13', b(find (mod (4 * b, 13) == 1, 1)), 10);
end

function check_domains (R)
  section ('Integral domains');
  [G, F] = meshgrid (1:63);
  sn ('pairs of nonzero GF(2)[x] polynomials of degree below 6', numel (F), 3969);
  P = clmul (F, G);
  require ('deg(fg) = deg f + deg g for each of them', all (all (degree (P) == degree (F) + degree (G))));
  f = F(P == 1);
  g = G(P == 1);
  parts = arrayfun (@(i) lst ([f(i), g(i)]), 1:numel (f), 'UniformOutput', false);
  same ('pairs among them with fg = 1', ['[' strjoin(parts, ', ') ']'], '[[1, 1]]');
  sl ('(2x)(2x) in Z_4[x], constant first', mod (conv ([0 2], [0 2]), 4), [0, 0, 0]);
  sl ('(2x + 1)(2x + 1) in Z_4[x], constant first', mod (conv ([1 2], [1 2]), 4), [1, 0, 0]);
  ok = true;
  for i = 1:numel (R)
    r = R{i};
    cancel = true;
    for a = 2:r.n
      cancel = cancel && numel (unique (r.mul(a, :))) == r.n;
    end
    ok = ok && cancel == isempty (zero_divisors (r));
  end
  require ('cancellation holds exactly when there are no zero divisors', ok);
  sn ('2 * 3 mod 6', mod (2 * 3, 6), 0);
end

function check_fields (R, stated)
  section ('Fields');
  dom = true; fld = true;
  for n = 2:64
    dom = dom && isempty (zn_zero_divisors (n)) == isprime (n);
    fld = fld && (numel (zn_units (n)) == n - 1) == isprime (n);
  end
  require ('Z_n is an integral domain exactly when n is prime, n = 2..64', dom);
  require ('Z_n is a field exactly when n is prime, n = 2..64', fld);
  finite_domain_ok = true;
  for i = 1:numel (R)
    r = R{i};
    comm = commutative (r.mul);
    domain = comm && isempty (zero_divisors (r));
    field = comm && numel (units (r)) == r.n - 1;
    finite_domain_ok = finite_domain_ok && (~domain || field);
    same (sprintf ('%s: integral domain, field', r.name), sprintf ('%s %s', yes (domain), yes (field)), stated{i});
  end
  require ('every finite integral domain among them is a field', finite_domain_ok);
  z4 = ring_zn (4);
  w2 = ring_words (2);
  f4 = ring_f4 ();
  F = f4_mul ();
  names = '01ab';
  sn ('Z_4: 2 * 2', z4.mul(3, 3), 0);
  same ('Z_4: units, zero divisors', [show_all(z4, units (z4)) '; ' show_all(z4, zero_divisors (z4))], '1 3; 2');
  same ('2-bit words: 01 AND 10', show (w2, w2.mul(2, 3)), '00');
  same ('F_4: a*a, a*b, b*b', [names(F(3, 3) + 1) ' ' names(F(3, 4) + 1) ' ' names(F(4, 4) + 1)], 'b 1 a');
  same ('F_4: 1+a, a+b, 1+b', [names(bitxor(1, 2) + 1) ' ' names(bitxor(2, 3) + 1) ' ' names(bitxor(1, 3) + 1)], 'b 1 a');
  same ('F_4: inverses of 1, a, b', show_all (f4, arrayfun (@(u) inverse_in (f4, u), 1:3)), '1 b a');
  same ('F_4: units', show_all (f4, units (f4)), '1 a b');
end

function check_characteristic (R, stated)
  section ('Characteristic');
  all_ok = true;
  for i = 1:numel (R)
    r = R{i};
    c = characteristic (r);
    sn (sprintf ('characteristic of %s', r.name), c, stated(i));
    s = zeros (1, r.n);
    for t = 1:c
      s = r.add(sub2ind (size (r.add), s + 1, 1:r.n));
    end
    all_ok = all_ok && all (s == 0);
  end
  require ('char * a = 0 for every element a, all eight rings', all_ok);
  f = 0:63;
  require ('f + f = 0, so -f = f, every f in GF(2)[x] of degree below 6', all (bitxor (f, f) == 0));
  is_n = true;
  ok = true;
  for n = 2:64
    c = characteristic (ring_zn (n));
    is_n = is_n && c == n;
    if (isempty (zn_zero_divisors (n)))
      ok = ok && isprime (c);
    end
  end
  require ('characteristic of Z_n is n, n = 2..64', is_n);
  require ('characteristic of Z_n prime whenever Z_n is a domain, n = 2..64', ok);
  same ('7 + 1 as integers, bits', dec2bin (7 + 1, 4), '1000');
  same ('0111 XOR 0001', dec2bin (bitxor (7, 1), 4), '0110');
  carries = [];
  c = 0;
  for i = 0:3
    x = bitand (bitshift (7, -i), 1);
    y = bitand (bitshift (1, -i), 1);
    c = bitor (bitor (bitand (x, y), bitand (x, c)), bitand (y, c));
    carries(end + 1) = c;
  end
  sl ('carries out of bits 0..3 for 0111 + 0001', carries, [1, 1, 1, 0]);
end

function check_freshman ()
  section ('The binomial theorem in characteristic p');
  primes_ = primes (49);
  sl ('primes below 50', primes_, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]);
  ok = true;
  for p = primes_
    row = pascal_row (p);
    ok = ok && all (mod (row(2:p), p) == 0);
  end
  require ('p divides C(p, i) for 0 < i < p, every prime p below 50', ok);
  ok = true;
  for n = 4:49
    if (~isprime (n))
      row = pascal_row (n);
      ok = ok && any (mod (row(2:n), n) ~= 0);
    end
  end
  require ('n fails to divide some C(n, i), every composite n from 4 to 49', ok);
  sl ('C(5, i) for i = 0..5', pascal_row (5), [1, 5, 10, 10, 5, 1]);
  row4 = pascal_row (4);
  sn ('C(4, 2)', row4(3), 6);
  ok = true;
  for p = primes_
    [B, A] = meshgrid (0:p - 1);
    ok = ok && all (all (powmod (A + B, p, p) == mod (powmod (A, p, p) + powmod (B, p, p), p)));
  end
  require ('(a + b)^p = a^p + b^p in Z_p for all a, b, every prime p below 50', ok);
  sl ('(1 + 1)^2 and 1^2 + 1^2 in Z_4', [mod((1 + 1) ^ 2, 4), mod(1 + 1, 4)], [0, 2]);
  [G, F] = meshgrid (0:63);
  S = bitxor (F, G);
  require ('(f + g)^2 = f^2 + g^2 in GF(2)[x], all f, g of degree below 6', ...
           all (all (clmul (S, S) == bitxor (clmul (F, F), clmul (G, G)))));
  f = 0:255;
  require ('squaring spreads bit i to bit 2i, every f of degree below 8', all (clmul (f, f) == spread (f)));
  same ('(x^3 + x + 1)^2 in GF(2)[x], bits', [bits(11) ' -> ' bits(clmul (11, 11))], '1011 -> 1000101');
  sl ('(x^3 + x + 1)^2 in Z[x], constant first', conv ([1 1 0 1], [1 1 0 1]), [1, 2, 1, 2, 2, 0, 1]);
  M = f4_mul ();
  [V, U] = meshgrid (0:3);
  S = bitxor (U, V);
  lhs = M(sub2ind ([4 4], S + 1, S + 1));
  rhs = bitxor (M(sub2ind ([4 4], U + 1, U + 1)), M(sub2ind ([4 4], V + 1, V + 1)));
  require ('(u + v)^2 = u^2 + v^2 in F_4 for all 16 pairs', all (all (lhs == rhs)));
end

function check_sizes ()
  section ('How big can a finite field be');
  f4 = ring_f4 ();
  z13 = ring_zn (13);
  M = f4_mul ();
  same ('multiples of 1 in F_4', show_all (f4, multiples_of_one (f4)), '0 1');
  pairs = [0 0; 0 1; 1 0; 1 1];
  vals = arrayfun (@(i) bitxor (M(pairs(i, 1) + 1, 3), M(pairs(i, 2) + 1, 2)), 1:4);
  same ('c1 * a + c0 * 1 in F_4 for c1 c0 = 00, 01, 10, 11', show_all (f4, vals), '0 1 a b');
  sn ('multiples of 1 in Z_13', numel (multiples_of_one (z13)), 13);
  n = 2:32;
  pp = arrayfun (@is_prime_power, n);
  sl ('sizes from 2 to 32 that are not prime powers', n(~pp), [6, 10, 12, 14, 15, 18, 20, 21, 22, 24, 26, 28, 30]);
  n = 2:64;
  pp = arrayfun (@is_prime_power, n);
  pr = isprime (n);
  sl ('sizes from 2 to 64: primes, higher prime powers, excluded', [sum(pr), sum(pp & ~pr), sum(~pp)], [18, 9, 36]);
  sl ('higher prime powers up to 64', n(pp & ~pr), [4, 8, 9, 16, 25, 27, 32, 49, 64]);
end

function check_examples ()
  section ('Worked examples');
  z15 = ring_zn (15);
  u15 = units (z15);
  sl ('units of Z_15', u15, [1, 2, 4, 7, 8, 11, 13, 14]);
  sl ('their inverses', arrayfun (@(u) inverse_in (z15, u), u15), [1, 8, 4, 13, 2, 11, 7, 14]);
  sl ('zero divisors of Z_15', zero_divisors (z15), [3, 5, 6, 9, 10, 12]);
  sl ('3 * 5, 6 * 5, 9 * 5, 10 * 3, 12 * 5 mod 15', mod ([3 * 5, 6 * 5, 9 * 5, 10 * 3, 12 * 5], 15), [0, 0, 0, 0, 0]);
  sl ('6x = 10 in Z_256', solutions (6, 10, 256), [87, 215]);
  sl ('6x = 7 in Z_256', solutions (6, 7, 256), []);
  b = 0:127;
  sn ('inverse of 3 modulo 128', b(find (mod (3 * b, 128) == 1, 1)), 43);
  sn ('5 * 43 mod 128', mod (5 * 43, 128), 87);
  sl ('6 * 87 and 6 * 215 mod 256', [mod(6 * 87, 256), mod(6 * 215, 256)], [10, 10]);
  m = ring_matrices ();
  N = 15;
  I = 9;
  same ('N * N for N = 11|11', show (m, mat_mul (N, N)), '00|00');
  same ('I + N', show (m, bitxor (I, N)), '01|10');
  same ('(I + N)^2', show (m, mat_mul (bitxor (I, N), bitxor (I, N))), '10|01');
  f = 211;
  same ('a byte squared by wiring, bits', [bits(f) ' -> ' bits(spread (f))], '11010011 -> 101000100000101');
  same ('the same square by multiplying', bits (clmul (f, f)), '101000100000101');
  sn ('degree of the square', degree (spread (f)), 14);
end

function check_exercises ()
  section ('Exercises');
  z20 = ring_zn (20);
  sl ('units of Z_20', units (z20), [1, 3, 7, 9, 11, 13, 17, 19]);
  sl ('zero divisors of Z_20', zero_divisors (z20), [2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 18]);
  sl ('12x = 6 in Z_18', solutions (12, 6, 18), [2, 5, 8, 11, 14, 17]);
  sl ('12x = 4 in Z_18', solutions (12, 4, 18), []);
  sn ('gcd(12, 18)', gcd (12, 18), 6);
  w3 = ring_words (3);
  same ('3-bit words: units, zero divisors, characteristic', ...
        sprintf ('%s; %d; %d', show_all (w3, units (w3)), numel (zero_divisors (w3)), characteristic (w3)), '111; 6; 2');
  sl ('C(7, i) for i = 0..7', pascal_row (7), [1, 7, 21, 35, 35, 21, 7, 1]);
  sl ('C(6, i) for i = 0..6', pascal_row (6), [1, 6, 15, 20, 15, 6, 1]);
  sl ('C(6, 2) mod 6, C(6, 3) mod 6', [mod(15, 6), mod(20, 6)], [3, 2]);
  same ('49, 50, 64, 100, 121 prime powers', ...
        strjoin (arrayfun (@(n) yes (is_prime_power (n)), [49 50 64 100 121], 'UniformOutput', false), ' '), ...
        'yes no yes no yes');
  sl ('(3 + 2x)(3 + 2x) in Z_4[x], constant first', mod (conv ([3 2], [3 2]), 4), [1, 0, 0]);
  M = f4_mul ();
  names = '01ab';
  x = find (M(3, :) == bitxor (1, 2), 1) - 1;
  same ('F_4: solution of a x = 1 + a', names(x + 1), 'a');
end


R = chapter_rings ();
check_opening ();
check_axioms (R, {'4 yes', '12 yes', '13 yes', '256 yes', '4 yes', '256 yes', '16 no', '4 yes'});
check_gallery ();
check_consequences (R);
check_homomorphism ();
check_units (R);
check_solving ();
check_domains (R);
check_fields (R, {'no no', 'no no', 'yes yes', 'no no', 'no no', 'no no', 'no no', 'yes yes'});
check_characteristic (R, [4, 12, 13, 256, 2, 2, 2, 2]);
check_freshman ();
check_sizes ();
check_examples ();
check_exercises ();
printf ('\nAll chapter 32 numbers reproduced.\n');
