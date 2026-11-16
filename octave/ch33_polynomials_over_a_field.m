% Reproduce every number in chapter 33 of Error Correction Coding,
% "Polynomials over a Field", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch33_polynomials_over_a_field.m
%
% It needs only core Octave 7 or later, with no pkg load and no Forge
% packages. It is the Octave twin of python/ch33_polynomials_over_a_field.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. All the arithmetic is on whole numbers small enough for a double
% to hold exactly, so every comparison is exact. Polynomials over GF(2) are
% held as whole numbers whose bit i is the coefficient of x^i, so addition is
% bitxor and multiplication is the carry-less product. Polynomials over GF(3),
% GF(5) and the integers modulo n are row vectors of coefficients, constant
% term first, with the zero polynomial an empty row.

1;


% Helpers ----------------------------------------------------------------------

function s = num (v)
  s = sprintf ('%d', v);
end

function s = jn (v)
  % Whole numbers joined by single spaces.
  if (isempty (v))
    s = '';
  else
    s = strjoin (arrayfun (@(x) sprintf ('%d', x), v(:)', 'UniformOutput', false), ' ');
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


% The coefficient rings -------------------------------------------------------

function ok = ring_axioms (n)
  e = 0:n-1;
  add = mod (e' + e, n);
  mul = mod (e' * e, n);
  ok = false;
  for a = 0:n-1
    if (add(a+1, 1) ~= a || mul(a+1, 2) ~= a)
      return;
    end
    if (~any (add(a+1, :) == 0))
      return;
    end
    for b = 0:n-1
      if (add(a+1, b+1) ~= add(b+1, a+1) || mul(a+1, b+1) ~= mul(b+1, a+1))
        return;
      end
      for c = 0:n-1
        if (add(add(a+1, b+1)+1, c+1) ~= add(a+1, add(b+1, c+1)+1))
          return;
        end
        if (mul(mul(a+1, b+1)+1, c+1) ~= mul(a+1, mul(b+1, c+1)+1))
          return;
        end
        if (mul(a+1, add(b+1, c+1)+1) ~= add(mul(a+1, b+1)+1, mul(a+1, c+1)+1))
          return;
        end
      end
    end
  end
  ok = true;
end

function u = units (n)
  u = [];
  for a = 1:n-1
    if (any (mod (a * (0:n-1), n) == 1))
      u(end+1) = a;
    end
  end
end

function z = zero_divisors (n)
  z = [];
  for a = 1:n-1
    if (any (mod (a * (1:n-1), n) == 0))
      z(end+1) = a;
    end
  end
end

function b = inverse (a, p)
  for b = 1:p-1
    if (mod (a * b, p) == 1)
      return;
    end
  end
  error ('%d has no inverse modulo %d', a, p);
end


% Polynomials over GF(2) as whole numbers -------------------------------------

function d = deg2 (a)
  % The degree, with -1 standing for the zero polynomial.
  if (a == 0)
    d = -1;
  else
    [~, e] = log2 (a);
    d = e - 1;
  end
end

function r = clmul (a, b)
  % The carry-less product: shift and XOR, so no carry ever crosses a column.
  r = 0;
  while (b > 0)
    if (bitand (b, 1))
      r = bitxor (r, a);
    end
    a = a * 2;
    b = floor (b / 2);
  end
end

function [q, f] = divmod2 (f, g)
  q = 0;
  dg = deg2 (g);
  while (f > 0 && deg2 (f) >= dg)
    s = deg2 (f) - dg;
    q = bitor (q, 2 ^ s);
    f = bitxor (f, g * 2 ^ s);
  end
end

function [shifts, rems] = division_steps2 (f, g)
  shifts = [];
  rems = {};
  dg = deg2 (g);
  while (f > 0 && deg2 (f) >= dg)
    s = deg2 (f) - dg;
    f = bitxor (f, g * 2 ^ s);
    shifts(end+1) = s;
    rems{end+1} = bits (f);
  end
end

function r = mod2 (f, g)
  [~, r] = divmod2 (f, g);
end

function t = divides2 (d, f)
  t = mod2 (f, d) == 0;
end

function s = bits (a)
  s = dec2bin (a);
end

function w = weight (a)
  w = sum (dec2bin (a) == '1');
end

function s = ptxt2 (a)
  if (a == 0)
    s = '0';
    return;
  end
  terms = {};
  for k = deg2 (a):-1:0
    if (bitget (a, k + 1))
      if (k == 0)
        terms{end+1} = '1';
      elseif (k == 1)
        terms{end+1} = 'x';
      else
        terms{end+1} = sprintf ('x^%d', k);
      end
    end
  end
  s = strjoin (terms, '+');
end

function s = jp2 (v)
  % Polynomials over GF(2) written out and joined by single spaces.
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = ptxt2 (v(i));
  end
  s = strjoin (parts, ' ');
end

function s = jb (v)
  % Bit strings joined by single spaces.
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = bits (v(i));
  end
  s = strjoin (parts, ' ');
end

function v = eval2 (f, a)
  if (a == 1)
    v = mod (weight (f), 2);
  else
    v = bitand (f, 1);
  end
end

function rows = xgcd2 (a, b)
  % Rows [r, q, s, t] of the extended Euclidean table. Row i carries the
  % quotient q_i of r_(i-1) by r_i, so r_(i+1) = r_(i-1) - q_i r_i, and the
  % last row holds the remainder 0 and q = -1.
  rows = [a, -1, 1, 0; b, -1, 0, 1];
  while (rows(end, 1) ~= 0)
    r0 = rows(end-1, 1); s0 = rows(end-1, 3); t0 = rows(end-1, 4);
    r1 = rows(end, 1); s1 = rows(end, 3); t1 = rows(end, 4);
    [q, r] = divmod2 (r0, r1);
    rows(end, 2) = q;
    rows(end+1, :) = [r, -1, bitxor(s0, clmul (q, s1)), bitxor(t0, clmul (q, t1))];
  end
end

function out = factor_list2 (f)
  % The irreducible factors of f by trial division, smallest first, with repeats.
  out = [];
  d = 2;
  while (deg2 (f) > 0)
    if (deg2 (d) * 2 > deg2 (f))
      out(end+1) = f;
      break;
    end
    [q, r] = divmod2 (f, d);
    if (r == 0)
      out(end+1) = d;
      f = q;
    else
      d = d + 1;
    end
  end
end

function s = factor_text2 (f)
  fs = factor_list2 (f);
  s = '';
  i = 1;
  while (i <= numel (fs))
    j = i;
    while (j <= numel (fs) && fs(j) == fs(i))
      j = j + 1;
    end
    s = [s '(' ptxt2(fs(i)) ')'];
    if (j - i > 1)
      s = [s sprintf('^%d', j - i)];
    end
    i = j;
  end
end

function t = irreducible_by_trial2 (f)
  t = false;
  if (deg2 (f) < 1)
    return;
  end
  for d = 2 : 2 ^ (floor (deg2 (f) / 2) + 1) - 1
    if (mod2 (f, d) == 0)
      return;
    end
  end
  t = true;
end

function irr = sieve2 (top)
  % Cross out every product of a known irreducible with a polynomial of
  % degree at least 1. What survives up to degree top is irreducible.
  limit = 2 ^ (top + 1);
  crossed = false (1, limit);
  for a = 2:limit-1
    if (crossed(a+1) || 2 * deg2 (a) > top)
      continue;
    end
    for b = 2 : 2 ^ (top - deg2 (a) + 1) - 1
      crossed(clmul (a, b) + 1) = true;
    end
  end
  irr = find (~crossed(3:end)) + 1;
end

function d = degs2 (v)
  d = zeros (1, numel (v));
  for i = 1:numel (v)
    d(i) = deg2 (v(i));
  end
end

function counts = irreducible_counts (top)
  irr = sieve2 (top);
  dg = degs2 (irr);
  counts = zeros (1, top);
  for m = 1:top
    counts(m) = sum (dg == m);
  end
end

function hits = build (start, product, irr6, hits)
  hits(product + 1) = hits(product + 1) + 1;
  for i = start:numel (irr6)
    if (deg2 (product) + deg2 (irr6(i)) <= 6)
      hits = build (i, clmul (product, irr6(i)), irr6, hits);
    end
  end
end


% Polynomials over GF(p) or Z_n as rows, constant term first -------------------

function f = trim (f)
  f = f(:)';
  while (~isempty (f) && f(end) == 0)
    f(end) = [];
  end
end

function d = pdeg (f)
  d = numel (trim (f)) - 1;
end

function h = padd (f, g, n)
  m = max (numel (f), numel (g));
  F = [f(:)', zeros(1, m - numel (f))];
  G = [g(:)', zeros(1, m - numel (g))];
  h = trim (mod (F + G, n));
end

function h = pneg (f, n)
  h = trim (mod (-f, n));
end

function c = pconv (f, g)
  % The coefficients of the product before any reduction: c_k = sum of f_i g_(k-i).
  if (isempty (f) || isempty (g))
    c = [];
    return;
  end
  c = zeros (1, numel (f) + numel (g) - 1);
  for i = 1:numel (f)
    for j = 1:numel (g)
      c(i + j - 1) = c(i + j - 1) + f(i) * g(j);
    end
  end
end

function h = pmul (f, g, n)
  h = trim (mod (pconv (f, g), n));
end

function [q, f] = pdivmod (f, g, p)
  f = trim (f);
  g = trim (g);
  iv = inverse (g(end), p);
  q = zeros (1, max (numel (f) - numel (g) + 1, 0));
  while (~isempty (f) && numel (f) >= numel (g))
    s = numel (f) - numel (g);
    c = mod (f(end) * iv, p);
    q(s + 1) = c;
    f = padd (f, pneg ([zeros(1, s), mod(c * g, p)], p), p);
  end
  q = trim (q);
end

function [terms, rems] = division_steps (f, g, p)
  f = trim (f);
  g = trim (g);
  iv = inverse (g(end), p);
  terms = {};
  rems = {};
  while (~isempty (f) && numel (f) >= numel (g))
    s = numel (f) - numel (g);
    c = mod (f(end) * iv, p);
    f = padd (f, pneg ([zeros(1, s), mod(c * g, p)], p), p);
    terms{end+1} = ptxt ([zeros(1, s), c]);
    rems{end+1} = ptxt (f);
  end
end

function v = peval (f, a, n)
  v = 0;
  for k = numel (f):-1:1
    v = mod (v * a + f(k), n);
  end
end

function v = pevals (f, as, n)
  v = zeros (1, numel (as));
  for i = 1:numel (as)
    v(i) = peval (f, as(i), n);
  end
end

function r = proots (f, n)
  r = [];
  for a = 0:n-1
    if (peval (f, a, n) == 0)
      r(end+1) = a;
    end
  end
end

function s = ptxt (f)
  f = trim (f);
  if (isempty (f))
    s = '0';
    return;
  end
  terms = {};
  for k = numel (f) - 1:-1:0
    c = f(k + 1);
    if (c == 0)
      continue;
    end
    if (k == 0)
      terms{end+1} = sprintf ('%d', c);
    else
      if (c == 1)
        coef = '';
      else
        coef = sprintf ('%d', c);
      end
      if (k == 1)
        terms{end+1} = [coef 'x'];
      else
        terms{end+1} = sprintf ('%sx^%d', coef, k);
      end
    end
  end
  s = strjoin (terms, '+');
end

function s = jpc (rows, col)
  % One column of a table of polynomials over GF(p), written out and joined.
  parts = cell (1, size (rows, 1));
  for i = 1:size (rows, 1)
    parts{i} = ptxt (rows{i, col});
  end
  s = strjoin (parts, ' ');
end

function rows = pxgcd (a, b, p)
  rows = {trim(a), [], 1, []; trim(b), [], [], 1};
  while (~isempty (rows{end, 1}))
    r0 = rows{end-1, 1}; s0 = rows{end-1, 3}; t0 = rows{end-1, 4};
    r1 = rows{end, 1}; s1 = rows{end, 3}; t1 = rows{end, 4};
    [q, r] = pdivmod (r0, r1, p);
    rows{end, 2} = q;
    rows(end+1, :) = {r, [], padd(s0, pneg (pmul (q, s1, p), p), p), padd(t0, pneg (pmul (q, t1, p), p), p)};
  end
end

function h = scale (f, c, p)
  h = trim (mod (c * f, p));
end

function out = from_index (i, p, len)
  out = zeros (1, len);
  for k = 1:len
    out(k) = mod (i, p);
    i = floor (i / p);
  end
  out = trim (out);
end

function m = monics (p, d)
  m = cell (1, p ^ d);
  for i = 0:p^d - 1
    m{i + 1} = [mod(floor (i ./ p .^ (0:d-1)), p), 1];
  end
end

function k = pkey (f, p)
  k = sum (f .* p .^ (0:numel (f) - 1));
end

function show_table2 (name, rows, st)
  same ([name ': remainders r_i'], jp2 (rows(:, 1)), st{1});
  same ([name ': quotients q_i'], jp2 (rows(rows(:, 2) >= 0, 2)), st{2});
  same ([name ': cofactors s_i'], jp2 (rows(:, 3)), st{3});
  same ([name ': cofactors t_i'], jp2 (rows(:, 4)), st{4});
end


% The checks, in the order the chapter states them ---------------------------

section ('The coefficient rings');
ns = [4, 8, 2, 3, 5];
names = {'Z_4', 'Z_8', 'GF(2)', 'GF(3)', 'GF(5)'};
for i = 1:5
  require (sprintf ('ring axioms on all %d triples of %s', ns(i) ^ 3, names{i}), ring_axioms (ns(i)));
end
for p = [2, 3, 5]
  sn (sprintf ('nonzero elements of GF(%d) with an inverse', p), numel (units (p)), p - 1);
end
same ('zero divisors of Z_4', jn (zero_divisors (4)), '2');
same ('zero divisors of Z_8', jn (zero_divisors (8)), '2 4 6');
invs = zeros (1, 4);
for a = 1:4
  invs(a) = inverse (a, 5);
end
same ('inverses of 1 2 3 4 in GF(5)', jn (invs), '1 3 2 4');

section ('The opening product');
a = 11; b = 7;
same ('1011 and 0111 as whole numbers', sprintf ('%d %d', a, b), '11 7');
sn ('their integer product', a * b, 77);
same ('77 in binary', bits (a * b), '1001101');
rows = [];
for i = 0:3
  if (bitget (b, i + 1))
    rows(end+1) = a * 2 ^ i;
  end
end
same ('the shifted copies of 1011', jb (rows), '1011 10110 101100');
sums = zeros (1, 6);
for k = 5:-1:0
  sums(6 - k) = sum (bitget (rows, k + 1));
end
same ('column sums, from x^5 down to the constant', jn (sums), '1 1 2 2 2 1');
carry = 0;
carries = zeros (1, 6);
for k = 0:5
  total = sum (bitget (rows, k + 1)) + carry;
  carry = floor (total / 2);
  carries(k + 1) = carry;
end
same ('carries out of the six columns, constant column first', jn (carries), '0 1 1 1 1 1');
p = clmul (a, b);
same ('carry-less product of 1011 and 0111', bits (p), '110001');
sn ('the carry-less product as a whole number', p, 49);
same ('(x^3+x+1)(x^2+x+1) over GF(2)', ptxt2 (p), 'x^5+x^4+1');
same ('column sums reduced modulo 2', jn (mod (sums, 2)), '1 1 0 0 0 1');

section ('Addition and multiplication');
same ('(x^3+x+1) + (x^2+x+1) over GF(2)', ptxt2 (bitxor (11, 7)), 'x^3+x^2');
same ('1011 XOR 0111', dec2bin (bitxor (11, 7), 4), '1100');
f3 = [1, 1, 2]; g3 = [2, 1];
same ('(2x^2+x+1) + (x+2) over GF(3)', ptxt (padd (f3, g3, 3)), '2x^2+2x');
same ('convolution sums of 2x^2+x+1 and x+2, top first', jn (fliplr (pconv (f3, g3))), '2 5 3 2');
same ('(2x^2+x+1)(x+2) over GF(3)', ptxt (pmul (f3, g3, 3)), '2x^3+2x^2+2');
f5 = [2, 4, 3]; g5 = [3, 2];
same ('(3x^2+4x+2) + (2x+3) over GF(5)', ptxt (padd (f5, g5, 5)), '3x^2+x');
same ('convolution sums of 3x^2+4x+2 and 2x+3, top first', jn (fliplr (pconv (f5, g5))), '6 17 16 6');
same ('(3x^2+4x+2)(2x+3) over GF(5)', ptxt (pmul (f5, g5, 5)), 'x^3+2x^2+x+1');
pairs = 0; ok = true;
for i = 1:26
  for j = 1:26
    f = from_index (i, 3, 3); g = from_index (j, 3, 3);
    pairs = pairs + 1;
    ok = ok && pdeg (pmul (f, g, 3)) == pdeg (f) + pdeg (g);
  end
end
sn ('pairs of nonzero GF(3) polynomials of degree below 3', pairs, 676);
require ('deg(fg) = deg f + deg g on every one of those pairs', ok);
same ('(2x+1)^2 over Z_4', ptxt (pmul ([1, 2], [1, 2], 4)), '1');
sq = {};
for a1 = 1:3
  for a0 = 0:3
    if (isequal (pmul ([a0, a1], [a0, a1], 4), 1))
      sq{end+1} = ptxt ([a0, a1]);
    end
  end
end
same ('degree-1 polynomials over Z_4 whose square is 1', strjoin (sq, ' '), '2x+1 2x+3');
same ('(2x)(2x) over Z_4', ptxt (pmul ([0, 2], [0, 2], 4)), '0');

section ('A polynomial is not its function');
f = 6;
same ('x^2+x over GF(2) at 0 and at 1', sprintf ('%d %d', eval2 (f, 0), eval2 (f, 1)), '0 0');
vanish = [];
for g = 0:15
  if (eval2 (g, 0) == 0 && eval2 (g, 1) == 0)
    vanish(end+1) = g;
  end
end
same ('polynomials of degree below 4 that vanish at 0 and 1', jp2 (vanish), '0 x^2+x x^3+x x^3+x^2');
sn ('how many of the 16 that is', numel (vanish), 4);
sn ('functions from GF(2) to GF(2)', 2 ^ 2, 4);
ok = true;
for v = vanish
  ok = ok && divides2 (6, v);
end
require ('each of those four is divisible by x^2+x', ok);
same ('x^3-x over GF(3) at 0, 1 and 2', jn (pevals ([0, 2, 0, 1], 0:2, 3)), '0 0 0');
same ('x^5-x over GF(5) at 0 through 4', jn (pevals ([0, 4, 0, 0, 0, 1], 0:4, 5)), '0 0 0 0 0');

section ('The division algorithm');
pairs = 0; ok = true; uniq = true;
for f = 0:63
  for g = 1:15
    [q, r] = divmod2 (f, g);
    pairs = pairs + 1;
    ok = ok && bitxor (clmul (q, g), r) == f && deg2 (r) < deg2 (g);
    count = 0;
    for rr = 0 : 2 ^ deg2 (g) - 1
      if (divides2 (g, bitxor (f, rr)))
        count = count + 1;
      end
    end
    uniq = uniq && count == 1;
  end
end
sn ('pairs f of degree below 6, g nonzero of degree below 4', pairs, 960);
require ('f = qg + r with deg r < deg g on every pair', ok);
require ('exactly one remainder of degree below deg g on every pair', uniq);
require ('2c = 1 has no integer solution, since 2c is even', all (mod (2 * (-50:50), 2) == 0));
f = 147; g = 11;
same ('x^7+x^4+x+1 as a bit string', bits (f), '10010011');
[shifts, rems] = division_steps2 (f, g);
same ('shifts of 1011 at each step', jn (shifts), '4 2 0');
same ('what is left after each step', strjoin (rems, ' '), '100011 1111 100');
[q, r] = divmod2 (f, g);
same ('quotient of x^7+x^4+x+1 by x^3+x+1', ptxt2 (q), 'x^4+x^2+1');
same ('remainder', ptxt2 (r), 'x^2');
same ('(x^4+x^2+1)(x^3+x+1) over GF(2)', ptxt2 (clmul (q, g)), 'x^7+x^4+x^2+x+1');
f5 = [2, 4, 1, 3]; g5 = [1, 3, 2];
sn ('inverse of the leading coefficient 2 in GF(5)', inverse (2, 5), 3);
[terms, rems5] = division_steps (f5, g5, 5);
same ('quotient terms of 3x^3+x^2+4x+2 by 2x^2+3x+1 over GF(5)', strjoin (terms, ' '), '4x 2');
same ('what is left after each step', strjoin (rems5, ' '), '4x^2+2 4x');
[q5, r5] = pdivmod (f5, g5, 5);
same ('quotient', ptxt (q5), '4x+2');
same ('remainder', ptxt (r5), '4x');
same ('(4x+2)(2x^2+3x+1) + 4x over GF(5)', ptxt (padd (pmul (q5, g5, 5), r5, 5)), '3x^3+x^2+4x+2');

section ('Remainders respect sums and products');
g = 11; a = 16; b = 12;
ra = mod2 (a, g); rb = mod2 (b, g);
same ('x^4 mod (x^3+x+1)', ptxt2 (ra), 'x^2+x');
same ('(x^3+x^2) mod (x^3+x+1)', ptxt2 (rb), 'x^2+x+1');
same ('(x^4)(x^3+x^2) = x^7+x^6, reduced', ptxt2 (mod2 (clmul (a, b), g)), 'x^2');
same ('(x^2+x)(x^2+x+1), before reduction', ptxt2 (clmul (ra, rb)), 'x^4+x');
same ('(x^2+x)(x^2+x+1), reduced', ptxt2 (mod2 (clmul (ra, rb), g)), 'x^2');
same ('x^4 + x^3 + x^2, reduced', ptxt2 (mod2 (bitxor (a, b), g)), '1');
same ('(x^2+x) + (x^2+x+1)', ptxt2 (bitxor (ra, rb)), '1');
pairs = 0; ok = true;
for a = 0:31
  for b = 0:31
    pairs = pairs + 1;
    ok = ok && mod2 (bitxor (a, b), g) == bitxor (mod2 (a, g), mod2 (b, g));
    ok = ok && mod2 (clmul (a, b), g) == mod2 (clmul (mod2 (a, g), mod2 (b, g)), g);
  end
end
sn ('pairs of polynomials of degree below 5', pairs, 1024);
require ('both rules hold modulo x^3+x+1 on every pair', ok);

section ('Greatest common divisors and Euclid''s algorithm');
a = 49; b = 21;
same ('x^5+x^4+1, factored', factor_text2 (a), '(x^2+x+1)(x^3+x+1)');
same ('x^4+x^2+1, factored', factor_text2 (b), '(x^2+x+1)^2');
rows = xgcd2 (a, b);
show_table2 ('GF(2) table', rows, {'x^5+x^4+1 x^4+x^2+1 x^3+x^2+x x^2+x+1 0', 'x+1 x+1 x', ...
                                  '1 0 1 x+1 x^2+x+1', '0 1 x+1 x^2 x^3+x+1'});
g = rows(end-1, 1); s = rows(end-1, 3); t = rows(end-1, 4);
same ('gcd(x^5+x^4+1, x^4+x^2+1)', ptxt2 (g), 'x^2+x+1');
same ('(x+1)(x^5+x^4+1) + x^2(x^4+x^2+1)', ptxt2 (bitxor (clmul (s, a), clmul (t, b))), 'x^2+x+1');
same ('degrees of the remainders', jn (degs2 (rows(1:end-1, 1))), '5 4 3 2');
divides = false (32, 32);
for d = 1:31
  for f = 0:31
    divides(d + 1, f + 1) = divides2 (d, f);
  end
end
dg = degs2 (0:31);
pairs = 0; ok = true;
for x = 1:31
  for y = 1:31
    pairs = pairs + 1;
    common = find (divides(:, x + 1) & divides(:, y + 1))' - 1;
    top = max (dg(common + 1));
    best = common(dg(common + 1) == top);
    rows = xgcd2 (x, y);
    gg = rows(end-1, 1); ss = rows(end-1, 3); tt = rows(end-1, 4);
    ok = ok && isequal (best, gg) && bitxor (clmul (ss, x), clmul (tt, y)) == gg;
    for d = common
      ok = ok && divides2 (d, gg);
    end
  end
end
sn ('pairs of nonzero GF(2) polynomials of degree below 5', pairs, 961);
require ('Euclid''s answer is the one common divisor of top degree', ok);
require ('and every common divisor divides it, with s a + t b equal to it', ok);
a3 = [1, 1, 1, 1]; b3 = [1, 0, 2];
same ('x^3+x^2+x+1 over GF(3), factored by hand', ptxt (pmul ([1, 1], [1, 0, 1], 3)), 'x^3+x^2+x+1');
same ('(x+1)(2x+1) over GF(3)', ptxt (pmul ([1, 1], [1, 2], 3)), '2x^2+1');
rows3 = pxgcd (a3, b3, 3);
same ('GF(3) table: remainders r_i', jpc (rows3, 1), 'x^3+x^2+x+1 2x^2+1 2x+2 0');
same ('GF(3) table: quotients q_i', jpc (rows3(2:end-1, :), 2), '2x+2 x+2');
same ('GF(3) table: cofactors s_i', jpc (rows3, 3), '1 0 1 2x+1');
same ('GF(3) table: cofactors t_i', jpc (rows3, 4), '0 1 x+1 2x^2+2');
g3 = rows3{end-1, 1}; s3 = rows3{end-1, 3}; t3 = rows3{end-1, 4};
same ('1(x^3+x^2+x+1) + (x+1)(2x^2+1) over GF(3)', ptxt (padd (pmul (s3, a3, 3), pmul (t3, b3, 3), 3)), '2x+2');
c = inverse (g3(end), 3);
sn ('inverse of the leading coefficient 2 in GF(3)', c, 2);
same ('the monic gcd', ptxt (scale (g3, c, 3)), 'x+1');
same ('the monic cofactors s and t', [ptxt(scale (s3, c, 3)) ' ' ptxt(scale (t3, c, 3))], '2 2x+2');
same ('2(x^3+x^2+x+1) + (2x+2)(2x^2+1) over GF(3)', ...
      ptxt (padd (pmul (scale (s3, c, 3), a3, 3), pmul (scale (t3, c, 3), b3, 3), 3)), 'x+1');
rows = xgcd2 (11, 5);
show_table2 ('inverse of x^2+1', rows, {'x^3+x+1 x^2+1 1 0', 'x x^2+1', '1 0 1 x^2+1', '0 1 x x^3+x+1'});
same ('inverse of x^2+1 modulo x^3+x+1', ptxt2 (rows(end-1, 4)), 'x');
same ('x(x^2+1) = x^3+x, reduced modulo x^3+x+1', ptxt2 (mod2 (clmul (2, 5), 11)), '1');
m3 = [1, 0, 1]; u3 = [1, 1];
rows = pxgcd (m3, u3, 3);
same ('x^2+1 and x+1 over GF(3): remainders r_i', jpc (rows, 1), 'x^2+1 x+1 2 0');
same ('x^2+1 and x+1 over GF(3): q_1, s_2 and t_2', ...
      [ptxt(rows{2, 2}) ' ' ptxt(rows{3, 3}) ' ' ptxt(rows{3, 4})], 'x+2 1 2x+1');
inv3 = scale (rows{3, 4}, inverse (rows{3, 1}(end), 3), 3);
same ('inverse of x+1 modulo x^2+1 over GF(3), 2 times t_2', ptxt (inv3), 'x+2');
[~, r3] = pdivmod (pmul (u3, inv3, 3), m3, 3);
same ('(x+1)(x+2) modulo x^2+1 over GF(3)', ptxt (r3), '1');

section ('Roots and factors');
f = [1, 2, 0, 1];
sn ('x^3+2x+1 at 2 over GF(5)', peval (f, 2, 5), 3);
[q, r] = pdivmod (f, [3, 1], 5);
same ('x - 2 written over GF(5)', ptxt ([3, 1]), 'x+3');
same ('quotient of x^3+2x+1 by x+3 over GF(5)', ptxt (q), 'x^2+2x+1');
same ('remainder', ptxt (r), '3');
polys = 0; ok = true; bound = true;
for i = 0:624
  g = from_index (i, 5, 4);
  polys = polys + 1;
  for a = 0:4
    [~, rr] = pdivmod (g, [mod(-a, 5), 1], 5);
    if (isempty (rr))
      r0 = 0;
    else
      r0 = rr(1);
    end
    ok = ok && r0 == peval (g, a, 5);
  end
  if (~isempty (g))
    bound = bound && numel (proots (g, 5)) <= pdeg (g);
  end
end
sn ('polynomials over GF(5) of degree below 4', polys, 625);
require ('f mod (x - a) = f(a) for every one of them and every a', ok);
require ('no nonzero one has more roots than its degree', bound);
same ('x^4-1 over GF(5): its roots', jn (proots ([4, 0, 0, 0, 1], 5)), '1 2 3 4');
same ('x^2-1 over Z_8: its roots', jn (proots ([7, 0, 1], 8)), '1 3 5 7');
same ('(x-3)(x+3) over Z_8', ptxt (pmul ([5, 1], [3, 1], 8)), 'x^2+7');
same ('(x-1)(x+1) over Z_8', ptxt (pmul ([7, 1], [1, 1], 8)), 'x^2+7');
root1 = true; root0 = true;
for g = 0:255
  root1 = root1 && ((eval2 (g, 1) == 0) == (mod (weight (g), 2) == 0));
  root0 = root0 && ((eval2 (g, 0) == 0) == (bitand (g, 1) == 0));
  if (g ~= 0 && eval2 (g, 1) == 0)
    root1 = root1 && divides2 (3, g);
  end
end
require ('over GF(2): 1 is a root exactly when the weight is even', root1);
require ('over GF(2): 0 is a root exactly when the constant is 0', root0);
same ('x^3+x^2+x+1 over GF(2): weight and root 1', sprintf ('%d %s', weight (15), yes (eval2 (15, 1) == 0)), '4 yes');

section ('Irreducible polynomials');
for p = [2, 3, 5]
  ok = true;
  for d = [2, 3]
    seen = false (1, p ^ (d + 1));
    G = monics (p, 1);
    H = monics (p, d - 1);
    for i = 1:numel (G)
      for j = 1:numel (H)
        seen(pkey (pmul (G{i}, H{j}, p), p) + 1) = true;
      end
    end
    F = monics (p, d);
    for i = 1:numel (F)
      ok = ok && ((~seen(pkey (F{i}, p) + 1)) == isempty (proots (F{i}, p)));
    end
  end
  require (sprintf ('GF(%d), degrees 2 and 3: irreducible exactly when rootless', p), ok);
end
nroots = 0;
for a = 0:1
  if (eval2 (21, a) == 0)
    nroots = nroots + 1;
  end
end
sn ('x^4+x^2+1 over GF(2): roots', nroots, 0);
same ('(x^2+x+1)^2 over GF(2)', ptxt2 (clmul (7, 7)), 'x^4+x^2+1');
irr = sieve2 (10);
irrdeg = degs2 (irr);
stated = {'10 11', '111', '1011 1101', '10011 11001 11111', '100101 101001 101111 110111 111011 111101'};
for d = 1:5
  same (sprintf ('irreducibles of degree %d over GF(2)', d), jb (irr(irrdeg == d)), stated{d});
end
counts = irreducible_counts (10);
same ('counts for degrees 1 to 5', jn (counts(1:5)), '2 1 2 3 6');
same ('counts for degrees 6 to 10', jn (counts(6:end)), '9 18 30 56 99');
trial = zeros (1, 8);
for m = 1:8
  for f = 2 ^ m : 2 ^ (m + 1) - 1
    trial(m) = trial(m) + irreducible_by_trial2 (f);
  end
end
same ('counts for degrees 1 to 8 by trial division', jn (trial), '2 1 2 3 6 9 18 30');
by_trial = [];
for f = 2 : 2 ^ 9 - 1
  if (irreducible_by_trial2 (f))
    by_trial(end+1) = f;
  end
end
upto8 = irr(irrdeg <= 8);
require ('trial division and the sieve list the same irreducibles through degree 8', isequal (by_trial(:)', upto8(:)'));
sn ('polynomials of degree 10 over GF(2)', 2 ^ 10, 1024);
same ('99 out of 1024, times 10, to three places', sprintf ('%.3f', 10 * 99 / 1024), '0.967');
stated = {'2 1 0 1', '4 2 0 2', '8 4 1 3', '16 8 2 6'};
for d = 2:5
  by_x = 0; by_x1 = 0; by_q = 0; kept = 0;
  for f = 2 ^ d : 2 ^ (d + 1) - 1
    odd_weight = mod (weight (f), 2) == 1;
    if (bitand (f, 1) == 0)
      by_x = by_x + 1;
    end
    if (bitand (f, 1) && ~odd_weight)
      by_x1 = by_x1 + 1;
    end
    if (bitand (f, 1) && odd_weight && f ~= 7 && divides2 (7, f))
      by_q = by_q + 1;
    end
    if (any (irr == f))
      kept = kept + 1;
    end
  end
  same (sprintf ('degree %d: crossed by x, by x+1, by x^2+x+1, kept', d), jn ([by_x, by_x1, by_q, kept]), stated{d - 1});
end
same ('x^2+x+1 times x^3+x^2+1', ptxt2 (clmul (7, 13)), 'x^5+x+1');
F = monics (3, 2);
irr3 = {};
for i = 1:numel (F)
  if (isempty (proots (F{i}, 3)))
    irr3{end+1} = ptxt (F{i});
  end
end
same ('monic irreducible quadratics over GF(3)', strjoin (irr3, ' '), 'x^2+1 x^2+x+2 x^2+2x+2');
sn ('monic quadratics over GF(3)', numel (F), 9);
split3 = 0;
for a = 0:2
  for b = a:2
    split3 = split3 + 1;
  end
end
sn ('products (x-a)(x-b) with a <= b over GF(3)', split3, 6);
same ('x^2+x+2 over GF(3) at 0, 1 and 2', jn (pevals ([2, 1, 1], 0:2, 3)), '2 1 2');

section ('Unique factorization');
irr = sieve2 (3);
checks = 0; ok = true;
for p = irr
  for a = 1:15
    for b = 1:15
      checks = checks + 1;
      if (divides2 (p, clmul (a, b)))
        ok = ok && (divides2 (p, a) || divides2 (p, b));
      end
    end
  end
end
sn ('irreducibles p of degree at most 3', numel (irr), 5);
sn ('triples p, a, b with a, b nonzero of degree below 4', checks, 1125);
require ('p | ab forces p | a or p | b on every triple', ok);
same ('x^2+1 = (x+1)^2 divides (x+1)(x+1) but not x+1', yes (divides2 (5, clmul (3, 3)) && ~divides2 (5, 3)), 'yes');
irr6 = sieve2 (6);
hits = build (1, 1, irr6, zeros (1, 128));
sn ('products of irreducibles of total degree at most 6', sum (hits) - 1, 126);
require ('every polynomial of degree 1 to 6 arises exactly once', all (hits(3:end) == 1));
same ('x^4+1 over GF(2), factored', factor_text2 (17), '(x+1)^4');
same ('(x+1)^4 over GF(2)', ptxt2 (clmul (clmul (3, 3), clmul (3, 3))), 'x^4+1');
same ('x^5+1 over GF(2), factored', factor_text2 (33), '(x+1)(x^4+x^3+x^2+x+1)');
require ('x^4+x^3+x^2+x+1 is irreducible', irreducible_by_trial2 (31));

section ('Polynomial arithmetic in hardware');
ms = [3, 8];
stated = {'3 9 4', '8 64 49'};
for i = 1:2
  m = ms(i);
  xors = 0;
  for k = 0 : 2 * m - 2
    cnt = 0;
    for j = 0:m-1
      if (k - j >= 0 && k - j < m)
        cnt = cnt + 1;
      end
    end
    xors = xors + cnt - 1;
  end
  same (sprintf ('%d-bit operands: adder XORs, multiplier ANDs and XORs', m), sprintf ('%d %d %d', m, m * m, xors), stated{i});
end
same ('x^3 mod (x^3+x+1)', ptxt2 (mod2 (8, 11)), 'x+1');
same ('x^4 mod (x^3+x+1)', ptxt2 (mod2 (16, 11)), 'x^2+x');
ok = true;
for c = 0:31
  bt = bitget (c, 1:5);
  wired = bitxor (bt(1), bt(4)) + 2 * bitxor (bitxor (bt(2), bt(4)), bt(5)) + 4 * bitxor (bt(3), bt(5));
  ok = ok && wired == mod2 (c, 11);
end
require ('r0 = c0+c3, r1 = c1+c3+c4, r2 = c2+c4 on all 32 inputs', ok);
sn ('XOR gates in that reduction, one per tap bit', weight (mod2 (8, 11)) + weight (mod2 (16, 11)), 4);
% A shift-register divider for 1011: s2' = s1, s1' = s0 + s2, s0' = in + s2.
state = 0; quotient = '';
for k = 7:-1:0
  fb = bitand (bitshift (state, -2), 1);
  state = bitor (bitand (bitshift (state, 1), 7), bitand (bitshift (147, -k), 1));
  if (fb)
    state = bitxor (state, 3);
  end
  quotient(end+1) = char ('0' + fb);
end
same ('divider for 1011 on 10010011: quotient and remainder', [quotient(4:end) ' ' bits(state)], '10101 100');
sn ('XOR gates in that divider, ones of 1011 below the top', weight (bitxor (11, 8)), 2);
same ('x^3 mod (x^3+x^2+1) and x^4 mod (x^3+x^2+1)', [ptxt2(mod2 (8, 13)) ' ' ptxt2(mod2 (16, 13))], 'x^2+1 x^2+x+1');

section ('Worked examples');
f = [2, 1, 0, 2, 1]; g = [1, 1, 2];
[terms, rems] = division_steps (f, g, 3);
same ('x^4+2x^3+x+2 by 2x^2+x+1 over GF(3): quotient terms', strjoin (terms, ' '), '2x^2 2');
same ('what is left after each step', strjoin (rems, ' '), 'x^2+x+2 2x');
[q, r] = pdivmod (f, g, 3);
same ('quotient and remainder', [ptxt(q) ' ' ptxt(r)], '2x^2+2 2x');
rows = xgcd2 (19, 4);
show_table2 ('inverse of x^2', rows, {'x^4+x+1 x^2 x+1 1 0', 'x^2 x+1 x+1', '1 0 1 x+1 x^2', ...
                                     '0 1 x^2 x^3+x^2+1 x^4+x+1'});
iv = rows(end-1, 4);
same ('inverse of x^2 modulo x^4+x+1', ptxt2 (iv), 'x^3+x^2+1');
same ('x^2(x^3+x^2+1) before reduction', ptxt2 (clmul (4, iv)), 'x^5+x^4+x^2');
same ('reduced modulo x^4+x+1', ptxt2 (mod2 (clmul (4, iv), 19)), '1');
f = 83;
same ('x^6+x^4+x+1: weight and constant term', sprintf ('%d %d', weight (f), bitand (f, 1)), '4 1');
[q1, r1] = divmod2 (f, 3);
same ('divided by x+1', [ptxt2(q1) ' ' ptxt2(r1)], 'x^5+x^4+1 0');
sn ('x^5+x^4+1 at 1', eval2 (q1, 1), 1);
[q2, r2] = divmod2 (q1, 7);
same ('x^5+x^4+1 divided by x^2+x+1', [ptxt2(q2) ' ' ptxt2(r2)], 'x^3+x+1 0');
same ('x^6+x^4+x+1, factored', factor_text2 (f), '(x+1)(x^2+x+1)(x^3+x+1)');
f = [1, 1, 0, 1];
same ('x^3+x+1 over GF(5) at 0 to 4', jn (pevals (f, 0:4, 5)), '1 3 1 1 4');
same ('x^3+x+1 over GF(3) at 0 to 2', jn (pevals (f, 0:2, 3)), '1 0 2');
[q, r] = pdivmod (f, [2, 1], 3);
same ('x^3+x+1 divided by x-1 = x+2 over GF(3)', [ptxt(q) ' ' ptxt(r)], 'x^2+x+2 0');
same ('x^3+x+1 over GF(3), the quadratic factor at 0 to 2', jn (pevals (q, 0:2, 3)), '2 1 2');

section ('Exercises');
same ('101 times 111, carry-less', bits (clmul (5, 7)), '11011');
same ('as a polynomial', ptxt2 (clmul (5, 7)), 'x^4+x^3+x+1');
same ('5 times 7 as integers, in binary', bits (5 * 7), '100011');
[q, r] = divmod2 (67, 7);
same ('x^6+x+1 by x^2+x+1: quotient and remainder', [ptxt2(q) ' ' ptxt2(r)], 'x^4+x^3+x+1 x');
show_table2 ('exercise gcd', xgcd2 (35, 23), {'x^5+x+1 x^4+x^2+x+1 x^3+x^2+1 0', 'x x+1', '1 0 1 x+1', ...
                                            '0 1 x x^2+x+1'});
same ('x^5+x+1, factored', factor_text2 (35), '(x^2+x+1)(x^3+x^2+1)');
same ('x^4+x^2+x+1, factored', factor_text2 (23), '(x+1)(x^3+x^2+1)');
same ('x^3+2x+1 over GF(3) at 0, 1 and 2', jn (pevals ([1, 2, 0, 1], 0:2, 3)), '1 1 1');
same ('roots of x^2-1 in Z_15', jn (proots ([14, 0, 1], 15)), '1 4 11 14');
rows = xgcd2 (11, 2);
same ('inverse of x modulo x^3+x+1', ptxt2 (rows(end-1, 4)), 'x^2+1');
same ('x^4+x^3+x^2+1 over GF(2), factored', factor_text2 (29), '(x+1)(x^3+x+1)');
same ('(2x+1)(3x+1) over Z_6', ptxt (pmul ([1, 2], [1, 3], 6)), '5x+1');
F = monics (5, 2);
n_irr = 0; split5 = 0;
for i = 1:numel (F)
  if (isempty (proots (F{i}, 5)))
    n_irr = n_irr + 1;
  else
    split5 = split5 + 1;
  end
end
sn ('monic irreducible quadratics over GF(5)', n_irr, 10);
same ('monic quadratics over GF(5), and the split ones', sprintf ('%d %d', numel (F), split5), '25 15');

printf ('\nAll chapter 33 numbers reproduced.\n');
