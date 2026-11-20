% Reproduce every number in chapter 34 of Error Correction Coding,
% "Ideals and Quotient Rings", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch34_ideals_and_quotient_rings.m
%
% It needs only core Octave 7 or later, with no pkg load and no Forge
% packages. It is the Octave twin of python/ch34_ideals_and_quotient_rings.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. All the arithmetic is on whole numbers small enough for a double
% to hold exactly, so every comparison is exact. Every finite ring is held as
% an addition table and a multiplication table on the elements 0 to n - 1,
% stored one place to the right because Octave counts from 1. Polynomials over
% GF(2) are whole numbers whose bit i is the coefficient of x^i, so addition is
% bitxor and multiplication is the carry-less product followed by reduction.
% The quotients of GF(3)[x] and GF(5)[x] by x^2 + k hold c + d*alpha as the
% number c + q*d, and the output writes alpha as the letter a. A subset of a
% ring with n elements is a whole number below 2^n whose bit e says whether
% the element e belongs to it.

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


% Polynomials over GF(2) as whole numbers --------------------------------------

function d = deg2 (a)
  if (a == 0)
    d = -1;
  else
    d = numel (dec2bin (a)) - 1;
  end
end

function r = clmul (a, b)
  r = 0;
  while (b > 0)
    if (bitand (b, 1))
      r = bitxor (r, a);
    end
    a = a * 2;
    b = floor (b / 2);
  end
end

function f = mod2 (f, g)
  dg = deg2 (g);
  df = deg2 (f);
  while (f > 0 && df >= dg)
    f = bitxor (f, g * 2 ^ (df - dg));
    df = deg2 (f);
  end
end

function q = div2 (f, g)
  q = 0;
  dg = deg2 (g);
  df = deg2 (f);
  while (f > 0 && df >= dg)
    q = q + 2 ^ (df - dg);
    f = bitxor (f, g * 2 ^ (df - dg));
    df = deg2 (f);
  end
end

function a = gcd2 (a, b)
  while (b > 0)
    t = mod2 (a, b);
    a = b;
    b = t;
  end
end

function r = mulmod (a, b, p)
  r = mod2 (clmul (a, b), p);
end

function s = bits (a, w)
  s = dec2bin (a, w);
end

function s = words (v, w)
  if (isempty (v))
    s = 'none';
  else
    s = strjoin (arrayfun (@(x) dec2bin (x, w), v(:)', 'UniformOutput', false), ' ');
  end
end

function s = binlist (v)
  s = strjoin (arrayfun (@(x) dec2bin (x), v(:)', 'UniformOutput', false), ' ');
end

function w = weight (a)
  w = sum (dec2bin (a) == '1');
end

function ok = is_irreducible2 (p)
  m = deg2 (p);
  ok = true;
  for d = 2:(2 ^ (floor (m / 2) + 1) - 1)
    if (mod2 (p, d) == 0)
      ok = false;
      return;
    end
  end
end

function out = euclid_remainders (a, b)
  out = [];
  while (b > 0)
    t = mod2 (a, b);
    a = b;
    b = t;
    if (b > 0)
      out(end+1) = b;
    end
  end
end


% Finite rings as tables -------------------------------------------------------

function [add, mul, one] = ring_zn (n)
  e = 0:n-1;
  add = mod (e' + e, n);
  mul = mod (e' * e, n);
  one = mod (1, n);
end

function [add, mul, one] = ring_mod2 (p)
  n = 2 ^ deg2 (p);
  add = zeros (n);
  mul = zeros (n);
  for a = 0:n-1
    for b = 0:n-1
      add(a+1, b+1) = bitxor (a, b);
      mul(a+1, b+1) = mulmod (a, b, p);
    end
  end
  one = 1;
end

function r = quad_mul (e, f, q, k)
  c = mod (e, q);
  d = floor (e / q);
  g = mod (f, q);
  h = floor (f / q);
  r = mod (c * g - k * d * h, q) + q * mod (c * h + d * g, q);
end

function [add, mul, one] = ring_quad (q, k)
  n = q * q;
  add = zeros (n);
  mul = zeros (n);
  for a = 0:n-1
    for b = 0:n-1
      add(a+1, b+1) = mod (mod (a, q) + mod (b, q), q) + q * mod (floor (a / q) + floor (b / q), q);
      mul(a+1, b+1) = quad_mul (a, b, q, k);
    end
  end
  one = 1;
end

function s = qtxt (e, q)
  c = mod (e, q);
  d = floor (e / q);
  if (d == 0)
    s = num (c);
    return;
  end
  if (d == 1)
    dpart = 'a';
  else
    dpart = sprintf ('%da', d);
  end
  if (c == 0)
    s = dpart;
  else
    s = sprintf ('%d+%s', c, dpart);
  end
end

function [add, mul, one] = bitwise_ring ()
  e = 0:3;
  add = bitxor (repmat (e', 1, 4), repmat (e, 4, 1));
  mul = bitand (repmat (e', 1, 4), repmat (e, 4, 1));
  one = 3;
end

function [add, mul, one] = f4_ring ()
  % F_4 exactly as Rings and Fields gives it, on the symbols 0, 1, a, b held as 0, 1, 2, 3.
  add = [0 1 2 3; 1 0 3 2; 2 3 0 1; 3 2 1 0];
  mul = [0 0 0 0; 0 1 2 3; 0 2 3 1; 0 3 1 2];
  one = 1;
end

function ok = is_comm_ring (add, mul, one)
  n = rows (add);
  ok = false;
  for a = 0:n-1
    if (add(a+1, 1) ~= a || mul(a+1, one+1) ~= a)
      return;
    end
    if (~any (add(a+1, :) == 0))
      return;
    end
    for b = 0:n-1
      if (add(a+1, b+1) ~= add(b+1, a+1) || mul(a+1, b+1) ~= mul(b+1, a+1))
        return;
      end
      ab = add(a+1, b+1);
      if (any (add(ab+1, :) ~= add(a+1, add(b+1, :)+1)))
        return;
      end
      mb = mul(a+1, b+1);
      if (any (mul(mb+1, :) ~= mul(a+1, mul(b+1, :)+1)))
        return;
      end
      if (any (mul(a+1, add(b+1, :)+1) ~= add(mb+1, mul(a+1, :)+1)))
        return;
      end
    end
  end
  ok = true;
end

function u = unit_list (mul, one)
  u = (find (any (mul == one, 2)) - 1)';
end

function z = zero_divisor_list (mul)
  n = rows (mul);
  z = [];
  for a = 1:n-1
    if (any (mul(a+1, 2:n) == 0))
      z(end+1) = a;
    end
  end
end

function out = nilpotent_list (mul)
  n = rows (mul);
  out = [];
  for a = 1:n-1
    p = a;
    for t = 1:n
      if (p == 0)
        break;
      end
      p = mul(p+1, a+1);
    end
    if (p == 0)
      out(end+1) = a;
    end
  end
end

function b = inverse_in (mul, one, a)
  b = find (mul(a+1, :) == one, 1) - 1;
  if (isempty (b))
    b = -1;
  end
end

function ok = is_field (mul, one)
  ok = numel (unit_list (mul, one)) == rows (mul) - 1;
end

function k = characteristic (add, one)
  k = 1;
  s = one;
  while (s ~= 0)
    s = add(s+1, one+1);
    k = k + 1;
  end
end

function found = ideal_masks (add, mul)
  n = rows (add);
  neg = zeros (1, n);
  for a = 0:n-1
    neg(a+1) = find (add(a+1, :) == 0, 1) - 1;
  end
  found = [];
  for mask = 1:2^n-1
    inside = logical (bitget (mask, 1:n));
    members = find (inside);
    s = add(members, members);
    if (~all (inside(s(:)+1)))
      continue;
    end
    if (~all (inside(neg(members)+1)))
      continue;
    end
    t = mul(:, members);
    if (all (inside(t(:)+1)))
      found(end+1) = mask;
    end
  end
end

function found = subgroup_masks (add)
  n = rows (add);
  neg = zeros (1, n);
  for a = 0:n-1
    neg(a+1) = find (add(a+1, :) == 0, 1) - 1;
  end
  found = [];
  for mask = 1:2^n-1
    inside = logical (bitget (mask, 1:n));
    members = find (inside);
    s = add(members, members);
    if (all (inside(s(:)+1)) && all (inside(neg(members)+1)))
      found(end+1) = mask;
    end
  end
end

function m = members_of (mask)
  m = find (bitget (mask, 1:53)) - 1;
end

function out = maximal_masks (masks, n)
  full = 2 ^ n - 1;
  out = [];
  for m = masks
    if (m == full)
      continue;
    end
    between = false;
    for j = masks
      if (j ~= m && j ~= full && bitand (j, m) == m)
        between = true;
      end
    end
    if (~between)
      out(end+1) = m;
    end
  end
end

function out = prime_masks (masks, mul)
  n = rows (mul);
  full = 2 ^ n - 1;
  out = [];
  for m = masks
    if (m == full)
      continue;
    end
    inside = logical (bitget (m, 1:n));
    prod_in = inside(mul + 1);
    either = inside' | inside;
    if (all (either(prod_in)))
      out(end+1) = m;
    end
  end
end

function d = generator_zn (mask, n)
  e = members_of (mask);
  e = e(e > 0);
  if (isempty (e))
    d = 0;
  else
    d = min (e);
  end
end

function v = generators_zn (masks, n)
  v = sort (arrayfun (@(m) generator_zn (m, n), masks));
end

function v = sizes_of (masks)
  v = sort (arrayfun (@(m) numel (members_of (m)), masks));
end

function [qadd, qmul, coset_of, ok] = quotient (add, mul, mask)
  n = rows (add);
  ideal = members_of (mask);
  coset_of = -ones (1, n);
  reps = [];
  for a = 0:n-1
    if (coset_of(a+1) < 0)
      for i = ideal
        coset_of(add(a+1, i+1)+1) = numel (reps);
      end
      reps(end+1) = a;
    end
  end
  k = numel (reps);
  qadd = -ones (k);
  qmul = -ones (k);
  ok = true;
  for a = 0:n-1
    for b = 0:n-1
      s = coset_of(add(a+1, b+1)+1);
      p = coset_of(mul(a+1, b+1)+1);
      ca = coset_of(a+1);
      cb = coset_of(b+1);
      if (qadd(ca+1, cb+1) < 0)
        qadd(ca+1, cb+1) = s;
      elseif (qadd(ca+1, cb+1) ~= s)
        ok = false;
      end
      if (qmul(ca+1, cb+1) < 0)
        qmul(ca+1, cb+1) = p;
      elseif (qmul(ca+1, cb+1) ~= p)
        ok = false;
      end
    end
  end
end

function found = matchings (add_a, mul_a, add_b, mul_b)
  found = zeros (0, 4);
  for m0 = 0:3
    for m1 = 0:3
      for m2 = 0:3
        for m3 = 0:3
          m = [m0 m1 m2 m3];
          if (numel (unique (m)) < 4)
            continue;
          end
          if (isequal (m(add_a + 1), add_b(m+1, m+1)) && isequal (m(mul_a + 1), mul_b(m+1, m+1)))
            found(end+1, :) = m;
          end
        end
      end
    end
  end
end

function c = bijection_count ()
  c = 0;
  for m0 = 0:3
    for m1 = 0:3
      for m2 = 0:3
        for m3 = 0:3
          if (numel (unique ([m0 m1 m2 m3])) == 4)
            c = c + 1;
          end
        end
      end
    end
  end
end


% Unit counts of every quotient GF(2)[x]/(p(x)) with p of degree 1 to 6 ---------

function c = unit_count2 (p)
  c = 0;
  for f = 1:2^deg2(p)-1
    if (gcd2 (p, f) == 1)
      c = c + 1;
    end
  end
end

function counts = unit_counts (top)
  counts = cell (1, top);
  for m = 1:top
    counts{m} = arrayfun (@unit_count2, 2^m:2^(m+1)-1);
  end
end


% The checks -------------------------------------------------------------------

function check_opening ()
  section ('The opening: two rules for multiplying 3-bit words');
  p = 11;
  q = 9;
  same ('raw product of 010 and 101', bits (clmul (2, 5), 4), '1010');
  same ('010 * 101 with x^3 = x + 1', bits (mulmod (2, 5, p), 3), '001');
  [add_p, mul_p, one_p] = ring_mod2 (p);
  same ('inverses of 001 to 111 with x^3 = x + 1', words (arrayfun (@(a) inverse_in (mul_p, 1, a), 1:7), 3), ...
        '001 101 110 111 010 011 100');
  sn ('nonzero words with an inverse, x^3 = x + 1', numel (unit_list (mul_p, 1)), 7);
  same ('x^3 + x + 1 at x = 0 and x = 1', jn ([bitand(p, 1), mod(weight (p), 2)]), '1 1');
  same ('raw product of 011 and 111', bits (clmul (3, 7), 4), '1001');
  same ('011 * 111 with x^3 = 1', bits (mulmod (3, 7, q), 3), '000');
  same ('(x + 1)(x^2 + x + 1)', bits (clmul (3, 7), 4), '1001');
  [add_q, mul_q, one_q] = ring_mod2 (q);
  u = unit_list (mul_q, 1);
  same ('units with x^3 = 1', words (u, 3), '001 010 100');
  same ('their inverses', words (arrayfun (@(a) inverse_in (mul_q, 1, a), u), 3), '001 100 010');
  same ('zero divisors with x^3 = 1', words (zero_divisor_list (mul_q), 3), '011 101 110 111');
  require ('each zero divisor times its partner is 000', ...
           mul_q(4, 8) == 0 && mul_q(6, 8) == 0 && mul_q(7, 8) == 0 && mul_q(8, 4) == 0);
  require ('both rules satisfy the ring axioms on all 512 triples', ...
           is_comm_ring (add_p, mul_p, one_p) && is_comm_ring (add_q, mul_q, one_q));
end

function check_ideals ()
  section ('Ideals');
  [add, mul, one] = ring_zn (12);
  masks = ideal_masks (add, mul);
  sn ('ideals of Z_12 among its 4095 nonempty subsets', numel (masks), 6);
  same ('their generators', jn (generators_zn (masks, 12)), '0 1 2 3 4 6');
  same ('their sizes', jn (sizes_of (masks)), '1 2 3 4 6 12');
  d = 1:12;
  same ('divisors of 12', jn (d(mod (12, d) == 0)), '1 2 3 4 6 12');
  ok = true;
  for m = masks
    g = generator_zn (m, 12);
    ok = ok && isequal (members_of (m), unique (mod (g * (0:11), 12)));
  end
  require ('each ideal of Z_12 is the multiples of its generator', ok);
  require ('the additive subgroups of Z_12 are exactly its ideals', isequal (subgroup_masks (add), masks));
  units12 = unit_list (mul, 1);
  same ('units of Z_12', jn (units12), '1 5 7 11');
  ok = true;
  for m = masks
    if (any (bitget (m, units12 + 1)) && m ~= 4095)
      ok = false;
    end
  end
  require ('an ideal of Z_12 holding a unit is all of Z_12', ok);
  k = -36:36;
  same ('kernel of reduction mod 12 within -36..36', jn (k(mod (k, 12) == 0)), '-36 -24 -12 0 12 24 36');
  two = [];
  agree = true;
  for n = 2:12
    [a2, m2, o2] = ring_zn (n);
    exactly_two = numel (ideal_masks (a2, m2)) == 2;
    agree = agree && (exactly_two == is_field (m2, o2));
    if (exactly_two)
      two(end+1) = n;
    end
  end
  same ('n from 2 to 12 with exactly two ideals in Z_n', jn (two), '2 3 5 7 11');
  require ('two ideals exactly when Z_n is a field, n = 2..12', agree);
  mult = unique (arrayfun (@(u) clmul (5, u), 0:3));
  same ('multiples of x^2 + 1 of degree below 4', words (mult, 4), '0000 0101 1010 1111');
  ok = true;
  for a = mult
    for b = mult
      ok = ok && any (mult == bitxor (a, b));
    end
  end
  require ('they are closed under addition', ok);
end

function check_principal ()
  section ('Principal ideals');
  f = 9;
  g = 5;
  same ('gcd(x^3 + 1, x^2 + 1)', bits (gcd2 (f, g), 2), '11');
  same ('Euclid''s remainders', words (euclid_remainders (f, g), 2), '11');
  same ('1 (x^3 + 1) + x (x^2 + 1)', bits (bitxor (clmul (1, f), clmul (2, g)), 2), '11');
  combos = [];
  for u = 0:7
    for v = 0:15
      combos(end+1) = bitxor (clmul (u, f), clmul (v, g));
    end
  end
  combos = unique (combos);
  sn ('distinct u f + v g, deg u < 3, deg v < 4', numel (combos), 32);
  mult = [];
  for h = 0:63
    if (mod2 (h, 3) == 0)
      mult(end+1) = h;
    end
  end
  require ('they are the multiples of x + 1 of degree below 6', isequal (combos, mult));
  constants = [];
  for a0 = -2:2
    for a1 = -2:2
      for b0 = -2:2
        for b1 = -2:2
          constants(end+1) = 2 * a0;
        end
      end
    end
  end
  sn ('combinations 2a + xb tried in Z[x]', numel (constants), 625);
  require ('each has an even constant term, so none is 1', all (mod (constants, 2) == 0));
  sn ('remainder of the coefficient 1 of x on division by 2', mod (1, 2), 1);
end

function check_quotients ()
  section ('Cosets and quotient rings');
  p = 7;
  stated = {'0000 0111 1001 1110', '0001 0110 1000 1111', '0010 0101 1011 1100', '0011 0100 1010 1101'};
  for r = 0:3
    col = [];
    for h = 0:15
      if (mod2 (h, p) == r)
        col(end+1) = h;
      end
    end
    same (sprintf ('the coset %s + I among degree below 4', bits (r, 2)), words (col, 4), stated{r+1});
  end
  same ('x^3 * x^2 reduced mod x^2 + x + 1', bits (mulmod (8, 4, p), 2), '11');
  same ('1 * (x + 1) reduced mod x^2 + x + 1', bits (mulmod (1, 3, p), 2), '11');
  [add, mul, one] = ring_zn (12);
  [qadd, qmul, coset_of, ok] = quotient (add, mul, 1 + 16 + 256);
  sn ('cosets of (4) in Z_12', rows (qadd), 4);
  same ('size of each coset', jn (unique (arrayfun (@(c) sum (coset_of == c), 0:3))), '3');
  require ('sums and products of cosets of (4) are well defined', ok);
  [a4, m4, o4] = ring_zn (4);
  require ('a + (4) -> a mod 4 matches Z_12/(4) with Z_4', isequal (qadd, a4) && isequal (qmul, m4));
  [add_w, mul_w, one_w] = bitwise_ring ();
  s = [0 3];
  ok = any (s == one_w);
  for a = s
    for b = s
      ok = ok && any (s == add_w(a+1, b+1)) && any (s == mul_w(a+1, b+1));
    end
  end
  require ('{00, 11} is a subring of the 2-bit words', ok);
  same ('01 AND 11, which leaves {00, 11}', bits (mul_w(2, 4), 2), '01');
  cosets = {[0 3], [1 2]};
  which = [0 1 1 0];
  sum_ok = 0;
  mul_ok = 0;
  for ci = 1:2
    for di = 1:2
      sums = [];
      prods = [];
      for a = cosets{ci}
        for b = cosets{di}
          sums(end+1) = which(add_w(a+1, b+1)+1);
          prods(end+1) = which(mul_w(a+1, b+1)+1);
        end
      end
      sum_ok = sum_ok + (numel (unique (sums)) == 1);
      mul_ok = mul_ok + (numel (unique (prods)) == 1);
    end
  end
  sn ('coset pairs with a well-defined sum', sum_ok, 4);
  sn ('coset pairs with a well-defined product', mul_ok, 1);
  same ('01 AND 01 and 01 AND 10', [bits(mul_w(2, 2), 2) ' ' bits(mul_w(2, 3), 2)], '01 00');
  [a9, m9, o9] = ring_quad (3, 1);
  [a2, m2, o2] = ring_mod2 (7);
  sn ('elements of GF(2)[x]/(x^2 + x + 1)', rows (a2), 4);
  [a3, m3, o3] = ring_mod2 (11);
  sn ('elements of GF(2)[x]/(x^3 + x + 1)', rows (a3), 8);
  sn ('elements of GF(3)[x]/(x^2 + 1)', rows (a9), 9);
end

function check_four ()
  section ('The four quotients with four elements');
  moduli = [7 6 5 4];
  names = {'x^2 + x + 1', 'x^2 + x', 'x^2 + 1', 'x^2'};
  stated = {{'00 10 11 01', '00 11 01 10', '01 10 11', 'none', 'none', 2}, ...
            {'00 10 10 00', '00 11 00 11', '01', '10 11', 'none', 4}, ...
            {'00 10 01 11', '00 11 11 00', '01 10', '11', '11', 3}, ...
            {'00 10 00 10', '00 11 10 01', '01 11', '10', '10', 3}};
  adds = cell (1, 4);
  muls = cell (1, 4);
  for i = 1:4
    [add, mul, one] = ring_mod2 (moduli(i));
    adds{i} = add;
    muls{i} = mul;
    st = stated{i};
    nm = names{i};
    require (sprintf ('mod %s: the ring axioms on all 64 triples', nm), is_comm_ring (add, mul, one));
    same (sprintf ('mod %s: the row of x', nm), words (mul(3, :), 2), st{1});
    same (sprintf ('mod %s: the row of x + 1', nm), words (mul(4, :), 2), st{2});
    same (sprintf ('mod %s: units', nm), words (unit_list (mul, one), 2), st{3});
    same (sprintf ('mod %s: zero divisors', nm), words (zero_divisor_list (mul), 2), st{4});
    same (sprintf ('mod %s: nonzero nilpotents', nm), words (nilpotent_list (mul), 2), st{5});
    sn (sprintf ('mod %s: ideals', nm), numel (ideal_masks (add, mul)), st{6});
  end
  sn ('one-to-one maps of four elements onto four', bijection_count (), 24);
  sym = {'0', '1', 'a', 'b'};
  [af, mf, of] = f4_ring ();
  found = matchings (adds{1}, muls{1}, af, mf);
  sn ('matchings of mod x^2 + x + 1 with F_4', rows (found), 2);
  same ('images of x and x + 1, first matching', [sym{found(1, 3)+1} ' ' sym{found(1, 4)+1}], 'a b');
  same ('images of x and x + 1, second matching', [sym{found(2, 3)+1} ' ' sym{found(2, 4)+1}], 'b a');
  [aw, mw, ow] = bitwise_ring ();
  found = matchings (adds{2}, muls{2}, aw, mw);
  sn ('matchings of mod x^2 + x with the 2-bit words', rows (found), 2);
  same ('image of 1 in both', words (unique (found(:, 2)), 2), '11');
  same ('images of x and x + 1, first matching', words (found(1, 3:4), 2), '01 10');
  same ('images of x and x + 1, second matching', words (found(2, 3:4), 2), '10 01');
  counts = [];
  for i = 1:4
    for j = i+1:4
      counts(end+1) = rows (matchings (adds{i}, muls{i}, adds{j}, muls{j}));
    end
  end
  same ('matchings between each pair of the four', jn (counts), '0 0 0 0 0 1');
  found = matchings (adds{3}, muls{3}, adds{4}, muls{4});
  same ('mod x^2 + 1 to mod x^2: images of x, x + 1', words (found(1, 3:4), 2), '11 10');
  [a4, m4, o4] = ring_zn (4);
  same ('matchings of Z_4 with each of the four', ...
        jn (arrayfun (@(i) rows (matchings (a4, m4, adds{i}, muls{i})), 1:4)), '0 0 0 0');
  sn ('characteristic of Z_4', characteristic (a4, 1), 4);
  same ('characteristic of each of the four', jn (arrayfun (@(i) characteristic (adds{i}, 1), 1:4)), '2 2 2 2');
  sn ('1 + 1 in Z_4', a4(2, 2), 2);
end

function v = ev3 (f, a)
  v = mod (mod (f, 3) + mod (floor (f / 3), 3) * a + floor (f / 9) * a * a, 3);
end

function e = mul3_linear (u)
  u0 = mod (u, 3);
  u1 = floor (u / 3);
  e = u0 + 3 * mod (u0 + u1, 3) + 9 * u1;
end

function v = pair_map (f)
  v = 2 * mod (weight (f), 2) + bitand (f, 1);
end

function check_first_iso ()
  section ('The first isomorphism theorem');
  all27 = 0:26;
  values = arrayfun (@(f) ev3 (f, 2), all27);
  kernel = all27(values == 0);
  sn ('kernel of evaluation at 2, degree below 3 over GF(3)', numel (kernel), 9);
  require ('it is the multiples of x + 1 = x - 2 of degree below 3', ...
           isequal (kernel, unique (arrayfun (@mul3_linear, 0:8))));
  sizes = arrayfun (@(v) sum (values == v), 0:2);
  same ('fiber sizes over the values 0, 1, 2', jn (sizes), '9 9 9');
  ok = true;
  for v = 0:2
    coset = unique (mod (mod (kernel, 3) + v, 3) + 3 * floor (kernel / 3));
    ok = ok && isequal (all27(values == v), coset);
  end
  require ('each fiber is the coset v + kernel', ok);
  all16 = 0:15;
  pv = arrayfun (@pair_map, all16);
  kern = all16(pv == 0);
  same ('kernel of f -> (f(1), f(0)), degree below 4', words (kern, 4), '0000 0110 1010 1100');
  require ('it is the multiples of x^2 + x of degree below 4', ...
           isequal (kern, unique (arrayfun (@(u) clmul (6, u), 0:3))));
  same ('values taken by f -> (f(1), f(0))', words (unique (pv), 2), '00 01 10 11');
  induced = arrayfun (@pair_map, 0:3);
  same ('images of 1, x and x + 1', words (induced(2:4), 2), '11 10 01');
  [a6, m6, o6] = ring_mod2 (6);
  [aw, mw, ow] = bitwise_ring ();
  found = matchings (a6, m6, aw, mw);
  require ('the induced map is the second matching found', isequal (found(2, :), induced));
end

function mask = ideal_zn (d, n)
  mask = sum (2 .^ unique (mod (d * (0:n-1), n)));
end

function check_maximal ()
  section ('Maximal and prime ideals');
  [add, mul, one] = ring_zn (12);
  masks = ideal_masks (add, mul);
  same ('generators of the maximal ideals of Z_12', jn (generators_zn (maximal_masks (masks, 12), 12)), '2 3');
  same ('generators of the prime ideals of Z_12', jn (generators_zn (prime_masks (masks, mul), 12)), '2 3');
  verdicts = {};
  for d = [2 3 4 6]
    [qadd, qmul, coset_of, ok] = quotient (add, mul, ideal_zn (d, 12));
    verdicts{end+1} = sprintf ('%d:%s', rows (qadd), yes (ok && is_field (qmul, 1)));
  end
  same ('size and field verdict of Z_12/(d), d = 2, 3, 4, 6', strjoin (verdicts, ' '), '2:yes 3:yes 4:no 6:no');
  above = masks(bitand (masks, 273) == 273);
  [a4, m4, o4] = ring_zn (4);
  sn ('ideals of Z_12 containing (4)', numel (above), 3);
  sn ('ideals of Z_12/(4), which is Z_4', numel (ideal_masks (a4, m4)), 3);
  [a8, m8r, o8] = ring_mod2 (9);
  m8 = ideal_masks (a8, m8r);
  sn ('ideals of GF(2)[x]/(x^3 + 1)', numel (m8), 4);
  same ('their sizes', jn (sizes_of (m8)), '1 2 4 8');
  same ('sizes of the maximal ones', jn (sizes_of (maximal_masks (m8, 8))), '2 4');
  same ('the ideal (x + 1)', words (unique (arrayfun (@(u) mulmod (3, u, 9), 0:7)), 3), '000 011 101 110');
  divs = [];
  for d = 1:15
    if (mod2 (9, d) == 0)
      divs(end+1) = d;
    end
  end
  same ('monic divisors of x^3 + 1', binlist (divs), '1 11 111 1001');
  agree = true;
  irreducibles = [];
  for m = 1:4
    for p = 2^m:2^(m+1)-1
      [ap, mp, op] = ring_mod2 (p);
      irr = is_irreducible2 (p);
      agree = agree && (is_field (mp, op) == irr);
      if (irr)
        irreducibles(end+1) = p;
      end
    end
  end
  require ('field exactly when p is irreducible, all 30 of degree 1 to 4', agree);
  same ('the irreducible p of degree 1 to 4', binlist (irreducibles), '10 11 111 1011 1101 10011 11001 11111');
  same ('(x^3 + x + 1) + x (x^2 + 1)', bits (bitxor (11, clmul (2, 5)), 1), '1');
  [ap, mul_p, op] = ring_mod2 (11);
  same ('inverse of 101 modulo x^3 + x + 1', bits (inverse_in (mul_p, 1, 5), 3), '010');
  same ('inverse of 010 modulo x^3 + x + 1', bits (inverse_in (mul_p, 1, 2), 3), '101');
  [aq, mul_q, oq] = ring_mod2 (9);
  same ('units and zero divisors modulo x^3 + 1', ...
        jn ([numel(unit_list (mul_q, 1)), numel(zero_divisor_list (mul_q))]), '3 4');
  by3 = arrayfun (@(h) mod2 (h, 3) == 0, 0:7);
  by7 = arrayfun (@(h) mod2 (h, 7) == 0, 0:7);
  same ('words divisible by x + 1, by x^2 + x + 1, by both', jn ([sum(by3), sum(by7), sum(by3 & by7)]), '4 2 1');
  sn ('8 - 4 - 2 + 1', 8 - 4 - 2 + 1, 3);
  counts = unit_counts (6);
  same ('moduli of each degree 1 to 6', jn (cellfun (@numel, counts)), '2 4 8 16 32 64');
  reach = arrayfun (@(m) sum (counts{m} == 2^m - 1), 1:6);
  same ('moduli reaching 2^m - 1 units, each degree', jn (reach), '2 1 2 3 6 9');
  ok = true;
  most = 0;
  for m = 1:6
    for p = 2^m:2^(m+1)-1
      irr = is_irreducible2 (p);
      ok = ok && ((counts{m}(p - 2^m + 1) == 2^m - 1) == irr);
      if (m == 6 && ~irr)
        most = max (most, counts{m}(p - 2^m + 1));
      end
    end
  end
  require ('2^m - 1 units exactly when p is irreducible, degree 1 to 6', ok);
  sn ('most units of a reducible p of degree 6', most, 56);
  sn ('fewest units of any p of degree 6', min (counts{6}), 12);
  ok = true;
  for m = 1:4
    for p = 2^m:2^(m+1)-1
      [ap, mp, op] = ring_mod2 (p);
      ok = ok && (counts{m}(p - 2^m + 1) == numel (unit_list (mp, 1)));
    end
  end
  require ('the gcd count agrees with the pair count, degree 1 to 4', ok);
end

function check_root ()
  section ('A root by construction');
  same ('alpha^2 in GF(2)[x]/(x^2 + x + 1)', bits (mulmod (2, 2, 7), 2), '11');
  same ('alpha^2 + alpha + 1 there', bits (bitxor (bitxor (mulmod (2, 2, 7), 2), 1), 2), '00');
  a = 0:2;
  same ('x^2 + 1 at 0, 1, 2 over GF(3)', jn (mod (a .* a + 1, 3)), '1 2 2');
  same ('squares in GF(3)', jn (unique (mod (a .* a, 3))), '0 1');
  [a9, m9, o9] = ring_quad (3, 1);
  require ('GF(3)[x]/(x^2 + 1): the ring axioms on all 729 triples', is_comm_ring (a9, m9, o9));
  sn ('its units, by all 81 pairs', numel (unit_list (m9, 1)), 8);
  same ('alpha^2 there', qtxt (m9(4, 4), 3), '2');
  same ('(1 + 2 alpha)(2 + alpha)', qtxt (m9(1 + 3 * 2 + 1, 2 + 3 + 1), 3), '2a');
  same ('inverse of 1 + alpha', qtxt (inverse_in (m9, 1, 4), 3), '2+a');
  f = clmul (7, 11);
  same ('(x^2 + x + 1)(x^3 + x + 1)', dec2bin (f), '110001');
  same ('x^5 + x^4 + 1 at x = 0 and x = 1', jn ([bitand(f, 1), mod(weight (f), 2)]), '1 1');
  sn ('x^5 + x^4 + 1 reduced mod x^2 + x + 1', mod2 (f, 7), 0);
  sn ('x^5 + x^4 + 1 reduced mod x^3 + x + 1', mod2 (f, 11), 0);
end

function [images, taps] = network (p)
  m = deg2 (p);
  images = arrayfun (@(k) mod2 (2 ^ k, p), 0:2*m-2);
  taps = cell (1, m);
  for j = 0:m-1
    taps{j+1} = find (bitget (images, j + 1)) - 1;
  end
end

function r = run_network (c, taps)
  r = 0;
  for j = 1:numel (taps)
    b = 0;
    for k = taps{j}
      b = bitxor (b, bitget (c, k + 1));
    end
    r = r + b * 2 ^ (j - 1);
  end
end

function check_hardware ()
  section ('Hardware');
  sn ('XOR gates to add two 3-bit words', 3, 3);
  sn ('AND gates in the 3-bit product', 3 * 3, 9);
  columns = zeros (1, 5);
  for i = 0:2
    for j = 0:2
      columns(i+j+1) = columns(i+j+1) + 1;
    end
  end
  same ('partial products in columns 0 to 4', jn (columns), '1 2 3 2 1');
  sn ('XOR gates to sum the columns', sum (columns - 1), 4);
  moduli = [11 9];
  names = {'x^3 + x + 1', 'x^3 + 1'};
  stated = {{'011 110', 'c0 c3|c1 c3 c4|c2 c4', 4}, {'001 010', 'c0 c3|c1 c4|c2', 2}};
  for i = 1:2
    p = moduli(i);
    nm = names{i};
    st = stated{i};
    [images, taps] = network (p);
    same (sprintf ('x^3 and x^4 modulo %s', nm), words (images(4:5), 3), st{1});
    parts = cellfun (@(t) strjoin (arrayfun (@(k) sprintf ('c%d', k), t, 'UniformOutput', false), ' '), ...
                     taps, 'UniformOutput', false);
    same (sprintf ('taps of r0|r1|r2 modulo %s', nm), strjoin (parts, '|'), st{2});
    sn (sprintf ('XOR gates in the reduction modulo %s', nm), sum (cellfun (@numel, taps) - 1), st{3});
    ok = true;
    for c = 0:31
      ok = ok && (run_network (c, taps) == mod2 (c, p));
    end
    require (sprintf ('the network matches reduction on all 32 products, %s', nm), ok);
    ok = true;
    for a = 0:7
      for b = 0:7
        ok = ok && (run_network (clmul (a, b), taps) == mulmod (a, b, p));
      end
    end
    require (sprintf ('AND array, XOR trees and network on all 64 pairs, %s', nm), ok);
  end
  sn ('XOR gates in the whole multiplier, x^3 + x + 1', 4 + 4, 8);
  sn ('XOR gates in the whole multiplier, x^3 + 1', 4 + 2, 6);
end

function check_worked ()
  section ('Worked examples');
  [a8, m8, o8] = ring_zn (8);
  masks = ideal_masks (a8, m8);
  same ('generators of the ideals of Z_8', jn (generators_zn (masks, 8)), '0 1 2 4');
  same ('generators of its maximal ideals', jn (generators_zn (maximal_masks (masks, 8), 8)), '2');
  [qadd, qmul, coset_of, ok] = quotient (a8, m8, 85);
  same ('Z_8/(2): size and field', [num(rows (qadd)) ' ' yes(ok && is_field (qmul, 1))], '2 yes');
  [qadd, qmul, coset_of, ok] = quotient (a8, m8, 17);
  same ('Z_8/(4): size and field', [num(rows (qadd)) ' ' yes(ok && is_field (qmul, 1))], '4 no');
  p = 19;
  same ('Euclid''s remainders for x^4 + x + 1 and x^2 + 1', binlist (euclid_remainders (p, 5)), '10 1');
  same ('quotients at the two steps', binlist ([div2(p, 5), div2(5, 2)]), '101 10');
  sn ('(x^2 + 1)(x^3 + x + 1) + x (x^4 + x + 1)', bitxor (clmul (5, 11), clmul (2, p)), 1);
  [ap, mp, op] = ring_mod2 (p);
  same ('inverse of x^2 + 1 modulo x^4 + x + 1', bits (inverse_in (mp, 1, 5), 4), '1011');
  q = 21;
  same ('(x^2 + x + 1)^2', dec2bin (clmul (7, 7)), '10101');
  [aq, mul_q, oq] = ring_mod2 (q);
  same ('units, zero divisors modulo x^4 + x^2 + 1', ...
        jn ([numel(unit_list (mul_q, 1)), numel(zero_divisor_list (mul_q))]), '12 3');
  same ('its zero divisors', words (zero_divisor_list (mul_q), 4), '0111 1001 1110');
  same ('its nonzero nilpotents', words (nilpotent_list (mul_q), 4), '0111 1001 1110');
  a = 0:4;
  same ('squares in GF(5)', jn (unique (mod (a .* a, 5))), '0 1 4');
  sn ('-2 in GF(5)', mod (-2, 5), 3);
  [a25, m25, o25] = ring_quad (5, 2);
  require ('GF(5)[x]/(x^2 + 2): the ring axioms on all 15625 triples', is_comm_ring (a25, m25, o25));
  sn ('its units, by all 625 pairs', numel (unit_list (m25, 1)), 24);
  same ('alpha^2 there', qtxt (m25(6, 6), 5), '3');
  same ('inverse of alpha', qtxt (inverse_in (m25, 1, 5), 5), '2a');
  same ('(1 + alpha)^2', qtxt (m25(7, 7), 5), '4+2a');
  [a12, m12, o12] = ring_zn (12);
  [a4, m4, o4] = ring_zn (4);
  r = mod (0:11, 4);
  require ('a -> a mod 4 keeps every sum and product of Z_12', ...
           isequal (mod (a12, 4), a4(r+1, r+1)) && isequal (mod (m12, 4), m4(r+1, r+1)));
  e = 0:11;
  same ('its kernel', jn (e(mod (e, 4) == 0)), '0 4 8');
  same ('its fiber sizes', jn (arrayfun (@(v) sum (mod (e, 4) == v), 0:3)), '3 3 3 3');
end

function check_exercises ()
  section ('Exercises');
  [a10, m10, o10] = ring_zn (10);
  masks = ideal_masks (a10, m10);
  same ('generators of the ideals of Z_10', jn (generators_zn (masks, 10)), '0 1 2 5');
  same ('generators of its maximal ideals', jn (generators_zn (maximal_masks (masks, 10), 10)), '2 5');
  ok = true;
  for f = 0:63
    ok = ok && ((mod (weight (f), 2) == 0) == (mod2 (f, 3) == 0));
  end
  require ('f(1) = 0 exactly for the multiples of x + 1, degree below 6', ok);
  [ap, mp, op] = ring_mod2 (19);
  same ('inverse of x modulo x^4 + x + 1', bits (inverse_in (mp, 1, 2), 4), '1001');
  same ('x^4 + x^3 + x^2 + x + 1 irreducible', yes (is_irreducible2 (31)), 'yes');
  [ap, mp, op] = ring_mod2 (31);
  sn ('units modulo x^4 + x^3 + x^2 + x + 1', numel (unit_list (mp, 1)), 15);
  same ('(x + 1)^4', dec2bin (clmul (clmul (3, 3), clmul (3, 3))), '10001');
  [ar, mul_r, orr] = ring_mod2 (17);
  same ('units and nonzero nilpotents modulo x^4 + 1', ...
        jn ([numel(unit_list (mul_r, 1)), numel(nilpotent_list (mul_r))]), '8 7');
  a = 0:4;
  same ('roots of x^2 + 1 in GF(5)', jn (a(mod (a .* a + 1, 5) == 0)), '2 3');
  [a5, m5, o5] = ring_quad (5, 1);
  same ('(x + 2)(x + 3) modulo x^2 + 1 over GF(5)', qtxt (m5(2 + 5 + 1, 3 + 5 + 1), 5), '0');
  sn ('units of GF(5)[x]/(x^2 + 1)', numel (unit_list (m5, 1)), 16);
  [a12, m12, o12] = ring_zn (12);
  masks12 = ideal_masks (a12, m12);
  same ('(4) and (6) in Z_12 share', jn (members_of (bitand (ideal_zn (4, 12), ideal_zn (6, 12)))), '0');
  same ('(2) and (3) in Z_12 share', jn (members_of (bitand (ideal_zn (2, 12), ideal_zn (3, 12)))), '0 6');
  require ('both intersections are ideals', ...
           any (masks12 == bitand (ideal_zn (4, 12), ideal_zn (6, 12))) && ...
           any (masks12 == bitand (ideal_zn (2, 12), ideal_zn (3, 12))));
end


check_opening ();
check_ideals ();
check_principal ();
check_quotients ();
check_four ();
check_first_iso ();
check_maximal ();
check_root ();
check_hardware ();
check_worked ();
check_exercises ();
printf ('\nAll chapter 34 numbers reproduced.\n');
