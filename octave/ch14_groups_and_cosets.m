% Reproduce every number in chapter 14 of Error Correction Coding,
% "Groups and Cosets", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch14_groups_and_cosets.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch14_groups_and_cosets.py and prints the same
% lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in this chapter is a whole number, a bit string or
% a list of them, so every comparison is exact. The small groups are checked
% exhaustively: the group axioms on every triple of elements, the subgroup
% tests on every subset of Z_n for n up to 12 and of the 3-bit strings, the
% coset rules on every pair of elements, and the syndrome map of the code
% C_5 on all 32 words of five bits. A run that reaches the final line has
% reproduced every number in the chapter.
%
% A finite group is held as a struct with the row vector e of its elements
% and a function handle op that combines two arrays element by element, so
% a whole Cayley table comes from one call on the grids that ndgrid makes. A
% word of bits is held in a double, which stores these small integers
% exactly, and XOR is bitxor.

1;


% Helpers ----------------------------------------------------------------------

function s = as_text (x)
  % The text Python's str() prints for a str, a bool or an int.
  if (ischar (x))
    s = x;
  elseif (islogical (x))
    if (x)
      s = 'True';
    else
      s = 'False';
    end
  else
    s = sprintf ('%d', x);
  end
end

function s = as_repr (x)
  if (ischar (x))
    s = ['''' x ''''];
  else
    s = as_text (x);
  end
end

function same (label, computed, stated)
  c = as_text (computed);
  s = as_text (stated);
  printf ('  %-62s %-16s %s\n', label, c, s);
  if (~strcmp (c, s))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, as_repr (computed), as_repr (stated));
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

function s = show (x, n)
  s = dec2bin (x, n);
end

function w = weight (x)
  w = zeros (size (x));
  while (any (x(:)))
    w = w + bitand (x, 1);
    x = bitshift (x, -1);
  end
end

function parts = texts (xs, n)
  parts = cell (1, numel (xs));
  for i = 1:numel (xs)
    if (n)
      parts{i} = dec2bin (xs(i), n);
    else
      parts{i} = sprintf ('%d', xs(i));
    end
  end
end

function s = as_list (xs, n)
  if (nargin < 2)
    n = 0;
  end
  s = strjoin (texts (xs, n), ', ');
end

function s = as_set (xs, n)
  if (nargin < 2)
    n = 0;
  end
  s = ['{' as_list(xs, n) '}'];
end

function s = as_row (xs, n)
  if (nargin < 2)
    n = 0;
  end
  s = strjoin (texts (xs, n), ' ');
end


% Groups -----------------------------------------------------------------------

function g = z_add (n)
  g.e = 0:n - 1;
  g.op = @(a, b) mod (a + b, n);
end

function g = z_units (n)
  e = 1:n - 1;
  g.e = e(gcd (e, n) == 1);
  g.op = @(a, b) mod (a .* b, n);
end

function g = bit_strings (n)
  g.e = 0:2 ^ n - 1;
  g.op = @(a, b) bitxor (a, b);
end

function ids = identities (e, op)
  ids = [];
  for x = e
    if (all (op (x, e) == e) && all (op (e, x) == e))
      ids(end + 1) = x;
    end
  end
end

function x = identity (e, op)
  ids = identities (e, op);
  x = ids(1);
end

function r = is_closed (h, op)
  [A, B] = ndgrid (h, h);
  T = op (A, B);
  r = all (ismember (T(:), h));
end

function r = is_group (h, op)
  r = false;
  if (isempty (h) || ~is_closed (h, op))
    return;
  end
  [A, B, C] = ndgrid (h, h, h);
  L = op (op (A, B), C);
  R = op (A, op (B, C));
  if (any (L(:) ~= R(:)))
    return;
  end
  ids = identities (h, op);
  if (isempty (ids))
    return;
  end
  x = ids(1);
  for a = h
    if (~any (op (a, h) == x & op (h, a) == x))
      return;
    end
  end
  r = true;
end

function b = inverses_of (e, op, a)
  x = identity (e, op);
  b = e(op (a, e) == x & op (e, a) == x);
end

function b = inverse (e, op, a)
  all_b = inverses_of (e, op, a);
  b = all_b(1);
end

function inv = inverse_map (e, op)
  % inv(x + 1) is the inverse of x, for every element x.
  inv = zeros (1, max (e) + 1);
  for a = e
    inv(a + 1) = inverse (e, op, a);
  end
end

function rows = table_rows (e, op, n)
  rows = cell (1, numel (e));
  for i = 1:numel (e)
    rows{i} = as_row (op (e(i), e), n);
  end
end

function r = is_latin (e, op)
  target = sort (e);
  r = true;
  for a = e
    r = r && isequal (sort (op (a, e)), target) && isequal (sort (op (e, a)), target);
  end
end

function out = powers (e, op, a)
  % a, a^2, a^3, ... up to and including the first power equal to the identity.
  x = identity (e, op);
  out = a;
  while (out(end) ~= x)
    out(end + 1) = op (out(end), a);
  end
end

function k = order (e, op, a)
  k = numel (powers (e, op, a));
end

function subs = all_subsets (e)
  % Every nonempty subset, smallest first. nchoosek reads a one-element
  % vector as the number n of n-choose-k, so that case is answered directly.
  if (numel (e) == 1)
    subs = {e};
    return;
  end
  subs = cell (1, 2 ^ numel (e) - 1);
  k = 0;
  for r = 1:numel (e)
    c = nchoosek (e, r);
    for i = 1:rows (c)
      k = k + 1;
      subs{k} = c(i, :);
    end
  end
end

function found = subgroups (e, op)
  found = {};
  subs = all_subsets (e);
  for i = 1:numel (subs)
    h = subs{i};
    if (is_closed (h, op) && is_group (h, op))
      found{end + 1} = h;
    end
  end
end

function r = subgroup_test (op, h, inv)
  [A, B] = ndgrid (h, h);
  T = op (A, inv(B + 1));
  r = ~isempty (h) && all (ismember (T(:), h));
end

function c = coset (op, a, h)
  c = unique (op (a, h));
end

function out = cosets (e, op, h)
  out = {};
  for a = e
    c = coset (op, a, h);
    known = false;
    for j = 1:numel (out)
      known = known || isequal (out{j}, c);
    end
    if (~known)
      out{end + 1} = c;
    end
  end
end

function r = is_hom (dom, op_dom, op_cod, f)
  [A, B] = ndgrid (dom, dom);
  L = f (op_dom (A, B));
  R = op_cod (f (A), f (B));
  r = all (L(:) == R(:));
end

function k = kernel (dom, f, e_cod)
  k = dom(f (dom) == e_cod);
end

function m = image_of (dom, f)
  m = unique (f (dom));
end


% The code C_5 and its syndromes -------------------------------------------------

function c = c5 ()
  c = [bin2dec('00000'), bin2dec('01101'), bin2dec('10110'), bin2dec('11011')];
end

function b = bit (x, i, n)
  % Bit i of an n-bit word, counting from 1 at the left.
  if (nargin < 3)
    n = 5;
  end
  b = bitand (bitshift (x, -(n - i)), 1);
end

function s = syndrome (x)
  % The three parity checks of C_5: x1+x2+x3, x1+x4 and x2+x5, each mod 2.
  s1 = bitxor (bitxor (bit (x, 1), bit (x, 2)), bit (x, 3));
  s2 = bitxor (bit (x, 1), bit (x, 4));
  s3 = bitxor (bit (x, 2), bit (x, 5));
  s = 4 * s1 + 2 * s2 + s3;
end

function ld = leader (row)
  % Least weight, and among equals the smaller word read as a binary number.
  w = weight (row);
  ld = min (row(w == min (w)));
end

function out = standard_array (code, n)
  % Each coset as bitxor(leader, code), rows by leader weight, then leftmost bit first.
  seen = false (1, 2 ^ n);
  leaders = [];
  for x = 0:2 ^ n - 1
    if (seen(x + 1))
      continue;
    end
    row = unique (bitxor (x, code));
    seen(row + 1) = true;
    leaders(end + 1) = leader (row);
  end
  [~, idx] = sortrows ([weight(leaders(:)), -leaders(:)]);
  leaders = leaders(idx);
  out = cell (1, numel (leaders));
  for i = 1:numel (leaders)
    out{i} = bitxor (leaders(i), code);
  end
end

function p = primes_below (n)
  p = primes (n - 1);
end

function k = order_of_two (p)
  g = z_units (p);
  k = order (g.e, g.op, 2);
end

function r = powmod (a, k, n)
  r = 1;
  for i = 1:k
    r = mod (r * a, n);
  end
end

function i = row_holding (rows, x)
  i = 0;
  for j = 1:numel (rows)
    if (any (rows{j} == x))
      i = j;
      return;
    end
  end
end


% The test groups --------------------------------------------------------------

function groups = test_groups ()
  groups = {};
  for n = 1:12
    groups{end + 1} = z_add (n);
  end
  for n = 1:3
    groups{end + 1} = bit_strings (n);
  end
  for n = [5 7 8 9 12]
    groups{end + 1} = z_units (n);
  end
end


% The chapter's numbers ----------------------------------------------------------

function check_tables ()
  section ('Three tables that behave alike');
  z4 = z_add (4);
  b2 = bit_strings (2);
  u5 = z_units (5);
  rows = table_rows (z4.e, z4.op, 0);
  stated = {'0 1 2 3', '1 2 3 0', '2 3 0 1', '3 0 1 2'};
  for i = 1:4
    same (sprintf ('Z_4 under addition, row %s', rows{i}(1)), rows{i}, stated{i});
  end
  rows = table_rows (b2.e, b2.op, 2);
  stated = {'00 01 10 11', '01 00 11 10', '10 11 00 01', '11 10 01 00'};
  for i = 1:4
    same (sprintf ('two-bit strings under XOR, row %s', rows{i}(1:2)), rows{i}, stated{i});
  end
  rows = table_rows (u5.e, u5.op, 0);
  stated = {'1 2 3 4', '2 4 1 3', '3 1 4 2', '4 3 2 1'};
  for i = 1:4
    same (sprintf ('nonzero residues mod 5 under multiplication, row %s', rows{i}(1)), rows{i}, stated{i});
  end
  rename = mod (2 .^ z4.e, 5);
  same ('2^k mod 5 for k = 0, 1, 2, 3', as_list (rename), '1, 2, 4, 3');
  [A, B] = ndgrid (z4.e, z4.e);
  L = rename(z4.op (A, B) + 1);
  R = u5.op (rename(A + 1), rename(B + 1));
  require ('2^(a + b mod 4) = 2^a 2^b mod 5 for all a, b', all (L(:) == R(:)));
  same ('1 + 1 in Z_4', z4.op (1, 1), 2);
  same ('elements a of Z_4 with a + a = 0', as_list (z4.e(z4.op (z4.e, z4.e) == 0)), '0, 2');
  same ('two-bit strings x with x XOR x = 00', sum (b2.op (b2.e, b2.e) == 0), 4);
  require ('every row and column of all three tables holds each element once', ...
           is_latin (z4.e, z4.op) && is_latin (b2.e, b2.op) && is_latin (u5.e, u5.op));
end

function check_axioms ()
  section ('The group axioms');
  groups = test_groups ();
  ok_group = true; ok_id = true; ok_inv = true; ok_invinv = true; ok_cancel = true; ok_latin = true; ok_prod = true;
  for i = 1:numel (groups)
    e = groups{i}.e;
    op = groups{i}.op;
    ok_group = ok_group && is_group (e, op);
    ok_id = ok_id && numel (identities (e, op)) == 1;
    for a = e
      ok_inv = ok_inv && numel (inverses_of (e, op, a)) == 1;
      ok_invinv = ok_invinv && inverse (e, op, inverse (e, op, a)) == a;
      ok_cancel = ok_cancel && numel (unique (op (a, e))) == numel (e);
      for b = e
        ok_prod = ok_prod && inverse (e, op, op (a, b)) == op (inverse (e, op, b), inverse (e, op, a));
      end
    end
    ok_latin = ok_latin && is_latin (e, op);
  end
  require ('Z_n, n = 1..12, bit strings of 1..3 bits, Z_5*, Z_7*, Z_8*, Z_9*, Z_12* are groups', ok_group);
  require ('each of them has exactly one identity', ok_id);
  require ('each element of each of them has exactly one inverse', ok_inv);
  require ('the inverse of the inverse of a is a', ok_invinv);
  require ('cancellation: a*b = a*c forces b = c', ok_cancel);
  require ('every Cayley table is a Latin square', ok_latin);
  require ('(a*b)^-1 = b^-1 * a^-1 in every test group', ok_prod);
  same ('2 x 3 mod 6', mod (2 * 3, 6), 0);
  require ('the nonzero residues mod 6 are not closed under multiplication', ~is_closed (1:5, @(a, b) mod (a .* b, 6)));
  same ('(5 - 3) - 1', (5 - 3) - 1, 1);
  same ('5 - (3 - 1)', 5 - (3 - 1), 3);
  u5 = z_units (5);
  same ('inverses of 1, 2, 3, 4 mod 5', as_list (arrayfun (@(a) inverse (u5.e, u5.op, a), u5.e)), '1, 3, 2, 4');
  b3 = bit_strings (3);
  require ('every 3-bit string is its own inverse under XOR', all (b3.op (b3.e, b3.e) == 0));
end

function check_subgroups ()
  section ('Subgroups');
  z12 = z_add (12);
  add12 = z12.op;
  require ('the even residues {0, 2, ..., 10} form a subgroup of Z_12', is_group ([0 2 4 6 8 10], add12));
  same ('-4 in Z_12', inverse (z12.e, add12, 4), 8);
  require ('{0, 4, 8} is a subgroup of Z_12', is_group ([0 4 8], add12));
  require ('{1, 5, 9} is not a subgroup of Z_12', ~is_group ([1 5 9], add12));
  same ('1 + 5 in Z_12', add12 (1, 5), 6);
  require ('{0, 1, 2, 3, 4, 5} is not a subgroup of Z_12', ~is_group ([0 1 2 3 4 5], add12));
  same ('3 + 4 in Z_12', add12 (3, 4), 7);
  ok_test = true;
  ok_finite = true;
  cases = {};
  for n = 1:12
    cases{end + 1} = z_add (n);
  end
  cases{end + 1} = bit_strings (3);
  for i = 1:numel (cases)
    e = cases{i}.e;
    op = cases{i}.op;
    inv = inverse_map (e, op);
    subs = all_subsets (e);
    for j = 1:numel (subs)
      h = subs{j};
      truth = is_group (h, op);
      ok_test = ok_test && (subgroup_test (op, h, inv) == truth);
      ok_finite = ok_finite && (is_closed (h, op) == truth);
    end
  end
  require ('subgroup test agrees with the axioms on every subset tried', ok_test);
  require ('closure alone decides every finite subset tried', ok_finite);
  same ('01101 XOR 10110', show (bitxor (bin2dec ('01101'), bin2dec ('10110')), 5), '11011');
  same ('01101 XOR 11011', show (bitxor (bin2dec ('01101'), bin2dec ('11011')), 5), '10110');
  same ('10110 XOR 11011', show (bitxor (bin2dec ('10110'), bin2dec ('11011')), 5), '01101');
  require ('C_5 is a subgroup of the 5-bit strings', is_group (c5 (), @(a, b) bitxor (a, b)));
  require ('the even-weight words {000, 011, 101, 110} form a subgroup', is_group ([0 3 5 6], @(a, b) bitxor (a, b)));
end

function check_cyclic ()
  section ('Cyclic groups and the order of an element');
  u7 = z_units (7);
  same ('powers of 3 mod 7', as_list (powers (u7.e, u7.op, 3)), '3, 2, 6, 4, 5, 1');
  same ('powers of 2 mod 7', as_list (powers (u7.e, u7.op, 2)), '2, 4, 1');
  same ('order of 6 mod 7', order (u7.e, u7.op, 6), 2);
  z12 = z_add (12);
  same ('multiples of 9 in Z_12', as_list (powers (z12.e, z12.op, 9)), '9, 6, 3, 0');
  orders = arrayfun (@(a) order (z12.e, z12.op, a), z12.e);
  same ('orders of 0, 1, ..., 11 in Z_12', as_list (orders), '1, 12, 6, 4, 3, 12, 2, 12, 3, 4, 6, 12');
  same ('generators of Z_12', as_list (z12.e(orders == 12)), '1, 5, 7, 11');
  ok = true;
  for n = 1:40
    g = z_add (n);
    for a = 0:n - 1
      ok = ok && order (g.e, g.op, a) == n / gcd (a, n);
    end
  end
  require ('the order of a in Z_n is n / gcd(a, n), n = 1..40', ok);
  groups = test_groups ();
  ok_size = true;
  ok_div = true;
  for i = 1:numel (groups)
    e = groups{i}.e;
    op = groups{i}.op;
    x0 = identity (e, op);
    for a = e
      p = powers (e, op, a);
      m = numel (p);
      ok_size = ok_size && numel (unique (p)) == m;
      x = x0;
      for k = 0:36
        ok_div = ok_div && ((x == x0) == (mod (k, m) == 0));
        x = op (x, a);
      end
    end
  end
  require ('<a> has exactly ord(a) distinct elements in every test group', ok_size);
  require ('a^k = e exactly when ord(a) divides k, k = 0..36', ok_div);
  ok = true;
  for n = 1:12
    g = z_add (n);
    found = subgroups (g.e, g.op);
    got = sort (cellfun (@(h) mat2str (h), found, 'UniformOutput', false));
    expected = {};
    for d = 1:n
      if (mod (n, d) == 0)
        expected{end + 1} = mat2str (unique (powers (g.e, g.op, mod (n / d, n))));
      end
    end
    ok = ok && isequal (got, sort (expected));
  end
  require ('the subgroups of Z_n are the <n/d> for d dividing n, n = 1..12', ok);
  found = subgroups (z12.e, z12.op);
  same ('number of subgroups of Z_12', numel (found), 6);
  [~, idx] = sort (cellfun (@numel, found));
  stated = {'{0}', '{0, 6}', '{0, 4, 8}', '{0, 3, 6, 9}', '', '{0, 2, 4, 6, 8, 10}'};
  for i = idx
    h = found{i};
    if (numel (h) <= 6)
      same (sprintf ('the subgroup of Z_12 of order %d', numel (h)), as_set (h), stated{numel (h)});
    else
      same ('the subgroup of Z_12 of order 12', numel (h), 12);
    end
  end
end

function check_cosets ()
  section ('Cosets');
  z12 = z_add (12);
  add12 = z12.op;
  found = cosets (z12.e, add12, [0 4 8]);
  stated = {'{0, 4, 8}', '{1, 5, 9}', '{2, 6, 10}', '{3, 7, 11}'};
  for i = 1:numel (found)
    same (sprintf ('coset %d + {0, 4, 8} in Z_12', found{i}(1)), as_set (found{i}), stated{i});
  end
  same ('number of cosets of {0, 4, 8}', numel (found), 4);
  same ('4 + {0, 4, 8}', as_set (coset (add12, 4, [0 4 8])), '{0, 4, 8}');
  found = cosets (z12.e, add12, [0 3 6 9]);
  stated = {'{0, 3, 6, 9}', '{1, 4, 7, 10}', '{2, 5, 8, 11}'};
  for i = 1:numel (found)
    same (sprintf ('coset %d + {0, 3, 6, 9} in Z_12', found{i}(1)), as_set (found{i}), stated{i});
  end
  b3 = bit_strings (3);
  r3 = [0 7];
  reps = [0 1 2 4];
  stated = {'{000, 111}', '{001, 110}', '{010, 101}', '{011, 100}'};
  for i = 1:4
    same (sprintf ('%s XOR R_3', show (reps(i), 3)), as_set (coset (b3.op, reps(i), r3), 3), stated{i});
  end
  same ('number of cosets of R_3', numel (cosets (b3.e, b3.op, r3)), 4);
  u7 = z_units (7);
  families = {z12, b3, u7};
  ok_crit = true; ok_part = true; ok_size = true; ok_own = true;
  for f = 1:numel (families)
    e = families{f}.e;
    op = families{f}.op;
    hs = subgroups (e, op);
    for j = 1:numel (hs)
      h = hs{j};
      for a = e
        ca = coset (op, a, h);
        ok_size = ok_size && numel (ca) == numel (h);
        ok_own = ok_own && any (ca == a);
        ia = inverse (e, op, a);
        for b = e
          cb = coset (op, b, h);
          ok_crit = ok_crit && (isequal (ca, cb) == any (h == op (ia, b)));
          ok_part = ok_part && (isequal (ca, cb) || isempty (intersect (ca, cb)));
        end
      end
    end
  end
  require ('aH = bH exactly when a^-1 b lies in H, every case', ok_crit);
  require ('two cosets are equal or share no element, every case', ok_part);
  require ('every coset has as many elements as H, every case', ok_size);
  require ('a lies in its own coset aH, every case', ok_own);
end

function check_lagrange ()
  section ('Lagrange''s theorem');
  z12 = z_add (12);
  found = subgroups (z12.e, z12.op);
  [~, idx] = sort (cellfun (@numel, found));
  found = found(idx);
  same ('orders of the subgroups of Z_12', as_list (cellfun (@numel, found)), '1, 2, 3, 4, 6, 12');
  counts = cellfun (@(h) numel (cosets (z12.e, z12.op, h)), found);
  same ('number of cosets of each, 12 / |H|', as_list (counts), '12, 6, 4, 3, 2, 1');
  u7 = z_units (7);
  same ('orders of 1, 2, 3, 4, 5, 6 mod 7', as_list (arrayfun (@(a) order (u7.e, u7.op, a), u7.e)), '1, 3, 6, 3, 6, 2');
  same ('powers of 5 mod 7', as_list (powers (u7.e, u7.op, 5)), '5, 4, 6, 2, 3, 1');
  same ('powers of 1 mod 7', as_list (powers (u7.e, u7.op, 1)), '1');
  same ('powers of 4 mod 7', as_list (powers (u7.e, u7.op, 4)), '4, 2, 1');
  same ('powers of 6 mod 7', as_list (powers (u7.e, u7.op, 6)), '6, 1');
  same ('coset 3{1, 2, 4} mod 7', as_set (coset (u7.op, 3, powers (u7.e, u7.op, 2))), '{3, 5, 6}');
  same ('elements of order 4 mod 7', sum (arrayfun (@(a) order (u7.e, u7.op, a), u7.e) == 4), 0);
  groups = test_groups ();
  ok_div = true;
  ok_pow = true;
  for i = 1:numel (groups)
    e = groups{i}.e;
    op = groups{i}.op;
    x0 = identity (e, op);
    for a = e
      ok_div = ok_div && mod (numel (e), order (e, op, a)) == 0;
      x = x0;
      for k = 1:numel (e)
        x = op (x, a);
      end
      ok_pow = ok_pow && x == x0;
    end
  end
  require ('every element order divides the group order, every test group', ok_div);
  require ('a^|G| = e in every test group', ok_pow);
  ok = true;
  for p = primes_below (60)
    for a = 1:p - 1
      ok = ok && powmod (a, p - 1, p) == 1;
    end
  end
  require ('a^(p - 1) = 1 mod p for every prime p < 60 and 0 < a < p', ok);
  ok = true;
  for n = 2:40
    u = z_units (n);
    for a = u.e
      ok = ok && powmod (a, numel (u.e), n) == 1;
    end
  end
  require ('a^phi(n) = 1 mod n for every n from 2 to 40 and a coprime to n', ok);
  ok = true;
  for p = primes_below (14)
    g = z_add (p);
    hs = subgroups (g.e, g.op);
    ok = ok && isequal (cellfun (@numel, hs), [1 p]);
    for a = 1:p - 1
      ok = ok && order (g.e, g.op, a) == p;
    end
  end
  require ('Z_p has only {0} and Z_p, and every a != 0 generates it, p < 14', ok);
  b3 = bit_strings (3);
  same ('sizes of the subgroups of the 3-bit strings', as_list (unique (cellfun (@numel, subgroups (b3.e, b3.op)))), '1, 2, 4, 8');
  same ('64 mod 12, so no 12-word linear code of length 6', mod (64, 12), 4);
  same ('32 mod 6, so no 6-word linear code of length 5', mod (32, 6), 2);
  same ('12 mod 5, so no subgroup of order 5 when |G| = 12', mod (12, 5), 2);
  same ('cosets of a 16-word linear code of length 7', 2 ^ 7 / 16, 8);
  odd = primes_below (100);
  odd = odd(odd > 2);
  orders = arrayfun (@order_of_two, odd);
  same ('odd primes below 100', numel (odd), 24);
  require ('the order of 2 mod p divides p - 1 for each of them', all (mod (odd - 1, orders) == 0));
  full = odd(orders == odd - 1);
  half = odd(2 * orders == odd - 1);
  rest = ~(orders == odd - 1 | 2 * orders == odd - 1);
  same ('primes where 2 has order p - 1', as_list (full), '3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83');
  same ('how many', numel (full), 12);
  same ('primes where 2 has order (p - 1) / 2', as_list (half), '7, 17, 23, 41, 47, 71, 79, 97');
  same ('how many', numel (half), 8);
  pairs = arrayfun (@(p, k) sprintf ('%d: %d', p, k), odd(rest), orders(rest), 'UniformOutput', false);
  same ('the other primes and the order of 2', strjoin (pairs, ', '), '31: 5, 43: 14, 73: 9, 89: 11');
end

function check_homomorphisms ()
  section ('Homomorphisms and kernels');
  z12 = z_add (12);
  z4 = z_add (4);
  mod4 = @(x) mod (x, 4);
  require ('x -> x mod 4 from Z_12 to Z_4 is a homomorphism', is_hom (z12.e, z12.op, z4.op, mod4));
  same ('its kernel', as_set (kernel (z12.e, mod4, 0)), '{0, 4, 8}');
  fibers = arrayfun (@(y) as_set (z12.e(mod4 (z12.e) == y)), z4.e, 'UniformOutput', false);
  same ('the words sent to 1', fibers{2}, '{1, 5, 9}');
  cs = cosets (z12.e, z12.op, [0 4 8]);
  ok = true;
  for i = 1:4
    ok = ok && strcmp (fibers{i}, as_set (cs{i}));
  end
  require ('the words sent to each value are a coset of the kernel', ok);
  b3 = bit_strings (3);
  parity = @(x) mod (weight (x), 2);
  require ('parity from 3-bit strings to Z_2 is a homomorphism', is_hom (b3.e, b3.op, @(a, b) mod (a + b, 2), parity));
  same ('kernel of parity on 3-bit strings', as_set (kernel (b3.e, parity, 0), 3), '{000, 011, 101, 110}');
  same ('the 3-bit strings of odd parity', as_set (b3.e(parity (b3.e) == 1), 3), '{001, 010, 100, 111}');
  ok_hom = true; ok_e = true; ok_sub = true; ok_fib = true; ok_inj = true;
  for k = 0:11
    f = @(x) mod (k * x, 12);
    ok_hom = ok_hom && is_hom (z12.e, z12.op, z12.op, f);
    ok_e = ok_e && f (0) == 0;
    for a = z12.e
      ok_e = ok_e && f (inverse (z12.e, z12.op, a)) == inverse (z12.e, z12.op, f (a));
    end
    ker = kernel (z12.e, f, 0);
    ok_sub = ok_sub && is_group (ker, z12.op) && is_group (image_of (z12.e, f), z12.op);
    for a = z12.e
      for b = z12.e
        ok_fib = ok_fib && ((f (a) == f (b)) == any (ker == z12.op (b, inverse (z12.e, z12.op, a))));
      end
    end
    ok_inj = ok_inj && ((numel (image_of (z12.e, f)) == 12) == isequal (ker, 0));
  end
  require ('every x -> kx mod 12 is a homomorphism of Z_12', ok_hom);
  require ('f(0) = 0 and f(-a) = -f(a) for each of them', ok_e);
  require ('kernel and image are subgroups for each of them', ok_sub);
  require ('f(a) = f(b) exactly when b - a is in the kernel, each map', ok_fib);
  require ('one-to-one exactly when the kernel is {0}, each map', ok_inj);
  u5 = z_units (5);
  two_to = @(k) mod (2 .^ k, 5);
  require ('k -> 2^k mod 5 is a homomorphism from Z_4', is_hom (z4.e, z4.op, u5.op, two_to));
  same ('2^4 mod 5', mod (2 ^ 4, 5), 1);
  require ('k -> 2^k mod 5 is one-to-one and onto', isequal (sort (two_to (z4.e)), u5.e));
  b2 = bit_strings (2);
  P = perms (b2.e);
  homs = 0;
  for i = 1:rows (P)
    p = P(i, :);
    if (is_hom (z4.e, z4.op, b2.op, @(k) p(k + 1)))
      homs = homs + 1;
    end
  end
  same ('bijections from Z_4 to 2-bit strings that are homomorphisms', homs, 0);
  same ('bijections tried', rows (P), 24);
end

function check_quotients ()
  section ('Quotient groups');
  z12 = z_add (12);
  add12 = z12.op;
  h = [0 4 8];
  cs = cosets (z12.e, add12, h);
  reps = cellfun (@(c) c(1), cs);
  name = @(x) min (coset (add12, x, h));
  stated = {'0 1 2 3', '1 2 3 0', '2 3 0 1', '3 0 1 2'};
  for i = 1:numel (reps)
    row = as_row (arrayfun (@(b) name (add12 (reps(i), b)), reps));
    same (sprintf ('cosets of {0, 4, 8} added, row %s', row(1)), row, stated{i});
  end
  same ('1 + 2 in Z_12', add12 (1, 2), 3);
  same ('5 + 10 in Z_12', add12 (5, 10), 3);
  same ('{1, 5, 9} + {2, 6, 10}', as_set (coset (add12, 3, h)), '{3, 7, 11}');
  b3 = bit_strings (3);
  families = {z12, b3};
  ok_well = true;
  ok_group = true;
  for f = 1:numel (families)
    e = families{f}.e;
    op = families{f}.op;
    ks = subgroups (e, op);
    for j = 1:numel (ks)
      k = ks{j};
      C = cell (1, max (e) + 1);
      for x = e
        C{x + 1} = coset (op, x, k);
      end
      for a = e
        for a2 = C{a + 1}
          for b = e
            for b2 = C{b + 1}
              ok_well = ok_well && isequal (C{op(a, b) + 1}, C{op(a2, b2) + 1});
            end
          end
        end
      end
      names = cellfun (@(c) c(1), C);
      qs = cosets (e, op, k);
      ok_group = ok_group && is_group (cellfun (@(c) c(1), qs), @(a, b) names(op (a, b) + 1));
    end
  end
  require ('the coset sum does not depend on the representatives', ok_well);
  require ('the cosets form a group under that sum, every case', ok_group);
  ok = true;
  for k = 0:11
    f = @(x) mod (k * x, 12);
    ok = ok && 12 / numel (kernel (z12.e, f, 0)) == numel (image_of (z12.e, f));
  end
  require ('|G/ker f| = |image f| for every x -> kx mod 12', ok);
  parity = @(x) mod (weight (x), 2);
  same ('cosets of the even-weight words in 3-bit strings', numel (cosets (b3.e, b3.op, kernel (b3.e, parity, 0))), 2);
end

function check_syndromes ()
  section ('Cosets and syndromes');
  C5 = c5 ();
  b5 = bit_strings (5);
  require ('all three parity checks are 0 on every codeword of C_5', all (syndrome (C5) == 0));
  require ('the syndrome of x XOR y is the XOR of their syndromes', is_hom (b5.e, b5.op, b5.op, @syndrome));
  same ('kernel of the syndrome map', as_set (kernel (b5.e, @syndrome, 0), 5), '{00000, 01101, 10110, 11011}');
  same ('number of syndromes reached', numel (image_of (b5.e, @syndrome)), 8);
  same ('number of cosets of C_5, 32 / 4', numel (cosets (b5.e, b5.op, C5)), 8);
  [X, Y] = ndgrid (b5.e, b5.e);
  L = syndrome (X) == syndrome (Y);
  R = ismember (bitxor (X, Y), C5);
  require ('same syndrome exactly when x XOR y lies in C_5', all (L(:) == R(:)));
  rows = standard_array (C5, 5);
  stated = {'00000 01101 10110 11011', '000'; '10000 11101 00110 01011', '110'; ...
            '01000 00101 11110 10011', '101'; '00100 01001 10010 11111', '100'; ...
            '00010 01111 10100 11001', '010'; '00001 01100 10111 11010', '001'; ...
            '01010 00111 11100 10001', '111'; '00011 01110 10101 11000', '011'};
  for i = 1:numel (rows)
    row = rows{i};
    same (sprintf ('coset with leader %s', show (row(1), 5)), as_row (row, 5), stated{i, 1});
    same ('its syndrome', show (syndrome (row(1)), 3), stated{i, 2});
    require ('every word of that coset has that syndrome', all (syndrome (row) == syndrome (row(1))));
  end
  same ('leader weights', as_list (cellfun (@(r) weight (r(1)), rows)), '0, 1, 1, 1, 1, 1, 2, 2');
  same ('8 cosets of 4 words', sum (cellfun (@numel, rows)), 32);
  same ('syndrome of 11111', show (syndrome (bin2dec ('11111')), 3), '100');
  same ('11111 XOR its leader 00100', show (bitxor (bin2dec ('11111'), bin2dec ('00100')), 5), '11011');
  ties = cell (1, 2);
  for i = 7:8
    r = rows{i};
    ties{i - 6} = as_list (sort (r(weight (r) == 2), 'descend'), 5);
  end
  same ('least-weight words in the coset of syndrome 111', ties{1}, '10001, 01010');
  same ('least-weight words in the coset of syndrome 011', ties{2}, '11000, 00011');
  doubles = b5.e(weight (b5.e) == 2);
  in_low = 0;
  for x = doubles
    if (weight (leader (rows{row_holding (rows, x)})) == 1)
      in_low = in_low + 1;
    end
  end
  same ('weight-2 words that lie in a coset with a weight-1 leader', in_low, 6);
  same ('weight-2 words in all', numel (doubles), 10);
  ok = true;
  for x = b5.e
    ok = ok && weight (leader (rows{row_holding (rows, x)})) == min (weight (bitxor (x, C5)));
  end
  require ('the leader of r''s coset is r XOR a nearest codeword, all r', ok);
end

function check_examples ()
  section ('Worked examples');
  u8 = z_units (8);
  rows = table_rows (u8.e, u8.op, 0);
  stated = {'1 3 5 7', '3 1 7 5', '5 7 1 3', '7 5 3 1'};
  for i = 1:4
    same (sprintf ('residues 1, 3, 5, 7 mod 8 multiplied, row %s', rows{i}(1)), rows{i}, stated{i});
  end
  require ('Z_8* is a group', is_group (u8.e, u8.op));
  same ('squares of 1, 3, 5, 7 mod 8', as_list (u8.op (u8.e, u8.e)), '1, 1, 1, 1');
  b2 = bit_strings (2);
  rename = zeros (1, 8);
  rename(u8.e + 1) = b2.e;
  require ('1, 3, 5, 7 -> 00, 01, 10, 11 is an isomorphism onto 2-bit strings', is_hom (u8.e, u8.op, b2.op, @(a) rename(a + 1)));
  u9 = z_units (9);
  same ('residues mod 9 that share no factor with 9', as_list (u9.e), '1, 2, 4, 5, 7, 8');
  same ('powers of 2 mod 9', as_list (powers (u9.e, u9.op, 2)), '2, 4, 8, 7, 5, 1');
  same ('powers of 4 mod 9', as_list (powers (u9.e, u9.op, 4)), '4, 7, 1');
  same ('powers of 8 mod 9', as_list (powers (u9.e, u9.op, 8)), '8, 1');
  same ('orders of 1, 2, 4, 5, 7, 8 mod 9', as_list (arrayfun (@(a) order (u9.e, u9.op, a), u9.e)), '1, 6, 3, 6, 3, 2');
  hs = subgroups (u9.e, u9.op);
  [~, idx] = sort (cellfun (@numel, hs));
  same ('subgroups of Z_9*', strjoin (cellfun (@as_set, hs(idx), 'UniformOutput', false), ', '), '{1}, {1, 8}, {1, 4, 7}, {1, 2, 4, 5, 7, 8}');
  same ('coset 2{1, 4, 7} mod 9', as_set (coset (u9.op, 2, [1, 4, 7])), '{2, 5, 8}');
  r = bin2dec ('10011');
  same ('syndrome of 10011', show (syndrome (r), 3), '101');
  sa = standard_array (c5 (), 5);
  for i = 1:numel (sa)
    if (syndrome (sa{i}(1)) == syndrome (r))
      ld = sa{i}(1);
      break;
    end
  end
  same ('its coset leader', show (ld, 5), '01000');
  same ('10011 XOR 01000', show (bitxor (r, ld), 5), '11011');
  same ('distances from 10011 to the four codewords', as_list (weight (bitxor (r, c5 ()))), '3, 4, 2, 1');
end

function check_exercises ()
  section ('Exercises');
  u5 = z_units (5);
  same ('order of 2 mod 5', order (u5.e, u5.op, 2), 4);
  same ('order of 4 mod 5', order (u5.e, u5.op, 4), 2);
  z12 = z_add (12);
  hs = {[0 6], [0 5], [0 3 6 9], [0 2 4 8]};
  stated = {true, false, true, false};
  for i = 1:4
    same (sprintf ('%s is a subgroup of Z_12', as_set (hs{i})), is_group (hs{i}, z12.op), stated{i});
  end
  same ('5 + 5 in Z_12', z12.op (5, 5), 10);
  same ('2 + 4 in Z_12', z12.op (2, 4), 6);
  same ('order of 8 in Z_12', order (z12.e, z12.op, 8), 3);
  same ('order of 10 in Z_12', order (z12.e, z12.op, 10), 6);
  z20 = z_add (20);
  same ('order of 6 in Z_20', order (z20.e, z20.op, 6), 10);
  b5 = bit_strings (5);
  r5 = [0 bin2dec('11111')];
  found = cosets (b5.e, b5.op, r5);
  same ('number of cosets of R_5', numel (found), 16);
  lw = cellfun (@(c) weight (leader (c)), found);
  same ('leaders of weight 0, 1, 2', as_list ([sum(lw == 0), sum(lw == 1), sum(lw == 2)]), '1, 5, 10');
  require ('each coset of R_5 holds exactly one word of weight at most 2', all (cellfun (@(c) sum (weight (c) <= 2) == 1, found)));
  triple = @(x) mod (3 * x, 12);
  same ('kernel of x -> 3x mod 12', as_set (kernel (z12.e, triple, 0)), '{0, 4, 8}');
  same ('image of x -> 3x mod 12', as_set (image_of (z12.e, triple)), '{0, 3, 6, 9}');
  r = bin2dec ('10101');
  same ('syndrome of 10101', show (syndrome (r), 3), '011');
  sa = standard_array (c5 (), 5);
  row = sa{row_holding (sa, r)};
  same ('the coset of 10101', as_set (sort (row), 5), '{00011, 01110, 10101, 11000}');
  same ('distances from 10101 to the four codewords', as_list (weight (bitxor (r, c5 ()))), '3, 2, 2, 3');
end


check_tables ();
check_axioms ();
check_subgroups ();
check_cyclic ();
check_cosets ();
check_lagrange ();
check_homomorphisms ();
check_quotients ();
check_syndromes ();
check_examples ();
check_exercises ();
printf ('\nAll chapter 14 numbers reproduced.\n');
