% Reproduce every number in chapter 15 of Error Correction Coding,
% "Vector Spaces over Finite Fields", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch15_vector_spaces_over_finite_fields.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch15_vector_spaces_over_finite_fields.py and
% prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. The small spaces are checked exhaustively: the vector space
% axioms on every triple of vectors, the subspace test on every subset of
% GF(2)^3, every subspace of GF(2)^4, GF(3)^3 and GF(5)^2 built and counted
% by dimension, and every dual found by testing each vector of the space
% against the code. The only numbers in the chapter that are not whole are
% the fractions of lists of vectors that form a basis, which are printed to
% three decimal places. A run that reaches the final line has reproduced
% every number in the chapter.
%
% A vector of GF(q)^n is held as its index, the number whose base-q digits
% are its entries with the first entry most significant, so sorting indices
% sorts vectors in the same dictionary order Python uses for tuples. A space
% is a struct holding its digit matrix and three tables built from the digits
% by arithmetic mod q: add for the sum of two vectors, scl for a scalar
% times a vector, and dot for the dot product. A set of vectors is a sorted
% row of indices.
%
% Doubles hold whole numbers exactly only up to 2^53, and the products that
% count ordered bases pass that for q = 3 and 5 at n = 8. So the number of
% subspaces of each dimension is computed here by the q-Pascal rule, which
% keeps every intermediate value below the answer, where Python divides one
% exact product by another. Both give the same whole numbers.

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

function s = as_ints (xs)
  s = strjoin (arrayfun (@(x) sprintf ('%d', x), xs, 'UniformOutput', false), ' ');
end

function s = as_tuple (d)
  % A row of integers written the way show prints a vector over GF(q), such as (3,6).
  s = ['(' strjoin(arrayfun (@(a) sprintf ('%d', a), d, 'UniformOutput', false), ',') ')'];
end


% Vectors over GF(q) -----------------------------------------------------------

function sp = make_space (q, n)
  sp.q = q;
  sp.n = n;
  sp.N = q ^ n;
  sp.w = q .^ (n - 1:-1:0);
  r = (0:sp.N - 1)';
  V = zeros (sp.N, n);
  for t = n:-1:1
    V(:, t) = mod (r, q);
    r = floor (r / q);
  end
  sp.V = V;
  [I, J] = ndgrid (1:sp.N, 1:sp.N);
  sp.add = reshape (mod (V(I(:), :) + V(J(:), :), q) * sp.w', sp.N, sp.N);
  sp.dot = mod (V * V', q);
  sp.scl = zeros (q, sp.N);
  for c = 0:q - 1
    sp.scl(c + 1, :) = (mod (c * V, q) * sp.w')';
  end
end

function x = ix (sp, digits)
  % The index of the vector written as a string of digits, such as '0110'.
  x = (digits - '0') * sp.w';
end

function s = show (sp, x)
  d = sp.V(x + 1, :);
  if (sp.q == 2)
    s = char (d + '0');
  else
    s = ['(' strjoin(arrayfun (@(a) sprintf ('%d', a), d, 'UniformOutput', false), ',') ')'];
  end
end

function s = as_set (sp, xs)
  xs = sort (xs);
  parts = arrayfun (@(x) show (sp, x), xs, 'UniformOutput', false);
  s = ['{' strjoin(parts, ', ') '}'];
end

function s = as_row (sp, xs)
  parts = arrayfun (@(x) show (sp, x), xs, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function w = weight (sp, x)
  w = sum (sp.V(x + 1, :) ~= 0, 2)';
end

function z = vadd (sp, u, v)
  z = sp.add(u + 1, v + 1);
end

function z = vscale (sp, c, v)
  z = sp.scl(c + 1, v + 1);
end

function S = span (sp, vs)
  % Every linear combination of the vectors vs, sorted.
  S = 0;
  for v = vs
    T = sp.add(S + 1, sp.scl(:, v + 1) + 1);
    S = unique (T(:))';
  end
end

function C = combos (sp, vs)
  % All q^k combinations of the k vectors vs, one per choice of coefficients.
  C = 0;
  for v = vs
    T = sp.add(C(:) + 1, sp.scl(:, v + 1) + 1);
    C = T(:)';
  end
end

function ok = independent (sp, vs)
  % True when only the all-zero choice of coefficients gives the zero vector.
  ok = sum (combos (sp, vs) == 0) == 1;
end

function L = all_lists (N, k)
  % Every ordered list of k indices from 0 to N - 1, one list per row.
  r = (0:N ^ k - 1)';
  L = zeros (N ^ k, k);
  for t = 1:k
    L(:, t) = mod (floor (r / N ^ (k - t)), N);
  end
end

function C = combos_rows (sp, L)
  % combos for every row of L at once.
  M = rows (L);
  C = zeros (M, 1);
  for j = 1:columns (L)
    m = columns (C);
    D = zeros (M, m * sp.q);
    for c = 0:sp.q - 1
      s = sp.scl(c + 1, L(:, j) + 1)';
      D(:, c * m + (1:m)) = sp.add(sub2ind ([sp.N, sp.N], C + 1, repmat (s + 1, 1, m)));
    end
    C = D;
  end
end

function n = lists_independent (sp, k)
  % Ordered lists of k vectors that are linearly independent, by trying every list.
  C = combos_rows (sp, all_lists (sp.N, k));
  n = sum (sum (C == 0, 2) == 1);
end

function ok = is_subspace (sp, S)
  ok = false;
  if (isempty (S))
    return;
  end
  if (~any (S == 0))
    return;
  end
  sums = sp.add(S + 1, S + 1);
  if (~all (ismember (sums(:), S)))
    return;
  end
  sc = sp.scl(:, S + 1);
  ok = all (ismember (sc(:), S));
end

function k = dimension (S, q)
  k = 0;
  sz = 1;
  while (sz < numel (S))
    k = k + 1;
    sz = sz * q;
  end
  if (sz ~= numel (S))
    k = -1;
  end
end

function D = dual (sp, S)
  % Every vector whose dot product with each member of S is 0.
  D = find (all (sp.dot(:, S + 1) == 0, 2))' - 1;
end

function levels = subspaces (sp)
  % All subspaces grouped by dimension: levels{k + 1} holds one sorted subspace per row.
  levels = cell (1, sp.n + 1);
  levels{1} = 0;
  for d = 1:sp.n
    prev = levels{d};
    found = zeros (0, sp.q ^ d);
    for r = 1:rows (prev)
      S = prev(r, :);
      for v = 0:sp.N - 1
        if (~any (S == v))
          T = sp.add(S + 1, sp.scl(:, v + 1) + 1);
          found(end + 1, :) = unique (T(:))';
        end
      end
    end
    levels{d + 1} = unique (found, 'rows');
  end
end


% Counting formulas ------------------------------------------------------------

function out = ordered_independent (q, n, k)
  % (q^n - 1)(q^n - q)...(q^n - q^(k-1)).
  out = prod (q ^ n - q .^ (0:k - 1));
end

function g = gaussian (n, k, q)
  % The number of k-dimensional subspaces of GF(q)^n, by the q-Pascal rule.
  if (k == 0 || k == n)
    g = 1;
  else
    g = gaussian (n - 1, k - 1, q) + q ^ k * gaussian (n - 1, k, q);
  end
end

function p = basis_fraction (q, n)
  p = 1.0;
  for i = 0:n - 1
    p = p * (1 - 1 / q ^ (n - i));
  end
end

function s = fraction_text (num, den)
  g = gcd (num, den);
  s = sprintf ('%d/%d', num / g, den / g);
end


% The codes of the chapter -----------------------------------------------------

function S = c5 (sp)
  S = sort ([ix(sp, '00000'), ix(sp, '01101'), ix(sp, '10110'), ix(sp, '11011')]);
end

function H = checks5 (sp)
  H = [ix(sp, '11100'), ix(sp, '10010'), ix(sp, '01001')];
end

function s = syndrome (x)
  % The syndrome of Groups and Cosets: (x1 + x2 + x3, x1 + x4, x2 + x5) over GF(2).
  s = [mod(x(1) + x(2) + x(3), 2), mod(x(1) + x(4), 2), mod(x(2) + x(5), 2)];
end

function E = even_words (sp)
  E = find (mod (sum (sp.V, 2), 2) == 0)' - 1;
end

function s = weight_counts (sp, S)
  w = weight (sp, S);
  u = unique (w);
  parts = arrayfun (@(a) sprintf ('%d:%d', a, sum (w == a)), u, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function ok = axioms_hold (sp)
  N = sp.N;
  A = sp.add;
  Sc = sp.scl;
  e = 0:N - 1;
  ok = all (A(:, 1)' == e) && all (Sc(2, :) == e);
  ok = ok && all (any (A == 0, 2));
  ok = ok && isequal (A, A');
  [U, V, W] = ndgrid (1:N, 1:N, 1:N);
  lhs = A(sub2ind ([N, N], A(sub2ind ([N, N], U(:), V(:))) + 1, W(:)));
  rhs = A(sub2ind ([N, N], U(:), A(sub2ind ([N, N], V(:), W(:))) + 1));
  ok = ok && isequal (lhs, rhs);
  [U2, V2] = ndgrid (1:N, 1:N);
  for c = 0:sp.q - 1
    left = Sc(c + 1, A(sub2ind ([N, N], U2(:), V2(:))) + 1);
    right = A(sub2ind ([N, N], Sc(c + 1, U2(:)) + 1, Sc(c + 1, V2(:)) + 1));
    ok = ok && isequal (left(:), right(:));
    for d = 0:sp.q - 1
      ok = ok && isequal (Sc(mod (c + d, sp.q) + 1, :), A(sub2ind ([N, N], Sc(c + 1, :) + 1, Sc(d + 1, :) + 1)));
      ok = ok && isequal (Sc(mod (c * d, sp.q) + 1, :), Sc(c + 1, Sc(d + 1, :) + 1));
    end
  end
end

function total = extend_count (sp, cands, chosen, k)
  if (numel (chosen) == k)
    total = 1;
    return;
  end
  members = span (sp, chosen);
  ok = ~ismember (cands, members);
  for u = chosen
    ok = ok & (sp.dot(cands + 1, u + 1)' == 0);
  end
  total = 0;
  for v = cands(ok)
    total = total + extend_count (sp, cands, [chosen, v], k);
  end
end

function n_codes = self_dual_count (n)
  % Ordered bases of n/2 pairwise orthogonal even words, over the bases per code.
  sp = make_space (2, n);
  k = n / 2;
  cands = find (diag (sp.dot) == 0)' - 1;
  cands = cands(cands ~= 0);
  n_codes = extend_count (sp, cands, [], k) / ordered_independent (2, k, k);
end

function ok = has_information_set (sp, S)
  % True when some k = dim S positions take every pattern exactly once across S.
  k = dimension (S, 2);
  if (k == 0)
    ok = numel (S) == 1;
    return;
  end
  J = nchoosek (1:sp.n, k);
  ok = false;
  for r = 1:rows (J)
    if (rows (unique (sp.V(S + 1, J(r, :)), 'rows')) == numel (S))
      ok = true;
      return;
    end
  end
end


% The checks -------------------------------------------------------------------

function check_scaling ()
  section ('Bit strings that can be scaled');
  s23 = make_space (2, 3);
  s52 = make_space (5, 2);
  s33 = make_space (3, 3);
  v = ix (s23, '101');
  same ('vectors in GF(2)^3', s23.N, 8);
  same ('0 times 101', show (s23, vscale (s23, 0, v)), '000');
  same ('1 times 101', show (s23, vscale (s23, 1, v)), '101');
  same ('101 + 011 over GF(2)', show (s23, vadd (s23, v, ix (s23, '011'))), '110');
  e = 0:s23.N - 1;
  require ('v + v = 0 for every v in GF(2)^3', all (s23.add(sub2ind ([8, 8], e + 1, e + 1)) == 0));
  same ('vectors in GF(5)^2', s52.N, 25);
  same ('vectors in GF(3)^3', s33.N, 27);
  a = ix (s52, '12');
  same ('(1,2) + (3,4) over GF(5)', show (s52, vadd (s52, a, ix (s52, '34'))), '(4,1)');
  same ('3 times (1,2) over GF(5)', show (s52, vscale (s52, 3, a)), '(3,1)');
  same ('(1,2) + (1,2) over GF(5)', show (s52, vadd (s52, a, a)), '(2,4)');
  same ('c times (1,2) for c = 0, 1, 2, 3, 4', as_row (s52, s52.scl(:, a + 1)'), '(0,0) (1,2) (2,4) (3,1) (4,3)');
  same ('3 times (1,2) over the integers', as_tuple (3 * [1, 2]), '(3,6)');
  same ('(1,2) + (3,4) over the integers', as_tuple ([1, 2] + [3, 4]), '(4,6)');
  ok = true;
  for v = 1:s52.N - 1
    m = s52.scl(2:5, v + 1);
    ok = ok && numel (unique (m)) == 4 && all (m ~= 0);
  end
  require ('each nonzero word of GF(5)^2 has four different nonzero multiples', ok);
  ok = true;
  for v = 1:s23.N - 1
    m = s23.scl(:, v + 1);
    ok = ok && isequal (m(m ~= 0), v);
  end
  require ('a nonzero binary word has one nonzero multiple, itself', ok);
  E = even_words (s23);
  same ('the even-weight words of length 3', as_set (s23, E), '{000, 011, 101, 110}');
  sums = s23.add(E + 1, E + 1);
  require ('they are closed under addition', all (ismember (sums(:), E)));
  same ('the coset 001 + E', as_set (s23, vadd (s23, ix (s23, '001'), E)), '{001, 010, 100, 111}');
  ok = true;
  for u = 0:7
    for v = 0:7
      if (weight (s23, vadd (s23, u, v)) == 1)
        ok = ok && mod (weight (s23, u), 2) ~= mod (weight (s23, v), 2);
      end
    end
  end
  require ('every edge of the cube joins an even-weight word to an odd one', ok);
end

function check_fields ()
  section ('Fields as scalars');
  ok = true;
  for n = 2:12
    has_inverses = true;
    for a = 1:n - 1
      has_inverses = has_inverses && any (mod (a * (0:n - 1), n) == 1);
    end
    ok = ok && (has_inverses == isprime (n));
  end
  require ('Z_n is a field exactly when n is prime, n = 2 to 12', ok);
  same ('an inverse of 2 modulo 4 exists', any (mod (2 * (0:3), 4) == 1), false);
  same ('2 times 2 modulo 4', mod (2 * 2, 4), 0);
  ok = true;
  for q = [2, 3, 5]
    [A, B, C] = ndgrid (0:q - 1, 0:q - 1, 0:q - 1);
    ok = ok && isequal (mod (A .* (B + C), q), mod (A .* B + A .* C, q)) && isequal (mod ((A + B) .* C, q), mod (A .* C + B .* C, q));
  end
  require ('both distributive laws in GF(2), GF(3) and GF(5)', ok);
  ok = true;
  for q = [2, 3, 5]
    ok = ok && all (mod (0 * (0:q - 1), q) == 0);
  end
  require ('0 a = 0 for every a in GF(2), GF(3) and GF(5)', ok);
  same ('the negative of 1 in GF(2)', mod (-1, 2), 1);
  same ('the negative of 2 in GF(5)', mod (-2, 5), 3);
end

function check_axioms ()
  section ('The vector space axioms');
  require ('the eight axioms on every triple of GF(2)^3', axioms_hold (make_space (2, 3)));
  require ('the eight axioms on every triple of GF(3)^2', axioms_hold (make_space (3, 2)));
  s52 = make_space (5, 2);
  require ('the eight axioms on every triple of GF(5)^2', axioms_hold (s52));
  e = 0:s52.N - 1;
  require ('0 v = 0 for every v in GF(5)^2', all (s52.scl(1, :) == 0));
  require ('c 0 = 0 for every c in GF(5)', all (s52.scl(:, 1) == 0));
  require ('(-1) v is the negative of v for every v in GF(5)^2', all (s52.add(sub2ind ([25, 25], e + 1, s52.scl(5, :) + 1)) == 0));
  same ('the negative of (1,2) in GF(5)^2', show (s52, vscale (s52, 4, ix (s52, '12'))), '(4,3)');
  same ('4 times (1,2) over the integers', as_tuple (4 * [1, 2]), '(4,8)');
  same ('(1,2) + (4,3) over the integers', as_tuple ([1, 2] + [4, 3]), '(5,5)');
  same ('(1,2) + (4,3) over GF(5)', show (s52, vadd (s52, ix (s52, '12'), ix (s52, '43'))), '(0,0)');
end

function check_subspaces ()
  section ('Subspaces');
  s23 = make_space (2, 3);
  n_sub = 0;
  n_closed = 0;
  for mask = 0:255
    S = find (bitand (mask, 2 .^ (0:7))) - 1;
    n_sub = n_sub + is_subspace (s23, S);
    if (~isempty (S))
      sums = s23.add(S + 1, S + 1);
      n_closed = n_closed + all (ismember (sums(:), S));
    end
  end
  same ('subsets of GF(2)^3', 2 ^ 8, 256);
  same ('subsets that pass the subspace test', n_sub, 16);
  same ('nonempty subsets closed under addition alone', n_closed, 16);
  s25 = make_space (2, 5);
  require ('C5 contains 00000 and is closed under addition', is_subspace (s25, c5 (s25)));
  odd3 = find (mod (sum (s23.V, 2), 2) == 1)' - 1;
  same ('the odd-weight words of length 3 contain 000', any (odd3 == 0), false);
  same ('100 + 010, two words of weight 1', show (s23, vadd (s23, ix (s23, '100'), ix (s23, '010'))), '110');
  s52 = make_space (5, 2);
  ok = true;
  for c = 0:4
    for v = 0:s52.N - 1
      t = 0;
      for i = 1:c
        t = vadd (s52, t, v);
      end
      ok = ok && vscale (s52, c, v) == t;
    end
  end
  require ('over GF(5), c v is v added to itself c times', ok);
  mul4 = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
  [A, B, C] = ndgrid (0:3, 0:3, 0:3);
  at4 = @(x, y) sub2ind ([4, 4], x + 1, y + 1);
  AB = mul4(at4 (A, B));
  field4 = isequal (mul4, mul4') ...
    && isequal (mul4(at4 (A, mul4(at4 (B, C)))), mul4(at4 (AB, C))) ...
    && isequal (mul4(at4 (A, bitxor (B, C))), bitxor (AB, mul4(at4 (A, C)))) ...
    && isequal (mul4(2, :), 0:3) && all (any (mul4(2:4, :) == 1, 2));
  require ('GF(4) with XOR and the product table x^2 = x + 1 is a field', field4);
  [a, b] = ndgrid (0:1, 0:1);
  same ('{0, 1} closed under addition in GF(4)', all (ismember (bitxor (a(:), b(:)), [0, 1])), true);
  same ('the scalar 2 of GF(4) times the vector 1', mul4(3, 2), 2);
  sc = mul4(:, 1:2);
  same ('{0, 1} closed under scaling by GF(4)', all (ismember (sc(:), [0, 1])), false);
end

function check_span ()
  section ('Span');
  s24 = make_space (2, 4);
  gens = [ix(s24, '1100'), ix(s24, '0110'), ix(s24, '0011')];
  sizes = arrayfun (@(j) numel (span (s24, gens(1:j))), 0:3);
  same ('sizes of the spans of the first 0, 1, 2, 3 vectors', as_ints (sizes), '1 2 4 8');
  same ('span of 1100', as_set (s24, span (s24, gens(1:1))), '{0000, 1100}');
  same ('span of 1100, 0110', as_set (s24, span (s24, gens(1:2))), '{0000, 0110, 1010, 1100}');
  same ('span of 1100, 0110, 0011', as_set (s24, span (s24, gens)), '{0000, 0011, 0101, 0110, 1001, 1010, 1100, 1111}');
  require ('that span is the set of even-weight words of length 4', isequal (span (s24, gens), even_words (s24)));
  s52 = make_space (5, 2);
  same ('span of (1,2) in GF(5)^2', as_set (s52, span (s52, ix (s52, '12'))), '{(0,0), (1,2), (2,4), (3,1), (4,3)}');
  same ('size of the span of (1,2), (0,1)', numel (span (s52, [ix(s52, '12'), ix(s52, '01')])), 25);
  same ('size of the span of (1,2), (2,4)', numel (span (s52, [ix(s52, '12'), ix(s52, '24')])), 5);
  ok = true;
  for k = 0:3
    L = all_lists (s24.N, k);
    for r = 1:rows (L)
      ok = ok && is_subspace (s24, span (s24, L(r, :)));
    end
  end
  require ('every span in GF(2)^4 of up to 3 vectors is a subspace', ok);
  ok = true;
  for qn = [2, 4; 3, 2; 5, 2]'
    sp = make_space (qn(1), qn(2));
    levels = subspaces (sp);
    for d = 1:numel (levels)
      for r = 1:rows (levels{d})
        S = levels{d}(r, :);
        for w = setdiff (0:sp.N - 1, S)
          T = sp.add(S + 1, sp.scl(:, w + 1)' + 1);
          cosets = unique (T(:))';
          ok = ok && isequal (span (sp, [S, w]), cosets) && numel (cosets) == sp.q * numel (S);
        end
      end
    end
  end
  require ('a vector outside a subspace multiplies its size by q', ok);
  same ('y - 2x mod 5 up the column x = 0', as_ints (mod ((0:4) - 2 * 0, 5)), '0 1 2 3 4');
  same ('y - 2x mod 5 up the column x = 1', as_ints (mod ((0:4) - 2 * 1, 5)), '3 4 0 1 2');
  L = span (s52, ix (s52, '12'));
  digit = mod (s52.V(:, 2) - 2 * s52.V(:, 1), 5);
  same ('points of GF(5)^2 with digit 0 form the span of (1,2)', isequal (find (digit == 0)' - 1, L), true);
  P = s52.V(vadd (s52, ix (s52, '01'), L) + 1, :);
  same ('the digit on the coset (0,1) + span of (1,2)', as_ints (unique (mod (P(:, 2) - 2 * P(:, 1), 5))'), '1');
end

function check_independence ()
  section ('Linear independence');
  s24 = make_space (2, 4);
  a = ix (s24, '1100');
  b = ix (s24, '0110');
  c = ix (s24, '0011');
  same ('1100, 0110, 0011 independent', independent (s24, [a, b, c]), true);
  same ('1100, 0110, 1010 independent', independent (s24, [a, b, ix(s24, '1010')]), false);
  same ('1100 + 0110', show (s24, vadd (s24, a, b)), '1010');
  same ('1100, 0110, 0011, 1111 independent', independent (s24, [a, b, c, ix(s24, '1111')]), false);
  same ('1100 + 0011', show (s24, vadd (s24, a, c)), '1111');
  s52 = make_space (5, 2);
  same ('(1,2), (2,4) independent in GF(5)^2', independent (s52, [ix(s52, '12'), ix(s52, '24')]), false);
  same ('(1,2), (0,1) independent in GF(5)^2', independent (s52, [ix(s52, '12'), ix(s52, '01')]), true);
  s23 = make_space (2, 3);
  ok = true;
  for k = 1:3
    C = combos_rows (s23, all_lists (s23.N, k));
    ok = ok && isequal (sum (C == 0, 2) == 1, all (C(:, 2:end) ~= 0, 2));
  end
  require ('over GF(2), independent means no nonempty sublist sums to 0', ok);
  ok = true;
  for qn = [2, 3; 3, 2]'
    sp = make_space (qn(1), qn(2));
    for k = 0:sp.n - 1
      L = all_lists (sp.N, k);
      for r = 1:rows (L)
        vs = L(r, :);
        if (independent (sp, vs))
          members = span (sp, vs);
          for v = setdiff (0:sp.N - 1, members)
            ok = ok && independent (sp, [vs, v]);
          end
        end
      end
    end
  end
  require ('adding a vector outside the span keeps a list independent', ok);
  same ('1100 + 0110 + 1010 over GF(2)', show (s24, vadd (s24, vadd (s24, ix (s24, '1100'), ix (s24, '0110')), ix (s24, '1010'))), '0000');
  same ('2 (1,2) + 4 (2,4) over the integers', as_tuple (2 * [1, 2] + 4 * [2, 4]), '(10,20)');
  same ('2 (1,2) + 4 (2,4) over GF(5)', show (s52, vadd (s52, vscale (s52, 2, ix (s52, '12')), vscale (s52, 4, ix (s52, '24')))), '(0,0)');
end

function check_basis ()
  section ('Basis and dimension');
  s23 = make_space (2, 3);
  b = [ix(s23, '110'), ix(s23, '011')];
  for target = {'101', '011'}
    found = {};
    for c1 = 0:1
      for c2 = 0:1
        if (vadd (s23, vscale (s23, c1, b(1)), vscale (s23, c2, b(2))) == ix (s23, target{1}))
          found{end + 1} = sprintf ('(%d, %d)', c1, c2);
        end
      end
    end
    if (strcmp (target{1}, '101'))
      stated = '(1, 1)';
    else
      stated = '(0, 1)';
    end
    same (sprintf ('coordinates of %s in the basis 110, 011', target{1}), strjoin (found, ' '), stated);
  end
  ok = true;
  for qnt = [2, 4, 3; 3, 3, 2; 5, 2, 2]'
    sp = make_space (qnt(1), qnt(2));
    for k = 0:qnt(3)
      C = combos_rows (sp, all_lists (sp.N, k));
      ind = sum (C == 0, 2) == 1;
      distinct = sum (diff (sort (C, 2), 1, 2) ~= 0, 2) + 1;
      ok = ok && all (distinct(ind) == sp.q ^ k);
    end
  end
  require ('k independent vectors give q^k distinct combinations', ok);
  s25 = make_space (2, 5);
  same ('dimension of C5', dimension (c5 (s25), 2), 2);
  same ('dimension of R5', dimension ([0, ix(s25, '11111')], 2), 1);
  same ('dimension of the even-weight words of length 4', dimension (even_words (make_space (2, 4)), 2), 3);
  s53 = make_space (5, 3);
  S = span (s53, [ix(s53, '120'), ix(s53, '013')]);
  same ('size of the span of (1,2,0), (0,1,3) in GF(5)^3', numel (S), 25);
  same ('its dimension', dimension (S, 5), 2);
  levels = subspaces (s23);
  planes = levels{3};
  sizes = [];
  for p = 1:rows (planes)
    P = planes(p, :);
    for k = 1:3
      L = P(all_lists (numel (P), k) + 1);
      if (k == 1)
        L = L(:);
      end
      for r = 1:rows (L)
        if (independent (s23, L(r, :)) && isequal (span (s23, L(r, :)), P))
          sizes(end + 1) = k;
        end
      end
    end
  end
  same ('sizes of the bases of the planes of GF(2)^3', as_ints (unique (sizes)), '2');
  g = [ix(s25, '10110'), ix(s25, '01101')];
  same ('01101 + 10110 over GF(2)', show (s25, vadd (s25, g(2), g(1))), '11011');
  same ('the message 11 with the basis 10110, 01101', show (s25, vadd (s25, g(1), g(2))), '11011');
  same ('the 4 messages encode to the words of C5, each once', isequal (sort (combos (s25, g)), c5 (s25)), true);
end

function check_counting ()
  section ('Counting bases and subspaces');
  s22 = make_space (2, 2);
  s23 = make_space (2, 3);
  same ('ordered bases of GF(2)^2, by trying every list', lists_independent (s22, 2), 6);
  same ('ordered bases of GF(2)^3, by trying every list', lists_independent (s23, 3), 168);
  same ('ordered bases of GF(2)^3 by the formula, 7 x 6 x 4', ordered_independent (2, 3, 3), 168);
  same ('ordered independent pairs in GF(2)^3', lists_independent (s23, 2), 42);
  same ('ordered bases of GF(2)^4, by trying every list', lists_independent (make_space (2, 4), 4), 20160);
  same ('ordered bases of GF(3)^2, by trying every list', lists_independent (make_space (3, 2), 2), 48);
  same ('ordered bases of GF(5)^2, by trying every list', lists_independent (make_space (5, 2), 2), 480);
  same ('ordered bases of GF(3)^3, by trying every list', lists_independent (make_space (3, 3), 3), 11232);
  same ('ordered bases of GF(3)^3 by the formula, 26 x 24 x 18', ordered_independent (3, 3, 3), 11232);
  same ('ordered bases of each plane of GF(2)^3', ordered_independent (2, 2, 2), 6);
  qn = [2, 2; 2, 3; 2, 4; 3, 2; 3, 3; 5, 2];
  stated = {'1 3 1', '1 7 7 1', '1 15 35 15 1', '1 4 1', '1 13 13 1', '1 6 1'};
  totals = [5, 16, 67, 6, 28, 8];
  for i = 1:rows (qn)
    q = qn(i, 1);
    n = qn(i, 2);
    levels = subspaces (make_space (q, n));
    counts = cellfun (@rows, levels);
    same (sprintf ('subspaces of GF(%d)^%d by dimension, built one by one', q, n), as_ints (counts), stated{i});
    same (sprintf ('subspaces of GF(%d)^%d by dimension, by the formula', q, n), as_ints (arrayfun (@(k) gaussian (n, k, q), 0:n)), stated{i});
    same (sprintf ('subspaces of GF(%d)^%d in all', q, n), sum (counts), totals(i));
  end
  g7 = arrayfun (@(k) gaussian (7, k, 2), 0:7);
  same ('subspaces of GF(2)^7 by dimension', as_ints (g7), '1 127 2667 11811 11811 2667 127 1');
  same ('subspaces of GF(2)^7 in all', sum (g7), 29212);
  same ('sets of k of the 7 positions, k = 0 to 7', as_ints (arrayfun (@(k) nchoosek (7, k), 0:7)), '1 7 21 35 35 21 7 1');
  ok = true;
  for q = [2, 3, 5]
    for n = 0:8
      for k = 0:n
        ok = ok && gaussian (n, k, q) == gaussian (n, n - k, q);
      end
    end
  end
  require ('the count for k equals the count for n - k, n up to 8, q = 2, 3, 5', ok);
  same ('fraction of lists of 2 vectors of GF(2)^2 that are bases', fraction_text (ordered_independent (2, 2, 2), 2 ^ 4), '3/8');
  same ('fraction of lists of 3 vectors of GF(2)^3 that are bases', fraction_text (ordered_independent (2, 3, 3), 2 ^ 9), '21/64');
  same ('that fraction from the product, GF(2)^2', sprintf ('%.3f', basis_fraction (2, 2)), '0.375');
  same ('fraction of lists of 10 vectors of GF(2)^10 that are bases', sprintf ('%.3f', basis_fraction (2, 10)), '0.289');
  same ('fraction of lists of 10 vectors of GF(3)^10 that are bases', sprintf ('%.3f', basis_fraction (3, 10)), '0.560');
  same ('fraction of lists of 10 vectors of GF(5)^10 that are bases', sprintf ('%.3f', basis_fraction (5, 10)), '0.760');
  same ('the factors 2^3 - 2^i for i = 0, 1, 2', as_ints (2 ^ 3 - 2 .^ (0:2)), '7 6 4');
  fq = [2, 2; 3, 2; 5, 2; 3, 3];
  ftext = {'3 2', '8 6', '24 20', '26 24 18'};
  for i = 1:rows (fq)
    q = fq(i, 1);
    n = fq(i, 2);
    same (sprintf ('the factors q^n - q^i for GF(%d)^%d', q, n), as_ints (q ^ n - q .^ (0:n - 1)), ftext{i});
  end
  same ('bases over lists, GF(2)^2, before reducing', sprintf ('%d/%d', ordered_independent (2, 2, 2), 2 ^ 4), '6/16');
  same ('bases over lists, GF(2)^3, before reducing', sprintf ('%d/%d', ordered_independent (2, 3, 3), 2 ^ 9), '168/512');
  p = 1.0;
  tail = 1.0;
  ok_low = true;
  ok_down = true;
  ok_tail = true;
  for n = 1:60
    nw = p * (1 - 2 ^ (-n));
    ok_down = ok_down && nw <= p;
    p = nw;
    ok_low = ok_low && p >= 0.25;
    if (n >= 2)
      tail = tail * (1 - 2 ^ (-n));
      ok_tail = ok_tail && tail >= 0.5;
    end
  end
  require ('over GF(2) the product never increases, n = 1 to 60', ok_down);
  require ('the factors after the first keep it at least 1/2, n = 2 to 60', ok_tail);
  require ('over GF(2) the product stays at or above 1/4, n = 1 to 60', ok_low);
end

function check_dot ()
  section ('Dot products that are not positive');
  s22 = make_space (2, 2);
  s23 = make_space (2, 3);
  s24 = make_space (2, 4);
  same ('11 . 11 over GF(2)', s22.dot(ix (s22, '11') + 1, ix (s22, '11') + 1), 0);
  same ('11 . 11 over the integers', sum ([1, 1] .* [1, 1]), 2);
  same ('1111 . 1111 over GF(2)', s24.dot(ix (s24, '1111') + 1, ix (s24, '1111') + 1), 0);
  same ('111 . 111 over GF(2)', s23.dot(ix (s23, '111') + 1, ix (s23, '111') + 1), 1);
  same ('1100 . 1010 over GF(2)', s24.dot(ix (s24, '1100') + 1, ix (s24, '1010') + 1), 1);
  same ('1011 . 1101 over GF(2)', s24.dot(ix (s24, '1011') + 1, ix (s24, '1101') + 1), 0);
  x = [1, 1, 1, 1];
  ands = x .* x;
  same ('the four AND outputs for 1111 and 1111', char (ands + '0'), '1111');
  same ('the two first-level XOR outputs', as_ints ([xor(ands(1), ands(2)), xor(ands(3), ands(4))]), '0 0');
  same ('the final XOR output', double (xor (xor (ands(1), ands(2)), xor (ands(3), ands(4)))), 0);
  ok = true;
  for n = 1:8
    V = make_space (2, n).V;
    ok = ok && isequal (mod (sum (V .* V, 2), 2) == 0, mod (sum (V ~= 0, 2), 2) == 0);
  end
  require ('over GF(2), v . v = 0 exactly when v has even weight, length 1 to 8', ok);
  same ('self-orthogonal vectors in GF(2)^4', sum (diag (s24.dot) == 0), 8);
  s33 = make_space (3, 3);
  same ('(1,1,1) . (1,1,1) over GF(3)', s33.dot(ix (s33, '111') + 1, ix (s33, '111') + 1), 0);
  s32 = make_space (3, 2);
  d32 = diag (s32.dot)';
  same ('nonzero self-orthogonal vectors in GF(3)^2', sum (d32(2:end) == 0), 0);
  s52 = make_space (5, 2);
  same ('(1,2) . (1,2) over GF(5)', s52.dot(ix (s52, '12') + 1, ix (s52, '12') + 1), 0);
  d52 = diag (s52.dot)';
  so5 = find (d52 == 0) - 1;
  so5 = so5(so5 ~= 0);
  same ('nonzero self-orthogonal vectors in GF(5)^2', numel (so5), 8);
  same ('those vectors', as_set (s52, so5), '{(1,2), (1,3), (2,1), (2,4), (3,1), (3,4), (4,2), (4,3)}');
  ok = true;
  for sp = {s23, s32, s52}
    s = sp{1};
    N = s.N;
    Dt = s.dot;
    ok = ok && isequal (Dt, Dt');
    [U, V, W] = ndgrid (1:N, 1:N, 1:N);
    lhs = Dt(sub2ind ([N, N], s.add(sub2ind ([N, N], U(:), V(:))) + 1, W(:)));
    rhs = mod (Dt(sub2ind ([N, N], U(:), W(:))) + Dt(sub2ind ([N, N], V(:), W(:))), s.q);
    ok = ok && isequal (lhs, rhs);
    for c = 0:s.q - 1
      ok = ok && isequal (Dt(s.scl(c + 1, :) + 1, :), mod (c * Dt, s.q));
    end
  end
  require ('symmetric and linear in each slot on GF(2)^3, GF(3)^2, GF(5)^2', ok);
  s25 = make_space (2, 5);
  H = checks5 (s25);
  ok = true;
  for w = 0:s25.N - 1
    ok = ok && isequal (syndrome (s25.V(w + 1, :)), s25.dot(w + 1, H + 1));
  end
  require ('the syndrome of C5 is x . 11100, x . 10010, x . 01001 for all 32 words', ok);
  same ('positions where 1100 and 1010 both have a 1', as_ints (find ([1, 1, 0, 0] & [1, 0, 1, 0])), '1');
  same ('1^2 and 2^2 in GF(3)', as_ints (mod ([1, 2] .^ 2, 3)), '1 1');
  require ('over GF(3)^3, v . v is the nonzero-entry count mod 3', all (diag (s33.dot) == mod (sum (s33.V ~= 0, 2), 3)));
  same ('(1,2) . (1,2) over the integers', sum ([1, 2] .* [1, 2]), 5);
  both = setdiff (union (span (s52, ix (s52, '12')), span (s52, ix (s52, '13'))), 0);
  same ('they fill the lines through (1,2) and (1,3), less the origin', isequal (both, so5), true);
  require ('(x, y) . (3, 1) = y - 2x mod 5 at every point of GF(5)^2', isequal (s52.dot(:, ix (s52, '31') + 1), mod (s52.V(:, 2) - 2 * s52.V(:, 1), 5)));
  same ('-2 modulo 5', mod (-2, 5), 3);
end

function check_dual ()
  section ('Orthogonal complements and dual codes');
  s25 = make_space (2, 5);
  C = c5 (s25);
  D = dual (s25, C);
  H = checks5 (s25);
  same ('the dual of C5', as_set (s25, D), '{00000, 00111, 01001, 01110, 10010, 10101, 11011, 11100}');
  same ('its size', numel (D), 8);
  same ('its dimension', dimension (D, 2), 3);
  same ('dim C5 + dim of its dual', dimension (C, 2) + dimension (D, 2), 5);
  same ('words in both C5 and its dual', as_set (s25, intersect (C, D)), '{00000, 11011}');
  require ('11100, 10010, 01001 lie in the dual of C5', all (ismember (H, D)));
  same ('11100, 10010, 01001 independent', independent (s25, H), true);
  same ('their span is the dual of C5', isequal (span (s25, H), D), true);
  same ('the dual of the dual of C5', as_set (s25, dual (s25, D)), '{00000, 01101, 10110, 11011}');
  heads = arrayfun (@(v) char (s25.V(v + 1, 1:2) + '0'), C, 'UniformOutput', false);
  same ('C5 on positions 1 and 2', strjoin (heads, ' '), '00 01 10 11');
  free = {'100', '010', '001'};
  words = cell (1, 3);
  for i = 1:3
    y = free{i} - '0';
    words{i} = char ([mod(y(1) + y(2), 2), mod(y(1) + y(3), 2), y] + '0');
  end
  same ('free choices 100, 010, 001 for y3 y4 y5 give', strjoin (words, ' '), '11100 10010 01001');
  s23 = make_space (2, 3);
  same ('the dual of R3', as_set (s23, dual (s23, [0, ix(s23, '111')])), '{000, 011, 101, 110}');
  R5 = [0, ix(s25, '11111')];
  same ('size of the dual of R5', numel (dual (s25, R5)), 16);
  require ('the dual of R5 is the set of even-weight words of length 5', isequal (dual (s25, R5), even_words (s25)));
  s24 = make_space (2, 4);
  Cd = span (s24, [ix(s24, '1100'), ix(s24, '0111')]);
  same ('the code spanned by 1100, 0111', as_set (s24, Cd), '{0000, 0111, 1011, 1100}');
  same ('its dual', as_set (s24, dual (s24, Cd)), '{0000, 0011, 1101, 1110}');
  same ('words in both', as_set (s24, intersect (Cd, dual (s24, Cd))), '{0000}');
  s52 = make_space (5, 2);
  L = span (s52, ix (s52, '12'));
  same ('the dual of the span of (1,2) in GF(5)^2', as_set (s52, dual (s52, L)), '{(0,0), (1,2), (2,4), (3,1), (4,3)}');
  same ('the coset (0,1) + span of (1,2)', as_set (s52, vadd (s52, ix (s52, '01'), L)), '{(0,1), (1,3), (2,0), (3,2), (4,4)}');
  ok = true;
  for c = 0:4
    P = s52.V(vadd (s52, c, L) + 1, :);
    ok = ok && numel (unique (mod (P(:, 2) - 2 * P(:, 1), 5))) == 1;
  end
  require ('y - 2x is the same on every point of each coset of the line', ok);
  ok = true;
  count = 0;
  for qn = [2, 4; 3, 3; 5, 2]'
    sp = make_space (qn(1), qn(2));
    levels = subspaces (sp);
    for d = 1:numel (levels)
      for r = 1:rows (levels{d})
        S = levels{d}(r, :);
        Dd = dual (sp, S);
        ok = ok && dimension (S, sp.q) + dimension (Dd, sp.q) == sp.n && isequal (dual (sp, Dd), S);
        count = count + 1;
      end
    end
  end
  same ('subspaces of GF(2)^4, GF(3)^3 and GF(5)^2 tested', count, 103);
  require ('dim C + dim of its dual = n, and the dual of the dual is C', ok);
  levels = subspaces (s24);
  ok = true;
  for d = 1:numel (levels)
    for r = 1:rows (levels{d})
      ok = ok && has_information_set (s24, levels{d}(r, :));
    end
  end
  require ('every subspace of GF(2)^4 has an information set', ok);
  Y = s25.V;
  ok = isequal (s25.dot(:, ix (s25, '10110') + 1), mod (Y(:, 1) + Y(:, 3) + Y(:, 4), 2)) ...
    && isequal (s25.dot(:, ix (s25, '01101') + 1), mod (Y(:, 2) + Y(:, 3) + Y(:, 5), 2));
  require ('y . 10110 = y1 + y3 + y4 and y . 01101 = y2 + y3 + y5, all y', ok);
  tails = arrayfun (@(h) char (s25.V(h + 1, 3:5) + '0'), H, 'UniformOutput', false);
  same ('11100, 10010, 01001 on positions 3, 4, 5', strjoin (tails, ' '), '100 010 001');
  D5 = dual (s25, c5 (s25));
  ok = true;
  for y = D5
    t = 0;
    for j = 1:3
      t = vadd (s25, t, vscale (s25, s25.V(y + 1, j + 2), H(j)));
    end
    ok = ok && t == y;
  end
  require ('every y in the dual is y3 11100 + y4 10010 + y5 01001', ok);
  X = s25.V(c5 (s25) + 1, :);
  require ('in C5, x3 = x1 + x2, x4 = x1 and x5 = x2', all (X(:, 3) == mod (X(:, 1) + X(:, 2), 2)) && all (X(:, 4) == X(:, 1)) && all (X(:, 5) == X(:, 2)));
end

function check_self_dual ()
  section ('Self-dual codes');
  s22 = make_space (2, 2);
  two = [0, ix(s22, '11')];
  same ('{00, 11} equals its dual', isequal (dual (s22, two), two), true);
  same ('binary self-dual codes of length 2, 4, 6', as_ints (arrayfun (@self_dual_count, [2, 4, 6])), '1 3 15');
  s24 = make_space (2, 4);
  levels = subspaces (s24);
  planes = levels{3};
  found = {};
  for r = 1:rows (planes)
    if (isequal (dual (s24, planes(r, :)), planes(r, :)))
      found{end + 1} = planes(r, :);
    end
  end
  same ('the self-dual codes of length 4', strjoin (cellfun (@(S) as_set (s24, S), found, 'UniformOutput', false), ' '), '{0000, 0011, 1100, 1111} {0000, 0101, 1010, 1111} {0000, 0110, 1001, 1111}');
  require ('1111 lies in each of them', all (cellfun (@(S) any (S == ix (s24, '1111')), found)));
  s28 = make_space (2, 8);
  G = [ix(s28, '11110000'), ix(s28, '00111100'), ix(s28, '00001111'), ix(s28, '01010101')];
  dots = [];
  for i = 1:4
    for j = i:4
      dots(end + 1) = s28.dot(G(i) + 1, G(j) + 1);
    end
  end
  same ('the ten dot products among the four rows', as_ints (dots), '0 0 0 0 0 0 0 0 0 0');
  same ('the four rows independent', independent (s28, G), true);
  H = span (s28, G);
  same ('words in the length-8 code', numel (H), 16);
  same ('weights of its words, weight:count', weight_counts (s28, H), '0:1 4:14 8:1');
  same ('it equals its dual', isequal (dual (s28, H), H), true);
  same ('it contains 11111111', any (H == ix (s28, '11111111')), true);
  s34 = make_space (3, 4);
  T = [ix(s34, '1110'), ix(s34, '0121')];
  tdots = [s34.dot(T(1) + 1, T(1) + 1), s34.dot(T(1) + 1, T(2) + 1), s34.dot(T(2) + 1, T(2) + 1)];
  same ('the three dot products among the two rows over GF(3)', as_ints (tdots), '0 0 0');
  tetra = span (s34, T);
  same ('words in the ternary code of length 4', numel (tetra), 9);
  same ('weights of its words, weight:count', weight_counts (s34, tetra), '0:1 3:8');
  same ('it equals its dual', isequal (dual (s34, tetra), tetra), true);
  s52 = make_space (5, 2);
  levels = subspaces (s52);
  lines = levels{2};
  pts = [];
  n_sd = 0;
  for r = 1:rows (lines)
    S = lines(r, :);
    if (isequal (dual (s52, S), S))
      n_sd = n_sd + 1;
      pts = [pts, S(s52.V(S + 1, 1)' == 1)];
    end
  end
  same ('self-dual lines of GF(5)^2', n_sd, 2);
  same ('their nonzero points with first entry 1', as_row (s52, pts), '(1,2) (1,3)');
  s32 = make_space (3, 2);
  levels = subspaces (s32);
  lines = levels{2};
  n_sd = 0;
  for r = 1:rows (lines)
    n_sd = n_sd + isequal (dual (s32, lines(r, :)), lines(r, :));
  end
  same ('self-dual lines of GF(3)^2', n_sd, 0);
  same ('prod of 2^i + 1 for i = 1 to n/2 - 1, n = 2, 4, 6', as_ints (arrayfun (@(n) prod (2 .^ (1:n / 2 - 1) + 1), [2, 4, 6])), '1 3 15');
  s23 = make_space (2, 3);
  E = even_words (s23);
  d23 = diag (s23.dot);
  require ('every word of E has v . v = 0', all (d23(E + 1) == 0));
  same ('011 . 101 over GF(2)', s23.dot(ix (s23, '011') + 1, ix (s23, '101') + 1), 1);
  Tv = s34.V(T + 1, :);
  same ('the ternary dot products over the integers', as_ints ([Tv(1, :) * Tv(1, :)', Tv(1, :) * Tv(2, :)', Tv(2, :) * Tv(2, :)']), '3 3 6');
  same ('1110, 0121 independent over GF(3)', independent (s34, T), true);
end

function check_worked_examples ()
  section ('Worked examples');
  s33 = make_space (3, 3);
  T = find (mod (sum (s33.V, 2), 3) == 0)' - 1;
  same ('vectors of GF(3)^3 whose entries sum to 0', numel (T), 9);
  same ('they pass the subspace test', is_subspace (s33, T), true);
  B = [ix(s33, '120'), ix(s33, '012')];
  same ('(1,2,0), (0,1,2) independent', independent (s33, B), true);
  same ('their span is that set', isequal (span (s33, B), T), true);
  line = span (s33, ix (s33, '111'));
  same ('that set is the dual of the span of (1,1,1)', isequal (dual (s33, line), T), true);
  same ('(1,1,1) lies in it', any (T == ix (s33, '111')), true);
  s53 = make_space (5, 3);
  v1 = ix (s53, '120');
  v2 = ix (s53, '013');
  v3 = ix (s53, '203');
  same ('2 (1,2,0) + (0,1,3) over GF(5)', show (s53, vadd (s53, vscale (s53, 2, v1), v2)), '(2,0,3)');
  same ('(1,2,0), (0,1,3), (2,0,3) independent', independent (s53, [v1, v2, v3]), false);
  same ('size of their span', numel (span (s53, [v1, v2, v3])), 25);
  same ('lines of GF(3)^3', gaussian (3, 1, 3), 13);
  same ('planes of GF(3)^3', gaussian (3, 2, 3), 13);
  same ('lines inside each plane of GF(3)^3', gaussian (2, 1, 3), 4);
  levels = subspaces (s33);
  planes = levels{3};
  n_in = 0;
  for r = 1:rows (planes)
    n_in = n_in + all (ismember (line, planes(r, :)));
  end
  same ('planes of GF(3)^3 that contain the line through (1,1,1)', n_in, 4);
  same ('those vectors', as_set (s33, T), '{(0,0,0), (0,1,2), (0,2,1), (1,0,2), (1,1,1), (1,2,0), (2,0,1), (2,1,0), (2,2,2)}');
  same ('their dimension', dimension (T, 3), 2);
  d1 = s53.V(v1 + 1, :);
  d2 = s53.V(v2 + 1, :);
  d3 = s53.V(v3 + 1, :);
  same ('2 (1,2,0) over GF(5)', show (s53, vscale (s53, 2, v1)), '(2,4,0)');
  same ('(2,4,0) + (0,1,3) over the integers', as_tuple (mod (2 * d1, 5) + d2), '(2,5,3)');
  same ('4 (2,0,3) over the integers', as_tuple (4 * d3), '(8,0,12)');
  same ('4 (2,0,3) over GF(5)', show (s53, vscale (s53, 4, v3)), '(3,0,2)');
  same ('2 v1 + v2 + 4 v3 over the integers', as_tuple (mod (2 * d1, 5) + d2 + mod (4 * d3, 5)), '(5,5,5)');
  same ('2 v1 + v2 + 4 v3 over GF(5)', show (s53, vadd (s53, vadd (s53, vscale (s53, 2, v1), v2), vscale (s53, 4, v3))), '(0,0,0)');
  same ('(1,2,0), (0,1,3) independent', independent (s53, [v1, v2]), true);
  found = {};
  for c1 = 0:4
    for c2 = 0:4
      if (vadd (s53, vscale (s53, c1, v1), vscale (s53, c2, v2)) == v3)
        found{end + 1} = sprintf ('(%d, %d)', c1, c2);
      end
    end
  end
  same ('coordinates of (2,0,3) in the basis (1,2,0), (0,1,3)', strjoin (found, ' '), '(2, 1)');
  same ('dimension of the span of all three', dimension (span (s53, [v1, v2, v3]), 5), 2);
  s24 = make_space (2, 4);
  C4 = [ix(s24, '0000'), ix(s24, '1100'), ix(s24, '0111'), ix(s24, '1011')];
  same ('0000, 1100, 0111, 1011 is the span of 1100, 0111', isequal (sort (C4), span (s24, [ix(s24, '1100'), ix(s24, '0111')])), true);
  heads = arrayfun (@(v) char (s24.V(v + 1, 1:2) + '0'), C4, 'UniformOutput', false);
  same ('positions 1, 2 of 0000, 1100, 0111, 1011', strjoin (heads, ' '), '00 11 01 10');
  free = [0, 0; 1, 0; 0, 1; 1, 1];
  words = zeros (1, 4);
  for i = 1:4
    s = mod (free(i, 1) + free(i, 2), 2);
    words(i) = [s, s, free(i, :)] * s24.w';
  end
  same ('free choices 00, 10, 01, 11 for y3 y4 give', as_row (s24, words), '0000 1110 1101 0011');
  same ('those words are the dual of that code', isequal (sort (words), dual (s24, C4)), true);
  same ('0111 . 0111 over GF(2)', s24.dot(ix (s24, '0111') + 1, ix (s24, '0111') + 1), 1);
  s28 = make_space (2, 8);
  R = [1, 1, 1, 1, 0, 0, 0, 0; 0, 0, 1, 1, 1, 1, 0, 0; 0, 0, 0, 0, 1, 1, 1, 1; 0, 1, 0, 1, 0, 1, 0, 1];
  pairs = [1, 2; 1, 4; 2, 3; 2, 4; 3, 4; 1, 3];
  shared = cell (1, rows (pairs));
  for t = 1:rows (pairs)
    s = sprintf ('%d', find (R(pairs(t, 1), :) & R(pairs(t, 2), :)));
    if (isempty (s))
      s = '-';
    end
    shared{t} = s;
  end
  same ('positions shared by r1r2 r1r4 r2r3 r2r4 r3r4 r1r3', strjoin (shared, ' '), '34 24 56 46 68 -');
  holders = arrayfun (@(p) sprintf ('%d', find (R(:, p))'), [1, 7, 8, 3], 'UniformOutput', false);
  same ('rows holding a 1 at positions 1, 7, 8 and 3', strjoin (holders, ' '), '1 3 34 12');
  g = (R * s28.w')';
  same ('r1 + r3', show (s28, vadd (s28, g(1), g(3))), '11111111');
  H = span (s28, g);
  X = s28.V(H + 1, :);
  M = X * X';
  W = sum (X, 2);
  S = s28.add(H + 1, H + 1);
  WS = reshape (weight (s28, S(:)), numel (H), numel (H));
  require ('wt(x + y) = wt(x) + wt(y) - 2m on every pair of the code', isequal (WS, W + W' - 2 * M));
  require ('m is even on every pair of the code', all (mod (M(:), 2) == 0));
  wH = weight (s28, H);
  same ('its words of weight 0 and 8', as_row (s28, H(wH == 0 | wH == 8)), '00000000 11111111');
end

function check_exercises ()
  section ('Exercises');
  s32 = make_space (3, 2);
  levels = subspaces (s32);
  lines = levels{2};
  gens = zeros (1, rows (lines));
  for r = 1:rows (lines)
    for v = lines(r, :)
      d = s32.V(v + 1, :);
      nz = d(d ~= 0);
      if (~isempty (nz) && nz(1) == 1)
        gens(r) = v;
        break;
      end
    end
  end
  same ('lines of GF(3)^2, each named by a point with leading entry 1', as_row (s32, sort (gens)), '(0,1) (1,0) (1,1) (1,2)');
  s24 = make_space (2, 4);
  same ('1101 + 0110', show (s24, vadd (s24, ix (s24, '1101'), ix (s24, '0110'))), '1011');
  same ('1101, 0110, 1011 independent', independent (s24, [ix(s24, '1101'), ix(s24, '0110'), ix(s24, '1011')]), false);
  same ('2-dimensional subspaces of GF(2)^4, 15 x 14 / (3 x 2)', gaussian (4, 2, 2), 35);
  s25 = make_space (2, 5);
  B = [ix(s25, '11000'), ix(s25, '01100'), ix(s25, '00110'), ix(s25, '00011')];
  same ('11000, 01100, 00110, 00011 independent', independent (s25, B), true);
  same ('their span is the dual of R5', isequal (span (s25, B), dual (s25, [0, ix(s25, '11111')])), true);
  same ('the squares in GF(3)', as_ints (unique (mod ((0:2) .^ 2, 3))), '0 1');
  same ('the solutions of b^2 = 4 in GF(5)', as_ints (find (mod ((0:4) .^ 2, 5) == 4) - 1), '2 3');
  s52 = make_space (5, 2);
  same ('(1,3) . (1,3) over GF(5)', s52.dot(ix (s52, '13') + 1, ix (s52, '13') + 1), 0);
  same ('points on each line of GF(3)^2', as_ints (arrayfun (@(r) numel (unique (lines(r, :))), 1:rows (lines))), '3 3 3 3');
  same ('nonzero points of GF(3)^2', sum (any (s32.V ~= 0, 2)), 8);
  ok = true;
  for v = 1:s32.N - 1
    ok = ok && sum (any (lines == v, 2)) == 1;
  end
  require ('each nonzero point of GF(3)^2 lies on exactly one line', ok);
  same ('the factors q^n - q^i for GF(2)^4', as_ints (2 ^ 4 - 2 .^ (0:3)), '15 14 12 8');
  same ('weights of 11000, 01100, 00110, 00011', as_ints (weight (s25, B)), '2 2 2 2');
  same ('0^2, 1^2, 2^2 in GF(3)', as_ints (mod ((0:2) .^ 2, 3)), '0 1 1');
  w2 = find (sum (s24.V, 2) == 2)' - 1;
  same ('words of weight 2 in GF(2)^4', numel (w2), 6);
  flips = {};
  for v = w2
    f = vadd (s24, v, ix (s24, '1111'));
    if (v < f)
      flips{end + 1} = as_set (s24, [v, f]);
    end
  end
  same ('each paired with its flip', strjoin (flips, ' '), '{0011, 1100} {0101, 1010} {0110, 1001}');
  same ('(0,1) . (0,1) over GF(5)', s52.dot(ix (s52, '01') + 1, ix (s52, '01') + 1), 1);
  same ('(1,3) . (1,3) over the integers', sum ([1, 3] .* [1, 3]), 10);
end


check_scaling ();
check_fields ();
check_axioms ();
check_subspaces ();
check_span ();
check_independence ();
check_basis ();
check_counting ();
check_dot ();
check_dual ();
check_self_dual ();
check_worked_examples ();
check_exercises ();
printf ('\nAll chapter 15 numbers reproduced.\n');
