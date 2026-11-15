% Reproduce every number in chapter 23 of Error Correction Coding,
% "Weight Distributions and Undetected Errors", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch23_weight_distributions_and_undetected_errors.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch23_weight_distributions_and_undetected_errors.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Whole numbers, bit patterns, weight distributions and polynomial
% coefficients are compared exactly. A probability is compared to half a unit
% in the last digit the chapter prints.
%
% A word is a row of bits, first bit on the left, and a code is a matrix
% whose rows are its codewords in counting order. The MacWilliams transform
% multiplies out (x + y) and (x - y) with integer coefficients. The longer
% codes of the last section are reached from their dual words through the
% same transform run backwards, as Krawtchouk sums, whose coefficients are
% read off the multiplied-out polynomials (the chapter proves the two agree,
% and the character-sum section checks it for n = 7).
%
% A probability curve is sampled at p = j/1000 for j = 0 to 500. Each grid
% point is formed by a scalar division and every power by repeated
% multiplication, in the same order as the Python route, so both routes print
% the same digits.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


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

function check (label, computed, stated)
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
  c = as_text (computed);
  s = as_text (stated);
  printf ('  %-62s %-16s %s\n', label, c, s);
  if (~strcmp (c, s))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, c, s);
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
    printf ('MISMATCH in ''%s'': the claim does not hold\n', label);
    exit (1);
  end
end

function section (title)
  printf ('\n%s\n', title);
end

function s = word (v)
  s = char (v + '0');
end

function s = listing (v)
  % Adding 0 turns a negative zero into 0, so it never prints as -0.
  s = strjoin (arrayfun (@(x) sprintf ('%d', x + 0), v, 'UniformOutput', false), ',');
end

function v = parse (text)
  v = double (text - '0');
end

function s = yesno (condition)
  if (condition)
    s = 'yes';
  else
    s = 'no';
  end
end

function s = thousandths (j)
  s = sprintf ('0.%03d', j);
end

function c = binom (n, k)
  if (k < 0 || k > n)
    c = 0;
    return;
  end
  c = 1;
  for t = 1:k
    c = c * (n - k + t) / t;
  end
end


% Codes ------------------------------------------------------------------------

function W = all_words (n)
  % Every word of length n, in counting order, first bit most significant.
  x = (0:2 ^ n - 1)';
  W = zeros (2 ^ n, n);
  for t = 1:n
    W(:, t) = bitget (x, n - t + 1);
  end
end

function s = syndrome (H, e)
  s = mod (H * e(:), 2)';
end

function C = code_from_checks (H)
  W = all_words (columns (H));
  C = W(~any (mod (W * H', 2), 2), :);
end

function C = span (G)
  % Every sum of rows of G, in counting order of the selecting message.
  C = mod (all_words (rows (G)) * G, 2);
end

function A = distribution (C)
  A = accumarray (sum (C, 2) + 1, 1, [columns(C) + 1, 1])';
end

function H = shortened (H3, n)
  H = H3(:, 1:n);
end

function r = reach (e)
  % Positions from the first flip of a nonzero pattern to the last, inclusive.
  at = find (e);
  r = at(end) - at(1) + 1;
end

function out = polymul (a, b)
  out = zeros (1, numel (a) + numel (b) - 1);
  for i = 1:numel (a)
    for j = 1:numel (b)
      out(i + j - 1) = out(i + j - 1) + a(i) * b(j);
    end
  end
end

function out = term (n, i)
  % The coefficients of y^0 to y^n in (x + y)^(n - i) (x - y)^i.
  out = 1;
  for t = 1:n - i
    out = [out, 0] + [0, out];
  end
  for t = 1:i
    out = [out, 0] - [0, out];
  end
end

function out = transform (A)
  % The coefficients of W(x + y, x - y), multiplied out.
  n = numel (A) - 1;
  out = zeros (1, n + 1);
  for i = 0:n
    if (A(i + 1))
      out = out + A(i + 1) * term (n, i);
    end
  end
end

function B = macwilliams (A)
  sz = sum (A);
  raw = transform (A);
  if (any (mod (raw, sz)))
    printf ('MISMATCH: W(x + y, x - y) is not divisible by the size of the code\n');
    exit (1);
  end
  B = raw / sz;
end

function k = krawtchouk (n, j, i)
  k = 0;
  for s = 0:j
    if (j - s <= n - i)
      k = k + (-1) ^ s * binom (i, s) * binom (n - i, j - s);
    end
  end
end

function Km = kmatrix (n)
  % Row i + 1 holds K_0(i) to K_n(i), the coefficients of term (n, i).
  Km = zeros (n + 1, n + 1);
  for i = 0:n
    Km(i + 1, :) = term (n, i);
  end
end

function B = columns_dual (cols, r)
  % The weight distribution of the dual code of the code whose check columns are cols.
  n = numel (cols);
  a = (0:2 ^ r - 1)';
  M = bitand (repmat (a, 1, n), repmat (cols(:)', 2 ^ r, 1));
  count = zeros (size (M));
  for b = 1:r
    count = count + bitget (M, b);
  end
  B = accumarray (sum (mod (count, 2), 2) + 1, 1, [n + 1, 1])';
end

function A = from_dual (B, r, Km)
  n = numel (B) - 1;
  if (nargin < 3)
    Km = kmatrix (n);
  end
  total = B * Km;
  if (any (mod (total, 2 ^ r)))
    printf ('MISMATCH: a Krawtchouk sum is not divisible by the size of the dual code\n');
    exit (1);
  end
  A = total / 2 ^ r;
end

function rk = gf2_rank (values)
  basis = [];
  for v = values(:)'
    for b = basis
      v = min (v, bitxor (v, b));
    end
    if (v)
      basis(end + 1) = v;
    end
  end
  rk = numel (basis);
end


% Probabilities ----------------------------------------------------------------

function out = powers (base, n)
  out = zeros (1, n + 1);
  out(1) = 1;
  for t = 1:n
    out(t + 1) = out(t) * base;
  end
end

function total = pud (A, p)
  n = numel (A) - 1;
  pw = powers (p, n);
  qw = powers (1 - p, n);
  total = 0;
  for i = 1:n
    total = total + A(i + 1) * pw(i + 1) * qw(n - i + 1);
  end
end

function [first, last] = pud_dual_parts (B, r, p)
  n = numel (B) - 1;
  tw = powers (1 - 2 * p, n);
  total = 0;
  for j = 0:n
    total = total + B(j + 1) * tw(j + 1);
  end
  first = total / 2 ^ r;
  qw = powers (1 - p, n);
  last = qw(n + 1);
end

function P = grid_points ()
  P = zeros (1, 501);
  for j = 0:500
    P(j + 1) = j / 1000;
  end
end

function v = curve (A)
  % pud (A, j/1000) for j = 0 to 500, every grid point in the same order of operations.
  n = numel (A) - 1;
  P = grid_points ();
  Q = 1 - P;
  pw = ones (n + 1, numel (P));
  qw = ones (n + 1, numel (P));
  for t = 1:n
    pw(t + 1, :) = pw(t, :) .* P;
    qw(t + 1, :) = qw(t, :) .* Q;
  end
  v = zeros (1, numel (P));
  for i = 1:n
    v = v + A(i + 1) * pw(i + 1, :) .* qw(n - i + 1, :);
  end
end

function [top, j] = peak (A)
  [top, index] = max (curve (A));
  j = index - 1;
end

function r = rising (v)
  r = all (v(2:end) >= v(1:end - 1));
end


% The checks -------------------------------------------------------------------

function check_silent_check (K)
  section ('When a check stays silent: the (7,4) Hamming code as a detector');
  H3 = K.H3;
  same ('row 1 of H_3', word (H3(1, :)), '0001111');
  same ('row 2 of H_3', word (H3(2, :)), '0110011');
  same ('row 3 of H_3', word (H3(3, :)), '1010101');
  examples = {'0010000', '011'; '0100100', '111'; '0000111', '100'; ...
              '1110000', '000'; '1001100', '000'; '0001111', '000'; ...
              '1111111', '000'};
  for t = 1:rows (examples)
    e = examples{t, 1};
    same (sprintf ('syndrome of the error pattern %s', e), word (syndrome (H3, parse (e))), examples{t, 2});
  end
  W = all_words (7);
  patterns = W(any (W, 2), :);
  missed = patterns(~any (mod (patterns * H3', 2), 2), :);
  same ('nonzero error patterns of length 7', rows (patterns), 127);
  same ('patterns with syndrome 000', rows (missed), 15);
  same ('patterns caught', rows (patterns) - rows (missed), 112);
  C = K.HAMMING;
  require ('the missed patterns are the nonzero codewords', isequal (missed, C(any (C, 2), :)));
  A = distribution (C);
  stated_total = [7 21 35 35 21 7 1];
  stated_missed = [0 0 7 7 0 0 1];
  for w = 1:7
    same (sprintf ('patterns of weight %d', w), binom (7, w), stated_total(w));
    same (sprintf ('patterns of weight %d that are missed', w), A(w + 1), stated_missed(w));
  end
  same ('weight-3 patterns caught', binom (7, 3) - A(4), 28);
  require ('every pattern of weight 1 or 2 is caught', A(2) == 0 && A(3) == 0);
end

function check_distributions (K)
  section ('Counting codewords by weight');
  same ('weight distribution of the repetition code of length 3', listing (distribution (K.REP3)), '1,0,0,1');
  same ('weight distribution of the even-weight code of length 3', listing (distribution (K.EVEN3)), '1,0,3,0');
  same ('weight distribution of the even-weight code of length 4', listing (distribution (K.EVEN4)), '1,0,6,0,1');
  same ('weight distribution of the repetition code of length 4', listing (distribution (K.REP4)), '1,0,0,0,1');
  same ('weight distribution of the (7,4) Hamming code', listing (distribution (K.HAMMING)), '1,0,0,7,7,0,0,1');
  same ('weight distribution of the (7,3) simplex code', listing (distribution (K.SIMPLEX)), '1,0,0,0,7,0,0,0');
  same ('weight distribution of the extended (8,4) Hamming code', listing (distribution (K.EXTENDED)), '1,0,0,0,14,0,0,0,1');
  same ('W(1,1) for the (7,4) Hamming code', sum (distribution (K.HAMMING)), 16);
  w = sum (K.HAMMING, 2);
  same ('minimum distance of the (7,4) code, first nonzero weight', min (w(w > 0)), 3);
  names = {'(7,4) Hamming', 'extended (8,4)', 'even-weight length 4'};
  sets = {K.HAMMING, K.EXTENDED, K.EVEN4};
  for t = 1:3
    C = sets{t};
    A = distribution (C);
    ok = true;
    for r = 1:rows (C)
      ok = ok && isequal (distribution (mod (C + C(r, :), 2)), A);
    end
    require (sprintf ('distances from every word of the %s code match A', names{t}), ok);
  end
  c = parse ('0110011');
  C = K.HAMMING;
  require ('0110011 is a codeword of the (7,4) code', ismember (c, C, 'rows'));
  D = mod (C + c, 2);
  same ('distances from 0110011, counted by value', listing (distribution (D)), '1,0,0,7,7,0,0,1');
  near = C(sum (D, 2) == 3, :);
  same ('the first codeword at distance 3 from 0110011', word (near(1, :)), '0010110');
  same ('their sum, a codeword of weight 3', word (mod (c + parse ('0010110'), 2)), '0100101');
  names = {'repetition length 3', 'even-weight length 4', '(7,4) Hamming', 'extended (8,4)'};
  sets = {K.REP3, K.EVEN4, K.HAMMING, K.EXTENDED};
  for t = 1:4
    A = distribution (sets{t});
    require (sprintf ('A_i = A_(n-i) for the %s code', names{t}), isequal (A, fliplr (A)));
  end
  A = distribution (K.EVEN3);
  require ('the even-weight code of length 3 lacks 111 and is not symmetric', ~isequal (A, fliplr (A)));
end

function check_probability (K)
  section ('The probability of an undetected error');
  A = distribution (K.HAMMING);
  p = 0.01;
  q = 1 - p;
  pw = powers (p, 7);
  qw = powers (q, 7);
  check ('7 p^3 (1 - p)^4 at p = 0.01', 7 * pw(4) * qw(5), '6.72417e-6');
  check ('7 p^4 (1 - p)^3 at p = 0.01', 7 * pw(5) * qw(4), '6.79209e-8');
  check ('p^7 at p = 0.01', pw(8), '1e-14');
  check ('P_ud of the (7,4) code at p = 0.01', pud (A, 0.01), '6.79209e-6');
  p = 0.1;
  q = 1 - p;
  pw = powers (p, 7);
  qw = powers (q, 7);
  check ('7 p^3 (1 - p)^4 at p = 0.1', 7 * pw(4) * qw(5), '0.0045927');
  check ('7 p^4 (1 - p)^3 at p = 0.1', 7 * pw(5) * qw(4), '0.0005103');
  check ('p^7 at p = 0.1', pw(8), '1e-7');
  check ('P_ud of the (7,4) code at p = 0.1', pud (A, 0.1), '0.0051031');
  check ('P_ud of the (7,4) code at p = 1/2', pud (A, 0.5), '0.1171875');
  same ('P_ud of the (7,4) code at p = 1/2, as a fraction', sprintf ('%d/%d', sum (A) - 1, 2 ^ 7), '15/128');
  w = powers (0.99, 7);
  anyerr = 1 - w(8);
  check ('probability of any error in 7 bits at p = 0.01', anyerr, '0.0679347');
  check ('share of erroneous words that are missed at p = 0.01', pud (A, 0.01) / anyerr, '9.998e-5');
  check ('share of erroneous words that are missed at p = 1/2', pud (A, 0.5) / (1 - 0.5 ^ 7), '0.11811');
  check ('the rule of thumb 2^-3', 2 ^ -3, '0.125');
  names = {'repetition length 3', 'even-weight length 4', '(7,4) Hamming', 'extended (8,4)'};
  sets = {K.REP3, K.EVEN4, K.HAMMING, K.EXTENDED};
  stated = {'1e-6', '0.001', '0.125', '1/8'; ...
            '5.8807e-4', '0.0487', '0.4375', '7/16'; ...
            '6.79209e-6', '0.0051031', '0.1171875', '15/128'; ...
            '1.34483e-7', '9.1855e-4', '0.05859375', '15/256'};
  for t = 1:4
    D = distribution (sets{t});
    n = numel (D) - 1;
    k = round (log2 (rows (sets{t})));
    check (sprintf ('P_ud at p = 0.01, %s', names{t}), pud (D, 0.01), stated{t, 1});
    check (sprintf ('P_ud at p = 0.1, %s', names{t}), pud (D, 0.1), stated{t, 2});
    check (sprintf ('P_ud at p = 1/2, %s', names{t}), pud (D, 0.5), stated{t, 3});
    same (sprintf ('(2^k - 1)/2^n for the %s code', names{t}), sprintf ('%d/%d', 2 ^ k - 1, 2 ^ n), stated{t, 4});
    require (sprintf ('P_ud(1/2) equals that fraction for the %s code', names{t}), pud (D, 0.5) == (2 ^ k - 1) / 2 ^ n);
    require (sprintf ('P_ud(1/2) = 2^-(n-k) - 2^-n for the %s code', names{t}), pud (D, 0.5) == 2 ^ -(n - k) - 2 ^ -n);
  end

  section ('Short bursts on the (7,4) code');
  W = all_words (7);
  W = W(2:end, :);
  keep = false (rows (W), 1);
  for t = 1:rows (W)
    keep(t) = reach (W(t, :)) <= 4;
  end
  B = W(keep, :);
  same ('nonzero patterns whose flips lie within four neighboring positions', rows (B), 39);
  stated = [7 0; 15 0; 13 2; 4 2];
  wt = sum (B, 2);
  missed = ~any (mod (B * K.H3', 2), 2);
  found = zeros (4, 2);
  for w = 1:4
    found(w, :) = [sum(wt == w), sum(wt == w & missed)];
    same (sprintf ('bursts of weight %d', w), found(w, 1), stated(w, 1));
    same (sprintf ('codewords among them, weight %d', w), found(w, 2), stated(w, 2));
  end
  require ('no burst has weight 5 or more', all (wt <= 4));
  runs = true;
  for t = find (wt == 4)'
    runs = runs && reach (B(t, :)) == 4;
  end
  require ('every burst of weight 4 is a run of four neighboring flips', runs);
  lists = {'0010110 1110000', '0001111 0111100'};
  for w = 3:4
    hits = sort (cellfun (@word, num2cell (B(wt == w & missed, :), 2), 'UniformOutput', false));
    same (sprintf ('the missed bursts of weight %d', w), strjoin (hits', ' '), lists{w - 2});
  end
  check ('share of weight-3 bursts missed, 2/13', found(3, 2) / found(3, 1), '0.154');
  check ('share of all weight-3 patterns missed, 7/35', A(4) / binom (7, 3), '0.2');
  check ('share of weight-4 bursts missed, 2/4', found(4, 2) / found(4, 1), '0.5');
  check ('share of all weight-4 patterns missed, 7/35', A(5) / binom (7, 4), '0.2');
  require ('weight-3 bursts are missed less often than weight-3 patterns overall', ...
           found(3, 2) / found(3, 1) < A(4) / binom (7, 3));
  require ('runs of four are missed more than twice as often as weight-4 patterns overall', ...
           found(4, 2) / found(4, 1) > 2 * A(5) / binom (7, 4));
  missed_all = sum (found(:, 2));
  same ('short bursts that are codewords', missed_all, 4);
  check ('share of the 39 short bursts that are missed, 4/39', missed_all / rows (B), '0.103');
  check ('share of all nonzero patterns that are missed, 15/127', (sum (A) - 1) / (2 ^ 7 - 1), '0.118');
  require ('short bursts are missed less often than nonzero patterns overall', ...
           missed_all / rows (B) < (sum (A) - 1) / (2 ^ 7 - 1));
end

function check_character_sum (K)
  section ('The MacWilliams identity: the character sum and the double sum');
  names = {'repetition length 3', '(7,4) Hamming'};
  codes = {K.REP3, K.HAMMING};
  duals = {K.EVEN3, K.SIMPLEX};
  for t = 1:2
    C = codes{t};
    U = all_words (columns (C));
    ok = true;
    for r = 1:rows (U)
      u = U(r, :);
      total = sum (1 - 2 * mod (C * u', 2));
      if (ismember (u, duals{t}, 'rows'))
        want = rows (C);
      else
        want = 0;
      end
      ok = ok && total == want;
    end
    require (sprintf ('sum of (-1)^(u.c) over the %s code, every u', names{t}), ok);
  end
  U = all_words (3);
  sums = zeros (1, 8);
  for r = 1:8
    sums(r) = sum (1 - 2 * mod (K.REP3 * U(r, :)', 2));
  end
  same ('row sums of the sign grid for the repetition code', listing (sums), '2,0,0,2,0,2,2,0');
  hits = arrayfun (@(r) word (U(r, :)), find (sums == 2), 'UniformOutput', false);
  same ('words u whose row sum is 2', strjoin (hits, ' '), '000 011 101 110');
  for t = 1:2
    C = codes{t};
    n = columns (C);
    U = all_words (n);
    wu = sum (U, 2);
    double_sum = zeros (1, n + 1);
    factored = true;
    for r = 1:rows (C)
      c = C(r, :);
      inner = accumarray (wu + 1, 1 - 2 * mod (U * c', 2), [n + 1, 1])';
      factored = factored && isequal (inner, term (n, sum (c)));
      double_sum = double_sum + inner;
    end
    require (sprintf ('inner sum over u factors per coordinate, %s', names{t}), factored);
    dual_dist = macwilliams (distribution (C));
    require (sprintf ('double sum equals |C| times the dual enumerator, %s', names{t}), isequal (double_sum, rows (C) * dual_dist));
    require (sprintf ('double sum equals W(x + y, x - y), %s', names{t}), isequal (double_sum, transform (distribution (C))));
  end
  n = 7;
  same ('K_4(0) for n = 7', krawtchouk (n, 4, 0), 35);
  same ('K_4(3) for n = 7', krawtchouk (n, 4, 3), 3);
  same ('K_4(4) for n = 7', krawtchouk (n, 4, 4), 3);
  same ('K_4(7) for n = 7', krawtchouk (n, 4, 7), 35);
  same ('B_4 = (35 + 7*3 + 7*3 + 35)/16', floor ((35 + 7 * 3 + 7 * 3 + 35) / 16), 7);
  A = distribution (K.HAMMING);
  B = zeros (1, n + 1);
  for j = 0:n
    total = 0;
    for i = 0:n
      total = total + A(i + 1) * krawtchouk (n, j, i);
    end
    B(j + 1) = floor (total / 16);
  end
  same ('dual distribution of the (7,4) code by Krawtchouk sums', listing (B), '1,0,0,0,7,0,0,0');
  ok = true;
  for i = 0:7
    coefficients = term (7, i);
    for j = 0:7
      ok = ok && coefficients(j + 1) == krawtchouk (7, j, i);
    end
  end
  require ('K_j(i) is the coefficient of x^(n-j) y^j in (x+y)^(n-i) (x-y)^i, n = 7', ok);
end

function check_small_identities (K)
  section ('Checking the identity on small codes');
  t = term (3, 0) + term (3, 3);
  same ('(x+y)^3 + (x-y)^3, coefficients of y^0..y^3', listing (t), '2,0,6,0');
  same ('divided by |C| = 2', listing (macwilliams (distribution (K.REP3))), '1,0,3,0');
  same ('dual distribution read off directly', listing (distribution (K.EVEN3)), '1,0,3,0');
  t = term (3, 0) + 3 * term (3, 2);
  same ('(x+y)^3 + 3(x+y)(x-y)^2', listing (t), '4,0,0,4');
  same ('divided by |C| = 4', listing (macwilliams (distribution (K.EVEN3))), '1,0,0,1');
  outer = term (7, 0) + term (7, 7);
  same ('(x+y)^7 + (x-y)^7', listing (outer), '2,0,42,0,70,0,14,0');
  inner = 7 * term (7, 3) + 7 * term (7, 4);
  same ('7(x+y)^4(x-y)^3 + 7(x+y)^3(x-y)^4', listing (inner), '14,0,-42,0,42,0,-14,0');
  cube = polymul ([1, 0, -1], polymul ([1, 0, -1], [1, 0, -1]));
  same ('14x(x^2 - y^2)^3, multiplied out', listing ([14 * cube, 0]), '14,0,-42,0,42,0,-14,0');
  same ('the whole of W(x+y, x-y) for the (7,4) code', listing (outer + inner), '16,0,0,0,112,0,0,0');
  same ('divided by 16', listing (macwilliams (distribution (K.HAMMING))), '1,0,0,0,7,0,0,0');
  same ('the simplex code counted directly', listing (distribution (K.SIMPLEX)), '1,0,0,0,7,0,0,0');
  w = sum (K.SIMPLEX, 2);
  require ('every nonzero simplex word has weight 4', all (w(w > 0) == 4));
  same ('back from the simplex code, divided by 8', listing (macwilliams (distribution (K.SIMPLEX))), '1,0,0,7,7,0,0,1');
  same ('the extended (8,4) code transformed, divided by 16', listing (macwilliams (distribution (K.EXTENDED))), '1,0,0,0,14,0,0,0,1');
  require ('the extended (8,4) code equals its own dual', isequal (sortrows (K.EXTENDED), sortrows (span (K.HEXT))));
  same ('(x+y)^4 + 6(x+y)^2(x-y)^2 + (x-y)^4', listing (transform (distribution (K.EVEN4))), '8,0,0,0,8');
  same ('even-weight length 4 transformed, divided by 8', listing (macwilliams (distribution (K.EVEN4))), '1,0,0,0,1');
  same ('the repetition code of length 4 counted directly', listing (distribution (K.REP4)), '1,0,0,0,1');
end

function check_dual_form (K)
  section ('Undetected errors from the dual code');
  B = distribution (K.SIMPLEX);
  [first, last] = pud_dual_parts (B, 3, 0.01);
  w = powers (0.98, 4);
  check ('(1 - 2p)^4 at p = 0.01', w(5), '0.92236816');
  check ('(1/8)(1 + 7(0.98)^4)', first, '0.93207214');
  check ('(0.99)^7', last, '0.93206535');
  check ('their difference, P_ud of the (7,4) code at p = 0.01', first - last, '6.79209e-6');
  [first, last] = pud_dual_parts (B, 3, 0.1);
  check ('(1/8)(1 + 7(0.8)^4)', first, '0.4834');
  check ('(0.9)^7', last, '0.4782969');
  check ('their difference at p = 0.1', first - last, '0.0051031');
  A = distribution (K.HAMMING);
  for j = 0:500
    p = j / 1000;
    [a, b] = pud_dual_parts (B, 3, p);
    if (abs ((a - b) - pud (A, p)) > 1e-15)
      require ('dual form agrees with the direct sum on the grid', false);
    end
  end
  require ('dual form agrees with the direct sum on the grid', true);
  ok = true;
  for j = 0:500
    u = powers (1 - j / 1000, 6);
    v = powers (1 - 2 * (j / 1000), 3);
    ok = ok && 7 * u(7) - 7 * v(4) >= 0;
  end
  require ('7(1-p)^6 - 7(1-2p)^3 >= 0 on the grid', ok);
  require ('P_ud of the (7,4) code rises along the grid', rising (curve (A)));
  E = distribution (K.EXTENDED);
  check ('P_ud of the extended (8,4) code at 0.01, direct sum', pud (E, 0.01), '1.34483e-7');
  [first, last] = pud_dual_parts (E, 4, 0.01);
  check ('(1/16)(1 + 14(0.98)^4 + (0.98)^8)', first, '0.922744829');
  check ('(0.99)^8', last, '0.922744694');
  check ('their difference', first - last, '1.34483e-7');
  check ('relative gap between the two routes, powers by repeated multiplication', abs ((first - last) - pud (E, 0.01)) / pud (E, 0.01), '8e-10');
end

function check_proper_and_good (K)
  section ('Proper codes and good codes');
  names = {'(7,4) Hamming', 'extended (8,4)'};
  sets = {K.HAMMING, K.EXTENDED};
  checks = [3 4];
  for t = 1:2
    [top, j] = peak (distribution (sets{t}));
    same (sprintf ('grid point of the largest P_ud, %s', names{t}), thousandths (j), '0.500');
    require (sprintf ('largest P_ud on the grid is below 2^-%d, %s', checks(t), names{t}), top < 2 ^ -checks(t));
  end
  names = {'(7,4) Hamming', 'extended (8,4)', '(6,3) shortened Hamming', 'triplicated parity (8,5)'};
  duals = {K.SIMPLEX, K.EXTENDED, span(shortened (K.H3, 6)), span(K.TRIPLE_H)};
  lengths = [7 8 6 8];
  good = [true true true false];
  lowest = [4 4 3 2];
  for t = 1:4
    w = sum (duals{t}, 2);
    low = min (w(w > 0));
    n = lengths(t);
    same (sprintf ('smallest nonzero dual weight, %s', names{t}), low, lowest(t));
    same (sprintf ('is it at least n/2 = %g? %s', n / 2, names{t}), yesno (2 * low >= n), yesno (good(t)));
  end
  TH = K.TRIPLE_H;
  same ('row 1 of the triplicated-parity H', word (TH(1, :)), '00000011');
  same ('row 2 of the triplicated-parity H', word (TH(2, :)), '00000101');
  same ('row 3 of the triplicated-parity H', word (TH(3, :)), '11111001');
  T = K.TRIPLE;
  same ('number of codewords', rows (T), 32);
  parity = mod (sum (T(:, 1:5), 2), 2);
  require ('each codeword is five bits and their parity three times', ...
           all (T(:, 6) == T(:, 7) & T(:, 7) == T(:, 8) & T(:, 8) == parity));
  A = distribution (T);
  same ('weight distribution of the triplicated-parity code', listing (A), '1,0,10,0,10,0,10,0,1');
  same ('weight distribution of its dual', listing (distribution (span (TH))), '1,0,3,0,0,0,3,0,1');
  same ('the same dual by MacWilliams', listing (macwilliams (A)), '1,0,3,0,0,0,3,0,1');
  check ('P_ud of the triplicated-parity code at p = 0.3', pud (A, 0.3), '0.12897');
  [top, j] = peak (A);
  check ('its largest value on the grid', top, '0.12899');
  same ('where the largest value sits', thousandths (j), '0.304');
  check ('its value at p = 1/2', pud (A, 0.5), '0.12109375');
  same ('as a fraction', sprintf ('%d/%d', rows (T) - 1, 2 ^ 8), '31/256');
  check ('share of weight-2 patterns that are codewords, 10/28', A(3) / binom (8, 2), '0.357');
  same ('number of weight-2 patterns of length 8', binom (8, 2), 28);
  require ('P_ud exceeds 1/8 at p = 0.3', pud (A, 0.3) > 0.125);
  stated_count = [1 3 7 14 24];
  stated_bad = [0 0 0 0 2];
  for n = 4:8
    % Each multiset of n nonzero 3-bit columns, as n + 6 choose n shifted down.
    combos = nchoosek (1:n + 6, n) - (0:n - 1);
    Km = kmatrix (n);
    found = zeros (0, n + 1);
    for t = 1:rows (combos)
      cols = combos(t, :);
      if (gf2_rank (cols) < 3)
        continue;
      end
      B = columns_dual (cols, 3);
      if (B(2))
        continue;
      end
      found(end + 1, :) = from_dual (B, 3, Km);
    end
    seen = unique (found, 'rows');
    bad = cell (0, 3);
    for t = 1:rows (seen)
      [top, j] = peak (seen(t, :));
      if (top > 0.125)
        bad(end + 1, :) = {seen(t, :), top, j};
      end
    end
    same (sprintf ('length %d: distinct weight distributions', n), rows (seen), stated_count(n - 3));
    same (sprintf ('length %d: distributions rising above 1/8 on the grid', n), rows (bad), stated_bad(n - 3));
  end
  same ('first not-good distribution at length 8', listing (bad{1, 1}), '1,0,10,0,10,0,10,0,1');
  same ('second not-good distribution at length 8', listing (bad{2, 1}), '1,0,10,1,5,10,0,5,0');
  check ('its largest value on the grid', bad{2, 2}, '0.12925');
  same ('where it sits', thousandths (bad{2, 3}), '0.302');
end

function check_lengths (K)
  section ('Why the length matters: the (7,4) code shortened to n = 4, 5, 6, 7');
  stated = {4, '1,0,0,1,0', '9.9e-7', '1/16', '0.0625'; ...
            5, '1,0,0,2,1,0', '1.9701e-6', '3/32', '0.09375'; ...
            6, '1,0,0,4,3,0,0', '3.9106e-6', '7/64', '0.109375'; ...
            7, '1,0,0,7,7,0,0,1', '6.79209e-6', '15/128', '0.1171875'};
  for t = 1:rows (stated)
    n = stated{t, 1};
    C = code_from_checks (shortened (K.H3, n));
    D = distribution (C);
    same (sprintf ('weight distribution at n = %d', n), listing (D), stated{t, 2});
    check (sprintf ('P_ud at p = 0.01, n = %d', n), pud (D, 0.01), stated{t, 3});
    same (sprintf ('P_ud at p = 1/2 as a fraction, n = %d', n), sprintf ('%d/%d', rows (C) - 1, 2 ^ n), stated{t, 4});
    check (sprintf ('P_ud at p = 1/2, n = %d', n), pud (D, 0.5), stated{t, 5});
    require (sprintf ('P_ud rises along the grid, n = %d', n), rising (curve (D)));
  end
  same ('dual distribution at n = 5', listing (distribution (span (shortened (K.H3, 5)))), '1,0,2,4,1,0');
  same ('dual distribution at n = 6', listing (distribution (span (shortened (K.H3, 6)))), '1,0,0,4,3,0,0');
  same ('dual distribution at n = 4', listing (distribution (span (shortened (K.H3, 4)))), '1,1,3,3,0');
  C4 = code_from_checks (shortened (K.H3, 4));
  require ('at n = 4 the fourth bit is 0 in every codeword', all (C4(:, 4) == 0));
  % C_5 of Distance, Detection and Correction.
  C5 = [parse('00000'); parse('01101'); parse('10110'); parse('11011')];
  five = code_from_checks (shortened (K.H3, 5));
  same ('weight distribution of C_5 of Distance, Detection and Correction', listing (distribution (C5)), '1,0,0,2,1,0');
  moved = sort (cellfun (@word, num2cell (C5(:, [3 1 4 2 5]), 2), 'UniformOutput', false));
  same ('C_5 with its positions read in the order 3, 1, 4, 2, 5', strjoin (moved', ' '), '00000 01111 10011 11100');
  listed = sort (cellfun (@word, num2cell (five, 2), 'UniformOutput', false));
  same ('the n = 5 code', strjoin (listed', ' '), '00000 01111 10011 11100');
  shared = 0;
  for t = 1:rows (C5)
    if (any (C5(t, :)) && ismember (C5(t, :), five, 'rows'))
      shared = shared + 1;
    end
  end
  same ('nonzero codewords that C_5 shares with the n = 5 code', shared, 0);
  check ('P_ud of C_5 at p = 0.01', pud (distribution (C5), 0.01), '1.9701e-6');
  P = grid_points ();
  require ('p^2 (3 - 4p) >= 0 on the grid', all (P .* P .* (3 - 4 * P) >= 0));
  require ('p^2 (6(1-p)^2 - p^2) >= 0 on the grid', all (P .* P .* (6 * (1 - P) .^ 2 - P .* P) >= 0));
  require ('(1-p)^4 - (1-2p)^2 >= 0 on the grid', all ((1 - P) .^ 4 - (1 - 2 * P) .^ 2 >= 0));
  check ('6(1-p)^2 = p^2 where p = sqrt(6)/(1 + sqrt(6))', 6 ^ 0.5 / (1 + 6 ^ 0.5), '0.71');

  section ('Why the length matters: single parity at even and odd length');
  stated = {8, 28, '0.0026368', '0.0026361', '0.99974'; ...
            9, 36, '0.0033566', '0.0033554', '0.99964'; ...
            10, 45, '0.0041543', '0.0041524', '0.99952'};
  values = zeros (1, 10);
  for t = 1:rows (stated)
    n = stated{t, 1};
    D = distribution (code_from_checks (ones (1, n)));
    total = pud (D, 0.01);
    a = powers (0.01, 2);
    b = powers (0.99, n - 2);
    w2 = D(3) * a(3) * b(n - 1);
    values(n) = total;
    same (sprintf ('weight-2 codewords at n = %d', n), D(3), stated{t, 2});
    check (sprintf ('P_ud at p = 0.01, n = %d', n), total, stated{t, 3});
    check (sprintf ('the weight-2 term, n = %d', n), w2, stated{t, 4});
    check (sprintf ('its share of P_ud, n = %d', n), w2 / total, stated{t, 5});
    same (sprintf ('all-ones word is a codeword at n = %d', n), yesno (D(n + 1) > 0), yesno (mod (n, 2) == 0));
    require (sprintf ('no odd-weight pattern is a codeword, n = %d', n), all (D(2:2:end) == 0));
  end
  check ('P_ud at n = 9 over P_ud at n = 8', values(9) / values(8), '1.273');
  check ('36/28 times 0.99', 36 / 28 * 0.99, '1.273');

  section ('Why the length matters: five parity checks, n = 16 to 31');
  B16 = columns_dual (1:16, 5);
  require ('n = 16 has a coordinate that is 0 in every codeword', B16(2) > 0);
  require ('n = 15 has only four independent checks', gf2_rank (1:15) == 4);
  % n, A_3, A_4, where the maximum sits in thousandths, and whether the code stays at or below 2^-5.
  stated = [17 36 112 290 0; 18 38 126 309 0; 19 41 148 500 1; 20 45 176 500 1; ...
            21 50 213 500 1; 22 56 259 500 1; 23 63 315 500 1; 24 71 371 500 1; ...
            25 80 438 500 1; 26 90 516 500 1; 27 101 606 500 1; 28 113 706 500 1; ...
            29 126 819 500 1; 30 140 945 500 1; 31 155 1085 500 1];
  shares = zeros (1, rows (stated));
  for t = 1:rows (stated)
    n = stated(t, 1);
    A = from_dual (columns_dual (1:n, 5), 5);
    shares(t) = A(4) / binom (n, 3);
    [top, j] = peak (A);
    same (sprintf ('A_3 at n = %d', n), A(4), stated(t, 2));
    same (sprintf ('A_4 at n = %d', n), A(5), stated(t, 3));
    same (sprintf ('where the largest P_ud on the grid sits, n = %d', n), thousandths (j), thousandths (stated(t, 4)));
    same (sprintf ('stays at or below 2^-5 on the grid, n = %d', n), yesno (top <= 2 ^ -5), yesno (stated(t, 5)));
    if (n == 17)
      check (sprintf ('2^5 times the largest P_ud on the grid, n = %d', n), 32 * top, '1.0961');
    elseif (n == 18)
      check (sprintf ('2^5 times the largest P_ud on the grid, n = %d', n), 32 * top, '1.0188');
    else
      require (sprintf ('the largest value is 2^-5 - 2^-n, at p = 1/2, n = %d', n), ...
               top == pud (A, 0.5) && pud (A, 0.5) == 2 ^ -5 - 2 ^ -n);
    end
  end
  require ('A_3 / C(n, 3) is above 2^-5 at every n from 17 to 31', all (shares > 2 ^ -5));
  A = from_dual (columns_dual (1:17, 5), 5);
  check ('largest P_ud on the grid at n = 17', peak (A), '0.034254');
  check ('2^-5', 2 ^ -5, '0.03125');
  same ('weight-3 patterns of length 17', binom (17, 3), 680);
  check ('share of them that are codewords, 36/680', A(4) / binom (17, 3), '0.0529');
  A15 = from_dual (columns_dual (1:15, 4), 4);
  same ('weight-3 codewords inside the first 15 positions', A15(4), 35);
  check ('their share of the weight-3 patterns there, 35/455', 35 / binom (15, 3), '0.0769');
  A = from_dual (columns_dual (1:31, 5), 5);
  same ('A_3 at n = 31 equals the Hamming count n(n-1)/6', A(4), 31 * 30 / 6);
  same ('weight-3 patterns of length 31', binom (31, 3), 4495);
  check ('share of them that are codewords, 155/4495', A(4) / binom (31, 3), '0.0345');
end

function check_examples (K)
  section ('Worked examples');
  D = distribution (code_from_checks (shortened (K.H3, 6)));
  a = powers (0.1, 4);
  b = powers (0.9, 3);
  check ('(6,3) code: 4 p^3 (1-p)^3 at p = 0.1', 4 * a(4) * b(4), '0.002916');
  check ('(6,3) code: 3 p^4 (1-p)^2 at p = 0.1', 3 * a(5) * b(3), '0.000243');
  check ('(6,3) code: P_ud at p = 0.1', pud (D, 0.1), '0.003159');
  a = powers (0.01, 8);
  b = powers (0.99, 4);
  check ('extended (8,4): 14 p^4 (1-p)^4 at p = 0.01', 14 * a(5) * b(5), '1.34483e-7');
  check ('extended (8,4): p^8 at p = 0.01', a(9), '1e-16');
  [first, last] = pud_dual_parts (distribution (K.EXTENDED), 4, 0.01);
  sa = sprintf ('%.9f', first);
  sb = sprintf ('%.9f', last);
  agree = 0;
  while (sa(3 + agree) == sb(3 + agree))
    agree = agree + 1;
  end
  same (sprintf ('leading digits shared by %s and %s', sa, sb), agree, 6);
end

function check_exercises (K)
  section ('Exercises');
  C = code_from_checks (shortened (K.H3, 5));
  texts = arrayfun (@(r) word (C(r, :)), 1:rows (C), 'UniformOutput', false);
  same ('(5,2) code: its four codewords', strjoin (texts, ' '), '00000 01111 10011 11100');
  same ('(5,2) code: dual distribution by MacWilliams', listing (macwilliams (distribution (C))), '1,0,2,4,1,0');
  check ('(5,2) code: P_ud at p = 1/2', pud (distribution (C), 0.5), '0.09375');
  D = distribution (code_from_checks (shortened (K.H3, 6)));
  same ('(6,3) code: nonzero codewords', sum (D) - 1, 7);
  check ('(6,3) code: P_ud at p = 1/2, 4/64 + 3/64', pud (D, 0.5), '0.109375');
  [first, last] = pud_dual_parts (distribution (K.EVEN3), 2, 0.1);
  check ('repetition length 3 at p = 0.1, direct p^3', pud (distribution (K.REP3), 0.1), '0.001');
  check ('repetition length 3 at p = 0.1, (1/4)(1 + 3(0.8)^2)', first, '0.73');
  check ('repetition length 3 at p = 0.1, (0.9)^3', last, '0.729');
  check ('repetition length 3 at p = 0.1, by the dual', first - last, '0.001');
  P = grid_points ();
  require ('12p(1-p)(1-2p) + 4p^3 >= 0 on the grid', all (12 * P .* (1 - P) .* (1 - 2 * P) + 4 * P .^ 3 >= 0));
  require ('P_ud of the even-weight code of length 4 rises along the grid', rising (curve (distribution (K.EVEN4))));
  check ('even-weight length 4: P_ud(1/2) = 7/16 against 2^-3', pud (distribution (K.EVEN4), 0.5) / 0.125, '3.5');
  same ('single parity, n = 10: weight-2 codewords', binom (10, 2), 45);
  D10 = distribution (code_from_checks (ones (1, 10)));
  D9 = distribution (code_from_checks (ones (1, 9)));
  check ('single parity, n = 10: P_ud at p = 0.01', pud (D10, 0.01), '0.0041543');
  require ('an inverted bus of 10 bits is missed, of 9 bits caught', D10(11) == 1 && D9(10) == 0);
end


% The codes of the chapter -----------------------------------------------------

H3 = zeros (3, 7);
for t = 0:2
  for c = 1:7
    H3(t + 1, c) = bitget (c, 3 - t);
  end
end
K.H3 = H3;
K.HAMMING = code_from_checks (H3);
K.REP3 = code_from_checks ([parse('110'); parse('011')]);
K.EVEN3 = code_from_checks (parse ('111'));
K.EVEN4 = code_from_checks (parse ('1111'));
K.REP4 = code_from_checks ([parse('1100'); parse('0110'); parse('0011')]);
K.SIMPLEX = span (H3);
K.HEXT = [H3, zeros(3, 1); ones(1, 8)];
K.EXTENDED = code_from_checks (K.HEXT);
K.TRIPLE_H = [parse('00000011'); parse('00000101'); parse('11111001')];
K.TRIPLE = code_from_checks (K.TRIPLE_H);

printf ('Chapter 23: Weight Distributions and Undetected Errors\n');
check_silent_check (K);
check_distributions (K);
check_probability (K);
check_character_sum (K);
check_small_identities (K);
check_dual_form (K);
check_proper_and_good (K);
check_lengths (K);
check_examples (K);
check_exercises (K);
printf ('\nAll chapter 23 numbers reproduced.\n');
