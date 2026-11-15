% Reproduce every number in chapter 4 of Error Correction Coding, "Information and Entropy".
%
% This is the GNU Octave twin of python/ch04_information_and_entropy.py. It needs core
% Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch04_information_and_entropy.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. A run that reaches
% the final line has reproduced every number in the chapter. The two routes print
% the same lines, and tools/compare-routes.sh checks that against the Python.
%
% The Python module also holds the entropy and mutual-information functions that
% ch04_charts.py plots. Charts stay Python-only, so DISTRIBUTIONS, which only the
% charts read, is left out here.
%
% Differences from the Python that a side-by-side reader meets:
%   - Octave indexes from 1, so a Python J[(x, y)] is J(x + 1, y + 1) and the
%     Python's counts[s] for a dice sum s is counts(s - 1).
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave functions cannot see script variables, so the module constants of the
%     Python (P_FLIP, LOADED_DIE and the rest) are functions that return a value.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's, and
%     every place the Python calls sum on floats calls py_sum here.
%   - Python's round (x, n) rounds the exact binary value to n decimals. Octave's
%     round takes no second argument, so rnd rounds through sprintf, which also
%     rounds the exact binary value.
%   - A Python dict of joint probabilities becomes a 2-by-2 matrix with the X
%     values down the rows and the Y values across the columns. A dict iterates in
%     insertion order, which for these tables is row by row, so row_major lists a
%     matrix in that same order wherever the order of a sum matters.
%   - Octave already owns the name xor, so the Python's xor table is xor3 here, and
%     the Python's I for mutual information is I_XY, since I is the imaginary unit.
%   - Python's max (range (1, 1000), key=f) returns the first k with the largest f.
%     Octave's max returns the first position of the largest value, which is the
%     same thing when f is evaluated on 1:999.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function check (label, computed, stated)
  value = str2double (strrep (stated, ',', ''));
  computed = double (computed);
  printf ('  %-60s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function require (label, condition)
  if (condition)
    verdict = 'holds';
  else
    verdict = 'FAILS';
  end
  printf ('  %-60s %s\n', label, verdict);
  if (~condition)
    printf ('MISMATCH in ''%s''\n', label);
    exit (1);
  end
end

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

function r = rnd (value, places)
  % The value as the chapter prints it, so sums of printed values can be checked too.
  % Python's round (value, places) and this sprintf both round the exact binary value.
  r = str2double (sprintf ('%.*f', places, value));
end

function v = row_major (M)
  % The entries of a matrix row by row, the order a Python dict of (x, y) keys iterates in.
  M = M.';
  v = M(:).';
end


% Information measures --------------------------------------------------------------

function r = P_FLIP ()
  r = 0.1;          % the flip probability of the chapter's binary symmetric channel
end

function r = RUNNING_Q ()
  r = 0.25;         % P(X = 1) in the running example
end

function s = surprise (p)
  % The surprise of an outcome of probability p, in bits.
  s = log2 (1 / p);
end

function H = entropy (probs)
  % H = sum of p log2(1/p), with 0 log 1/0 taken as 0.
  p = probs(probs > 0);
  H = py_sum (p .* log2 (1 ./ p));
end

function r = h (p)
  % The binary entropy function.
  r = entropy ([p, 1 - p]);
end

function r = h_slope (p)
  % dh/dp = log2((1 - p) / p).
  r = log2 ((1 - p) / p);
end

function total = relative_entropy (P, Q)
  % D(P || Q) in bits, infinite when Q gives zero to an outcome P does not.
  total = 0.0;
  for i = 1:numel (P)
    p = P(i);
    q = Q(i);
    if (p > 0)
      if (q == 0)
        total = Inf;
        return;
      end
      total += p * log2 (p / q);       % a plain running sum, as in the Python
    end
  end
end

function J = bsc_joint (q, p)
  % Joint distribution of (X, Y) when P(X = 1) = q and X crosses BSC(p).
  % J(x + 1, y + 1) is P(X = x, Y = y).
  J = [(1 - q) * (1 - p), (1 - q) * p;
       q * p,             q * (1 - p)];
end

function r = bsc_output_one (q, p)
  % P(Y = 1) for the binary symmetric channel, r = q(1 - p) + (1 - q)p.
  r = q * (1 - p) + (1 - q) * p;
end

function r = bsc_mi (q, p)
  % I(X; Y) = h(r) - h(p) for the binary symmetric channel.
  r = h (bsc_output_one (q, p)) - h (p);
end

function r = z_mi (q, p)
  % I(X; Y) = h(q(1 - p)) - q h(p) for the Z-channel, where only a 1 can flip, to 0.
  r = h (q * (1 - p)) - q * h (p);
end

function r = bec_mi (q, eps_)
  % I(X; Y) = (1 - eps) h(q) for the binary erasure channel.
  % The argument is eps_ because eps is Octave's machine epsilon.
  r = (1 - eps_) * h (q);
end

function r = mi_from_joint (J)
  % I = sum of p(x, y) log2(p(x, y) / (p(x) p(y))), straight from a joint table.
  px = sum (J, 2);          % px(x + 1) is P(X = x), the row sums
  py = sum (J, 1);          % py(y + 1) is P(Y = y), the column sums
  terms = [];
  for x = 1:rows (J)
    for y = 1:columns (J)
      v = J(x, y);
      if (v > 0)
        terms(end + 1) = v * log2 (v / (px(x) * py(y)));
      end
    end
  end
  r = py_sum (terms);
end

function counts = dice_sum_counts ()
  % How many of the 36 rolls of two dice give each sum. counts(s - 1) is the number
  % of rolls that sum to s, for s = 2 to 12, so the entries run in the order the
  % Python's Counter first sees the sums.
  [A, B] = ndgrid (1:6, 1:6);
  S = A(:) + B(:);
  counts = zeros (1, 11);
  for s = 2:12
    counts(s - 1) = sum (S == s);
  end
end

function d = LOADED_DIE ()
  d = [0.5, 0.1, 0.1, 0.1, 0.1, 0.1];
end

function d = DYADIC_DIE ()
  d = [1/2, 1/4, 1/8, 1/16, 1/32, 1/32];
end

function d = LOPSIDED_DIE ()
  d = [0.3, 0.2, 0.2, 0.1, 0.1, 0.1];
end


function main ()
  printf ('Measuring Surprise\n');
  check ('fair coin heads', surprise (1 / 2), '1');
  check ('one of four equally likely outcomes', surprise (1 / 4), '2');
  check ('one face of a fair die', surprise (1 / 6), '2.585');
  check ('an ace, probability 1/13', surprise (1 / 13), '3.700');
  check ('a heart, probability 1/4', surprise (1 / 4), '2');
  check ('the ace of hearts, probability 1/52', surprise (1 / 52), '5.700');
  check ('rank plus suit', surprise (1 / 13) + surprise (1 / 4), '5.700');
  check ('biased coin, the common outcome 0.9', surprise (0.9), '0.152');
  check ('biased coin, the rare outcome 0.1', surprise (0.1), '3.322');
  check ('two dice both showing six, 1/36', surprise (1 / 36), '5.170');
  check ('die face, rounded, doubled', 2 * rnd (surprise (1 / 6), 3), '5.170');
  check ('bits per nat, log2 e', log2 (e), '1.4427');
  require ('surprise is additive: s(1/52) = s(1/13) + s(1/4)', ...
           abs (surprise (1 / 52) - surprise (1 / 13) - surprise (1 / 4)) < 1e-12);

  printf ('\nEntropy as Average Surprise\n');
  check ('biased coin: 0.9 times its surprise', 0.9 * rnd (surprise (0.9), 3), '0.1368');
  check ('biased coin: 0.1 times its surprise', 0.1 * rnd (surprise (0.1), 3), '0.3322');
  check ('biased coin entropy', entropy ([0.9, 0.1]), '0.4690');
  check ('four-sided source 1/2, 1/4, 1/8, 1/8', entropy ([0.5, 0.25, 0.125, 0.125]), '1.75');
  check ('fair die', entropy (repmat (1 / 6, 1, 6)), '2.5850');
  check ('loaded die', entropy (LOADED_DIE ()), '2.1610');
  check ('loaded die, rounded surprise 3.3219', 0.5 + 5 * 0.1 * 3.3219, '2.1610');
  check ('figure: fair coin', entropy ([0.5, 0.5]), '1');
  check ('figure: biased coin', entropy ([0.9, 0.1]), '0.469');
  check ('figure: fair die', entropy (repmat (1 / 6, 1, 6)), '2.585');
  check ('figure: loaded die', entropy (LOADED_DIE ()), '2.161');
  for p = [1e-2, 1e-4, 1e-6]
    require (sprintf ('p log2(1/p) at p = %g is below 2 sqrt(p) / ln 2', p), ...
             p * log2 (1 / p) < 2 * p ^ 0.5 / log (2));
  end
  check ('p log2(1/p) at p = 0.001', 0.001 * log2 (1000), '0.00997');

  printf ('\nThe Binary Entropy Function\n');
  check ('h(0.01)', h (0.01), '0.0808');
  check ('h(0.01) / 0.01, about eight', h (0.01) / 0.01, '8.08');
  check ('h(0.1)', h (0.1), '0.4690');
  check ('h(0.25)', h (0.25), '0.8113');
  check ('h(0.5)', h (0.5), '1');
  check ('slope at p = 0.01, log2 99', h_slope (0.01), '6.629');
  check ('slope at p = 0.1, log2 9', h_slope (0.1), '3.170');
  check ('slope at p = 0.25, log2 3', h_slope (0.25), '1.585');
  ok = true;
  for k = 1:99
    ok = ok && abs (h (k / 100) - h (1 - k / 100)) < 1e-12;
  end
  require ('h(p) = h(1 - p) on a grid', ok);
  hk = zeros (1, 999);                   % hk(k) is h(k / 1000)
  for k = 1:999
    hk(k) = h (k / 1000);
  end
  require ('h rises below 1/2 and falls above it', ...
           all (hk(1:499) < hk(2:500)) && all (hk(500:998) > hk(501:999)));
  step = 1e-6;
  require ('finite-difference slope matches log2((1-p)/p) at p = 0.2', ...
           abs ((h (0.2 + step) - h (0.2 - step)) / (2 * step) - h_slope (0.2)) < 1e-6);
  check ('1 - h(0.1), chapter 1''s capacity', 1 - h (0.1), '0.531');

  printf ('\nConvexity and Jensen''s Inequality\n');
  check ('X in {1, 4}: E[log2 X]', (log2 (1) + log2 (4)) / 2, '1');
  check ('X in {1, 4}: log2 E[X]', log2 (2.5), '1.3219');
  check ('X in {1, 4}: E[X^2]', (1 + 16) / 2, '8.5');
  check ('X in {1, 4}: (E[X])^2', 2.5 ^ 2, '6.25');
  check ('ln 2', log (2), '0.6931');
  check ('ln 0.5', log (0.5), '-0.6931');
  k = 1:1000;
  require ('ln x <= x - 1 on a grid from 0.01 to 10', all (log (k / 100) <= k / 100 - 1 + 1e-15));
  xs = [0.5, 2.0, 7.0];
  ws = [0.2, 0.5, 0.3];
  require ('Jensen for log2 on three points', py_sum (ws .* log2 (xs)) <= log2 (py_sum (ws .* xs)));

  printf ('\nRelative Entropy and Gibbs'' Inequality\n');
  check ('log2(0.9 / 0.5)', log2 (0.9 / 0.5), '0.8480');
  check ('0.9 log2(0.9 / 0.5)', 0.9 * log2 (0.9 / 0.5), '0.7632');
  check ('0.1 log2(0.1 / 0.5)', 0.1 * log2 (0.1 / 0.5), '-0.2322');
  check ('D(biased || fair)', relative_entropy ([0.9, 0.1], [0.5, 0.5]), '0.5310');
  check ('cross-entropy of the biased coin under the fair model', 0.9 * 1 + 0.1 * 1, '1');
  check ('that minus the true entropy', 1 - entropy ([0.9, 0.1]), '0.5310');
  check ('0.5 log2(0.5 / 0.9)', 0.5 * log2 (0.5 / 0.9), '-0.4240');
  check ('0.5 log2(0.5 / 0.1)', 0.5 * log2 (0.5 / 0.1), '1.1610');
  check ('D(fair || biased)', relative_entropy ([0.5, 0.5], [0.9, 0.1]), '0.7370');
  require ('D is infinite when Q forbids an outcome P allows', ...
           relative_entropy ([0.5, 0.5], [1.0, 0.0]) == Inf);

  printf ('\nHow Large Entropy Can Be\n');
  check ('log2 6 minus the loaded die''s entropy', log2 (6) - entropy (LOADED_DIE ()), '0.4240');
  check ('D(loaded || uniform)', relative_entropy (LOADED_DIE (), repmat (1 / 6, 1, 6)), '0.4240');
  require ('D(P || U) = log2|X| - H(P) for the loaded die', ...
           abs (relative_entropy (LOADED_DIE (), repmat (1 / 6, 1, 6)) - (log2 (6) - entropy (LOADED_DIE ()))) < 1e-12);

  printf ('\nJoint and Conditional Entropy (running example, q = 0.25, BSC(0.1))\n');
  q = RUNNING_Q ();
  p = P_FLIP ();
  J = bsc_joint (q, p);
  check ('p(0, 0)', J(1, 1), '0.675');
  check ('p(0, 1)', J(1, 2), '0.075');
  check ('p(1, 0)', J(2, 1), '0.025');
  check ('p(1, 1)', J(2, 2), '0.225');
  check ('P(X = 0)', J(1, 1) + J(1, 2), '0.75');
  check ('P(Y = 0)', J(1, 1) + J(2, 1), '0.7');
  check ('P(Y = 1)', J(1, 2) + J(2, 2), '0.3');
  terms = J .* log2 (1 ./ J);            % v log2(1/v) for each entry of the table
  check ('term 0.675 log2(1/0.675)', terms(1, 1), '0.3828');
  check ('term 0.075 log2(1/0.075)', terms(1, 2), '0.2803');
  check ('term 0.025 log2(1/0.025)', terms(2, 1), '0.1330');
  check ('term 0.225 log2(1/0.225)', terms(2, 2), '0.4842');
  check ('sum of the four rounded terms', py_sum (arrayfun (@(t) rnd (t, 4), row_major (terms))), '1.2803');
  HXY = entropy (row_major (J));
  check ('H(X, Y)', HXY, '1.2803');
  HX = h (q);
  check ('0.75 log2(4/3)', 0.75 * log2 (4 / 3), '0.3113');
  check ('H(X) = h(0.25)', HX, '0.8113');
  HYgX = (1 - q) * h (p) + q * h (p);
  check ('H(Y | X = 0) = h(0.1)', h (p), '0.4690');
  check ('H(Y | X)', HYgX, '0.4690');
  check ('chain rule: H(X) + H(Y | X)', HX + HYgX, '1.2803');
  check ('chain rule, rounded parts', rnd (HX, 4) + rnd (HYgX, 4), '1.2803');
  HY = h (0.3);
  check ('0.3 log2(10/3)', 0.3 * log2 (10 / 3), '0.5211');
  check ('0.7 log2(10/7)', 0.7 * log2 (10 / 7), '0.3602');
  check ('H(Y) = h(0.3)', HY, '0.8813');
  HXgY = HXY - HY;
  check ('H(X | Y) = H(X, Y) - H(Y)', HXgY, '0.3990');
  check ('H(X | Y), rounded parts', rnd (HXY, 4) - rnd (HY, 4), '0.3990');
  post0 = J(2, 1) / 0.7;
  post1 = J(2, 2) / 0.3;
  check ('P(X = 1 | Y = 0) = 0.025 / 0.7', post0, '0.0357');
  check ('P(X = 1 | Y = 1) = 0.225 / 0.3', post1, '0.75');
  check ('h(0.0357)', h (post0), '0.2223');
  check ('h(0.75)', h (post1), '0.8113');
  check ('H(X | Y) = 0.7 h(1/28) + 0.3 h(0.75)', 0.7 * h (post0) + 0.3 * h (post1), '0.3990');
  check ('same from the rounded entropies', 0.7 * 0.2223 + 0.3 * 0.8113, '0.3990');
  check ('chain rule the other way: H(Y) + H(X | Y)', HY + HXgY, '1.2803');
  require ('H(X | Y = 1) equals H(X) exactly here', abs (h (post1) - HX) < 1e-12);
  % Three-variable chain rule on a small check: two fair bits and their XOR. Each of
  % the four triples (a, b, a xor b) has probability 0.25.
  xor3 = [0.25, 0.25, 0.25, 0.25];
  check ('H(X, Y, Z) for two fair bits and their XOR', entropy (xor3), '2');

  printf ('\nConditioning Reduces Entropy\n');
  require ('H(X | Y) <= H(X) in the running example', HXgY <= HX);
  qa = 0.05;
  ya1 = bsc_output_one (qa, p);
  pa1 = qa * (1 - p) / ya1;
  check ('alarm: P(X = 1)', qa, '0.05');
  check ('alarm: H(X) = h(0.05)', h (qa), '0.2864');
  check ('alarm: P(Y = 1)', ya1, '0.14');
  check ('alarm: P(X = 1 | Y = 1) = 0.045 / 0.14', pa1, '0.3214');
  check ('alarm: H(X | Y = 1)', h (pa1), '0.9059');
  require ('alarm: a received 1 raises the uncertainty', h (pa1) > h (qa));

  printf ('\nMutual Information\n');
  I_XY = HX - HXgY;
  check ('I(X; Y) = H(X) - H(X | Y)', I_XY, '0.4123');
  check ('rounded: 0.8113 - 0.3990', rnd (HX, 4) - rnd (HXgY, 4), '0.4123');
  check ('I(X; Y) = H(Y) - H(Y | X)', HY - HYgX, '0.4123');
  check ('rounded: 0.8813 - 0.4690', rnd (HY, 4) - rnd (HYgX, 4), '0.4123');
  check ('I(X; Y) = H(X) + H(Y) - H(X, Y)', HX + HY - HXY, '0.4123');
  check ('rounded: 0.8113 + 0.8813 - 1.2803', rnd (HX, 4) + rnd (HY, 4) - rnd (HXY, 4), '0.4123');
  check ('I from the joint table directly', mi_from_joint (J), '0.4123');
  check ('share of H(X) that crosses, I / H(X)', I_XY / HX, '0.508');
  % The independent pair: P(x, y) = P(x) P(y), with P(X) = (0.75, 0.25) and P(Y) = (0.6, 0.4).
  indep = [0.75 * 0.6, 0.75 * 0.4;
           0.25 * 0.6, 0.25 * 0.4];
  check ('independent pair: I', mi_from_joint (indep), '0');
  % Y = X puts all the probability on the diagonal.
  same_joint = [0.75, 0;
                0,    0.25];
  check ('Y = X: I(X; X) = H(X)', mi_from_joint (same_joint), '0.8113');

  printf ('\nThe Information Diagram\n');
  check ('left crescent H(X | Y)', HXgY, '0.3990');
  check ('overlap I(X; Y)', I_XY, '0.4123');
  check ('right crescent H(Y | X)', HYgX, '0.4690');
  check ('union H(X, Y)', HXgY + I_XY + HYgX, '1.2803');
  check ('XOR: I(X; Y), X and Y independent', 0, '0');
  % I(X; Y | Z) = H(X | Z) - H(X | Y, Z) for the XOR triple.
  hx_z = entropy (xor3) - 1;             % H(X, Z) = 2 and H(Z) = 1
  check ('XOR: H(X | Z)', hx_z, '1');
  check ('XOR: H(X | Y, Z)', 0, '0');
  check ('XOR: I(X; Y | Z)', hx_z - 0, '1');
  check ('XOR: the centre region I(X; Y) - I(X; Y | Z)', 0 - (hx_z - 0), '-1');

  printf ('\nMutual Information of the Binary Symmetric Channel\n');
  ok = true;
  for a = 0:20
    for b = 0:20
      ok = ok && abs (bsc_output_one (a / 20, b / 20) - 0.5 - (1 - 2 * b / 20) * (a / 20 - 0.5)) < 1e-12;
    end
  end
  require ('r - 1/2 = (1 - 2p)(q - 1/2) on a grid', ok);
  check ('uniform input: 1 - h(0.1)', bsc_mi (0.5, p), '0.5310');
  check ('uniform input: P(Y = 1)', bsc_output_one (0.5, p), '0.5');
  % Each row: q, then the stated P(Y = 1), H(Y) and I.
  tbl = {0.5,  '0.5',  '1',      '0.5310';
         0.4,  '0.42', '0.9815', '0.5125';
         0.3,  '0.34', '0.9248', '0.4558';
         0.25, '0.3',  '0.8813', '0.4123';
         0.1,  '0.18', '0.6801', '0.2111';
         0.05, '0.14', '0.5842', '0.1152'};
  for i = 1:rows (tbl)
    qq = tbl{i, 1};
    r_s = tbl{i, 2};
    hy_s = tbl{i, 3};
    i_s = tbl{i, 4};
    r = bsc_output_one (qq, p);
    check (sprintf ('table q = %s: P(Y = 1)', py_str (qq)), r, r_s);
    check (sprintf ('table q = %s: H(Y)', py_str (qq)), h (r), hy_s);
    check (sprintf ('table q = %s: I', py_str (qq)), bsc_mi (qq, p), i_s);
    check (sprintf ('table q = %s: rounded H(Y) - 0.4690', py_str (qq)), rnd (h (r), 4) - 0.4690, i_s);
  end
  mi_grid = zeros (1, 999);              % mi_grid(k) is the information at q = k / 1000
  for k = 1:999
    mi_grid(k) = bsc_mi (k / 1000, p);
  end
  [~, best_k] = max (mi_grid);
  require ('uniform input is best on a fine grid for BSC(0.1)', best_k == 500);
  ok = true;
  for k = 1:99
    ok = ok && abs (bsc_mi (0.5, k / 100) - bsc_mi (0.5, 1 - k / 100)) < 1e-12;
  end
  require ('uniform input: I(1 - p) = I(p)', ok);
  check ('BSC(0.5): I', bsc_mi (0.5, 0.5), '0');
  check ('figure: q = 0.25 at p = 0 carries all of H(X)', h (0.25), '0.811');
  check ('BEC(0.1), uniform input', bec_mi (0.5, 0.1), '0.9');
  check ('Z(0.1), uniform input: h(0.45)', h (0.45), '0.9928');
  check ('Z(0.1), uniform input: 0.5 h(0.1)', 0.5 * h (0.1), '0.2345');
  check ('Z(0.1), uniform input: I', z_mi (0.5, p), '0.7583');

  printf ('\nLooking Ahead\n');
  check ('four-sided source: average questions', 0.5 * 1 + 0.25 * 2 + 0.125 * 3 + 0.125 * 3, '1.75');

  printf ('\nWorked Examples\n');
  counts = dice_sum_counts ();
  probs = counts / 36;                   % P(S = s) for s = 2 to 12
  check ('dice sum: log2 36', log2 (36), '5.1699');
  check ('dice sum: log2 18', log2 (18), '4.1699');
  check ('dice sum: log2 12', log2 (12), '3.5850');
  check ('dice sum: log2 9', log2 (9), '3.1699');
  check ('dice sum: log2 7.2', log2 (7.2), '2.8480');
  check ('dice sum: log2 6', log2 (6), '2.5850');
  HS = entropy (probs);
  check ('dice sum: H(S)', HS, '3.2744');
  check ('dice sum: from the rounded surprises', ...
         py_sum (arrayfun (@(c) c / 36 * rnd (log2 (36 / c), 4), counts)), '3.2744');
  check ('dice sum: H(D1, D2)', log2 (36), '5.1699');
  check ('dice sum: H(D1, D2 | S)', log2 (36) - HS, '1.8955');
  check ('dice sum: H(S | D1)', log2 (6), '2.5850');
  check ('dice sum: I(D1; S)', HS - log2 (6), '0.6894');
  check ('dice sum: rounded I', 3.2744 - 2.5850, '0.6894');
  check ('dice sum: rounded H(D1, D2 | S)', 5.1699 - 3.2744, '1.8955');
  check ('dice sum: P(S = 7)', counts(7 - 1) / 36, '0.1667');

  ia = bsc_mi (qa, p);
  check ('alarm: h(0.14)', h (ya1), '0.5842');
  check ('alarm: I = h(0.14) - h(0.1)', ia, '0.1152');
  check ('alarm: share of H(X) that crosses', ia / h (qa), '0.40');
  pa0 = qa * p / (1 - ya1);
  check ('alarm: P(X = 1 | Y = 0) = 0.005 / 0.86', pa0, '0.0058');
  check ('alarm: h(0.0058...)', h (pa0), '0.05154');
  check ('alarm: h(0.3214...)', h (pa1), '0.90593');
  check ('alarm: H(X | Y) from five-place entropies', 0.86 * 0.05154 + 0.14 * 0.90593, '0.1712');
  check ('alarm: H(X | Y) exactly', (1 - ya1) * h (pa0) + ya1 * h (pa1), '0.1712');
  check ('alarm: H(X) - H(X | Y)', h (qa) - ((1 - ya1) * h (pa0) + ya1 * h (pa1)), '0.1152');
  check ('alarm: rounded 0.2864 - 0.1712', 0.2864 - 0.1712, '0.1152');

  check ('Z: uniform I', z_mi (0.5, p), '0.7583');
  check ('Z: q = 0.4', z_mi (0.4, p), '0.7551');
  z_grid = zeros (1, 999);               % z_grid(k) is the information at q = k / 1000
  for k = 1:999
    z_grid(k) = z_mi (k / 1000, p);
  end
  [~, best_k] = max (z_grid);
  best_q = best_k / 1000;
  check ('Z: best q on a grid of step 0.001', best_q, '0.456');
  check ('Z: I at that q', z_mi (best_q, p), '0.7628');
  check ('Z: P(Y = 1) at q = 0.456', 0.456 * 0.9, '0.4104');
  check ('Z: gain over the uniform input', z_mi (best_q, p) - z_mi (0.5, p), '0.0046');
  require ('Z: the gain is under 0.005 bits', z_mi (best_q, p) - z_mi (0.5, p) < 0.005);
  % The optimum from setting the derivative to zero: (1 - p) log2((1 - r)/r) = h(p), r = (1 - p) q.
  slope_target = h (p) / (1 - p);
  check ('Z optimum: h(0.1) / 0.9', slope_target, '0.5211');
  check ('Z optimum: 2 to that power', 2 ^ slope_target, '1.4351');
  r_opt = 1 / (1 + 2 ^ slope_target);
  check ('Z optimum: 1/r = 1 + 1.4351', 1 / r_opt, '2.4351');
  check ('Z optimum: r = 1 / (1 + 1.4351)', r_opt, '0.4107');
  check ('Z optimum: same from the rounded 1.4351', 1 / (1 + 1.4351), '0.4107');
  check ('Z optimum: q = r / 0.9', r_opt / (1 - p), '0.4563');
  check ('Z optimum: same from the rounded 0.4107', 0.4107 / 0.9, '0.4563');
  check ('Z optimum: I there', z_mi (r_opt / (1 - p), p), '0.7628');
  require ('Z optimum: the derivative vanishes there', ...
           abs ((1 - p) * log2 ((1 - r_opt) / r_opt) - h (p)) < 1e-12);
  require ('Z optimum: the grid''s best q is the optimum to three places', rnd (r_opt / (1 - p), 3) == best_q);
  check ('BSC(0.1) uniform, for comparison', bsc_mi (0.5, p), '0.5310');

  check ('dyadic die entropy', entropy (DYADIC_DIE ()), '1.9375');
  require ('dyadic die surprises are 1, 2, 3, 4, 5, 5', ...
           isequal (arrayfun (@(x) rnd (surprise (x), 12), DYADIC_DIE ()), [1, 2, 3, 4, 5, 5]));
  check ('dyadic die: average questions', py_sum (DYADIC_DIE () .* arrayfun (@surprise, DYADIC_DIE ())), '1.9375');
  check ('dyadic die: log2 6 - H', log2 (6) - entropy (DYADIC_DIE ()), '0.6475');
  check ('dyadic die: D(P || U)', relative_entropy (DYADIC_DIE (), repmat (1 / 6, 1, 6)), '0.6475');

  printf ('\nExercises\n');
  check ('face card, probability 3/13: log2(13/3)', surprise (3 / 13), '2.115');
  check ('lopsided die entropy', entropy (LOPSIDED_DIE ()), '2.4464');
  check ('lopsided die: below log2 6 by', log2 (6) - entropy (LOPSIDED_DIE ()), '0.1385');
  check ('lopsided die: 0.3 log2(10/3), five places', 0.3 * log2 (10 / 3), '0.52109');
  check ('lopsided die: 0.2 log2 5, five places', 0.2 * log2 (5), '0.46439');
  check ('lopsided die: 0.1 log2 10, five places', 0.1 * log2 (10), '0.33219');
  check ('lopsided die: rounded sum 0.52109 + 2(0.46439) + 3(0.33219)', 0.52109 + 2 * 0.46439 + 3 * 0.33219, ...
         '2.44644');
  check ('lopsided die entropy, five places', entropy (LOPSIDED_DIE ()), '2.44644');
  check ('log2 6, five places', log2 (6), '2.58496');
  check ('lopsided die: gap, five places', log2 (6) - entropy (LOPSIDED_DIE ()), '0.13852');
  check ('lopsided die: rounded 2.58496 - 2.44644', 2.58496 - 2.44644, '0.13852');
  check ('h(0.2): log2 5', log2 (5), '2.3219');
  check ('h(0.2): log2 1.25', log2 (1.25), '0.3219');
  check ('h(0.2)', h (0.2), '0.7219');
  check ('h(0.2) from the rounded logarithms', 0.2 * 2.3219 + 0.8 * 0.3219, '0.7219');
  check ('X in {1, 5}: E[X^2]', (1 + 25) / 2, '13');
  check ('X in {1, 5}: (E[X])^2', 3 ^ 2, '9');
  check ('X in {1, 5}: variance', (1 + 25) / 2 - 9, '4');
  check ('D(fair || biased) again', relative_entropy ([0.5, 0.5], [0.9, 0.1]), '0.7370');
  check ('D(uniform || loaded die)', relative_entropy (repmat (1 / 6, 1, 6), LOADED_DIE ()), '0.3500');
  check ('BSC(0.01), uniform input: 1 - h(0.01)', bsc_mi (0.5, 0.01), '0.9192');
  check ('BEC(0.1), q = 0.25', bec_mi (0.25, 0.1), '0.7302');
  check ('BEC(0.1), q = 0.25 from the rounded H(X)', 0.9 * 0.8113, '0.7302');
  check ('BSC(0.1), q = 0.25', bsc_mi (0.25, 0.1), '0.4123');
  check ('ln 6, the fair die in nats', log (6), '1.7918');
  check ('ln 6 in bits', log (6) * log2 (e), '2.5850');
  check ('rounded: 1.7918 times 1.4427', 1.7918 * 1.4427, '2.5850');

  printf ('\n');
  printf ('All chapter 4 numbers reproduced.\n');
end


main ();
