% Reproduce every number in chapter 31 of Error Correction Coding,
% "Lab: Measuring Miscorrection", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch31_lab_measuring_miscorrection.m
%
% It needs only core Octave 7 or later, with no pkg load and no Forge
% packages. It is the Octave twin of python/ch31_lab_measuring_miscorrection.py
% and prints the same lines, which tools/compare-routes.sh checks. The Python
% script's options, which change the experiment, are not repeated here.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Whole numbers are compared exactly, and every fractional value is
% printed with a fixed number of places, so its comparison is exact too.
%
% The three (72,64) check matrices are built as in "SECDED and Hsiao Codes"
% and compared with the columns of the model that generated the bench,
% tb/secded72/secded72_model.py, as the hardware run records them, read from
% data/ch31/miscorrection_counts.json with jsondecode. Every three-bit and
% four-bit error pattern is classified as one row of a matrix that nchoosek
% builds, 59,640 and 1,028,790 rows for each code, and the counts of weight 5
% and 6 come from a recursion over syndrome counts. The seeded campaign walks
% the draws of Python's random.Random(31) from common/PyRandom.m the way the
% Python script's draw() does, so its counts are the counts the Python script
% and the hardware find, and they are checked against the data file.
%
% A column of a check matrix with r rows is held as a whole number whose bit i
% is the entry in row i. An error pattern is a row of codeword positions,
% numbered from 0 as in Python, so position 0 is the first column.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

% Helpers ----------------------------------------------------------------------

function s = as_text (x)
  % The text Python's str() prints for a str or an int.
  if (ischar (x))
    s = x;
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

function s = pct (x)
  s = sprintf ('%.3f', 100 * x);
end

function s = joined (v)
  % Whole numbers separated by single spaces, as Python's " ".join prints them.
  s = strjoin (arrayfun (@(x) sprintf ('%d', x), v, 'UniformOutput', false), ' ');
end

function v = bc (n, k)
  % The binomial coefficient, 0 when k lies outside 0..n.
  if (k < 0 || k > n)
    v = 0;
  else
    k = min (k, n - k);
    v = 1;
    for t = 1:k
      v = v * (n - k + t) / t;
    end
  end
end

function n = ceil_whole (x)
  % The smallest whole number at least x, forgiving the rounding of a value that is whole on paper.
  n = ceil (x - 1e-9);
end

function v = col (rows)
  v = sum (2 .^ rows);
end

function v = cols_of (C)
  % One column for each row of C, a list of row sets.
  v = sum (2 .^ C, 2)';
end

function w = weight (v)
  w = zeros (size (v));
  while (any (v(:)))
    w = w + bitand (v, 1);
    v = bitshift (v, -1);
  end
end

% The codes --------------------------------------------------------------------

function code = make_code (name)
  N = 72;
  w3 = cols_of (nchoosek (0:7, 3));
  units = 2 .^ (0:7);
  if (strcmp (name, 'hsiao'))
    runs = zeros (1, 8);
    for i = 0:7
      runs(i + 1) = col (mod (i + (0:4), 8));
    end
    cols = [w3, runs, units];
    checks = 64:71;
  elseif (strcmp (name, 'hsiao_dictionary'))
    C5 = nchoosek (0:7, 5);
    cols = [w3, cols_of(C5(1:8, :)), units];
    checks = 64:71;
  else
    cols = bitor (0:N - 1, 128);
    checks = [0, 1, 2, 4, 8, 16, 32, 64];
  end
  code.name = name;
  code.cols = cols;
  code.checks = checks;
  code.isdata = true (1, N);
  code.isdata(checks + 1) = false;
  code.pos_of = -ones (1, 256);
  code.pos_of(cols + 1) = 0:N - 1;
end

function out = outcomes (code, P)
  % One outcome for each row of P: 1 correct, 2 flagged, 3 miscorrected, 4 undetected.
  % The syndrome is the XOR of the flipped columns. A zero syndrome leaves the
  % word alone, a syndrome equal to column j flips bit j, and any other
  % syndrome raises the flag. The data come out wrong when the flipped bits,
  % after the decoder's own flip, include a data position.
  n = rows (P);
  s = zeros (n, 1);
  for t = 1:columns (P)
    s = bitxor (s, reshape (code.cols(P(:, t) + 1), n, 1));
  end
  ndata = sum (reshape (code.isdata(P + 1), size (P)), 2);
  j = reshape (code.pos_of(s + 1), n, 1);
  matched = j >= 0;
  inside = any (P == j, 2);
  jd = zeros (n, 1);
  jd(matched) = code.isdata(j(matched) + 1);
  dw = ndata + jd .* (1 - 2 * inside);
  out = 2 * ones (n, 1);
  out(s == 0 & ndata == 0) = 1;
  out(s == 0 & ndata > 0) = 4;
  out(matched & dw == 0) = 1;
  out(matched & dw > 0) = 3;
end

function c = tally (out)
  % Counts of the four outcomes, in the order correct, flagged, miscorrected, undetected.
  c = accumarray (out(:), 1, [4, 1])';
end

function counts = syndrome_counts (cols, top)
  % counts(w + 1, s + 1), the number of weight-w error patterns whose syndrome is s.
  counts = zeros (top + 1, 256);
  counts(1, 1) = 1;
  for c = cols
    idx = bitxor (0:255, c) + 1;
    for w = top:-1:1
      counts(w + 1, idx) = counts(w + 1, idx) + counts(w, :);
    end
  end
end

function q = check_only (code)
  % q(w + 1), the number of codewords of weight w whose ones all sit in check positions.
  q = zeros (1, 73);
  k = numel (code.checks);
  for subset = 0:2 ^ k - 1
    pick = bitand (bitshift (subset, -(0:k - 1)), 1) == 1;
    s = 0;
    for p = code.checks(pick)
      s = bitxor (s, code.cols(p + 1));
    end
    if (s == 0)
      q(sum (pick) + 1) = q(sum (pick) + 1) + 1;
    end
  end
end

function c = exact_split (code, counts, w)
  % The four outcome counts of weight w, from the syndrome counts and the check-only codewords.
  N = 72;
  q = check_only (code);
  a_w = counts(w + 1, 1);
  m_w = sum (counts(w + 1, code.cols + 1));
  q_up = (w + 1) * q(w + 2);
  q_down = 0;
  if (w >= 1)
    q_down = (N - w + 1) * q(w);
  end
  c = [q(w + 1) + q_up + q_down, bc(N, w) - a_w - m_w, m_w - q_up - q_down, a_w - q(w + 1)];
end

function b = dual_weights (code)
  % b(j + 1), the number of the 256 codewords of the dual code with weight j.
  R = zeros (8, 72);
  S = zeros (256, 8);
  for i = 0:7
    R(i + 1, :) = bitand (bitshift (code.cols, -i), 1);
    S(:, i + 1) = bitand (bitshift ((0:255)', -i), 1);
  end
  wt = sum (mod (S * R, 2), 2);
  b = accumarray (wt + 1, 1, [73, 1])';
end

function a = macwilliams (b, i)
  % A_i of the code from the dual's weights, by the MacWilliams identity.
  N = 72;
  total = 0;
  for x = 0:N
    if (b(x + 1))
      k = 0;
      for t = 0:i
        k = k + (-1) ^ t * bc (x, t) * bc (N - x, i - t);
      end
      total = total + b(x + 1) * k;
    end
  end
  a = total / 256;
end

% Bursts -----------------------------------------------------------------------

function [order, steps] = search_layout (cols)
  % A cell order in which no three neighboring cells miscorrect, by depth-first search.
  N = 72;
  present = false (1, 256);
  present(cols + 1) = true;
  order = zeros (1, N);
  used = false (1, N);
  choice = zeros (1, N + 1);
  depth = 0;
  steps = 0;
  while (true)
    steps = steps + 1;
    if (depth == N)
      return;
    end
    j = choice(depth + 1);
    while (j < N && (used(j + 1) || (depth >= 2 && ...
           present(bitxor (bitxor (cols(order(depth - 1) + 1), cols(order(depth) + 1)), cols(j + 1)) + 1))))
      j = j + 1;
    end
    if (j < N)
      choice(depth + 1) = j + 1;
      depth = depth + 1;
      order(depth) = j;
      used(j + 1) = true;
      if (depth < N)
        choice(depth + 1) = 0;
      end
    else
      if (depth == 0)
        order = [];
        return;
      end
      used(order(depth) + 1) = false;
      depth = depth - 1;
    end
  end
end

function order = layout (kind, code)
  % order(k + 1) is the codeword position stored in cell k of the word.
  if (strcmp (kind, 'natural'))
    order = 0:71;
  elseif (strncmp (kind, 'stride', 6))
    order = mod (str2double (kind(7:end)) * (0:71), 72);
  elseif (strcmp (code.name, 'hamming'))
    % No order flags every window of extended Hamming, whose floor is 46
    % windows in any layout, so no search runs for it and stride 11, which
    % meets that floor, stands in.
    order = layout ('stride11', code);
  else
    order = search_layout (code.cols);
  end
end

function c = adjacent (code, order)
  % Outcome counts of the 70 triples of three neighboring cells.
  c = tally (outcomes (code, [order(1:70)', order(2:71)', order(3:72)']));
end

% The campaign -----------------------------------------------------------------

function P = campaign_patterns (seed, draws, weights)
  % The error patterns of the campaign, one matrix for each weight, drawn as
  % the Python script's draw() draws them: two random() values for the data
  % word, then one value for each position, scaled to 72 and truncated, drawn
  % again whenever it repeats a position already chosen.
  rng = PyRandom (seed);
  v = rng.random (12 * draws);
  idx = 1;
  P = cell (1, numel (weights));
  for t = 1:numel (weights)
    w = weights(t);
    M = zeros (draws, w);
    for d = 1:draws
      idx = idx + 2;
      k = 0;
      pos = -ones (1, w);
      while (k < w)
        if (idx > numel (v))
          v = [v; rng.random(draws)];
        end
        p = floor (v(idx) * 72);
        idx = idx + 1;
        if (~any (pos(1:k) == p))
          k = k + 1;
          pos(k) = p;
        end
      end
      M(d, :) = pos;
    end
    P{t} = M;
  end
end

function [lo, hi] = wilson (m, n)
  % The Wilson score interval for a proportion, m successes in n trials.
  z = 1.96;
  ph = m / n;
  center = ph + z * z / (2 * n);
  half = z * sqrt (ph * (1 - ph) / n + z * z / (4 * n * n));
  scale = 1 + z * z / n;
  lo = (center - half) / scale;
  hi = (center + half) / scale;
end

function v = silent_per_year (share, t_hours)
  % Silent corruptions a year when each word is decoded once per interval of t_hours.
  N = 72;
  mu = 1e-6 * t_hours;
  p3 = mu ^ 3 / 6 * exp (-mu);
  v = 4096 * (8760 / t_hours) * p3 * ((N - 1) * (N - 2) / N ^ 2) * share;
end

function v = counts_of (s)
  % A record's four outcome counts as a row.
  v = [s.correct, s.flagged, s.miscorrected, s.undetected];
end

% The checks -------------------------------------------------------------------

N = 72;
SEED = 31;
DRAWS = 100000;
Z = 1.96;
CHECKPOINTS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
TRIPLES = bc (N, 3);
QUADRUPLES = bc (N, 4);
NAMES = {'hsiao', 'hsiao_dictionary', 'hamming'};
LABELS = {'Hsiao balanced', 'Hsiao dictionary', 'Extended Hamming'};
TRACKED = {'', '', 'miscorrected', 'undetected'};

STATED_EXACT = [33568, 26072, 8392, 1020398; 34032, 25608, 8508, 1020282; 45304, 14336, 11326, 1017464];
STATED_SHARES = {'56.3', '56.284', '0.816'; '57.1', '57.062', '0.827'; '76.0', '75.962', '1.101'};
STATED_WEIGHT_6 = [1216136, 1215217, 1446144];
STATED_WEIGHT_5 = {7867472, '56.230'; 7869846, '56.247'; 9447032, '67.520'};
% For each code and weight 3, then 4: tracked count, flagged, percent, Wilson low, Wilson high.
STATED_CAMPAIGN = {56222, 43778, '56.222', '55.914', '56.529'; 864, 99136, '0.864', '0.808', '0.923';
                   56984, 43016, '56.984', '56.677', '57.291'; 829, 99171, '0.829', '0.775', '0.887';
                   75975, 24025, '75.975', '75.709', '76.239'; 1122, 98878, '1.122', '1.059', '1.189'};
ADJ_KINDS = {'natural', 'stride5', 'stride7', 'stride11', 'stride13', 'stride29', 'searched'};
STATED_ADJACENT = [35, 28, 70; 49, 48, 65; 43, 42, 53; 39, 42, 46; 39, 33, 46; 33, 29, 46; 0, 0, 46];
STATED_EXPECTED_WINDOWS = {'39.40', '39.94', '53.17'};
STATED_SEARCH_STEPS = [5039, 8291];
SCRUB_NAMES = {'1 day', '30 days', '1 year'};
SCRUB_HOURS = [24, 720, 8760];
STATED_SILENT = {'1.86e-09', '1.67e-06', '2.45e-04'; '2.51e-09', '2.26e-06', '3.31e-04'};
HARDWARE_LAYOUTS = {'natural', 'stride11', 'searched'};

lab = fileparts (fileparts (mfilename ('fullpath')));
data_path = fullfile (lab, 'data', 'ch31', 'miscorrection_counts.json');
if (~exist (data_path, 'file'))
  printf ('%s is missing; run tb/miscorrection/run_miscorrection.py to measure it\n', data_path);
  exit (1);
end
data = jsondecode (fileread (data_path));

codes = cellfun (@make_code, NAMES, 'UniformOutput', false);

section ('The three check matrices, against the model that generated the bench');
for c = 1:3
  require (sprintf ('%s: columns equal the model''s that generated the bench', LABELS{c}), ...
           isequal (codes{c}.cols(:), data.codes.(NAMES{c}).columns(:)));
end
for c = 1:2
  require (sprintf ('%s: every column has odd weight', LABELS{c}), all (mod (weight (codes{c}.cols), 2) == 1));
end
require ('Extended Hamming: every column has a 1 in the parity row', all (bitand (codes{3}.cols, 128) == 128));
for c = 1:3
  q = check_only (codes{c});
  same (sprintf ('%s: codewords in the check positions alone', LABELS{c}), sum (q(2:end)), 0);
end

section (sprintf ('Every triple and every quadruple, enumerated (%d and %d of them)', TRIPLES, QUADRUPLES));
T3 = nchoosek (0:N - 1, 3);
T4 = nchoosek (0:N - 1, 4);
exact = cell (3, 4);
for c = 1:3
  t = tally (outcomes (codes{c}, T3));
  f = tally (outcomes (codes{c}, T4));
  exact{c, 3} = t;
  exact{c, 4} = f;
  L = LABELS{c};
  same (sprintf ('%s: triples miscorrected', L), t(3), STATED_EXACT(c, 1));
  same (sprintf ('%s: triples flagged', L), t(2), STATED_EXACT(c, 2));
  same (sprintf ('%s: triples correct or undetected', L), t(1) + t(4), 0);
  same (sprintf ('%s: percent of triples miscorrected', L), sprintf ('%.1f', 100 * t(3) / TRIPLES), STATED_SHARES{c, 1});
  same (sprintf ('%s: the same percent to three places', L), pct (t(3) / TRIPLES), STATED_SHARES{c, 2});
  same (sprintf ('%s: quadruples undetected', L), f(4), STATED_EXACT(c, 3));
  same (sprintf ('%s: quadruples flagged', L), f(2), STATED_EXACT(c, 4));
  same (sprintf ('%s: quadruples correct or miscorrected', L), f(1) + f(3), 0);
  same (sprintf ('%s: percent of quadruples undetected', L), pct (f(4) / QUADRUPLES), STATED_SHARES{c, 3});
end

section ('Weight-4 codewords, the syndrome-count recursion and the MacWilliams identity');
for c = 1:3
  code = codes{c};
  L = LABELS{c};
  counts = syndrome_counts (code.cols, 6);
  b = dual_weights (code);
  a4 = counts(5, 1);
  a6 = counts(7, 1);
  same (sprintf ('%s: A_4 from the recursion', L), a4, STATED_EXACT(c, 3));
  same (sprintf ('%s: A_4 from the MacWilliams identity', L), macwilliams (b, 4), a4);
  require (sprintf ('%s: triples miscorrected equal 4 A_4', L), exact{c, 3}(3) == 4 * a4);
  require (sprintf ('%s: recursion split equals enumeration, weights 3, 4', L), ...
           isequal (exact_split (code, counts, 3), exact{c, 3}) && isequal (exact_split (code, counts, 4), exact{c, 4}));
  same (sprintf ('%s: A_6 from the recursion', L), a6, STATED_WEIGHT_6(c));
  same (sprintf ('%s: A_6 from the MacWilliams identity', L), macwilliams (b, 6), a6);
  five = exact_split (code, counts, 5);
  same (sprintf ('%s: weight-5 patterns miscorrected', L), five(3), STATED_WEIGHT_5{c, 1});
  require (sprintf ('%s: weight 5 miscorrected = 68 A_4 + 6 A_6', L), five(3) == 68 * a4 + 6 * a6);
  same (sprintf ('%s: weight 5 correct or undetected', L), five(1) + five(4), 0);
  same (sprintf ('%s: percent of weight 5 miscorrected', L), pct (five(3) / bc (N, 5)), STATED_WEIGHT_5{c, 2});
end

section (sprintf ('The seeded campaign: seed %d, %d draws of weight 3 and of weight 4', SEED, DRAWS));
P = campaign_patterns (SEED, DRAWS, [3, 4]);
camp = cell (3, 4);
running = cell (3, 4);
outs3 = cell (1, 3);
for c = 1:3
  for w = 3:4
    o = outcomes (codes{c}, P{w - 2});
    camp{c, w} = tally (o);
    hits = cumsum (o == w);
    running{c, w} = hits(CHECKPOINTS)';
    if (w == 3)
      outs3{c} = o;
    end
  end
end
for c = 1:3
  L = LABELS{c};
  for w = 3:4
    k = camp{c, w}(w);
    row = 2 * (c - 1) + w - 2;
    [lo, hi] = wilson (k, DRAWS);
    share = exact{c, w}(w) / bc (N, w);
    same (sprintf ('%s, weight %d: %s', L, w, TRACKED{w}), k, STATED_CAMPAIGN{row, 1});
    same (sprintf ('%s, weight %d: flagged', L, w), camp{c, w}(2), STATED_CAMPAIGN{row, 2});
    same (sprintf ('%s, weight %d: other outcomes', L, w), DRAWS - k - camp{c, w}(2), 0);
    same (sprintf ('%s, weight %d: percent %s', L, w, TRACKED{w}), pct (k / DRAWS), STATED_CAMPAIGN{row, 3});
    same (sprintf ('%s, weight %d: Wilson interval', L, w), [pct(lo) ' to ' pct(hi)], ...
          [STATED_CAMPAIGN{row, 4} ' to ' STATED_CAMPAIGN{row, 5}]);
    require (sprintf ('%s, weight %d: interval holds the exact share', L, w), lo <= share && share <= hi);
  end
end
same ('chance that six 95% intervals all hold their shares', sprintf ('%.3f', 0.95 ^ 6), '0.735');

section ('The hardware campaign recorded in data/ch31/miscorrection_counts.json');
same ('seed of the hardware campaign', data.seed, SEED);
same ('draws of each weight', data.draws, DRAWS);
same ('weights drawn', joined (data.weights'), '3 4');
for c = 1:3
  record = data.codes.(NAMES{c});
  L = LABELS{c};
  require (sprintf ('%s: hardware counts equal this program''s', L), ...
           isequal (counts_of (record.campaign.weight_3), camp{c, 3}) && isequal (counts_of (record.campaign.weight_4), camp{c, 4}));
  require (sprintf ('%s: hardware running counts equal this program''s', L), ...
           isequal (record.running.weight_3(:)', running{c, 3}) && isequal (record.running.weight_4(:)', running{c, 4}));
  same (sprintf ('%s: cases where the hardware differs from the model', L), record.model_mismatches, 0);
  same (sprintf ('%s: cases where it differs from the prediction', L), record.prediction_mismatches, 0);
  ok = true;
  for h = 1:numel (HARDWARE_LAYOUTS)
    kind = HARDWARE_LAYOUTS{h};
    ok = ok && isequal (counts_of (record.adjacent.(kind)), adjacent (codes{c}, layout (kind, codes{c})));
  end
  require (sprintf ('%s: hardware adjacent windows equal this program''s', L), ok);
end
only_h = sum (outs3{1} == 3 & outs3{2} ~= 3);
only_d = sum (outs3{2} == 3 & outs3{1} ~= 3);
same ('weight-3 draws miscorrected by the balanced matrix alone', data.paired.hsiao_only, only_h);
same ('weight-3 draws miscorrected by the dictionary matrix alone', data.paired.dictionary_only, only_d);

section ('Reading a rate');
p = STATED_EXACT(1, 1) / TRIPLES;
same ('standard error at the balanced share, 100000 draws, points', pct (sqrt (p * (1 - p) / DRAWS)), '0.157');
[lo, hi] = wilson (STATED_CAMPAIGN{1, 1}, DRAWS);
same ('half-width of the balanced Wilson interval, points', pct ((hi - lo) / 2), '0.307');
same ('draws for a half-width of 1 point at the balanced share', ceil_whole (Z * Z * p * (1 - p) / (0.01 * 0.01)), 9453);
same ('draws for a half-width of 1 point at p = 1/2', ceil_whole (Z * Z * 0.5 * (1 - 0.5) / (0.01 * 0.01)), 9604);
same ('draws for 1 point with 95% certainty, Chebyshev', ceil_whole (1 / (4 * 0.05 * 0.01 * 0.01)), 50000);
same ('draws for 1 point with 95% certainty, Hoeffding', ceil_whole (log (2 / 0.05) / (2 * 0.01 * 0.01)), 18445);
same ('Chebyshev half-width at 100000 draws, points', pct (sqrt (1 / (4 * 0.05 * DRAWS))), '0.707');
same ('Hoeffding half-width at 100000 draws, points', pct (sqrt (log (2 / 0.05) / (2 * DRAWS))), '0.429');
same ('normal half-width at p = 1/2 and 100000 draws, points', pct (Z * sqrt (0.25 / DRAWS)), '0.310');
o_h = outcomes (codes{1}, T3);
o_d = outcomes (codes{2}, T3);
ex_h = sum (o_h == 3 & o_d ~= 3);
ex_d = sum (o_d == 3 & o_h ~= 3);
p1 = STATED_EXACT(1, 1) / TRIPLES;
p2 = STATED_EXACT(2, 1) / TRIPLES;
diff = p2 - p1;
same ('triples miscorrected by the balanced matrix alone', ex_h, 5583);
same ('triples miscorrected by the dictionary matrix alone', ex_d, 6047);
same ('difference of the two shares, points', pct (diff), '0.778');
same ('draws to resolve it, separate patterns for each matrix', ...
      ceil_whole (Z * Z * (p1 * (1 - p1) + p2 * (1 - p2)) / diff ^ 2), 31167);
same ('draws to resolve it, the same patterns for both', ...
      ceil_whole (Z * Z * ((ex_h + ex_d) / TRIPLES - diff ^ 2) / diff ^ 2), 12373);
dh = (only_d - only_h) / DRAWS;
half = Z * sqrt (((only_d + only_h) / DRAWS - dh * dh) / DRAWS);
same ('campaign''s paired difference, points', pct (dh), '0.762');
same ('its 95% interval, points', [pct(dh - half) ' to ' pct(dh + half)], '0.488 to 1.036');

section ('Triples of adjacent cells, 70 windows in a 72-cell word');
for a = 1:numel (ADJ_KINDS)
  kind = ADJ_KINDS{a};
  got = zeros (1, 3);
  for c = 1:3
    counts = adjacent (codes{c}, layout (kind, codes{c}));
    got(c) = counts(3);
  end
  same (sprintf ('miscorrected windows, %s layout, three matrices', kind), joined (got), joined (STATED_ADJACENT(a, :)));
end
for c = 1:2
  [~, steps] = search_layout (codes{c}.cols);
  same (sprintf ('%s: steps the layout search took', LABELS{c}), steps, STATED_SEARCH_STEPS(c));
end
require ('Extended Hamming: no search runs, its searched layout repeats stride 11', ...
         isequal (layout ('searched', codes{3}), layout ('stride11', codes{3})));
require ('Extended Hamming: three low positions always miscorrect', all (outcomes (codes{3}, nchoosek (0:63, 3)) == 3));
ok = true;
for a = 1:6
  counts = adjacent (codes{3}, layout (ADJ_KINDS{a}, codes{3}));
  ok = ok && counts(2) <= 24;
end
require ('Extended Hamming: no layout tried flags more than 24 windows', ok);
for c = 1:3
  same (sprintf ('%s: expected miscorrected windows, random layout', LABELS{c}), ...
        sprintf ('%.2f', 70 * STATED_EXACT(c, 1) / TRIPLES), STATED_EXPECTED_WINDOWS{c});
end

section (sprintf ('A designer''s number: %d words at %g upsets per word-hour', 4096, 1e-6));
same ('chance three strikes land on three different cells', sprintf ('%.4f', (N - 1) * (N - 2) / N ^ 2), '0.9587');
rows_of = [1, 3];
for r = 1:2
  c = rows_of(r);
  share = STATED_EXACT(c, 1) / TRIPLES;
  for s = 1:3
    same (sprintf ('%s: silent corruptions a year, scrub every %s', LABELS{c}, SCRUB_NAMES{s}), ...
          sprintf ('%.2e', silent_per_year (share, SCRUB_HOURS(s))), STATED_SILENT{r, s});
  end
end
same ('ratio of the two, any interval', sprintf ('%.3f', STATED_EXACT(3, 1) / STATED_EXACT(1, 1)), '1.350');
share = STATED_EXACT(1, 1) / TRIPLES;
same ('factor a 30-day to 1-day scrub change buys, balanced', ...
      sprintf ('%.1f', silent_per_year (share, 720) / silent_per_year (share, 24)), '899.4');

printf ('\nAll chapter 31 numbers reproduced.\n');
