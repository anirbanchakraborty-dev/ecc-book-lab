% Reproduce every number in chapter 5 of Error Correction Coding, "Source Coding and Prefix Codes".
%
% This is the GNU Octave twin of python/ch05_source_coding_and_prefix_codes.py. It needs
% core Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch05_source_coding_and_prefix_codes.m
%
% Each check prints the value computed here beside the value the chapter states,
% and the program stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints. Codewords, parses
% and decoded strings must match the chapter exactly. A run that reaches the
% final line has reproduced every number in the chapter. The two routes print
% the same lines, and tools/compare-routes.sh checks that against the Python.
%
% The Python module also holds the Huffman, Shannon-code and block-coding functions
% that ch05_charts.py plots. Charts stay Python-only, and the functions are here
% because the checks use them too.
%
% Differences from the Python that a side-by-side reader meets:
%   - Python's fractions.Fraction is exact. Here a Fraction is the row
%     [numerator denominator] in lowest terms, built by fr and combined by fr_add,
%     fr_sub, fr_mul and fr_div. A list of Fractions is a matrix with one Fraction
%     per row. Every numerator and denominator in this chapter stays far below
%     2^53, the largest range in which a double holds every integer exactly.
%   - A list of floats is a column of doubles, so a function such as entropy tells
%     the two kinds of list apart by their number of columns.
%   - Octave has one number type where Python has two, and py_str prints a double
%     as a Python float. Python ints are converted to int64 before printing, and
%     ilist and flist turn a vector into the cell row that prints as a Python list,
%     so that a list of one item still prints with its brackets.
%   - A Fraction prints as 3/2 on its own and as Fraction(3, 2) inside a list, as
%     Python's str and repr do. py_text below adds that to py_str.
%   - A Python dict of codewords is a struct whose field names are the letters.
%     Octave keeps a struct's fields in the order they were written, which is the
%     order a Python dict iterates in.
%   - Octave indexes from 1, so the Python's bits[2] is bits(3) here, and bit
%     positions the Python reports from 0 subtract 1 from what find returns.
%   - Since Python 3.12 the built-in sum adds floats with compensated (Neumaier)
%     summation, and Octave's sum does not. py_sum below reproduces Python's.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%     Octave functions cannot see script variables, so the module constants of the
%     Python (FOUR, SIX, the codes and the trees) are functions that return a value.


1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% Helpers ------------------------------------------------------------------------

function check (label, computed, stated)
  % A Fraction arrives as [numerator denominator] and is divided once, which
  % rounds the exact quotient the way Python's float () of a Fraction does.
  computed = as_float (computed);
  value = str2double (strrep (stated, ',', ''));
  printf ('  %-60s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  % Python compares the two values with !=. Here they are compared through their
  % repr text, which for the strings, ints, floats and Fractions of this chapter
  % is equal exactly when the values are.
  printf ('  %-60s %-16s %s\n', label, py_text (computed), py_text (stated));
  if (~strcmp (py_text (computed, true), py_text (stated, true)))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_text (computed, true), py_text (stated, true));
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

function s = py_text (x, quoted)
  % py_str, extended to the Fractions made by frac: 3/2 for str and
  % Fraction(3, 2) for repr, which is what a list or a tuple prints its items with.
  if (nargin < 2)
    quoted = false;
  end
  if (isstruct (x) && isfield (x, 'frac'))
    if (quoted)
      s = sprintf ('Fraction(%d, %d)', x.frac(1), x.frac(2));
    elseif (x.frac(2) == 1)
      s = sprintf ('%d', x.frac(1));
    else
      s = sprintf ('%d/%d', x.frac(1), x.frac(2));
    end
  elseif (isstruct (x) && isfield (x, 'tuple'))
    items = cellfun (@(v) py_text (v, true), x.tuple, 'UniformOutput', false);
    s = ['(' strjoin(items, ', ') ')'];
  elseif (iscell (x))
    items = cellfun (@(v) py_text (v, true), x, 'UniformOutput', false);
    s = ['[' strjoin(items, ', ') ']'];
  else
    s = py_str (x, quoted);
  end
end

function t = tup (varargin)
  % A Python tuple, in the form py_str prints.
  t = struct ('tuple', {varargin});
end

function c = ilist (v)
  % A Python list of ints, as a cell row so that one item keeps its brackets.
  c = num2cell (int64 (v(:).'));
end

function c = flist (v)
  % A Python list of floats.
  c = num2cell (double (v(:).'));
end


% Fractions ----------------------------------------------------------------------

function f = fr (n, d)
  % The Fraction n/d as [numerator denominator] in lowest terms, with the
  % denominator positive, as fractions.Fraction keeps it.
  if (nargin < 2)
    d = 1;
  end
  g = gcd (n, d);
  if (d < 0)
    g = -g;
  end
  f = [n / g, d / g];
end

function f = fr_add (a, b)
  f = fr (a(1) * b(2) + b(1) * a(2), a(2) * b(2));
end

function f = fr_sub (a, b)
  f = fr (a(1) * b(2) - b(1) * a(2), a(2) * b(2));
end

function f = fr_mul (a, b)
  f = fr (a(1) * b(1), a(2) * b(2));
end

function f = fr_div (a, b)
  f = fr (a(1) * b(2), a(2) * b(1));
end

function s = fr_cmp (a, b)
  % -1, 0 or 1 as a is below, equal to or above b, by cross-multiplying.
  s = sign (a(1) * b(2) - b(1) * a(2));
end

function f = fr_abs (a)
  f = [abs(a(1)), a(2)];
end

function total = fr_sum (F)
  % Python's sum of a list of Fractions, one per row of F. An empty list sums to 0.
  total = fr (0);
  for i = 1:rows (F)
    total = fr_add (total, F(i, :));
  end
end

function q = fr_int (a)
  % Python's int () of a Fraction, which truncates toward zero.
  q = (a(1) - rem (a(1), a(2))) / a(2);
end

function F = fracs (numerators, den)
  % [Fraction (n, den) for n in numerators], one Fraction per row.
  F = zeros (numel (numerators), 2);
  for i = 1:numel (numerators)
    F(i, :) = fr (numerators(i), den);
  end
end

function x = as_float (v)
  % Python's float () of a Fraction [n d] or of a number.
  if (numel (v) == 2)
    x = v(1) / v(2);
  else
    x = double (v);
  end
end

function x = as_floats (probs)
  % [float (p) for p in probs] as a column, for a list of Fractions or of floats.
  if (columns (probs) == 2)
    x = probs(:, 1) ./ probs(:, 2);
  else
    x = probs(:);
  end
end

function s = frac (f)
  % A Fraction marked for printing, so that py_text writes it the Python way.
  s = struct ('frac', f);
end

function c = frac_list (F)
  % A list of Fractions marked for printing.
  c = cell (1, rows (F));
  for i = 1:rows (F)
    c{i} = frac (F(i, :));
  end
end

function [v, L] = common_den (F)
  % The Fractions of F written over their least common denominator L, so that
  % adding and comparing them is integer arithmetic on the numerators v.
  L = 1;
  for i = 1:rows (F)
    L = lcm (L, F(i, 2));
  end
  v = F(:, 1) .* (L ./ F(:, 2));
end

function F = fr_sort_desc (F)
  % Python's sorted (F, reverse=True) for a list of Fractions. Octave's sort is
  % stable in descending order too, as Python's is.
  [v, ~] = common_den (F);
  [~, order] = sort (v, 'descend');
  F = F(order, :);
end


% Sources, entropy and average length --------------------------------------------

function H = entropy (probs)
  % Entropy in bits of a distribution given as a list of probabilities.
  p = as_floats (probs);
  p = p(p > 0);
  H = py_sum (p .* log2 (1 ./ p));
end

function total = average_length (probs, lengths)
  if (columns (probs) == 2)
    total = fr (0);
    for i = 1:rows (probs)
      total = fr_add (total, fr_mul (probs(i, :), fr (lengths(i))));
    end
  else
    total = py_sum (probs(:) .* lengths(:));
  end
end

function v = length_variance (probs, lengths)
  mean = average_length (probs, lengths);
  second = fr (0);
  for i = 1:rows (probs)
    second = fr_add (second, fr_mul (probs(i, :), fr (lengths(i) * lengths(i))));
  end
  v = fr_sub (second, fr_mul (mean, mean));
end

function D = relative_entropy (p_dist, q_dist)
  p = as_floats (p_dist);
  q = as_floats (q_dist);
  keep = p > 0;
  D = py_sum (p(keep) .* log2 (p(keep) ./ q(keep)));
end

function k = kraft_sum (lengths)
  k = fr (0);
  for l = lengths
    k = fr_add (k, fr (1, 2 ^ l));
  end
end


% Codes as lists of codewords ----------------------------------------------------

function r = is_prefix_free (words)
  r = true;
  for i = 1:numel (words)
    for j = 1:numel (words)
      a = words{i};
      b = words{j};
      if (~strcmp (a, b) && strncmp (b, a, numel (a)))
        r = false;
        return;
      end
    end
  end
end

function found = parses (bits, words)
  % Every way to split a bit string into codewords, as strings of letters.
  if (isempty (bits))
    found = {''};
    return;
  end
  found = {};
  letters = fieldnames (words);
  for i = 1:numel (letters)
    word = words.(letters{i});
    % strncmp is false when bits is shorter than word, as startswith is.
    if (strncmp (bits, word, numel (word)))
      rest = parses (bits(numel (word) + 1:end), words);
      found = [found, cellfun(@(r) [letters{i} r], rest, 'UniformOutput', false)];
    end
  end
  % Octave sorts a cell of strings by character code, as Python's sorted does.
  found = sort (found);
end

function [letters, current] = decode_prefix (bits, words)
  % Decode left to right with a prefix code, returning the letters and any leftover bits.
  % The Python returns the pair as a tuple, and the callers here build tup from it.
  names = fieldnames (words);
  codewords = cellfun (@(c) words.(c), names, 'UniformOutput', false);
  letters = '';
  current = '';
  for b = bits
    current(end + 1) = b;
    k = find (strcmp (codewords, current), 1, 'last');
    if (~isempty (k))
      letters(end + 1) = names{k};
      current = '';
    end
  end
end

function bits = encode (message, words)
  bits = '';
  for c = message
    bits = [bits words.(c)];
  end
end

function out = dangling (a_set, b_set)
  % The inner function of the Python's sardinas_patterson. A set of strings is a
  % sorted cell row without repeats, which unique and union return.
  out = {};
  for i = 1:numel (a_set)
    for j = 1:numel (b_set)
      a = a_set{i};
      b = b_set{j};
      if (~strcmp (a, b) && strncmp (b, a, numel (a)))
        out{end + 1} = b(numel (a) + 1:end);
      end
    end
  end
  out = unique (out);
end

function [ud, history] = sardinas_patterson (words)
  % Return (uniquely decodable, list of dangling-suffix sets in order).
  code = unique (words);
  current = dangling (code, code);
  history = {};
  seen = {};
  while (~isempty (current))
    history{end + 1} = current;
    if (~isempty (intersect (current, code)))
      ud = false;
      return;
    end
    if (any (cellfun (@(s) isequal (s, current), seen)))
      ud = true;
      return;
    end
    seen{end + 1} = current;
    current = union (dangling (code, current), dangling (current, code));
  end
  ud = true;
end

function s = binary_digits (x, n)
  % The first n bits after the binary point of a number in [0, 1), with the
  % Fraction x carried as its numerator over a fixed denominator.
  num = x(1);
  den = x(2);
  s = '';
  for i = 1:n
    num = 2 * num;
    bit = num >= den;
    s(end + 1) = char ('0' + bit);
    num = num - bit * den;
  end
end

function words = kraft_code (lengths)
  % The prefix code built in the proof of Kraft's inequality, in the order the lengths are given.
  [~, order] = sort (lengths);
  words = cell (1, numel (lengths));
  point = fr (0);
  for i = order
    words{i} = binary_digits (point, lengths(i));
    point = fr_add (point, fr (1, 2 ^ lengths(i)));
  end
end


% Shannon's code, Fano's splitting and Huffman's algorithm -----------------------

function lengths = shannon_lengths (probs)
  % The smallest l with 2^-l at most p, which is the ceiling of log2(1/p), computed exactly.
  lengths = zeros (1, rows (probs));
  for i = 1:rows (probs)
    l = 0;
    while (fr_cmp (fr (1, 2 ^ l), probs(i, :)) > 0)
      l = l + 1;
    end
    lengths(i) = l;
  end
end

function words = shannon_code (probs)
  % Shannon's construction: probabilities in decreasing order, codeword i the first l_i bits of F_i.
  lengths = shannon_lengths (probs);
  words = {};
  cumulative = fr (0);
  for i = 1:rows (probs)
    words{end + 1} = binary_digits (cumulative, lengths(i));
    cumulative = fr_add (cumulative, probs(i, :));
  end
end

function [words, tied] = fano_split (probs, indices, words, tied)
  % The inner function split of the Python's fano_code. Octave has no closures
  % that change the caller's variables, so words and tied go in and come back out.
  if (numel (indices) == 1)
    return;
  end
  total = fr_sum (probs(indices, :));
  best = 0;
  best_gap = [];
  upper = fr (0);
  gaps = zeros (0, 2);
  for k = 1:numel (indices) - 1
    upper = fr_add (upper, probs(indices(k), :));
    gap = fr_abs (fr_sub (fr_mul (fr (2), upper), total));
    gaps(end + 1, :) = gap;
    if (isempty (best_gap) || fr_cmp (gap, best_gap) < 0)
      best = k;
      best_gap = gap;
    end
  end
  if (sum (all (gaps == best_gap, 2)) > 1)
    tied = true;
  end
  for i = indices(1:best)
    words{i}(end + 1) = '0';
  end
  for i = indices(best + 1:end)
    words{i}(end + 1) = '1';
  end
  [words, tied] = fano_split (probs, indices(1:best), words, tied);
  [words, tied] = fano_split (probs, indices(best + 1:end), words, tied);
end

function [words, tied] = fano_code (probs)
  % Fano's splitting of a list in decreasing order. Returns the codewords and whether a split was tied.
  n = rows (probs);
  [words, tied] = fano_split (probs, 1:n, repmat ({''}, 1, n), false);
end

function [lengths, merges] = huffman (probs, merged_high)
  % Huffman's algorithm by the sorted-list method.
  %
  % The list is kept in decreasing order of probability and the last two entries are
  % merged at each step. A merged entry is placed above every entry of equal
  % probability when merged_high is true and below them otherwise, which is the only
  % choice the algorithm leaves open. Returns the codeword lengths in the order of
  % probs and the list of merged probabilities in the order they were formed.
  %
  % A list of Fractions is first written over one common denominator, so that every
  % sum and comparison below is exact integer arithmetic on the numerators. A list
  % of floats is added and compared as the Python adds and compares it.
  if (nargin < 2)
    merged_high = true;
  end
  exact = columns (probs) == 2;
  if (exact)
    [v, L] = common_den (probs);
  else
    v = probs(:);
  end
  n = numel (v);
  [p, order] = sort (v, 'descend');
  members = num2cell (order(:).');
  lengths = zeros (1, n);
  merges = zeros (n - 1, 1);
  for step = 1:n - 1
    weight = p(end - 1) + p(end);
    group = [members{end - 1}, members{end}];
    p(end - 1:end) = [];
    members(end - 1:end) = [];
    lengths(group) = lengths(group) + 1;
    merges(step) = weight;
    if (merged_high)
      spot = find (p <= weight, 1);
    else
      spot = find (p < weight, 1);
    end
    if (isempty (spot))
      spot = numel (p) + 1;
    end
    p = [p(1:spot - 1); weight; p(spot:end)];
    members = [members(1:spot - 1), {group}, members(spot:end)];
  end
  if (exact)
    F = zeros (n - 1, 2);
    for i = 1:n - 1
      F(i, :) = fr (merges(i), L);
    end
    merges = F;
  end
end

function L = huffman_length (probs)
  % Average Huffman codeword length, as the sum of the merged probabilities.
  [~, merges] = huffman (as_floats (probs));
  L = py_sum (merges);
end


% Blocks of a binary source ------------------------------------------------------

function probs = block_probs (p_one, n)
  % Probabilities of the 2^n blocks of n independent letters, each 1 with probability p_one.
  % itertools.product lists the blocks as the binary numbers 0 to 2^n - 1, so
  % block y holds as many 1s as y has 1 bits.
  ones_in = sum (dec2bin ((0:2 ^ n - 1).', n) == '1', 2);
  if (numel (p_one) == 2)
    a = p_one(1);
    b = p_one(2);
    probs = zeros (2 ^ n, 2);
    for i = 1:2 ^ n
      probs(i, :) = fr (a ^ ones_in(i) * (b - a) ^ (n - ones_in(i)), b ^ n);
    end
  else
    % A scalar raised to a vector of powers calls pow for each one, as Python's
    % float ** int does.
    probs = (p_one .^ ones_in) .* ((1 - p_one) .^ (n - ones_in));
  end
end

function r = block_rate (p_one, n)
  % Huffman bits per letter when blocks of n letters are coded.
  r = huffman_length (block_probs (p_one, n)) / n;
end

function h = binary_entropy (p)
  h = entropy ([p; 1 - p]);
end


% McMillan -----------------------------------------------------------------------

function n = mcmillan_first_failure (lengths)
  % Smallest n with K^n above n(l_max - 1) + 1, the bound every uniquely decodable code obeys.
  % The Python raises the exact Fraction K to the power n. Here K is a dyadic
  % fraction, which a double holds exactly, and K^n is rounded once by pow. No K^n
  % in this chapter comes near n(l_max - 1) + 1, so the comparison comes out the same.
  k = kraft_sum (lengths);
  top = max (lengths);
  kf = k(1) / k(2);
  n = 1;
  while (kf ^ n <= n * (top - 1) + 1)
    n = n + 1;
  end
end


% The chapter --------------------------------------------------------------------

function F = FOUR ()
  F = [1 2; 1 4; 1 8; 1 8];
end

function c = CODE_FIXED ()
  c = struct ('a', '00', 'b', '01', 'c', '10', 'd', '11');
end

function c = CODE_PREFIX ()
  c = struct ('a', '0', 'b', '10', 'c', '110', 'd', '111');
end

function c = CODE_SUFFIX ()
  c = struct ('a', '0', 'b', '01', 'c', '011', 'd', '111');
end

function c = CODE_AMBIG ()
  c = struct ('a', '0', 'b', '01', 'c', '100', 'd', '101');
end

function c = CODE_SHORT ()
  c = struct ('a', '0', 'b', '1', 'c', '01', 'd', '10');
end

function F = SIX ()
  F = fracs ([45 20 10 10 10 5], 100);
end

function s = LETTERS6 ()
  s = 'abcdef';
end

function c = TREE_A ()
  c = struct ('a', '0', 'b', '10', 'c', '1100', 'd', '1101', 'e', '1110', 'f', '1111');
end

function c = TREE_B ()
  c = struct ('a', '0', 'b', '110', 'c', '100', 'd', '101', 'e', '1110', 'f', '1111');
end


function lengths = lengths_of (code)
  lengths = cellfun (@(c) numel (code.(c)), fieldnames (code)).';
end

function w = values_of (code)
  % The Python's code.values (): the codewords in the order of the letters.
  w = struct2cell (code).';
end

function lengths = lengths_in (code, letters)
  % [len (code[c]) for c in letters]
  lengths = arrayfun (@(c) numel (code.(c)), letters);
end

function r = siblings (u, v)
  r = numel (u) == numel (v) && strcmp (u(1:end - 1), v(1:end - 1)) && ~strcmp (u, v);
end

function points = starting_points (lengths)
  % [int (sum (Fraction (1, 2 ** l) for l in lengths[:i]) * 16) for i in range (len (lengths))]
  points = zeros (1, numel (lengths));
  for i = 1:numel (lengths)
    points(i) = fr_int (fr_mul (kraft_sum (lengths(1:i - 1)), fr (16)));
  end
end


function main ()
  printf ('Codes for a Source\n');
  check ('four-letter source entropy H', entropy (FOUR ()), '1.75');
  check ('fixed-length code 00 01 10 11, average length', average_length (FOUR (), lengths_of (CODE_FIXED ())), '2');
  check ('code 0 10 110 111, average length', average_length (FOUR (), lengths_of (CODE_PREFIX ())), '1.75');
  % A Fraction minus a float is a float in Python, so the Fraction is converted first.
  check ('fixed-length code, bits above entropy', as_float (average_length (FOUR (), [2 2 2 2])) - entropy (FOUR ()), '0.25');
  check ('code 0 1 01 10, average length', average_length (FOUR (), lengths_of (CODE_SHORT ())), '1.25');
  check ('code 0 1 01 10, average below entropy by', entropy (FOUR ()) - as_float (average_length (FOUR (), lengths_of (CODE_SHORT ()))), '0.5');
  same ('code 0 1 01 10, readings of 01', parses ('01', CODE_SHORT ()), {'ab', 'c'});
  same ('code 0 1 01 10, readings of 0110', parses ('0110', CODE_SHORT ()), {'abba', 'abd', 'cba', 'cd'});

  printf ('\nWhen a Code Can Be Decoded\n');
  names = {'0 10 110 111', '0 01 011 111', '0 01 100 101'};
  codes = {CODE_PREFIX(), CODE_SUFFIX(), CODE_AMBIG()};
  for i = 1:3
    check (sprintf ('code %s, Kraft sum', names{i}), kraft_sum (lengths_of (codes{i})), '1');
    check (sprintf ('code %s, average length', names{i}), average_length (FOUR (), lengths_of (codes{i})), '1.75');
  end
  same ('code 0 01 100 101, readings of 0101', parses ('0101', CODE_AMBIG ()), {'ad', 'bb'});
  [ud, hist] = sardinas_patterson (values_of (CODE_AMBIG ()));
  same ('code 0 01 100 101, dangling suffixes', hist, {{'1'}, {'00', '01'}});
  require ('code 0 01 100 101 is not uniquely decodable', ~ud);
  [ud, hist] = sardinas_patterson (values_of (CODE_SUFFIX ()));
  require ('code 0 01 011 111 is uniquely decodable', ud);
  require ('code 0 01 011 111 is not a prefix code', ~is_prefix_free (values_of (CODE_SUFFIX ())));
  require ('code 0 10 110 111 is a prefix code', is_prefix_free (values_of (CODE_PREFIX ())));
  same ('code 0 01 011 111, readings of 0 then six 1s', parses (['0' repmat('1', 1, 6)], CODE_SUFFIX ()), {'add'});
  same ('code 0 01 011 111, readings of 0 then seven 1s', parses (['0' repmat('1', 1, 7)], CODE_SUFFIX ()), {'bdd'});
  same ('code 0 01 011 111, readings of 0 then eight 1s', parses (['0' repmat('1', 1, 8)], CODE_SUFFIX ()), {'cdd'});

  printf ('\nPrefix Codes and Code Trees\n');
  [letters, rest] = decode_prefix ('1100100111', CODE_PREFIX ());
  same ('decode 1100100111 with 0 10 110 111', tup (letters, rest), tup ('cabad', ''));
  same ('encode cabad', encode ('cabad', CODE_PREFIX ()), '1100100111');

  printf ('\nThe Kraft Inequality\n');
  seven = [2, 2, 3, 3, 3, 4, 4];
  check ('lengths 2 2 3 3 3 4 4, Kraft sum', kraft_sum (seven), '1');
  same ('depth-4 leaves below each codeword', ilist (2 .^ (4 - seven)), ilist ([4, 4, 2, 2, 2, 1, 1]));
  same ('total depth-4 leaves', int64 (sum (2 .^ (4 - seven))), int64 (16));
  same ('Kraft construction for 2 2 3 3 3 4 4', kraft_code (seven), {'00', '01', '100', '101', '110', '1110', '1111'});
  same ('starting points in sixteenths', ilist (starting_points (seven)), ilist ([0, 4, 8, 10, 12, 14, 15]));
  four_lengths = [1, 3, 3, 4];
  check ('lengths 1 3 3 4, Kraft sum', kraft_sum (four_lengths), '0.8125');
  same ('Kraft construction for 1 3 3 4', kraft_code (four_lengths), {'0', '100', '101', '1100'});
  same ('unused depth-4 leaves for 1 3 3 4', int64 (16 - sum (2 .^ (4 - four_lengths))), int64 (3));
  check ('lengths 1 2 2 3, Kraft sum', kraft_sum ([1, 2, 2, 3]), '1.125');

  printf ('\nMcMillan''s Inequality\n');
  % counts(m) is the Python's counts[m], the number of pairs whose encoding has m bits.
  counts = zeros (1, 6);
  lens = lengths_of (CODE_PREFIX ());
  for x = lens
    for y = lens
      counts(x + y) = counts(x + y) + 1;
    end
  end
  same ('code 0 10 110 111, pairs by encoded length 2..6', ilist (counts(2:6)), ilist ([1, 2, 5, 4, 4]));
  total = fr (0);
  for m = find (counts > 0)
    total = fr_add (total, fr (counts(m), 2 ^ m));
  end
  check ('sum of A_m 2^-m for pairs', total, '1');
  k = kraft_sum ([1, 1, 2]);
  check ('lengths 1 1 2, Kraft sum', k, '1.25');
  % float () of the exact Fraction (5/4)^n equals pow (1.25, n) for these n.
  check ('1.25^11', as_float (k) ^ 11, '11.64');
  check ('1.25^12', as_float (k) ^ 12, '14.55');
  same ('first n with 1.25^n above n + 1', int64 (mcmillan_first_failure ([1, 1, 2])), int64 (12));

  printf ('\nThe Entropy Lower Bound\n');
  check ('four-letter source, D(P||U)', relative_entropy (FOUR (), repmat (fr (1, 4), 4, 1)), '0.25');
  h6 = entropy (SIX ());
  check ('six-letter source entropy', h6, '2.1955');
  sh_len = shannon_lengths (SIX ());
  k_sh = kraft_sum (sh_len);
  check ('Shannon code, Kraft sum', k_sh, '0.59375');
  same ('Shannon code, Kraft sum as a fraction', frac (k_sh), frac (fr (19, 32)));
  check ('Shannon code, log2(1/K)', log2 (1 / as_float (k_sh)), '0.7521');
  q_sh = zeros (numel (sh_len), 2);
  for i = 1:numel (sh_len)
    q_sh(i, :) = fr_div (fr (1, 2 ^ sh_len(i)), k_sh);
  end
  check ('Shannon code, D(P||Q)', relative_entropy (SIX (), q_sh), '0.0025');
  check ('Shannon code, L - H', as_float (average_length (SIX (), sh_len)) - h6, '0.7545');
  check ('Shannon code, D + log2(1/K)', relative_entropy (SIX (), q_sh) + log2 (1 / as_float (k_sh)), '0.7545');

  printf ('\nShannon''s Code and Fano''s Splitting\n');
  ps = [0.45, 0.20, 0.10, 0.05];
  stated = {'1.152', '2.322', '3.322', '4.322'};
  for i = 1:4
    check (sprintf ('surprise of a letter with p = %s', py_str (ps(i))), log2 (1 / ps(i)), stated{i});
  end
  same ('Shannon lengths', ilist (sh_len), ilist ([2, 3, 4, 4, 4, 5]));
  check ('Shannon code, average length', average_length (SIX (), sh_len), '2.95');
  check ('H + 1', h6 + 1, '3.1955');
  six = SIX ();
  hundredths = zeros (1, 6);
  for i = 1:6
    hundredths(i) = fr_int (fr_mul (fr_sum (six(1:i - 1, :)), fr (100)));
  end
  same ('cumulative probabilities F in hundredths', ilist (hundredths), ilist ([0, 45, 65, 75, 85, 95]));
  same ('binary expansion of 0.45, eight bits', binary_digits (fr (45, 100), 8), '01110011');
  same ('binary expansion of 0.65, eight bits', binary_digits (fr (65, 100), 8), '10100110');
  same ('binary expansion of 0.85, eight bits', binary_digits (fr (85, 100), 8), '11011001');
  same ('binary expansion of 0.95, eight bits', binary_digits (fr (95, 100), 8), '11110011');
  sh_words = shannon_code (SIX ());
  same ('Shannon codewords', sh_words, {'00', '011', '1010', '1100', '1101', '11110'});
  require ('Shannon codewords form a prefix code', is_prefix_free (sh_words));
  [fano_words, fano_tied] = fano_code (SIX ());
  same ('Fano codewords', fano_words, {'0', '100', '101', '110', '1110', '1111'});
  require ('no Fano split of the six-letter source is tied', ~fano_tied);
  check ('Fano code, average length', average_length (SIX (), cellfun (@numel, fano_words)), '2.25');
  check ('Fano first split, upper group a', six(1, :), '0.45');
  check ('Fano second split, b and c', fr_add (six(2, :), six(3, :)), '0.30');
  check ('Fano second split, d e f', fr_add (fr_add (six(4, :), six(5, :)), six(6, :)), '0.25');
  check ('Fano third split, e and f', fr_add (six(5, :), six(6, :)), '0.15');

  printf ('\nHuffman''s Algorithm\n');
  [len_b, merges_b] = huffman (SIX (), true);
  [len_a, merges_a] = huffman (SIX (), false);
  same ('merges, merged entry placed high', flist (as_floats (merges_b)), flist ([0.15, 0.2, 0.35, 0.55, 1.0]));
  same ('merges, merged entry placed low', flist (as_floats (merges_a)), flist ([0.15, 0.2, 0.35, 0.55, 1.0]));
  same ('lengths, merged entry placed low (tree A)', ilist (len_a), ilist ([1, 2, 4, 4, 4, 4]));
  same ('lengths, merged entry placed high (tree B)', ilist (len_b), ilist ([1, 3, 3, 3, 4, 4]));
  check ('tree A average length', average_length (SIX (), len_a), '2.25');
  check ('tree B average length', average_length (SIX (), len_b), '2.25');
  check ('sum of merged probabilities', fr_sum (merges_b), '2.25');
  check ('tree A, 0.45 + 2(0.20) + 4(0.35)', fr_add (fr_add (fr (45, 100), fr_mul (fr (2), fr (20, 100))), fr_mul (fr (4), fr (35, 100))), '2.25');
  check ('tree B, 0.45 + 3(0.40) + 4(0.15)', fr_add (fr_add (fr (45, 100), fr_mul (fr (3), fr (40, 100))), fr_mul (fr (4), fr (15, 100))), '2.25');
  check ('tree A variance of lengths', length_variance (SIX (), len_a), '1.7875');
  check ('tree B variance of lengths', length_variance (SIX (), len_b), '1.3875');
  same ('tree A codewords have tree A lengths', ilist (lengths_in (TREE_A (), LETTERS6 ())), ilist (len_a));
  same ('tree B codewords have tree B lengths', ilist (lengths_in (TREE_B (), LETTERS6 ())), ilist (len_b));
  require ('tree A and tree B codewords are prefix codes', is_prefix_free (values_of (TREE_A ())) && is_prefix_free (values_of (TREE_B ())));
  ta = TREE_A ();
  require ('tree A: c,d siblings, e,f siblings, cd and ef siblings, b and cdef siblings', ...
           siblings (ta.c, ta.d) && siblings (ta.e, ta.f) ...
           && siblings (ta.c(1:end - 1), ta.e(1:end - 1)) && siblings (ta.b, ta.c(1:end - 2)));
  tb = TREE_B ();
  require ('tree B: c,d siblings, e,f siblings, b and ef siblings, cd and bef siblings', ...
           siblings (tb.c, tb.d) && siblings (tb.e, tb.f) ...
           && siblings (tb.b, tb.e(1:end - 1)) && siblings (tb.c(1:end - 1), tb.b(1:end - 1)));
  check ('Huffman redundancy L - H', 2.25 - h6, '0.0545');
  check ('Huffman efficiency H / L', h6 / 2.25, '0.976');
  same ('Huffman on the four-letter source', ilist (huffman (FOUR ())), ilist ([1, 2, 3, 3]));

  printf ('\nShannon''s Source Coding Theorem\n');
  hb = binary_entropy (0.1);
  check ('h(0.1)', hb, '0.4690');
  pair = block_probs (fr (1, 10), 2);
  same ('pair probabilities', flist (sort (as_floats (pair), 'descend')), flist ([0.81, 0.09, 0.09, 0.01]));
  pair_sorted = fr_sort_desc (pair);
  [pair_lengths, pair_merges] = huffman (pair_sorted);
  same ('pair Huffman lengths', ilist (pair_lengths), ilist ([1, 2, 3, 3]));
  same ('pair merges', flist (as_floats (pair_merges)), flist ([0.1, 0.19, 1.0]));
  l2 = fr_sum (pair_merges);
  check ('L_2', l2, '1.29');
  check ('L_2 / 2', fr_div (l2, fr (2)), '0.645');
  triple = fr_sort_desc (block_probs (fr (1, 10), 3));
  [~, triple_merges] = huffman (triple);
  l3 = fr_sum (triple_merges);
  same ('triple merges', flist (as_floats (triple_merges)), flist ([0.01, 0.018, 0.028, 0.109, 0.162, 0.271, 1.0]));
  check ('L_3', l3, '1.598');
  check ('L_3 / 3', fr_div (l3, fr (3)), '0.5327');
  quad = fr_sort_desc (block_probs (fr (1, 10), 4));
  [~, quad_merges] = huffman (quad);
  l4 = fr_sum (quad_merges);
  check ('L_4', l4, '1.9702');
  check ('L_4 / 4', fr_div (l4, fr (4)), '0.4926');
  stated = {'1.4690', '0.9690', '0.8023', '0.7190'};
  for n = 1:4
    check (sprintf ('H + 1/%d', n), hb + 1 / n, stated{n});
  end
  stated = {'0.4802', '0.4702', '0.4743', '0.4758'};
  for n = 5:8
    check (sprintf ('L_%d / %d', n, n), block_rate (0.1, n), stated{n - 4});
  end
  require ('L_7 / 7 exceeds L_6 / 6, so the rate is not monotone in n', block_rate (0.1, 7) > block_rate (0.1, 6));
  block_totals = {'2.6667', '5.2222', '7.8148'};
  per_rolls = {'2.6667', '2.6111', '2.6049'};
  for n = 1:3
    die_block = repmat (fr (1, 6 ^ n), 6 ^ n, 1);
    [lengths_die, merges_die] = huffman (die_block);
    check (sprintf ('fair die, blocks of %d, L_%d', n, n), fr_sum (merges_die), block_totals{n});
    check (sprintf ('fair die, blocks of %d, per roll', n), fr_div (fr_sum (merges_die), fr (n)), per_rolls{n});
    % Python's int () truncates toward zero, which fix does.
    same (sprintf ('fair die, blocks of %d, lengths are the two nearest whole numbers to log2 6^%d', n, n), ...
          ilist (unique (lengths_die)), ilist ([fix(n * log2 (6)), fix(n * log2 (6)) + 1]));
  end
  same ('fair die, pairs: codewords of length 6', int64 (sum (huffman (repmat (fr (1, 36), 36, 1)) == 6)), int64 (8));
  same ('fair die, triples: codewords of length 8', int64 (sum (huffman (repmat (fr (1, 216), 216, 1)) == 8)), int64 (176));
  check ('fair die entropy log2 6', log2 (6), '2.5850');

  printf ('\nCompression and Its Fragility\n');
  message = 'fadedcab';
  bits = encode (message, TREE_B ());
  same ('fadedcab encoded with tree B', bits, '1111010111101011000110');
  same ('number of bits', int64 (numel (bits)), int64 (22));
  % The Python's bits[2], the third bit, is bits(3) here.
  if (bits(3) == '0')
    flip = '1';
  else
    flip = '0';
  end
  flipped = [bits(1:2), flip, bits(4:end)];
  same ('third bit flipped', flipped, '1101010111101011000110');
  [letters, rest] = decode_prefix (flipped, TREE_B ());
  same ('decoded after the flip', tup (letters, rest), tup ('bdafadcab', ''));
  same ('letters decoded', int64 (numel (letters)), int64 (9));

  printf ('\nLooking Ahead\n');
  check ('h(0.01)', binary_entropy (0.01), '0.0808');
  check ('1 bit against h(0.01)', 1 / binary_entropy (0.01), '12.4');
  check ('p = 0.01, pairs, bits per letter', block_rate (0.01, 2), '0.5150');
  check ('p = 0.01, triples, bits per letter', block_rate (0.01, 3), '0.3533');

  printf ('\nWorked Examples\n');
  five = fracs ([42, 15, 15, 15, 13], 100);
  check ('five-letter source entropy', entropy (five), '2.1399');
  [len5, merges5] = huffman (five);
  same ('five-letter Huffman lengths', ilist (len5), ilist ([1, 3, 3, 3, 3]));
  same ('five-letter merges', flist (as_floats (merges5)), flist ([0.28, 0.3, 0.58, 1.0]));
  check ('five-letter Huffman average', average_length (five, len5), '2.16');
  [f5, tied5] = fano_code (five);
  same ('five-letter Fano codewords', f5, {'00', '01', '10', '110', '111'});
  require ('no Fano split of the five-letter source is tied', ~tied5);
  check ('five-letter Fano average', average_length (five, cellfun (@numel, f5)), '2.28');
  check ('Fano first split gap with a alone, 0.58 - 0.42', fr_sub (fr (58, 100), fr (42, 100)), '0.16');
  check ('Fano first split gap with a and b, 0.57 - 0.43', fr_sub (fr (57, 100), fr (43, 100)), '0.14');
  check ('five-letter, Fano minus Huffman', fr_sub (fr (228, 100), fr (216, 100)), '0.12');

  lean = fracs ([40, 30, 25, 5], 100);
  check ('four-letter lean source entropy', entropy (lean), '1.7660');
  [len_lean, merges_lean] = huffman (lean);
  same ('lean source Huffman lengths', ilist (len_lean), ilist ([1, 2, 3, 3]));
  same ('lean source merges', flist (as_floats (merges_lean)), flist ([0.3, 0.6, 1.0]));
  check ('lean source Huffman average', average_length (lean, len_lean), '1.9');
  same ('lean source Shannon lengths', ilist (shannon_lengths (lean)), ilist ([2, 2, 2, 5]));
  check ('lean source Shannon average', average_length (lean, shannon_lengths (lean)), '2.15');
  check ('surprise of 0.25', log2 (4), '2');
  check ('surprise of 0.05', log2 (20), '4.322');

  steep = fracs ([40, 25, 15, 10, 5, 5], 100);
  [len_steep, merges_steep] = huffman (steep);
  same ('steep source Huffman lengths', ilist (len_steep), ilist ([1, 2, 3, 4, 5, 5]));
  same ('steep source merges', flist (as_floats (merges_steep)), flist ([0.1, 0.2, 0.35, 0.6, 1.0]));
  check ('steep source Huffman average', average_length (steep, len_steep), '2.25');
  check ('steep source entropy', entropy (steep), '2.2037');
  check ('steep source redundancy', 2.25 - entropy (steep), '0.0463');
  check ('surprise of 0.05 rounded up', shannon_lengths (fr (5, 100))(1), '5');

  uniform6 = repmat (fr (1, 6), 6, 1);
  check ('tree B used on a uniform source, average length', average_length (uniform6, len_b), '3');
  check ('log2 6', log2 (6), '2.585');
  check ('excess over log2 6', 3 - log2 (6), '0.415');
  dyadic_b = zeros (numel (len_b), 2);
  for i = 1:numel (len_b)
    dyadic_b(i, :) = fr (1, 2 ^ len_b(i));
  end
  check ('D(U || 2^-l) for tree B', relative_entropy (uniform6, dyadic_b), '0.415');
  len_u = huffman (uniform6);
  same ('Huffman lengths for the uniform source', ilist (sort (len_u)), ilist ([2, 2, 3, 3, 3, 3]));
  check ('Huffman average for the uniform source', average_length (uniform6, len_u), '2.667');
  check ('tree B on the six-letter source, cross-check', average_length (SIX (), len_b), '2.25');

  printf ('\nExercises\n');
  cases = {{'1', '01', '001', '000'}, true, true; ...
           {'1', '10', '00'}, true, false; ...
           {'0', '10', '01'}, false, false};
  for i = 1:rows (cases)
    words = cases{i, 1};
    ud = sardinas_patterson (words);
    same (sprintf ('code %s: uniquely decodable, prefix', strjoin (words, ' ')), ...
          tup (ud, is_prefix_free (words)), tup (cases{i, 2}, cases{i, 3}));
  end
  same ('code 1 10 00, readings of 1000', parses ('1000', struct ('x', '1', 'y', '10', 'z', '00')), {'yz'});
  same ('code 0 10 01, readings of 010', parses ('010', struct ('x', '0', 'y', '10', 'z', '01')), {'xy', 'zx'});
  check ('code 1 01 001 000, Kraft sum', kraft_sum ([1, 2, 3, 3]), '1');
  ex_lengths = [1, 2, 4, 4, 4];
  check ('lengths 1 2 4 4 4, Kraft sum', kraft_sum (ex_lengths), '0.9375');
  same ('Kraft construction for 1 2 4 4 4', kraft_code (ex_lengths), {'0', '10', '1100', '1101', '1110'});
  k = kraft_sum ([1, 2, 2, 3]);
  n_fail = mcmillan_first_failure ([1, 2, 2, 3]);
  same ('first n with 1.125^n above 2n + 1', int64 (n_fail), int64 (37));
  % The exact Fraction (9/8)^36 has a numerator past 2^53, so it cannot be held in
  % doubles. float () of it equals pow (1.125, 36), and likewise for 37.
  check ('1.125^36', as_float (k) ^ 36, '69.42');
  check ('1.125^37', as_float (k) ^ 37, '78.10');
  q4 = fracs ([4, 3, 2, 1], 10);
  check ('source 0.4 0.3 0.2 0.1, entropy', entropy (q4), '1.8464');
  same ('source 0.4 0.3 0.2 0.1, Huffman lengths', ilist (huffman (q4)), ilist ([1, 2, 3, 3]));
  check ('source 0.4 0.3 0.2 0.1, Huffman average', average_length (q4, huffman (q4)), '1.9');
  same ('source 0.4 0.3 0.2 0.1, Shannon lengths', ilist (shannon_lengths (q4)), ilist ([2, 2, 3, 4]));
  check ('source 0.4 0.3 0.2 0.1, Shannon average', average_length (q4, shannon_lengths (q4)), '2.4');
  check ('h(0.2)', binary_entropy (0.2), '0.7219');
  pair8 = fr_sort_desc (block_probs (fr (2, 10), 2));
  same ('pairs of the 0.8 coin', flist (as_floats (pair8)), flist ([0.64, 0.16, 0.16, 0.04]));
  [~, merges8] = huffman (pair8);
  same ('pairs of the 0.8 coin, merges', flist (as_floats (merges8)), flist ([0.2, 0.36, 1.0]));
  check ('pairs of the 0.8 coin, L_2', fr_sum (merges8), '1.56');
  check ('pairs of the 0.8 coin, per letter', fr_div (fr_sum (merges8), fr (2)), '0.78');
  uniform5 = repmat (fr (1, 5), 5, 1);
  len_u5 = huffman (uniform5);
  same ('five equal letters, Huffman lengths', ilist (sort (len_u5)), ilist ([2, 2, 2, 3, 3]));
  check ('five equal letters, Huffman average', average_length (uniform5, len_u5), '2.4');
  check ('log2 5', log2 (5), '2.3219');
  check ('five equal letters, redundancy', 2.4 - log2 (5), '0.0781');
  sent = encode ('decafd', TREE_B ());
  same ('decafd encoded with tree B', sent, '101111010001111101');
  received = '001111010001111101';
  % find counts from 1 and the Python's bit positions from 0.
  same ('received stream differs from sent in the first bit only', ilist (find (sent(1:18) ~= received(1:18)) - 1), ilist (0));
  [letters, rest] = decode_prefix (received, TREE_B ());
  same ('received tree-B stream decoded', tup (letters, rest), tup ('aafacafd', ''));
  [late, late_rest] = decode_prefix ([sent(1:16), '1', sent(18:end)], TREE_B ());
  same ('seventeenth bit flipped instead, decoded with leftover', tup (late, late_rest), tup ('decaf', '111'));
  same ('seventeenth bit of the sent stream, before the flip', sent(17), '0');
  decafd = 'decafd';
  same ('first bit flipped: the last four letters are the message''s', letters(end - 3:end), decafd(end - 3:end));

  further_checks ();

  printf ('\n');
  printf ('All chapter 5 numbers reproduced.\n');
end


function ends = codeword_ends (bits, words)
  % Bit positions, counted from 1, at which each decoded codeword ends.
  names = fieldnames (words);
  codewords = cellfun (@(c) words.(c), names, 'UniformOutput', false);
  ends = [];
  current = '';
  for position = 1:numel (bits)
    current(end + 1) = bits(position);
    if (any (strcmp (codewords, current)))
      ends(end + 1) = position;
      current = '';
    end
  end
end

function gaps = split_gaps (probs)
  % The gap between the two groups for every place a sorted list can be split.
  total = fr_sum (probs);
  gaps = zeros (rows (probs) - 1, 2);
  for k = 1:rows (probs) - 1
    gaps(k, :) = fr_abs (fr_sub (fr_mul (fr (2), fr_sum (probs(1:k, :))), total));
  end
end

function further_checks ()
  % Numbers stated in the prose and the figure texts that the sections above do not reach.
  printf ('\n');
  printf ('Further checks on the text\n');

  same ('code 0 1 01 10, Kraft sum', frac (kraft_sum (lengths_of (CODE_SHORT ()))), frac (fr (3, 2)));
  [len_four, merges_four] = huffman (FOUR ());
  same ('Huffman merges on the four-letter source', frac_list (merges_four), frac_list ([fr(1, 4); fr(1, 2); fr(1)]));
  check ('four-letter merges add to', fr_sum (merges_four), '1.75');
  same ('Huffman lengths on the four-letter source are code A''s', ilist (len_four), ilist (lengths_of (CODE_PREFIX ())));

  % The Python's 0b10 is the binary literal for 2, which bin2dec reads here.
  same ('interval named by the string 10', tup (frac (fr (bin2dec ('10'), 4)), frac (fr (bin2dec ('10') + 1, 4))), ...
        tup (frac (fr (1, 2)), frac (fr (3, 4))));
  same ('lengths 1 3 3 4, Kraft sum as a fraction', frac (kraft_sum ([1, 3, 3, 4])), frac (fr (13, 16)));
  same ('lengths 1 3 3 4, unused fraction of the line', frac (fr_sub (fr (1), kraft_sum ([1, 3, 3, 4]))), frac (fr (3, 16)));
  shortened = {'0', '100', '101', '11'};
  same ('shortened code 0 100 101 11, Kraft sum', frac (kraft_sum (cellfun (@numel, shortened))), frac (fr (1)));
  require ('shortened code 0 100 101 11 is a prefix code', is_prefix_free (shortened));

  code_a = CODE_PREFIX ();
  pair_lengths = [];
  for x = 'abcd'
    for y = 'abcd'
      pair_lengths(end + 1) = numel (code_a.(x)) + numel (code_a.(y));
    end
  end
  same ('code A pairs, A_m for m = 2 to 6', ilist (arrayfun (@(m) sum (pair_lengths == m), 2:6)), ilist ([1, 2, 5, 4, 4]));

  shannon_k = kraft_sum (shannon_lengths (SIX ()));
  same ('Shannon code, Kraft sum as 1/4 + 1/8 + 3/16 + 1/32', ...
        frac (fr_add (fr_add (fr_add (fr (1, 4), fr (1, 8)), fr (3, 16)), fr (1, 32))), frac (shannon_k));
  same ('Shannon code, unused fraction of the tree', frac (fr_sub (fr (1), shannon_k)), frac (fr (13, 32)));
  check ('Shannon code, unused fraction in percent', 100 * as_float (fr_sub (fr (1), shannon_k)), '41');
  same ('binary expansion of 0.75, four bits', binary_digits (fr (3, 4), 4), '1100');
  fano_words = fano_code (SIX ());
  check ('Shannon average minus Fano average', ...
         fr_sub (average_length (SIX (), shannon_lengths (SIX ())), average_length (SIX (), cellfun (@numel, fano_words))), '0.70');
  six = SIX ();
  first_gaps = split_gaps (six);
  check ('Fano first split, gap with a alone', first_gaps(1, :), '0.10');
  check ('Fano first split, gap with a and b', first_gaps(2, :), '0.30');
  rest_gaps = split_gaps (six(2:end, :));
  same ('Fano split of b to f, gaps after b, c, d, e', flist (as_floats (rest_gaps)), flist ([0.15, 0.05, 0.25, 0.45]));
  def_gaps = split_gaps (six(4:end, :));
  same ('Fano split of d e f, gaps after d and after e', flist (as_floats (def_gaps)), flist ([0.05, 0.15]));
  check ('tree A, 0.45 + 2(0.20) + 4(0.35)', average_length (SIX (), lengths_in (TREE_A (), LETTERS6 ())), '2.25');

  five = fracs ([42, 15, 15, 15, 13], 100);
  cde_gaps = split_gaps (five(3:end, :));
  check ('five-letter Fano split of c d e, gap with c alone', cde_gaps(1, :), '0.13');
  check ('five-letter Fano split of c d e, gap splitting off e alone', cde_gaps(2, :), '0.17');
  lean = fracs ([40, 30, 25, 5], 100);
  same ('lean source Shannon lengths, Kraft sum', frac (kraft_sum (shannon_lengths (lean))), frac (fr (25, 32)));
  uniform6 = repmat (fr (1, 6), 6, 1);
  check ('tree B on the uniform source minus the uniform Huffman code', ...
         fr_sub (average_length (uniform6, lengths_of (TREE_B ())), average_length (uniform6, huffman (uniform6))), '0.3333');

  h01 = binary_entropy (0.1);
  stated = {'0.6690', '0.6357', '0.6119', '0.5940'};
  for n = 5:8
    check (sprintf ('H + 1/%d', n), h01 + 1 / n, stated{n - 4});
  end
  % Python's all () stops at the first rate outside the interval, and so does this loop.
  inside = true;
  for n = 1:8
    rate = block_rate (0.1, n);
    if (~(h01 <= rate && rate < h01 + 1 / n))
      inside = false;
      break;
    end
  end
  require ('every rate from n = 1 to 8 lies in [h(0.1), h(0.1) + 1/n)', inside);
  counts = [5, 6, 36, 216];
  expected = [3, 2; 2, 4; 28, 8; 40, 176];
  for i = 1:4
    count = counts(i);
    % count.bit_length () - 1 is the exponent of the largest power of 2 not above count.
    k = numel (dec2bin (count)) - 1;
    lengths_n = huffman (repmat (fr (1, count), count, 1));
    want = tup (int64 (expected(i, 1)), int64 (expected(i, 2)));
    same (sprintf ('%d equal letters, Huffman leaves at depths %d and %d', count, k, k + 1), ...
          tup (int64 (sum (lengths_n == k)), int64 (sum (lengths_n == k + 1))), want);
    same (sprintf ('%d equal letters, 2^(k+1) - N and 2(N - 2^k)', count), ...
          tup (int64 (2 ^ (k + 1) - count), int64 (2 * (count - 2 ^ k))), want);
  end

  message = 'fadedcab';
  sent = encode (message, TREE_B ());
  if (sent(3) == '1')
    flip = '0';
  else
    flip = '1';
  end
  received = [sent(1:2), flip, sent(4:end)];
  sent_ends = codeword_ends (sent, TREE_B ());
  received_ends = codeword_ends (received, TREE_B ());
  common = intersect (sent_ends, received_ends);
  same ('fadedcab, first codeword end shared by sent and received', int64 (common(1)), int64 (12));
  require ('fadedcab, every end from bit 12 on is shared', isequal (sent_ends(sent_ends >= 12), received_ends(received_ends >= 12)));
  decoded = decode_prefix (received, TREE_B ());
  % zip stops at the shorter of the two strings.
  m = min (5, numel (decoded));
  require ('fadedcab, the first five decoded letters are all wrong', all (decoded(1:m) ~= message(1:m)));
  same ('fadedcab, the last four decoded letters, one place late', tup (decoded(6:end), message(5:end)), tup ('dcab', 'dcab'));

  same ('code 1 10 00 written backward', cellfun (@fliplr, {'1', '10', '00'}, 'UniformOutput', false), {'1', '01', '00'});
  require ('code 1 10 00 written backward is a prefix code', is_prefix_free ({'1', '01', '00'}));
  same ('starting points for lengths 1 2 4 4 4, in sixteenths', ilist (starting_points ([1, 2, 4, 4, 4])), ilist ([0, 8, 12, 13, 14]));
  same ('lengths 1 2 4 4 4 with 1110 shortened to 111, Kraft sum', frac (kraft_sum ([1, 2, 4, 4, 3])), frac (fr (1)));
  require ('codewords 0 10 1100 1101 111 form a prefix code', is_prefix_free ({'0', '10', '1100', '1101', '111'}));
  dyadic = [repmat(fr (1, 4), 3, 1); repmat(fr (1, 8), 2, 1)];
  dyadic_lengths = shannon_lengths (dyadic);
  same ('dyadic source, lengths log2(1/p)', ilist (dyadic_lengths), ilist ([2, 2, 2, 3, 3]));
  same ('dyadic source, Kraft construction', kraft_code (dyadic_lengths), {'00', '01', '10', '110', '111'});
  check ('dyadic source, average length', average_length (dyadic, dyadic_lengths), '2.25');
  check ('dyadic source, entropy', entropy (dyadic), '2.25');
end


main ();
