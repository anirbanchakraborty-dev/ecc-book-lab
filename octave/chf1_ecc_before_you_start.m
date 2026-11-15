% Reproduce every number in front-matter chapter f1 of Error Correction Coding, "Before You Start".
%
% This is the GNU Octave twin of python/chf1_ecc_before_you_start.py. It needs core
% Octave only (no pkg load, no Forge packages). From the repository root:
%
%     octave --no-window-system --norc --quiet octave/chf1_ecc_before_you_start.m
%
% The chapter is an orientation and a notation primer, so its numbers are of two
% kinds. The first describes the book itself: the chapters in each Part, the
% totals at each level, and the layout of Part I. Those counts are embedded below
% as a table copied from the book's table of contents, and the script checks that
% they add up to the 188 chapters the book plans. The second kind is the worked
% arithmetic of the primer: sums and products, factorials and binomial
% coefficients, powers and base-2 logarithms, floors and ceilings, XOR on bits and
% bit vectors, small sets, a few probabilities, and the cases the three proofs
% illustrate. Each of those is recomputed from scratch, and each theorem the
% chapter proves is also checked by brute force over a range of inputs.
%
% Every value computed here is printed beside the value the chapter states. The
% script stops with exit status 1 at the first pair that disagrees by more than
% half a unit in the last digit the chapter prints, and at the first count, list
% or bit pattern that differs at all. The two routes print the same lines, and
% tools/compare-routes.sh checks that against the Python.
%
% Differences from the Python that a side-by-side reader meets:
%   - Octave indexes from 1, so a Python index i becomes i + 1 wherever it touches
%     an array. Where the Python adds 1 to turn an index into a position,
%     Octave's find already returns the position.
%   - A script file may not start with a function, hence the "1;" below, and Octave
%     defines each function when execution reaches it, so main is called last.
%   - Octave functions cannot see script variables, so the module constants of the
%     Python (PARTS and PART_I) are functions that return the table.
%   - Octave has one default number type, the double. Python's ints are int64 here,
%     so integer results print without a decimal point and integer division is
%     idivide. Mixing an int64 with a double gives an int64, and Octave's integer
%     types saturate instead of growing, so every value below is far from the edge.
%   - Octave already owns the names xor and factorial, so the Python's xor is
%     xor_bits and its factorial is factorial_int, which stays exact in int64.
%   - Python's int() cuts the fraction off, which is Octave's fix. Octave's int64()
%     rounds to the nearest integer instead.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));


% The book's table of contents: label, title, level, chapters.
function p = PARTS ()
  p = {
    'Front matter', 'Before You Start', 'front', int64(3);
    'Part I', 'Errors, Information and Channels', 'undergraduate', int64(12);
    'Part II', 'Linear Block Codes', 'undergraduate', int64(19);
    'Part III', 'Finite Fields and Their Arithmetic', 'undergraduate', int64(15);
    'Part IV', 'Cyclic Codes, CRCs and Burst Errors', 'undergraduate', int64(12);
    'Part V', 'BCH and Reed-Solomon Codes', 'undergraduate', int64(13);
    'Part VI', 'Soft Decisions and Concatenation', 'graduate', int64(9);
    'Part VII', 'Convolutional Codes and Viterbi Decoding', 'graduate', int64(10);
    'Part VIII', 'Turbo Codes and Iterative Decoding', 'graduate', int64(10);
    'Part IX', 'LDPC Codes', 'graduate', int64(15);
    'Part X', 'Polar Codes', 'graduate', int64(8);
    'Part XI', 'ECC in Memory, Storage and Links', 'graduate', int64(14);
    'Part XII', 'Coded Modulation and Wireless', 'graduate', int64(9);
    'Part XIII', 'Advanced Algebraic Coding Theory', 'research', int64(19);
    'Part XIV', 'Codes, Cryptography and Security', 'research', int64(10);
    'Part XV', 'Frontiers', 'research', int64(7);
    'Reference', 'Reference', 'back', int64(3)
  };
end

% Part I as published: the nine teaching chapters, then the case study, the
% project and the lab, numbered 1 to 12.
function p = PART_I ()
  p = {
    'why-codes-exist', 'teaching';
    'probability-for-coding', 'teaching';
    'channels-signals-and-likelihoods', 'teaching';
    'information-and-entropy', 'teaching';
    'source-coding-and-prefix-codes', 'teaching';
    'compression-in-practice', 'teaching';
    'capacity-and-shannons-theorem', 'teaching';
    'distance-detection-and-correction', 'teaching';
    'error-control-strategies-and-performance', 'teaching';
    'case-study-codes-in-the-wild', 'case study';
    'project-channel-simulator', 'project';
    'lab-first-hardware-checker', 'lab'
  };
end


function check (label, computed, stated)
  value = str2double (regexprep (strrep (stated, ',', ''), '%+$', ''));
  if (stated(end) == '%')
    computed = 100 * double (computed);
  end
  computed = double (computed);
  printf ('  %-64s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed), stated);
    exit (1);
  end
end

function same (label, computed, stated)
  printf ('  %-64s %-16s %s\n', label, py_str (computed), py_str (stated));
  if (~isequal (computed, stated))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, py_str (computed, true), py_str (stated, true));
    exit (1);
  end
end

function t = as_tuple (varargin)
  % A Python tuple, which py_str prints with parentheses.
  t = struct ('tuple', {varargin});
end

function r = comb (n, k)
  % math.comb(n, k), exact. Multiplying before dividing keeps every step an integer,
  % and int64 holds each value this script reaches (the largest, comb (60, 30), is
  % about 1.2e17 against int64's 9.2e18). Octave's own nchoosek is exact only below 2^53.
  n = int64 (n);
  k = int64 (k);
  r = int64 (1);
  for i = int64 (1):k
    r = idivide (r * (n - k + i), i, 'floor');
  end
end

function r = factorial_int (n)
  % math.factorial(n) as an int64. The prod of an empty range is 1, as in Python.
  r = int64 (prod (int64 (1:n)));
end

function t = bit_tuples (n)
  % Every n-tuple of 0s and 1s, one per row, in the order itertools.product gives:
  % (0, 0), (0, 1), (1, 0), (1, 1) for n = 2. Counting in binary lists them so.
  if (n == 0)
    t = zeros (1, 0);              % one empty tuple
  else
    t = dec2bin (0:2 ^ n - 1, n) - '0';
  end
end

function s = strings (n)
  % Every n-bit string, in counting order, as a row of char rows.
  if (n == 0)
    s = {''};
  else
    s = cellstr (char (bit_tuples (n) + '0')).';
  end
end

function c = count_ones (s)
  % Python's s.count("1"), as an int64.
  c = int64 (sum (s == '1'));
end

function r = xor_bits (a, b)
  % Bitwise XOR of two bit strings of equal length.
  assert (numel (a) == numel (b));
  r = char ('0' + double (a ~= b));
end

function r = parity (bits)
  % The XOR of every bit of the string: 1 when it holds an odd number of 1s.
  r = mod (count_ones (bits), int64 (2));
end

function n = chapters (level)
  P = PARTS ();
  n = int64 (sum ([P{strcmp(P(:, 3), level).', 4}]));
end


function check_book ()
  printf ('The book\n');
  P = PARTS ();
  c = [P{:, 4}];                          % chapters in each section, as an int64 row
  is_part = strncmp (P(:, 1), 'Part ', 5).';
  same ('chapters in the book', int64 (sum (c)), int64 (188));
  same ('Parts, front and back matter aside', int64 (sum (is_part)), int64 (15));
  same ('chapters inside the fifteen Parts', int64 (sum (c(is_part))), int64 (182));
  same ('chapters of front matter', chapters ('front'), int64 (3));
  same ('chapters of back matter', chapters ('back'), int64 (3));
  same ('front and back matter together', chapters ('front') + chapters ('back'), int64 (6));
  same ('undergraduate chapters, Parts I to V', chapters ('undergraduate'), int64 (71));
  same ('graduate chapters, Parts VI to XII', chapters ('graduate'), int64 (75));
  same ('research chapters, Parts XIII to XV', chapters ('research'), int64 (36));
  same ('undergraduate Parts', P(strcmp (P(:, 3), 'undergraduate'), 1).',
        {'Part I', 'Part II', 'Part III', 'Part IV', 'Part V'});
  same ('graduate Parts', P(strcmp (P(:, 3), 'graduate'), 1).',
        {'Part VI', 'Part VII', 'Part VIII', 'Part IX', 'Part X', 'Part XI', 'Part XII'});
  same ('research Parts', P(strcmp (P(:, 3), 'research'), 1).',
        {'Part XIII', 'Part XIV', 'Part XV'});
  same ('chapters in each section, front matter to Reference', c,
        int64 ([3, 12, 19, 15, 12, 13, 9, 10, 10, 15, 8, 14, 9, 19, 10, 7, 3]));
  parts = P(is_part, [1 4]);              % label and chapters of the fifteen Parts
  counts = [parts{:, 2}];
  most = max (counts);
  same ('most chapters in one Part', most, int64 (19));
  same ('the Parts with that many', parts(counts == most, 1).', {'Part II', 'Part XIII'});
  fewest = min (counts);
  same ('fewest chapters in one Part', fewest, int64 (7));
  same ('the Part with that few', parts(counts == fewest, 1).', {'Part XV'});
  PI = PART_I ();
  same ('chapters in Part I', int64 (rows (PI)), P{strcmp (P(:, 1), 'Part I'), 4});
  same ('teaching chapters in Part I', int64 (sum (strcmp (PI(:, 2), 'teaching'))), int64 (9));
  same ('Part I closing chapters, in order', PI(10:end, 2).', {'case study', 'project', 'lab'});
  % find returns a 1-based position, which is the Python's index() + 1 already.
  same ('chapter number of the case study', ...
        int64 (find (strcmp (PI(:, 1), 'case-study-codes-in-the-wild') & strcmp (PI(:, 2), 'case study'))), int64 (10));
  same ('chapter number of the project', ...
        int64 (find (strcmp (PI(:, 1), 'project-channel-simulator') & strcmp (PI(:, 2), 'project'))), int64 (11));
  same ('chapter number of the lab', ...
        int64 (find (strcmp (PI(:, 1), 'lab-first-hardware-checker') & strcmp (PI(:, 2), 'lab'))), int64 (12));
  same ('reading order of this chapter', int64 (1), int64 (1));
  same ('reading order of chapter 1, after the three front chapters', chapters ('front') + 1, int64 (4));
end


function check_sums_and_products ()
  printf ('Sums and products\n');
  same ('1 + 2 + 3 + 4 + 5', int64 (sum (int64 (1:5))), int64 (15));
  same ('1 x 2 x 3 x 4', int64 (prod (int64 (1:4))), int64 (24));
  same ('a sum with no terms', int64 (sum ([])), int64 (0));
  same ('a product with no factors', int64 (prod ([])), int64 (1));
  p = 0.01;
  check ('all 8 bits survive at p = 0.01, (0.99)^8', (1 - p) ^ 8, '0.9227');
  check ('at least one of 8 bits flips at p = 0.01', 1 - (1 - p) ^ 8, '0.0773');
  check ('the same as a product of eight factors of 0.99', prod (repmat (1 - p, 1, 8)), '0.9227');
end


function check_binomials ()
  printf ('Factorials and binomial coefficients\n');
  same ('orders of three letters, 3!', factorial_int (3), int64 (6));
  same ('the orders listed', sort (cellstr (perms ('abc'))).',
        {'abc', 'acb', 'bac', 'bca', 'cab', 'cba'});
  same ('5!', factorial_int (5), int64 (120));
  same ('0!', factorial_int (0), int64 (1));
  s4 = strings (4);
  two_of_four = s4(cellfun (@count_ones, s4) == 2);
  same ('4-bit strings with exactly two 1s', sort (two_of_four, 'descend'),
        {'1100', '1010', '1001', '0110', '0101', '0011'});
  same ('their count, 4 choose 2', comb (4, 2), int64 (6));
  same ('4! / (2! 2!)', idivide (factorial_int (4), factorial_int (2) * factorial_int (2), 'floor'), int64 (6));
  same ('4! and 2! 2!', as_tuple (factorial_int (4), factorial_int (2) * factorial_int (2)),
        as_tuple (int64 (24), int64 (4)));
  same ('ordered choices of 2 positions from 4, 4 x 3', int64 (4) * 3, int64 (12));
  same ('each unordered pair counted 2! times, 12 / 2', idivide (int64 (12), factorial_int (2), 'floor'), int64 (6));
  same ('ways exactly 2 of 8 bits flip, 8 choose 2', comb (8, 2), int64 (28));
  same ('8 x 7 / 2', idivide (int64 (8) * 7, int64 (2), 'floor'), int64 (28));
  same ('Pascal: 7 choose 3', comb (7, 3), int64 (35));
  same ('Pascal: 6 choose 2 + 6 choose 3', comb (6, 2) + comb (6, 3), int64 (35));
  same ('6 choose 2 and 6 choose 3', as_tuple (comb (6, 2), comb (6, 3)), as_tuple (int64 (15), int64 (20)));
  row7 = int64 ([]);
  for i = 0:7
    row7(end + 1) = comb (7, i);
  end
  same ('row 7 of the binomial coefficients', row7, int64 ([1, 7, 21, 35, 35, 21, 7, 1]));
  same ('its sum', int64 (sum (row7)), int64 (128));
  same ('2^7', int64 (2) ^ 7, int64 (128));
  same ('row 7 reads the same backward', isequal (row7, fliplr (row7)), true);
  ok = true;
  for n = 2:40
    for k = 1:n - 1
      ok = ok && comb (n, k) == comb (n - 1, k - 1) + comb (n - 1, k);
    end
  end
  same ('Pascal''s rule holds for every 1 <= k <= n - 1, n <= 40', ok, true);
  ok = true;
  for n = 0:12
    ones_in = cellfun (@count_ones, strings (n));     % 1s in each n-bit string, counted directly
    for k = 0:n
      ok = ok && int64 (sum (ones_in == k)) == idivide (factorial_int (n), factorial_int (k) * factorial_int (n - k), 'floor');
    end
  end
  same ('formula n!/(k!(n-k)!) matches a direct count, n <= 12', ok, true);
  same ('ordered picks of 2 of 4 positions, 4!/2!', idivide (factorial_int (4), factorial_int (2), 'floor'), int64 (12));
  ok = true;
  for n = 1:12
    for k = 0:n
      ok = ok && int64 (prod (int64 (n - k + 1:n))) == idivide (factorial_int (n), factorial_int (n - k), 'floor');
    end
  end
  same ('n(n-1)...(n-k+1) = n!/(n-k)!, n up to 12', ok, true);
  first_rows = cell (1, 4);
  for n = 0:3
    row = int64 ([]);
    for k = 0:n
      row(end + 1) = comb (n, k);
    end
    first_rows{n + 1} = num2cell (row);       % a cell, so a one-entry row still prints as a list, [1]
  end
  same ('first rows of Pascal''s triangle', first_rows,
        {{int64(1)}, {int64(1), int64(1)}, {int64(1), int64(2), int64(1)}, {int64(1), int64(3), int64(3), int64(1)}});
end


function check_powers_and_logs ()
  printf ('Exponents and base-2 logarithms\n');
  same ('2^10', int64 (2) ^ 10, int64 (1024));
  same ('values of an 8-bit byte, 2^8', int64 (2) ^ 8, int64 (256));
  same ('2^3 x 2^4 = 2^7', int64 (2) ^ 3 * int64 (2) ^ 4, int64 (2) ^ 7);
  same ('2^7', int64 (2) ^ 7, int64 (128));
  same ('(2^3)^2 = 2^6', (int64 (2) ^ 3) ^ 2, int64 (64));
  same ('2^0', int64 (2) ^ 0, int64 (1));
  check ('2^-1', 2 ^ -1, '0.5');
  check ('2^-3', 2 ^ -3, '0.125');
  same ('log2 1024', log2 (1024), 10);          % a double, so both sides print as 10.0
  check ('log2 1000', log2 (1000), '9.966');
  check ('log2 10', log2 (10), '3.322');
  check ('3 x log2 10, the exponent that makes 2 reach 1000', 3 * log2 (10), '9.966');
  same ('log2 (8 x 32) = log2 8 + log2 32', as_tuple (log2 (8 * 32), log2 (8) + log2 (32)), as_tuple (8, 8));
  same ('log2 256', log2 (256), 8);
  same ('seconds per error at 10^-12 on a 10^10 bit/s link',
        idivide (int64 (10) ^ 12, int64 (10) ^ 10, 'floor'), int64 (100));
  check ('2^(1/2), the square root of 2', 2 ^ 0.5, '1.414');
end


function check_floor_and_ceiling ()
  printf ('Floor and ceiling\n');
  same ('floor 3.7', int64 (floor (3.7)), int64 (3));
  same ('ceiling 3.2', int64 (ceil (3.2)), int64 (4));
  same ('floor and ceiling of -3.7', as_tuple (int64 (floor (-3.7)), int64 (ceil (-3.7))),
        as_tuple (int64 (-4), int64 (-3)));
  same ('cutting the fraction off -3.7', int64 (fix (-3.7)), int64 (-3));
  same ('floor 4 and ceiling 4', as_tuple (int64 (floor (4)), int64 (ceil (4))), as_tuple (int64 (4), int64 (4)));
  check ('100 / 8', 100 / 8, '12.5');
  same ('bytes for 100 bits, ceiling of 100/8', int64 (ceil (100 / 8)), int64 (13));
  same ('bits to name 8 things, ceiling of log2 8', int64 (ceil (log2 (8))), int64 (3));
  same ('bits to name 9 things, ceiling of log2 9', int64 (ceil (log2 (9))), int64 (4));
  same ('whole 2s that fit in 5, floor of 5/2', int64 (floor (5 / 2)), int64 (2));
end


function check_mod2_and_xor ()
  printf ('Mod-2 arithmetic and XOR\n');
  same ('7 mod 2', mod (int64 (7), int64 (2)), int64 (1));
  same ('13 mod 2', mod (int64 (13), int64 (2)), int64 (1));
  same ('1 + 1', int64 (1) + int64 (1), int64 (2));
  same ('(1 + 1) mod 2', mod (int64 (1) + int64 (1), int64 (2)), int64 (0));
  addition = cell (1, 2);
  for a = 0:1
    addition{a + 1} = mod (int64 (a) + int64 ([0 1]), int64 (2));
  end
  same ('mod-2 addition table', addition, {int64([0 1]), int64([1 0])});
  truth = int64 ([]);
  for a = 0:1
    for b = 0:1
      truth(end + 1) = bitxor (int64 (a), int64 (b));
    end
  end
  same ('XOR truth table', truth, int64 ([0 1 1 0]));
  ok = true;
  for a = 0:1
    for b = 0:1
      ok = ok && bitxor (a, b) == mod (a + b, 2);
    end
  end
  same ('XOR agrees with addition mod 2 on every pair', ok, true);
  ok = true;
  for a = 0:1
    for b = 0:1
      ok = ok && mod (a * b, 2) == bitand (a, b);
    end
  end
  same ('mod-2 multiplication equals AND', ok, true);
  same ('1011 XOR 0110', xor_bits ('1011', '0110'), '1101');
  same ('1101 XOR 0110, the change undone', xor_bits ('1101', '0110'), '1011');
  same ('parity of 1011', parity ('1011'), int64 (1));
  same ('cases in the commutativity proof', int64 (rows (bit_tuples (2))), int64 (4));
  same ('cases in the associativity proof', int64 (rows (bit_tuples (3))), int64 (8));
  vs = strings (4);
  zero = '0000';
  ok = true;
  for i = 1:numel (vs)
    for j = 1:numel (vs)
      ok = ok && strcmp (xor_bits (vs{i}, vs{j}), xor_bits (vs{j}, vs{i}));
    end
  end
  same ('commutative on all 4-bit vectors', ok, true);
  ok = true;
  for i = 1:numel (vs)
    for j = 1:numel (vs)
      for k = 1:numel (vs)
        u = vs{i};
        v = vs{j};
        w = vs{k};
        ok = ok && strcmp (xor_bits (xor_bits (u, v), w), xor_bits (u, xor_bits (v, w)));
      end
    end
  end
  same ('associative on all 4-bit vectors', ok, true);
  ok = true;
  for i = 1:numel (vs)
    ok = ok && strcmp (xor_bits (vs{i}, zero), vs{i});
  end
  same ('zero is the identity', ok, true);
  ok = true;
  for i = 1:numel (vs)
    ok = ok && strcmp (xor_bits (vs{i}, vs{i}), zero);
  end
  same ('every vector is its own inverse', ok, true);
  ok = true;
  for i = 1:numel (vs)
    for j = 1:numel (vs)
      ok = ok && strcmp (xor_bits (xor_bits (vs{i}, vs{j}), vs{j}), vs{i});
    end
  end
  same ('r = m XOR e gives back m = r XOR e', ok, true);
end


function check_sets ()
  printf ('Sets\n');
  same ('the set {0,1}^3', strings (3), {'000', '001', '010', '011', '100', '101', '110', '111'});
  same ('its size', int64 (numel (strings (3))), int64 (8));
  s4 = strings (4);
  same ('{x in {0,1}^4 : exactly two 1s} has', int64 (sum (cellfun (@count_ones, s4) == 2)), int64 (6));
  ok = true;
  for n = 0:12
    ok = ok && int64 (numel (strings (n))) == int64 (2) ^ n;
  end
  same ('|{0,1}^n| = 2^n for n <= 12', ok, true);
  same ('{0,1}^1', strings (1), {'0', '1'});
end


function check_proofs ()
  printf ('The three kinds of proof\n');
  same ('7 x 8', int64 (7) * 8, int64 (56));
  ok = true;
  for n = -100:100
    ok = ok && mod (n * (n + 1), 2) == 0;
  end
  same ('n(n + 1) even for every n from -100 to 100', ok, true);
  same ('1 + 2 + 4 + ... + 128', int64 (sum (int64 (2) .^ int64 (0:7))), int64 (255));
  same ('2^8 - 1', int64 (2) ^ 8 - 1, int64 (255));
  same ('11111111 in binary', int64 (bin2dec ('11111111')), int64 (255));
  % The sums reach 2^64 - 1, which only a uint64 holds. 2^n - 1 is n ones in
  % binary, so it is also the all-ones uint64 shifted right by 64 - n places.
  ok = true;
  for n = 0:64
    s = uint64 (0);
    for i = 0:n - 1
      s = s + bitshift (uint64 (1), i);
    end
    if (n == 0)
      ones_n = uint64 (0);       % Octave's bitshift leaves a shift of 64 places alone
    else
      ones_n = bitshift (intmax ('uint64'), n - 64);
    end
    ok = ok && s == ones_n;
  end
  same ('sum of 2^i for i < n equals 2^n - 1, n <= 64', ok, true);
  same ('induction base, n = 1: 2^0 and 2^1 - 1', as_tuple (int64 (2) ^ 0, int64 (2) ^ 1 - 1),
        as_tuple (int64 (1), int64 (1)));
  row3 = int64 ([]);
  for i = 0:3
    row3(end + 1) = comb (3, i);
  end
  same ('row 3 of the binomials', row3, int64 ([1 3 3 1]));
  same ('its sum, 2^3', int64 (sum (row3)), int64 (8));
  s3 = strings (3);
  ones_in = cellfun (@count_ones, s3);
  grouped = cell (1, 4);
  for i = 0:3
    grouped{i + 1} = s3(ones_in == i);
  end
  same ('3-bit strings grouped by number of 1s', grouped,
        {{'000'}, {'001', '010', '100'}, {'011', '101', '110'}, {'111'}});
  % Row sums reach 2^60, and the entries comb (60, 30) near 1.2e17, both inside int64.
  ok = true;
  for n = 0:60
    total = int64 (0);
    for i = 0:n
      total = total + comb (n, i);
    end
    ok = ok && total == int64 (2) ^ n;
  end
  same ('rows of binomials sum to 2^n, n <= 60', ok, true);
  same ('1101 in binary', int64 (bin2dec ('1101')), int64 (13));
  same ('8 + 4 + 0 + 1', int64 (8 + 4 + 0 + 1), int64 (13));
end


function check_worked_examples ()
  printf ('Worked examples\n');
  same ('at most two flips in a byte, 1 + 8 + 28', comb (8, 0) + comb (8, 1) + comb (8, 2), int64 (37));
  same ('the three terms', int64 ([comb(8, 0), comb(8, 1), comb(8, 2)]), int64 ([1 8 28]));
  check ('share of the 256 patterns', 37 / 256, '14.45%');
  m = '10110010';
  e = '00010100';
  r = xor_bits (m, e);
  same ('10110010 XOR 00010100', r, '10100110');
  same ('10100110 XOR 00010100', xor_bits (r, e), m);
  same ('1s in the error pattern', count_ones (e), int64 (2));
  same ('their positions, counted from the left from 1', int64 (find (e == '1')), int64 ([4 6]));
  same ('1s in m and in r', as_tuple (count_ones (m), count_ones (r)), as_tuple (int64 (4), int64 (4)));
  same ('parity of m and of r', as_tuple (parity (m), parity (r)), as_tuple (int64 (0), int64 (0)));
  same ('address bits for 1000 sensors, ceiling of log2 1000', int64 (ceil (log2 (1000))), int64 (10));
  same ('2^9 and 2^10', as_tuple (int64 (2) ^ 9, int64 (2) ^ 10), as_tuple (int64 (512), int64 (1024)));
  same ('addresses left over', int64 (2) ^ 10 - 1000, int64 (24));
  same ('address bits for 1025 sensors', int64 (ceil (log2 (1025))), int64 (11));
  p = 1e-6;
  n = 8000;                % kept a double, since a double raised to an int64 would round to an integer
  same ('bits in a 1,000-byte packet', int64 (1000) * 8, int64 (n));
  check ('packet arrives with no flip at p = 10^-6', (1 - p) ^ n, '0.9920');
  check ('packet suffers at least one flip', 1 - (1 - p) ^ n, '0.0080');
  check ('damaged packets per 1,000', 1000 * (1 - (1 - p) ^ n), '8');
  check ('exponent of the approximation, -n p', -n * p, '-0.008');
  check ('e^(-np) as the approximation', 2.718281828459045 ^ (-n * p), '0.9920');
  bits = '11010110';
  same ('11010110 in decimal', int64 (sum ((fliplr (bits) - '0') .* 2 .^ (0:numel (bits) - 1))), int64 (214));
  same ('128 + 64 + 16 + 4 + 2', int64 (128 + 64 + 16 + 4 + 2), int64 (214));
  same ('214 mod 2, the last bit', mod (int64 (214), int64 (2)), int64 (0));
  same ('1s in 11010110', count_ones (bits), int64 (5));
  same ('parity of 11010110, unlike its last bit', parity (bits), int64 (1));
  same ('1011 + 0110 with no carries', xor_bits ('1011', '0110'), '1101');
  same ('1011 + 0110 in ordinary binary', dec2bin (bin2dec ('1011') + bin2dec ('0110')), '10001');
  check ('e, the base of the approximation', exp (1), '2.718');
  check ('(1 + 1/1000)^1000, on its way to e', (1 + 1 / 1000) ^ 1000, '2.717');
  check ('e^-1, where (1 - 1/m)^m heads', exp (-1), '0.368');
end


function check_exercises ()
  printf ('Exercise answers\n');
  same ('sum of i^2 for i = 1 to 6', int64 (sum (int64 (1:6) .* int64 (1:6))), int64 (91));
  same ('1 + 4 + 9 + 16 + 25 + 36', int64 (1 + 4 + 9 + 16 + 25 + 36), int64 (91));
  same ('product of 2 for i = 1 to 5', int64 (prod (repmat (int64 (2), 1, 5))), int64 (32));
  row4 = int64 ([]);
  for i = 0:4
    row4(end + 1) = comb (4, i);
  end
  same ('row 4, the starting row', row4, int64 ([1 4 6 4 1]));
  row5 = int64 ([]);
  for i = 0:5
    row5(end + 1) = comb (5, i);
  end
  same ('row 5', row5, int64 ([1 5 10 10 5 1]));
  row6 = int64 (1);
  for k = 1:5
    row6(end + 1) = comb (5, k - 1) + comb (5, k);
  end
  row6(end + 1) = 1;
  same ('row 6 by Pascal''s rule', row6, int64 ([1 6 15 20 15 6 1]));
  total = int64 (0);
  for i = 0:6
    total = total + comb (6, i);
  end
  same ('row 6 sum, 2^6', total, int64 (64));
  same ('11001010 XOR 01100110', xor_bits ('11001010', '01100110'), '10101100');
  same ('10101100 XOR 01100110', xor_bits ('10101100', '01100110'), '11001010');
  check ('1500 / 8', 1500 / 8, '187.5');
  same ('bytes for 1500 bits', int64 (ceil (1500 / 8)), int64 (188));
  same ('XOR tree depth over 64 inputs', int64 (ceil (log2 (64))), int64 (6));
  same ('XOR tree depth over 65 inputs', int64 (ceil (log2 (65))), int64 (7));
  same ('log2 4096', log2 (4096), 12);
  check ('log2 of one million', log2 (10 ^ 6), '19.93');
  same ('bits to label a million things', int64 (ceil (log2 (10 ^ 6))), int64 (20));
  same ('2^19 and 2^20', as_tuple (int64 (2) ^ 19, int64 (2) ^ 20), as_tuple (int64 (524288), int64 (1048576)));
  ok = true;
  for n = -100:100
    ok = ok && mod (n * n, 2) == mod (n, 2);
  end
  same ('n^2 is even exactly when n is, n from -100 to 100', ok, true);
  same ('1 + 2 + ... + 100', int64 (sum (int64 (1:100))), int64 (5050));
  same ('100 x 101 / 2', idivide (int64 (100) * 101, int64 (2), 'floor'), int64 (5050));
  ok = true;
  for n = 0:500
    ok = ok && int64 (sum (int64 (1:n))) == idivide (int64 (n) * (int64 (n) + 1), int64 (2), 'floor');
  end
  same ('1 + ... + n = n(n + 1)/2 for n <= 500', ok, true);
  same ('8 choose 3 and 8 choose 5', as_tuple (comb (8, 3), comb (8, 5)), as_tuple (int64 (56), int64 (56)));
  same ('8 x 7 x 6 / (3 x 2 x 1)', idivide (int64 (8) * 7 * 6, int64 (3 * 2 * 1), 'floor'), int64 (56));
  ok = true;
  for n = 0:40
    for k = 0:n
      ok = ok && comb (n, k) == comb (n, n - k);
    end
  end
  same ('n choose k = n choose n-k for n <= 40', ok, true);
  s3 = strings (3);
  odd3 = s3(cellfun (@parity, s3) == 1);
  same ('3-bit strings with an odd number of 1s', odd3, {'001', '010', '100', '111'});
  same ('their count', int64 (numel (odd3)), int64 (4));
  check ('1000 bits with no flip at p = 0.001', 0.999 ^ 1000, '0.3677');
  check ('e^-1, the approximation', 2.718281828459045 ^ -1, '0.3679');
  check ('log of 0.999 times 1000', 1000 * log (0.999), '-1.0005');
  check ('n p for 1000 bits at 0.001', 1000 * 0.001, '1');
  check ('gap between e^-1 and the exact value', 2.718281828459045 ^ -1 - 0.999 ^ 1000, '0.0002');
end


function check_looking_ahead ()
  printf ('Looking ahead\n');
  same ('01100100 in decimal', int64 (bin2dec ('01100100')), int64 (100));
  same ('its high bit flipped, 11100100', int64 (bin2dec ('11100100')), int64 (228));
  same ('the flip as an XOR with 10000000', xor_bits ('01100100', '10000000'), '11100100');
end


function main ()
  check_book ();
  check_sums_and_products ();
  check_binomials ();
  check_powers_and_logs ();
  check_floor_and_ceiling ();
  check_mod2_and_xor ();
  check_sets ();
  check_proofs ();
  check_worked_examples ();
  check_exercises ();
  check_looking_ahead ();
  printf ('\nAll chapter f1 numbers reproduced.\n');
end


main ();
