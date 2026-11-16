% Reproduce every number in chapter 28 of Error Correction Coding,
% "Verifying ECC Hardware", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch28_verifying_ecc_hardware.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch28_verifying_ecc_hardware.py and prints the
% same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Whole numbers, bit patterns and decoder outputs are compared
% exactly, and every fractional value is printed with a fixed number of
% places, so its comparison is exact too. Numbers too large for a double to
% hold as an integer in the usual way are powers of two, or a small number
% times one, which a double stores exactly, and they are printed with %.0f.
%
% The running design is the Hsiao (13,8) SEC-DED decoder. Its check matrix
% is built here the way python/ch20_secded_and_hsiao_codes.py builds it, as
% the 3-of-5 columns in combination order with two dropped, and compared with
% the columns the chapter prints. The decoder model plants the same four bugs
% the SystemVerilog decoder plants through its MUTANT parameter, and it works
% on a whole column of received words at once, so every count over the
% 23,552 envelope cases takes one call. The random campaigns use PyRandom
% from octave/common, which draws the numbers Python's random.Random(28)
% draws.
%
% A column or a syndrome of r rows is an integer whose bit i is row i,
% printed with row 0 first. A data word is printed as an 8-bit binary number
% with d7 on the left. Codeword bits 0 to 7 are d0 to d7 and bits 8 to 12 are
% c0 to c4. A case index k runs from 0 to 23,551 and stands for data word
% floor(k / 92) with envelope pattern mod(k, 92).

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

function s = big (x)
  % An exactly stored integer, printed in full.
  s = sprintf ('%.0f', x);
end

function s = joined (v)
  s = sprintf ('%d,', v);
  s = s(1:end-1);
end


% The code ---------------------------------------------------------------------

function w = weight (v, width)
  % The number of 1 bits in each entry of v, over the low width bits.
  w = zeros (size (v));
  for i = 1:width
    w = w + bitget (v, i);
  end
end

function s = bits (v, r)
  % A column or syndrome as a string, row 0 first.
  s = char ('0' + bitget (v, 1:r));
end

function s = word (v)
  % A data word as an 8-bit binary number, d7 on the left.
  s = dec2bin (v, 8);
end

function cols = hsiao_13_8 ()
  % The 3-of-5 columns in combination order, less two, then the unit columns.
  triples = nchoosek (0:4, 3);
  cols = [];
  for t = 1:rows (triples)
    c = sum (2 .^ triples(t, :));
    if (c ~= sum (2 .^ [0 1 2]) && c ~= sum (2 .^ [0 3 4]))
      cols(end+1) = c;
    end
  end
  cols = [cols, 2 .^ (0:4)];
end

function cw = encode (d, cols)
  % The codewords of a column of data words.
  check = zeros (size (d));
  for j = 1:8
    check = bitxor (check, bitget (d, j) .* cols(j));
  end
  cw = d + check * 256;
end

function p = envelope ()
  % The 92 patterns of weight at most 2, in the order the Python script lists them.
  pairs = nchoosek (0:12, 2);
  p = [0, 2 .^ (0:12), (2 .^ pairs(:, 1) + 2 .^ pairs(:, 2))'];
end

function [data, corr, unc] = decode (r, mutant, cols)
  % (data, corrected, uncorrectable) for a column of received words, one bug planted.
  dcols = cols;
  if (mutant == 1)
    dcols(7) = bitxor (dcols(7), 2);             % received bit 6 left out of row 1
  end
  s = zeros (size (r));
  for j = 1:13
    s = bitxor (s, bitget (r, j) .* dcols(j));
  end
  pos = -ones (size (r));
  for j = 1:13
    pos(s == cols(j)) = j - 1;
  end
  target = pos;
  target(pos >= 8) = -1;
  if (mutant == 3)
    target(pos == 5) = 4;                        % d5's syndrome flips d4
  end
  if (mutant == 4)
    target(pos == 0 & bitget (r, 8) == 1) = -1;  % d0's flip suppressed while bit 7 is 1
  end
  flip = zeros (size (r));
  hit = target >= 0;
  flip(hit) = 2 .^ target(hit);
  data = bitxor (bitand (r, 255), flip);
  corr = double (pos >= 0);
  unc = double (s ~= 0 & pos < 0);
  if (mutant == 2)
    unc(weight (s, 5) ~= 2) = 0;                 % raised only for two-1 syndromes
  end
end

function [wrong, silent] = judge (d, e, mutant, cols)
  % Columns of (wrong, silent) for columns of cases (d, e) inside the envelope.
  [data, corr, unc] = decode (bitxor (encode (d, cols), e), mutant, cols);
  w = weight (e, 13);
  want = d;
  want(w == 2) = bitxor (d(w == 2), bitand (e(w == 2), 255));
  wrong = data ~= want | corr ~= (w == 1) | unc ~= (w == 2);
  silent = wrong & data ~= d & unc == 0;
end

function [d, e] = all_cases ()
  pats = envelope ();
  k = (0:23551)';
  d = floor (k / 92);
  e = pats(mod (k, 92) + 1)';
end

function n = envelope_size (len)
  n = 1 + len + nchoosek (len, 2);
end

function h = harmonic (n)
  h = 0;
  for k = 1:n
    h = h + 1 / k;
  end
end

function t = campaign_size (m, n)
  t = ceil (log (1e-6) / log (1 - m / n));
end


% The checks -------------------------------------------------------------------

COLS = hsiao_13_8 ();
RATE = 1e9;                  % assumed checks per second
YEAR = 31557600;             % seconds in a year of 365.25 days
[D, E] = all_cases ();
PATS = envelope ();

section ('The Hsiao (13,8) code');
c = cell (1, 8);
for j = 1:8
  c{j} = bits (COLS(j), 5);
end
same ('data columns d0..d7, row 0 first', strjoin (c, ','), '11010,11001,10110,10101,01110,01101,01011,00111');
c = cell (1, 5);
for j = 1:5
  c{j} = bits (COLS(8 + j), 5);
end
same ('check columns c0..c4, row 0 first', strjoin (c, ','), '10000,01000,00100,00010,00001');
require ('every column has odd weight', all (mod (weight (COLS, 5), 2) == 1));
distance = 0;
for w = 1:4
  subsets = nchoosek (1:13, w);
  s = zeros (rows (subsets), 1);
  for t = 1:w
    s = bitxor (s, COLS(subsets(:, t))');
  end
  if (any (s == 0))
    distance = w;
    break;
  end
end
same ('minimum distance (fewest columns summing to zero)', distance, 4);
d = 105;
cw = encode (d, COLS);
same ('codeword of 01101001: check bits c0..c4', bits (floor (cw / 256), 5), '01001');
syn = @(r) bitxor_all (COLS(logical (bitget (r, 1:13))));
function s = bitxor_all (v)
  s = 0;
  for i = 1:numel (v)
    s = bitxor (s, v(i));
  end
end
[o1, o2, o3] = decode (bitxor (cw, 4), 0, COLS);
same ('single error in d2: syndrome', bits (syn (bitxor (cw, 4)), 5), '10110');
same ('single error in d2: (data, corrected, uncorrectable)', sprintf ('%s,%d,%d', word (o1), o2, o3), '01101001,1,0');
r = bitxor (bitxor (cw, 4), 512);
[o1, o2, o3] = decode (r, 0, COLS);
same ('double error in d2 and c1: syndrome', bits (syn (r), 5), '11110');
same ('double error in d2 and c1: (data, corrected, uncorr.)', sprintf ('%s,%d,%d', word (o1), o2, o3), '01101101,0,1');
[o1, o2, o3] = decode (bitxor (cw, 1024), 0, COLS);
same ('single error in c2: (data, corrected, uncorrectable)', sprintf ('%s,%d,%d', word (o1), o2, o3), '01101001,1,0');

section ('Guarantees and the design envelope');
same ('patterns of weight 0, 1, 2', sprintf ('%d,%d,%d', nchoosek (13, 0), nchoosek (13, 1), nchoosek (13, 2)), '1,13,78');
same ('patterns in the envelope', numel (PATS), 92);
same ('patterns outside the envelope', 2 ^ 13 - 92, 8100);
counts = zeros (1, 11);
for w = 3:13
  counts(w - 2) = nchoosek (13, w);
end
same ('patterns of weight 3 to 13', joined (counts), '286,715,1287,1716,1716,1287,715,286,78,13,1');
require ('the 13 single-error syndromes are distinct and odd', numel (unique (COLS)) == 13 && all (mod (weight (COLS, 5), 2) == 1));
pairs = nchoosek (1:13, 2);
dsyn = bitxor (COLS(pairs(:, 1)), COLS(pairs(:, 2)));
same ('double-error syndromes of weight 2 and of weight 4', sprintf ('%d,%d', sum (weight (dsyn, 5) == 2), sum (weight (dsyn, 5) == 4)), '52,26');
groups = {'check pairs', 9:13, 9:13, '10,0'; 'data-and-check pairs', 1:8, 9:13, '24,16'; 'data pairs', 1:8, 1:8, '18,10'};
for g = 1:rows (groups)
  [I, J] = ndgrid (groups{g, 2}, groups{g, 3});
  keep = I(:) < J(:);
  sp = bitxor (COLS(I(keep)), COLS(J(keep)));
  same (sprintf ('%s: syndromes of weight 2 and of weight 4', groups{g, 1}), sprintf ('%d,%d', sum (weight (sp, 5) == 2), sum (weight (sp, 5) == 4)), groups{g, 4});
end
left_out = zeros (1, 5);
for row = 1:5
  left_out(row) = sum (bitget (COLS(1:8), row) == 0);
end
same ('data columns leaving out rows 0 to 4', joined (left_out), '4,3,3,3,3');
same ('data pairs sharing two rows, sum of C(left out, 2)', sum (left_out .* (left_out - 1) / 2), 18);
[wr, si] = judge (D, E, 0, COLS);
require ('the correct decoder meets the guarantees on all 23,552 cases', sum (wr) == 0 && sum (si) == 0);
trip = nchoosek (0:12, 3);
[o1, o2, o3] = decode (sum (2 .^ trip, 2), 0, COLS);
same ('triples: total, miscorrected, flagged', sprintf ('%d,%d,%d', nchoosek (13, 3), sum (o2), sum (o3)), '286,220,66');
same ('triples: percent miscorrected', sprintf ('%.1f', 100 * sum (o2) / nchoosek (13, 3)), '76.9');

section ('Counting the cases');
same ('(13,8): data words', 256, 256);
same ('(13,8): error patterns', 2 ^ 13, 8192);
same ('(13,8): envelope cases, 256 x 92', 256 * 92, 23552);
same ('(13,8): all cases, 256 x 8192', 256 * 2 ^ 13, 2097152);
same ('(13,8): input combinations of the formal wrapper, 2^21', 2 ^ 21, 2 ^ 21);
same ('(72,64): envelope patterns 1 + 72 + 2556', sprintf ('%d=%d+%d+%d', 1 + 72 + nchoosek (72, 2), 1, 72, nchoosek (72, 2)), '2629=1+72+2556');
same ('(72,64): data words', big (2 ^ 64), '18446744073709551616');
same ('(72,64): error patterns', big (2 ^ 72), '4722366482869645213696');
same ('(72,64): envelope cases, 2629 x 2^64', big (2629 * 2 ^ 64), '48496490169782411198464');
same ('(72,64): all cases, 2^136', big (2 ^ 136), '87112285931760246646623899502532662132736');
section ('Times at an assumed 10^9 checks per second');
same ('(13,8) envelope, microseconds', sprintf ('%.1f', 256 * 92 * 1e6 / RATE), '23.6');
same ('(13,8) all cases, milliseconds', sprintf ('%.1f', 256 * 2 ^ 13 * 1e3 / RATE), '2.1');
same ('(72,64) envelope on one codeword, microseconds', sprintf ('%.1f', 2629 * 1e6 / RATE), '2.6');
same ('(72,64) all patterns on one codeword, whole years', floor (2 ^ 72 / (RATE * YEAR)), 149642);
same ('(72,64) envelope cases, millions of years', sprintf ('%.2f', 2629 * 2 ^ 64 / (RATE * YEAR) / 1e6), '1.54');
same ('(72,64) all cases, years', sprintf ('%.2e', 2 ^ 136 / (RATE * YEAR)), '2.76e+24');
section ('Envelope and all-pattern sweeps for four Hsiao codes');
codes = {[13 8], 92, '23552', '2097152'; ...
         [22 16], 254, '16646144', '274877906944'; ...
         [39 32], 781, '3354369458176', '2361183241434822606848'; ...
         [72 64], 2629, '48496490169782411198464', '87112285931760246646623899502532662132736'};
for i = 1:rows (codes)
  n = codes{i, 1}(1);
  k = codes{i, 1}(2);
  same (sprintf ('(%d,%d): envelope patterns', n, k), envelope_size (n), codes{i, 2});
  same (sprintf ('(%d,%d): envelope cases', n, k), big (2 ^ k * envelope_size (n)), codes{i, 3});
  same (sprintf ('(%d,%d): all cases', n, k), big (2 ^ k * 2 ^ n), codes{i, 4});
end
same ('(39,32) envelope sweep at 10^9 per second, minutes', sprintf ('%.1f', 2 ^ 32 * 781 / RATE / 60), '55.9');

section ('The reduction theorem, checked by brute force');
CW = encode (D, COLS);
pidx = mod ((0:23551)', 92) + 1;
broken = zeros (1, 5);
for mutant = 0:4
  [z1, z2, z3] = decode (PATS', mutant, COLS);
  [o1, o2, o3] = decode (bitxor (CW, E), mutant, COLS);
  broken(mutant + 1) = sum (o1 ~= bitxor (D, z1(pidx)) | o2 ~= z2(pidx) | o3 ~= z3(pidx));
end
for mutant = [0 2 3]
  require (sprintf ('mutant %d: output on c + e is d + output on e, every case', mutant), broken(mutant + 1) == 0);
end
same ('mutant 1: cases where that identity fails', broken(2), 9856);
same ('mutant 4: cases where that identity fails', broken(5), 128);
[o1, o2, o3] = decode (encode ((0:255)', COLS), 1, COLS);
require ('mutant 1: some codeword has a nonzero syndrome', any (o1 ~= (0:255)' | o2 ~= 0 | o3 ~= 0));
zero_stated = [0 12 26 1 0];
for mutant = 0:4
  wr = judge (zeros (92, 1), PATS', mutant, COLS);
  same (sprintf ('mutant %d: patterns wrong on the zero codeword', mutant), sum (wr), zero_stated(mutant + 1));
end
flips = PATS(bitget (PATS, 7) == 1);
wr = judge (zeros (numel (flips), 1), flips', 1, COLS);
same ('mutant 1, zero codeword: patterns flipping d6, wrong ones', sprintf ('%d,%d', numel (flips), sum (wr)), '13,12');
wr = judge (zeros (92, 1), PATS', 1, COLS);
require ('mutant 1, zero codeword: every wrong pattern flips d6', all (bitget (PATS(wr), 7) == 1));
MCOLS = COLS;
MCOLS(7) = bitxor (MCOLS(7), 2);
msyn = @(r) bitxor_all (MCOLS(logical (bitget (r, 1:13))));
e = 64 + 256;
[o1, o2, o3] = decode (e, 1, COLS);
same ('mutant 1, zero codeword, d6 and c0: mutant syndrome, output', sprintf ('%s,%s,%d,%d', bits (msyn (e), 5), word (o1), o2, o3), '10011,01000000,0,1');

section ('Coupon collector and miss probabilities');
same ('H_92', sprintf ('%.4f', harmonic (92)), '5.1044');
same ('92 H_92, expected draws to see all 92 patterns', sprintf ('%.1f', 92 * harmonic (92)), '469.6');
same ('H_2629', sprintf ('%.4f', harmonic (2629)), '8.4518');
same ('2629 H_2629, expected draws to see all 2629 patterns', sprintf ('%.1f', 2629 * harmonic (2629)), '22219.7');
same ('13 H_13, expected draws to see all 13 singles', sprintf ('%.1f', 13 * harmonic (13)), '41.3');
same ('T for miss <= 10^-6, m = 1, N = 92', campaign_size (1, 92), 1265);
same ('T for miss <= 10^-6, m = 1, N = 2629', campaign_size (1, 2629), 36315);
same ('T for miss <= 10^-6, m = 1, N = 23552', campaign_size (1, 23552), 325376);
same ('approximation N ln(10^6), N = 92', sprintf ('%.0f', 92 * log (1e6)), '1271');
same ('approximation N ln(10^6), N = 2629', sprintf ('%.0f', 2629 * log (1e6)), '36321');
same ('T for miss <= 10^-6, m = 26, N = 92', campaign_size (26, 92), 42);
same ('T for miss <= 10^-6, m = 1, N = 13', campaign_size (1, 13), 173);
same ('one-case campaign over the 23,552-case sweep', sprintf ('%.1f', campaign_size (1, 23552) / 23552), '13.8');

section ('The planted bugs over all 23,552 envelope cases');
stated = [0 0; 9856 6528; 6656 6656; 256 256; 128 128];
WRONG = false (23552, 5);
M = zeros (1, 5);
for mutant = 0:4
  [wr, si] = judge (D, E, mutant, COLS);
  WRONG(:, mutant + 1) = wr;
  M(mutant + 1) = sum (wr);
  same (sprintf ('mutant %d: wrong cases m, silent cases', mutant), sprintf ('%d,%d', sum (wr), sum (si)), sprintf ('%d,%d', stated(mutant + 1, 1), stated(mutant + 1, 2)));
end
w0 = judge (zeros (92, 1), PATS', 1, COLS);
w1 = judge (64 * ones (92, 1), PATS', 1, COLS);
same ('mutant 1: wrong patterns per word, d6 = 0 and d6 = 1', sprintf ('%d,%d', sum (w0), sum (w1)), '12,65');
rowsm = {1, '0.4185', '0.000', 26; 2, '0.2826', '0.000', 42; 3, '0.0109', '0.335', 1265; 4, '0.0054', '0.580', 2536};
for i = 1:rows (rowsm)
  mutant = rowsm{i, 1};
  m = M(mutant + 1);
  p = m / 23552;
  same (sprintf ('mutant %d: chance one draw hits it', mutant), sprintf ('%.4f', p), rowsm{i, 2});
  same (sprintf ('mutant %d: miss probability at T = 100', mutant), sprintf ('%.3f', (1 - p) ^ 100), rowsm{i, 3});
  same (sprintf ('mutant %d: T for miss <= 10^-6', mutant), campaign_size (m, 23552), rowsm{i, 4});
end
m4 = M(5);
require ('mutant 4: m divides N', mod (23552, m4) == 0);
same ('mutant 4: m/N as 1/x', sprintf ('1/%d', 23552 / m4), '1/184');
same ('mutant 4: sweep over campaign, N / T', sprintf ('%.2f', 23552 / campaign_size (m4, 23552)), '9.29');
narrow = weight (E, 13) < 2;
same ('mutant 2: wrong cases with wt(e) < 2', sum (WRONG(:, 3) & narrow), 0);
same ('mutant 3: wrong cases with wt(e) < 2', sum (WRONG(:, 4) & narrow), 256);
same ('cases with wt(e) < 2', 256 * 14, 3584);

section ('Seeded campaigns, seed 28');
rng = PyRandom (28);
first = floor (rng.random (100) * 23552);
same ('first draw as (data word, pattern index)', sprintf ('%s,%d', word (floor (first(1) / 92)), mod (first(1), 92)), '00011100,84');
hits = zeros (1, 4);
for mutant = 1:4
  hits(mutant) = sum (WRONG(first + 1, mutant + 1));
end
same ('one campaign of 100: hits on mutants 1, 2, 3, 4', joined (hits), '33,29,1,0');
same ('one campaign of 100: distinct patterns, distinct words', sprintf ('%d,%d', numel (unique (mod (first, 92))), numel (unique (floor (first / 92)))), '59,77');
rng = PyRandom (28);
all_draws = reshape (floor (rng.random (100 * 10000) * 23552), 100, 10000);
missed = zeros (1, 4);
for mutant = 1:4
  caught = any (reshape (WRONG(all_draws(:) + 1, mutant + 1), 100, 10000), 1);
  missed(mutant) = sum (~caught);
end
same ('10,000 campaigns of 100: missed mutants 1, 2, 3, 4', joined (missed), '0,0,3401,5883');
same ('mutant 3: fraction missed', sprintf ('%.4f', missed(3) / 10000), '0.3401');
same ('mutant 4: fraction missed', sprintf ('%.4f', missed(4) / 10000), '0.5883');
rng = PyRandom (28);
k = floor (rng.random (700) * 23552);
seen = false (1, 92);
all_at = 0;
for t = 1:700
  seen(mod (k(t), 92) + 1) = true;
  if (all_at == 0 && all (seen))
    all_at = t;
  end
end
same ('coupon campaign, seed 28: draw that completes all 92 patterns', all_at, 454);

section ('The formal counterexamples the chapter reads');
d = 188;
e = 32;
[o1, o2, o3] = decode (bitxor (encode (d, COLS), e), 3, COLS);
same ('mutant 3 at d = 10111100, e = d5: decoder output', sprintf ('%s,%d,%d', word (o1), o2, o3), '10001100,1,0');
[wr, si] = judge (d, e, 3, COLS);
same ('mutant 3 at that case: wrong, silent', sprintf ('%d,%d', wr, si), '1,1');
same ('syndrome of that received word', bits (syn (bitxor (encode (d, COLS), e)), 5), '01101');
cases = {1, 252, 3, '0000000000011', '00011', '10111111,1,0', '11111111,0,1'; ...
         2, 64, 4 + 512, '0001000000100', '11110', '01000100,0,0', '01000100,0,1'; ...
         4, 153, 1, '0000000000001', '11010', '10011000,1,0', '10011001,1,0'};
for i = 1:rows (cases)
  mutant = cases{i, 1};
  d = cases{i, 2};
  e = cases{i, 3};
  r = bitxor (encode (d, COLS), e);
  [o1, o2, o3] = decode (r, mutant, COLS);
  w = weight (e, 13);
  want = d;
  if (w == 2)
    want = bitxor (d, bitand (e, 255));
  end
  same (sprintf ('mutant %d at d = %s: pattern, position 12 first', mutant, word (d)), dec2bin (e, 13), cases{i, 4});
  same (sprintf ('mutant %d at that case: correct syndrome', mutant), bits (syn (r), 5), cases{i, 5});
  same (sprintf ('mutant %d at that case: decoder output', mutant), sprintf ('%s,%d,%d', word (o1), o2, o3), cases{i, 6});
  same (sprintf ('mutant %d at that case: the guarantee', mutant), sprintf ('%s,%d,%d', word (want), w == 1, w == 2), cases{i, 7});
  [wr, si] = judge (d, e, mutant, COLS);
  same (sprintf ('mutant %d at that case: wrong, silent', mutant), sprintf ('%d,%d', wr, si), '1,1');
end
r = bitxor (encode (252, COLS), 3);
same ('mutant 1 at d = 11111100: mutant syndrome', bits (msyn (r), 5), '01011');

printf ('\nAll chapter 28 numbers reproduced.\n');
