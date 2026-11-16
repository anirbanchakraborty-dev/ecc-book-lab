% Reproduce every number in chapter 29 of Error Correction Coding,
% "Case Study: ECC in Processor Caches", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch29_case_study_ecc_in_processor_caches.m
%
% It needs only core Octave, with no pkg load and no Forge packages, and it
% reads nothing from the network. It is the Octave twin of
% python/ch29_case_study_ecc_in_processor_caches.py and prints the same
% lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Whole numbers, bit patterns and hexadecimal words are compared
% exactly, and every fractional value is printed in one fixed format, exponent
% form with four decimals or with the one decimal the chapter states, or a
% decimal or a percentage to a stated number of places, so the comparison of
% the printed text is exact too.
%
% The open design the chapter studies is the SEC-DED generator of the
% OpenTitan project, read at commit 4ced702ca207ad420fbae5226ac8685a17ebef7b.
% The files read are listed, with their addresses, at the top of the Python
% script. The list of generated codes and the encoder masks below were
% transcribed by hand from util/design/data/secded_cfg.hjson and from
% hw/ip/prim/rtl/prim_secded_39_32_enc.sv and prim_secded_72_64_enc.sv at
% that commit, and every property the chapter claims for them is recomputed
% here from the transcription.
%
% A column of a check matrix with r rows is held as a whole number whose bit
% i is the entry in row i, and a bit pattern is printed with row 0 first. A
% 64-bit mask or data word is held as a pair of 32-bit halves, high half
% first, and printed as sixteen hexadecimal digits, as the Python script
% holds and prints it. Positions are counted from 0, as in the chapter, so
% position p is element p + 1 here. The book's own Hsiao matrices are built
% exactly as ch20_secded_and_hsiao_codes.m builds them, and Octave's nchoosek
% lists combinations in the same order as Python's itertools.combinations.

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

function w = weight (v)
  % The number of ones in each entry of v.
  w = zeros (size (v));
  while (any (v(:)))
    w = w + bitand (v, 1);
    v = bitshift (v, -1);
  end
end

function v = col (rows)
  v = sum (2 .^ rows);
end

function v = cols_of (C)
  v = sum (2 .^ C, 2)';
end

function s = bits (v, r)
  % A column or syndrome as text, row 0 first.
  s = char ('0' + bitand (bitshift (v, -(0:r - 1)), 1));
end

function s = joined (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ',');
end

function s = hex64 (h)
  % A 64-bit word held as [high low] halves, as sixteen hexadecimal digits.
  s = [dec2hex(h(1), 8), dec2hex(h(2), 8)];
end

function s = e4 (x)
  s = sprintf ('%.4e', x);
end

function d = levels (w)
  % Levels of two-input gates in a balanced tree over w signals.
  d = 0;
  reach = 1;
  while (reach < w)
    reach = reach * 2;
    d = d + 1;
  end
end

function loads = row_loads (cols, r)
  loads = zeros (1, r);
  for i = 0:r - 1
    loads(i + 1) = sum (bitand (bitshift (cols, -i), 1));
  end
end

function s = rows_of (v, r)
  idx = find (bitand (bitshift (v, -(0:r - 1)), 1)) - 1;
  s = ['(' joined(idx) ')'];
end

function b = word_bits (h)
  % The 64 bits of a word held as [high low], bit j at element j + 1.
  b = [bitand(bitshift(h(2), -(0:31)), 1), bitand(bitshift(h(1), -(0:31)), 1)];
end


% OpenTitan data ---------------------------------------------------------------
% The 22 entries of util/design/data/secded_cfg.hjson at commit 4ced702ca207,
% in the order the file lists them, as {code_type, k, m}.

function c = config ()
  c = {'hsiao', 16, 6; 'hsiao', 22, 6; 'hsiao', 32, 7; 'hsiao', 57, 7; ...
       'hsiao', 64, 8; 'hsiao', 128, 9; ...
       'hamming', 16, 6; 'hamming', 32, 7; 'hamming', 64, 8; 'hamming', 68, 8; ...
       'hamming', 128, 9; ...
       'inv_hsiao', 16, 6; 'inv_hsiao', 22, 6; 'inv_hsiao', 32, 7; 'inv_hsiao', 57, 7; ...
       'inv_hsiao', 64, 8; 'inv_hsiao', 128, 9; ...
       'inv_hamming', 16, 6; 'inv_hamming', 32, 7; 'inv_hamming', 64, 8; ...
       'inv_hamming', 68, 8; 'inv_hamming', 128, 9};
end

function s = infix (kind)
  % The module-name infix secded_gen.py gives each code type (CODE_OPTIONS).
  switch (kind)
    case 'hsiao'
      s = '';
    case 'inv_hsiao'
      s = '_inv';
    case 'hamming'
      s = '_hamming';
    otherwise
      s = '_inv_hamming';
  end
end

function m = masks_39 ()
  % Encoder masks of hw/ip/prim/rtl/prim_secded_39_32_enc.sv at commit
  % 4ced702ca207: check bit 32 + i is the XOR of the data bits mask i selects.
  % prim_secded_inv_39_32_enc.sv carries the same seven masks.
  m = hex2dec ({'2606BD25', 'DEBA8050', '413D89AA', '31234ED1', 'C2C1323B', ...
                '2DCC624C', '98505586'})';
end

function m = masks_72 ()
  % Encoder masks of hw/ip/prim/rtl/prim_secded_72_64_enc.sv at commit
  % 4ced702ca207, one row per mask as [high 32 bits, low 32 bits]: check bit
  % 64 + i is the XOR of the data bits mask i selects.
  % prim_secded_inv_72_64_enc.sv carries the same eight.
  m = [hex2dec({'B9000000'; '5E00000F'; '67003FF0'; 'CD0FC0F0'; 'B671C711'; ...
                'B5B65926'; 'CBDAAA4A'; '7AED348D'}), ...
       hex2dec({'001FFFFF'; 'FFE0003F'; '03E007C1'; '3C207842'; 'C4438884'; ...
                '488C9108'; '91152210'; '221A4420'})];
end

% The constants the inverted decoders XOR into the received word before the
% syndrome, from prim_secded_inv_39_32_dec.sv (39'h2A00000000) and
% prim_secded_inv_72_64_dec.sv (72'hAA0000000000000000), as check-bit patterns.
function v = inv_check_39 ()
  v = hex2dec ('2A');
end

function v = inv_check_72 ()
  v = hex2dec ('AA');
end

function cols = ot_columns_39 ()
  M = masks_39 ();
  r = numel (M);
  data = zeros (1, 32);
  for j = 0:31
    data(j + 1) = sum (bitand (bitshift (M, -j), 1) .* 2 .^ (0:r - 1));
  end
  cols = [data, 2 .^ (0:r - 1)];
end

function cols = ot_columns_72 ()
  M = masks_72 ();
  r = rows (M);
  data = zeros (1, 64);
  for j = 0:63
    if (j < 32)
      b = bitand (bitshift (M(:, 2)', -j), 1);
    else
      b = bitand (bitshift (M(:, 1)', -(j - 32)), 1);
    end
    data(j + 1) = sum (b .* 2 .^ (0:r - 1));
  end
  cols = [data, 2 .^ (0:r - 1)];
end

function h = hsiao_39_32 ()
  data = cols_of (nchoosek (0:6, 3));
  data = data(~ismember (data, [col([0 1 2]), col([0 3 4]), col([1 5 6])]));
  h = [data, 2 .^ (0:6)];
end

function h = hsiao_72_64 ()
  runs = zeros (1, 8);
  for i = 0:7
    runs(i + 1) = col (mod (i + (0:4), 8));
  end
  h = [cols_of(nchoosek (0:7, 3)), runs, 2 .^ (0:7)];
end

function f = ideal_fanin (k, m)
  % The largest data fan-in of a check bit that balanced odd columns allow,
  % computed as ideal_fanin() in secded_gen.py computes it.
  f = 0;
  needed = k;
  for w = 3:2:m
    avail = nchoosek (m, w);
    take = min (avail, needed);
    f = f + floor ((take * w + m - 1) / m);
    needed = needed - take;
    if (needed == 0)
      break;
    end
  end
end

function ok = distance_at_least_4 (cols)
  ok = false;
  if (any (cols == 0) || numel (unique (cols)) ~= numel (cols))
    return;
  end
  P = nchoosek (1:numel (cols), 2);
  x = bitxor (cols(P(:, 1)), cols(P(:, 2)));
  ok = ~any (ismember (x, cols));
end

function a4 = weight4_codewords (cols)
  % Codewords of weight 4: pairs of column pairs with the same sum, over 3.
  P = nchoosek (1:numel (cols), 2);
  x = bitxor (cols(P(:, 1)), cols(P(:, 2)));
  c = accumarray (x(:) + 1, 1);
  a4 = sum (c .* (c - 1) / 2) / 3;
end

function s = intersections (cols)
  % How many rows each pair of the given columns shares, as size:count.
  P = nchoosek (1:numel (cols), 2);
  w = weight (bitand (cols(P(:, 1)), cols(P(:, 2))));
  sizes = unique (w);
  parts = cell (1, numel (sizes));
  for i = 1:numel (sizes)
    parts{i} = sprintf ('%d:%d', sizes(i), sum (w == sizes(i)));
  end
  s = strjoin (parts, ' ');
end

function s = syndrome_class (v)
  if (v == 0)
    s = 'none';
  elseif (mod (weight (v), 2) == 1)
    s = 'single';
  else
    s = 'double';
  end
end

function c = check_of (cols, b)
  % Check bits for the data bits b, a row of 0s and 1s, bit j at element j + 1.
  c = 0;
  for j = find (b)
    c = bitxor (c, cols(j));
  end
end


% Checks -----------------------------------------------------------------------

function check_config ()
  section ('The codes secded_cfg.hjson asks the generator for, at commit 4ced702ca207');
  C = config ();
  for e = 1:rows (C)
    kind = C{e, 1};
    k = C{e, 2};
    m = C{e, 3};
    n = k + m;
    name = sprintf ('prim_secded%s_%d_%d', infix (kind), n, k);
    printf ('  %-34s %-12s (%d,%d)\n', name, kind, n, k);
  end
  kinds = C(:, 1);
  K = cell2mat (C(:, 2))';
  M = cell2mat (C(:, 3))';
  same ('entries in the configuration', rows (C), 22);
  same ('Hsiao codes', sum (strcmp (kinds, 'hsiao')), 6);
  same ('extended Hamming codes', sum (strcmp (kinds, 'hamming')), 5);
  same ('inverted Hsiao codes', sum (strcmp (kinds, 'inv_hsiao')), 6);
  same ('inverted extended Hamming codes', sum (strcmp (kinds, 'inv_hamming')), 5);
  h = strcmp (kinds, 'hsiao')';
  same ('Hsiao lengths n', joined (K(h) + M(h)), '22,28,39,64,72,137');
  g = strcmp (kinds, 'hamming')';
  same ('extended Hamming lengths n', joined (K(g) + M(g)), '22,39,72,76,137');
  fewest = all (2 .^ (M - 1) >= K + M & 2 .^ (M - 2) < K + M - 1);
  require ('every entry uses the fewest check bits the SEC-DED bound allows', fewest);
  same ('check bits of the (64,57) code', M(find (h & K == 57, 1)), 7);
  same ('odd-weight 7-bit columns there are', sum (mod (weight (1:127), 2) == 1), 64);
  same ('so the (64,57) code uses every one of them', 57 + 7, 64);
end

function check_matrix (label, cols, r, k, st)
  section (sprintf ('The OpenTitan Hsiao %s check matrix, from the transcribed masks', label));
  data = cols(1:k);
  require ('every column has odd weight', all (mod (weight (cols), 2) == 1));
  require ('the columns are distinct', numel (unique (cols)) == numel (cols));
  require ('no column is the sum of two others, so d >= 4', distance_at_least_4 (cols));
  same ('data column weights present', joined (unique (weight (data))), st.wts);
  loads = row_loads (cols, r);
  same ('ones in each row', joined (loads), st.rows);
  same ('ones in the whole matrix', sum (loads), st.ones);
  same ('XOR gates for the syndrome, one fewer than each row', sum (loads - 1), st.syn);
  same ('XOR gates for the encoder, data ones less one per row', sum (loads - 2), st.enc);
  same ('levels of two-input gates for the longest row', levels (max (loads)), st.levels);
  same ('ideal largest data fan-in, secded_gen.py''s rule', ideal_fanin (k, r), st.fanin);
  same ('largest data fan-in of the matrix', max (row_loads (data, r)), st.fanin);
end

function check_comparison (label, ot, book, r, k, st)
  section (sprintf ('OpenTitan''s %s matrix against the book''s, from SECDED and Hsiao Codes', label));
  od = ot(1:k);
  bd = book(1:k);
  same ('data columns the two matrices share', numel (intersect (od, bd)), st.common);
  same ('data positions holding the same column', sum (od == bd), st.agree);
  only_ot = od(~ismember (od, bd));
  parts = cell (1, numel (only_ot));
  for i = 1:numel (only_ot)
    parts{i} = rows_of (only_ot(i), r);
  end
  same ('columns only OpenTitan uses', strjoin (parts, ' '), st.ot_only);
  only_book = bd(~ismember (bd, od));
  parts = cell (1, numel (only_book));
  for i = 1:numel (only_book)
    parts{i} = rows_of (only_book(i), r);
  end
  same ('columns only the book uses', strjoin (parts, ' '), st.book_only);
  same ('book row weights', joined (row_loads (book, r)), st.book_rows);
  same ('book ones in the whole matrix', sum (row_loads (book, r)), st.book_ones);
  rows_differ = ~isequal (sort (row_loads (ot, r)), sort (row_loads (book, r)));
  same ('row-weight lists differ after sorting', rows_differ, st.rows_differ);
  if (isfield (st, 'ot_w5'))
    same ('OpenTitan weight-5 columns, shared rows per pair', ...
          intersections (od(weight (od) == 5)), st.ot_w5);
    same ('book weight-5 columns, shared rows per pair', ...
          intersections (bd(weight (bd) == 5)), st.book_w5);
  end
  same ('weight-4 codewords, OpenTitan', weight4_codewords (ot), st.ot_a4);
  same ('weight-4 codewords, book', weight4_codewords (book), st.book_a4);
end

function check_inverted (label, cols, r, k, inv, st)
  section (sprintf ('The inverted %s code, prim_secded_inv', label));
  zeros_k = repmat ('0', 1, k / 4);
  rule = sum (2 .^ (1:2:r - 1));
  same ('constant from the rule in secded_gen.py, hexadecimal', [dec2hex(rule, 2), zeros_k], st.hex);
  same ('constant copied from the inverted decoder, hexadecimal', [dec2hex(inv, 2), zeros_k], st.hex);
  same ('check part of the constant, row 0 first', bits (inv, r), st.check);
  h_const = inv;  % the constant is zero on the data, so its syndrome is its check part
  require ('the constant is not a codeword, so the code is not linear', h_const ~= 0);
  s_ones = 0;
  for c = cols
    s_ones = bitxor (s_ones, c);
  end
  same ('syndrome of the all-zero read', bits (h_const, r), st.s_zero);
  same ('its weight', weight (h_const), st.w_zero);
  same ('what the decoder flags for it', syndrome_class (h_const), st.c_zero);
  if (isfield (st, 'col_zero'))
    j = find (cols == h_const, 1) - 1;
    same ('data column it equals', j, st.col_zero);
    same ('data the inverted decoder returns, hexadecimal', dec2hex (2 ^ j, k / 4), st.d_zero);
  end
  same ('syndrome of the all-one word in the plain code', bits (s_ones, r), st.s_ones_plain);
  s_read = bitxor (s_ones, h_const);
  same ('syndrome of the all-one read in the inverted code', bits (s_read, r), st.s_ones);
  same ('its weight', weight (s_read), st.w_ones);
  same ('what the decoder flags for it', syndrome_class (s_read), st.c_ones);
  require ('the all-zero word is not in the inverted code', h_const ~= 0);
  require ('the all-one word is not in the inverted code', s_read ~= 0);
  % the two stored words the generator's formal assertion is about
  full = 2 ^ r - 1;
  chk_zero = bitxor (0, inv);
  chk_ones = bitxor (check_of (cols(1:k), ones (1, k)), inv);
  % all-zero data gives the all-zero word only if its check part is 0, and
  % all-one data gives the all-one word only if its check part is all ones
  require ('encoding all-zero or all-one data gives neither 0 nor all ones', ...
           chk_zero ~= 0 && chk_ones ~= full);
  same ('check part stored with all-zero data', bits (chk_zero, r), st.check);
  same ('check part stored with all-one data', bits (chk_ones, r), st.enc_ones);
end

function check_array ()
  section ('A 32 KiB data array, protected two ways');
  data_bits = 32 * 1024 * 8;
  words = data_bits / 64;
  same ('data bytes', 32 * 1024, 32768);
  same ('data bits', data_bits, 262144);
  same ('64-bit words', words, 4096);
  same ('extra bits for byte parity, one per byte', data_bits / 8, 32768);
  same ('extra bits for (72,64) SEC-DED, eight per word', words * 8, 32768);
  same ('overhead of either, percent', sprintf ('%.1f', 100 * 32768 / data_bits), '12.5');
  same ('bits stored either way', data_bits + 32768, 294912);
  pairs = 72 * 71 / 2;
  silent = 8 * (9 * 8 / 2);
  same ('two-bit errors in a 72-bit stored word', pairs, 2556);
  same ('of them inside one byte, so byte parity misses them', silent, 288);
  same ('share byte parity misses, percent', sprintf ('%.2f', 100 * silent / pairs), '11.27');
end

function check_double_errors (cols)
  section ('Every two-bit error in the OpenTitan (72,64) word');
  P = nchoosek (1:numel (cols), 2);
  x = bitxor (cols(P(:, 1)), cols(P(:, 2)));
  flagged = sum (x ~= 0 & mod (weight (x), 2) == 0);
  same ('two-bit errors flagged as double', flagged, 2556);
  same ('two-bit errors passing silently', 2556 - flagged, 0);
end

function check_partial_write (cols)
  section ('A one-byte write into a (72,64) word, OpenTitan masks');
  d0 = [hex2dec('01234567'), hex2dec('89ABCDEF')];
  old = bitand (bitshift (d0(2), -16), 255);
  new = hex2dec ('5A');
  delta = bitxor (old, new);
  d1 = [d0(1), d0(2) - old * 2 ^ 16 + new * 2 ^ 16];
  same ('old data word', hex64 (d0), '0123456789ABCDEF');
  same ('old byte 2', dec2hex (old, 2), 'AB');
  same ('new byte 2', dec2hex (new, 2), '5A');
  same ('change in the byte', dec2hex (delta, 2), 'F1');
  changed = 16 + find (bitand (bitshift (delta, -(0:7)), 1)) - 1;
  same ('data bits that change', joined (changed), '16,20,21,22,23');
  same ('new data word', hex64 (d1), '01234567895ACDEF');
  c0 = check_of (cols(1:64), word_bits (d0));
  c1 = check_of (cols(1:64), word_bits (d1));
  change = check_of (cols(17:24), bitand (bitshift (delta, -(0:7)), 1));
  same ('old check bits, row 0 first', bits (c0, 8), '01101010');
  same ('new check bits, re-encoded, row 0 first', bits (c1, 8), '00011111');
  same ('sum of columns 16, 20, 21, 22 and 23', bits (change, 8), '01110101');
  require ('old check bits plus that sum equal the re-encoded bits', bitxor (c0, change) == c1);
  stated = {'10001010', '10000011', '01110000', '01101000', '01100100'};
  js = [16 20 21 22 23];
  for i = 1:5
    same (sprintf ('column %d', js(i)), bits (cols(js(i) + 1), 8), stated{i});
  end
  % the same write when stored data bit 19 has flipped and the update skips the check
  stored = [d0(1), bitxor(d0(2), 2 ^ 19)];
  old_bad = bitand (bitshift (stored(2), -16), 255);
  delta_bad = bitxor (old_bad, new);
  c_bad = bitxor (c0, check_of (cols(17:24), bitand (bitshift (delta_bad, -(0:7)), 1)));
  s = bitxor (check_of (cols(1:64), word_bits (d1)), c_bad);
  same ('byte 2 as read with bit 19 flipped', dec2hex (old_bad, 2), 'A3');
  same ('change computed from that byte', dec2hex (delta_bad, 2), 'F9');
  same ('check bits written by the unchecked update', bits (c_bad, 8), '10011010');
  same ('syndrome of the next read', bits (s, 8), '10000101');
  same ('column 19', bits (cols(20), 8), '10000101');
  j = find (cols == s, 1) - 1;
  same ('data column the syndrome matches', j, 19);
  same ('data the decoder returns', hex64 ([d1(1), bitxor(d1(2), 2 ^ j)]), '012345678952CDEF');
  same ('bit 19 of the data written', bitand (bitshift (d1(2), -19), 1), 1);
  same ('even-parity bit stored with byte 5A', mod (weight (new), 2), 0);
end

function f = at_least_two (mu)
  % P(two or more upsets) for a Poisson count of mean mu, summed from j = 2.
  term = mu * mu / 2;
  total = term;
  for j = 3:30
    term = term * mu / j;
    total = total + term;
  end
  f = exp (-mu) * total;
end

function f = at_least_three (mu)
  % P(three or more upsets) for a Poisson count of mean mu, summed from j = 3.
  term = mu * mu * mu / 6;
  total = term;
  for j = 4:30
    term = term * mu / j;
    total = total + term;
  end
  f = exp (-mu) * total;
end

function check_scrubbing ()
  section ('Scrubbing a 4,096-word array, illustrative upset rate');
  lam = 1e-6;
  words = 4096;
  year = 8760;
  same ('upset rate per word per hour', e4 (lam), '1.0000e-06');
  same ('the same rate in FIT per word', round (lam * 1e9), 1000);
  T = [1 24 168 720];
  f_stated = {'5.0000e-13', '2.8800e-10', '1.4110e-08', '2.5908e-07'};
  lost_stated = {'1.7940e-05', '4.3056e-04', '3.0137e-03', '1.2911e-02'};
  lost = zeros (1, 4);
  mus = zeros (1, 4);
  for i = 1:4
    t = T(i);
    mu = lam * t;
    f = at_least_two (mu);
    lost(i) = words * (year / t) * f;
    mus(i) = mu;
    same (sprintf ('P(two or more in one word), T = %d h', t), e4 (f), f_stated{i});
    same (sprintf ('words lost a year, T = %d h', t), e4 (lost(i)), lost_stated{i});
    require (sprintf ('mu^2/2 - mu^3/3 <= P <= mu^2/2, T = %d h', t), ...
             mu * mu / 2 - mu ^ 3 / 3 <= f && f <= mu * mu / 2);
    require (sprintf ('mu^2/2 overstates P by at most 2mu/3 of itself, T = %d h', t), ...
             mu * mu / 2 - f <= (2 * mu / 3) * (mu * mu / 2));
  end
  same ('lambda T for a daily scrub', e4 (mus(2)), '2.4000e-05');
  same ('square law''s error bound 2 lambda T / 3, daily scrub', sprintf ('%.1e', 2 * mus(2) / 3), '1.6e-05');
  same ('linear estimate W Y lambda^2 T / 2, T = 24 h', e4 (words * year * lam * lam * 24 / 2), '4.3057e-04');
  share_T = [24 720 year];
  share_stated = {'8.0e-06', '2.4e-04', '2.9e-03'};
  for i = 1:3
    mu = lam * share_T(i);
    same (sprintf ('share of lost words with three or more upsets, T = %d h', share_T(i)), ...
          sprintf ('%.1e', at_least_three (mu) / at_least_two (mu)), share_stated{i});
  end
  mu_year = lam * year;
  lost_year = words * at_least_two (mu_year);
  same ('lambda T with no scrub for a year', e4 (mu_year), '8.7600e-03');
  same ('words holding two or more upsets after that year', e4 (lost_year), '1.5624e-01');
  same ('the fraction of the 4,096 words that is', e4 (at_least_two (mu_year)), '3.8145e-05');
  same ('ratio of a year unscrubbed to a daily scrub', sprintf ('%.0f', lost_year / lost(2)), '363');
  same ('factor 1 - 2 lambda T / 3 with no scrub for a year', sprintf ('%.3f', 1 - 2 * mu_year / 3), '0.994');
end

function m = worst (owner, start, len)
  w = owner (start:start + len - 1);
  m = max (accumarray (w(:) + 1, 1));
end

function check_interleaving ()
  section ('Bursts in a 288-cell row holding four 72-bit words');
  cells = 288;
  interleaved = @(c) mod (c, 4);
  side_by_side = @(c) floor (c / 72);
  groups = {[1 2 3 4], [5 6 7 8]};
  stated_count = [1146 1130];
  stated_worst = [1 2];
  for g = 1:2
    lengths = groups{g};
    count = 0;
    most = 0;
    for len = lengths
      for start = 0:cells - len
        count = count + 1;
        most = max (most, worst (interleaved, start, len));
      end
    end
    span = sprintf ('%d to %d', lengths(1), lengths(end));
    same (sprintf ('bursts of %s cells, interleaved', span), count, stated_count(g));
    same ('most flips one word takes from such a burst', most, stated_worst(g));
  end
  same ('word and bit held by cell 9', sprintf ('%d,%d', mod (9, 4), floor (9 / 4)), '1,2');
  bad = 0;
  for len = [2 3 4]
    for start = 0:cells - len
      if (worst (side_by_side, start, len) >= 2)
        bad = bad + 1;
      end
    end
  end
  same ('bursts of 2 to 4 cells, words side by side, uncorrectable', bad, 855);
  s = 0:cells - 2;
  same ('bursts of 2 cells across a word boundary', sum (side_by_side (s) ~= side_by_side (s + 1)), 3);
end


% Main -------------------------------------------------------------------------

o39 = ot_columns_39 ();
o72 = ot_columns_72 ();
section ('The transcribed encoder masks, hexadecimal');
M39 = masks_39 ();
parts = cell (1, 7);
for i = 1:7
  parts{i} = dec2hex (M39(i), 8);
end
same ('(39,32) masks', strjoin (parts, ' '), ...
      '2606BD25 DEBA8050 413D89AA 31234ED1 C2C1323B 2DCC624C 98505586');
M72 = masks_72 ();
for i = 1:8
  printf ('  (72,64) mask %d: %s\n', i - 1, hex64 (M72(i, :)));
end
check_config ();
check_matrix ('(39,32)', o39, 7, 32, struct ('wts', '3', 'rows', '15,15,15,15,15,15,13', ...
              'ones', 103, 'syn', 96, 'enc', 89, 'levels', 4, 'fanin', 14));
check_matrix ('(72,64)', o72, 8, 64, struct ('wts', '3,5', 'rows', '27,27,27,27,27,27,27,27', ...
              'ones', 216, 'syn', 208, 'enc', 200, 'levels', 5, 'fanin', 26));
check_comparison ('(39,32)', o39, hsiao_39_32 (), 7, 32, struct ( ...
  'common', 29, 'agree', 1, 'ot_only', '(0,3,4) (0,1,2) (1,5,6)', ...
  'book_only', '(0,1,6) (2,5,6) (3,4,6)', 'book_rows', '14,14,15,15,15,15,15', ...
  'book_ones', 103, 'rows_differ', true, 'ot_a4', 1366, 'book_a4', 1363));
check_comparison ('(72,64)', o72, hsiao_72_64 (), 8, 64, struct ( ...
  'common', 57, 'agree', 56, ...
  'ot_only', '(0,2,3,5,6) (1,2,4,6,7) (0,1,3,6,7) (0,1,4,5,7) (0,2,4,5,7) (1,2,3,6,7) (0,3,4,5,6)', ...
  'book_only', '(0,1,2,3,4) (2,3,4,5,6) (3,4,5,6,7) (0,4,5,6,7) (0,1,5,6,7) (0,1,2,6,7) (0,1,2,3,7)', ...
  'book_rows', '27,27,27,27,27,27,27,27', 'book_ones', 216, 'rows_differ', false, ...
  'ot_w5', '2:8 3:16 4:4', 'book_w5', '2:12 3:8 4:8', 'ot_a4', 8396, 'book_a4', 8392));
check_inverted ('(39,32)', o39, 7, 32, inv_check_39 (), struct ( ...
  'hex', '2A00000000', 'check', '0101010', 's_zero', '0101010', 'w_zero', 3, ...
  'c_zero', 'single', 'col_zero', 6, 'd_zero', '00000040', 's_ones_plain', '1111111', ...
  's_ones', '1010101', 'w_ones', 4, 'c_ones', 'double', 'enc_ones', '0101010'));
check_inverted ('(72,64)', o72, 8, 64, inv_check_72 (), struct ( ...
  'hex', 'AA0000000000000000', 'check', '01010101', 's_zero', '01010101', 'w_zero', 4, ...
  'c_zero', 'double', 's_ones_plain', '11111111', 's_ones', '10101010', 'w_ones', 4, ...
  'c_ones', 'double', 'enc_ones', '01010101'));
check_array ();
check_double_errors (o72);
check_partial_write (o72);
check_scrubbing ();
check_interleaving ();
printf ('\nAll chapter 29 numbers reproduced.\n');
