% Reproduce every number in chapter 30 of Error Correction Coding,
% "Project: A SEC-DED Encoder and Decoder", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch30_project_secded_encoder_and_decoder.m
%
% It needs Octave 7 or later for jsondecode, and only core Octave, with no pkg
% load and no Forge packages. It is the Octave twin of
% python/ch30_project_secded_encoder_and_decoder.py and prints the same lines,
% which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. The check matrices are rebuilt with the constructions of
% octave/ch20_secded_and_hsiao_codes.m and compared with the generated package
% hdl/secded72/secded72_pkg.sv, column for column. The flag counts come from
% every error pattern of weight at most two and from every syndrome. The
% measured gate counts and depths are read from data/ch30/secded72_gates.json,
% which synth/synth_secded72.py writes from a real Yosys netlist, and every
% percentage is computed here from that file. The Python route also checks the
% bench model in tb/secded72/secded72_model.py without printing; Octave cannot
% load that module, so this route leaves that check to Python.
%
% A column of a check matrix with r rows is held as a whole number whose bit i
% is the entry in row i, and a code is a row vector of such columns in
% codeword order. A bit pattern is printed with row 0 first. Positions are
% counted from 0, so position p of a code is element p + 1 here.

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
  s = char ('0' + bitand (bitshift (v, -(0:r - 1)), 1));
end

function s = joined (v)
  parts = cell (1, numel (v));
  for i = 1:numel (v)
    parts{i} = sprintf ('%d', v(i));
  end
  s = strjoin (parts, ',');
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

function s = pct (part, whole)
  s = sprintf ('%.1f', 100 * part / whole);
end

function loads = row_loads (cols, r)
  loads = zeros (1, r);
  for i = 0:r - 1
    loads(i + 1) = sum (bitand (bitshift (cols, -i), 1));
  end
end

function v = package_columns (text, prefix)
  % The column list a generated package holds, position 0 first. Every code
  % is stored in 72 entries of 8 bits, so a shorter code ends in zeros.
  block = regexp (text, ['localparam logic \[[^\]]*\] ' prefix '_COLUMNS = \{(.*?)\};'], 'tokens', 'once');
  words = regexp (block{1}, '8''b([01]+)', 'tokens');
  v = zeros (1, numel (words));
  for i = 1:numel (words)
    v(numel (words) + 1 - i) = bin2dec (words{i}{1});
  end
end

function fan = check_fan_in (cols, r, checks)
  % Data bits in each check equation, in the order of the check positions.
  % Every syndrome value is a sum of check columns in exactly one way, so a
  % table over all 2^r subsets of them solves each data column.
  table = zeros (1, 2 ^ r);
  for subset = 0:2 ^ r - 1
    v = 0;
    for j = 0:r - 1
      if (bitand (bitshift (subset, -j), 1))
        v = bitxor (v, cols(checks(j + 1) + 1));
      end
    end
    table(v + 1) = subset;
  end
  data = setdiff (0:numel (cols) - 1, checks);
  fan = zeros (1, r);
  for j = 0:r - 1
    fan(j + 1) = sum (bitand (bitshift (table(cols(data + 1) + 1), -j), 1));
  end
end

function ok = parity_ok (name, s, r)
  % The test the parity form of the flag applies: odd weight, or the parity row set.
  if (strcmp (name, 'hamming'))
    ok = bitand (bitshift (s, -(r - 1)), 1) == 1;
  else
    ok = mod (weight (s), 2) == 1;
  end
end


% The codes, as octave/ch20_secded_and_hsiao_codes.m builds them --------------

function h = hsiao_72_64 ()
  runs = zeros (1, 8);
  for i = 0:7
    runs(i + 1) = col (mod (i + (0:4), 8));
  end
  h = [cols_of(nchoosek (0:7, 3)), runs, 2 .^ (0:7)];
end

function h = hsiao_72_64_dictionary ()
  C5 = nchoosek (0:7, 5);
  h = [cols_of(nchoosek (0:7, 3)), cols_of(C5(1:8, :)), 2 .^ (0:7)];
end

function h = hsiao_39_32 ()
  data = cols_of (nchoosek (0:6, 3));
  data = data(~ismember (data, [col([0 1 2]), col([0 3 4]), col([1 5 6])]));
  h = [data, 2 .^ (0:6)];
end

function h = extended_hamming (n, rh)
  h = bitor (0:n - 1, 2 ^ rh);
end

function c = codes ()
  c = struct ('name', {}, 'prefix', {}, 'construction', {}, 'cols', {}, 'r', {}, 'checks', {}, 'label', {});
  c(1) = struct ('name', 'hsiao', 'prefix', 'HSIAO', 'construction', 'hsiao_72_64()', ...
                 'cols', hsiao_72_64 (), 'r', 8, 'checks', 64:71, 'label', 'Hsiao');
  c(2) = struct ('name', 'hamming', 'prefix', 'HAMMING', 'construction', 'extended_hamming(72, 7)', ...
                 'cols', extended_hamming (72, 7), 'r', 8, 'checks', [0, 2 .^ (0:6)], 'label', 'Hamming');
  c(3) = struct ('name', 'hsiao_dictionary', 'prefix', 'HSIAO_DICT', 'construction', 'hsiao_72_64_dictionary()', ...
                 'cols', hsiao_72_64_dictionary (), 'r', 8, 'checks', 64:71, 'label', 'Dictionary');
  c(4) = struct ('name', 'hsiao_39_32', 'prefix', 'HSIAO_39_32', 'construction', 'hsiao_39_32()', ...
                 'cols', hsiao_39_32 (), 'r', 7, 'checks', 32:38, 'label', 'Hsiao (39,32)');
end

function st = stated (name)
  switch (name)
    case 'hsiao'
      st = struct ('ones', 216, 'rows', '27,27,27,27,27,27,27,27', 'largest', 27, 'levels', 5, 'sep', 208, ...
                   'fan', '26,26,26,26,26,26,26,26', 'enc', 200, 'enc_levels', 5, ...
                   'matched', 72, 'nomatch', 183, 'parity', 127, 'only', 56, ...
                   'hand0', 855, 'hand1', 862, 'depth', '9,15,16,9', ...
                   'syn_whole', 169, 'syn_whole_depth', '5,5,5,5,5,5,5,5', 'enc_whole', 159, 'enc_whole_depth', 5, ...
                   'dec0', '550,9,15,16', 'dec1', '457,9,15,10', 'cells0', '139,154,257,39');
    case 'hamming'
      st = struct ('ones', 284, 'rows', '36,36,36,32,32,32,8,72', 'largest', 72, 'levels', 7, 'sep', 276, ...
                   'fan', '35,35,35,35,31,31,31,7', 'enc', 232, 'enc_levels', 6, ...
                   'matched', 72, 'nomatch', 183, 'parity', 127, 'only', 56, ...
                   'hand0', 923, 'hand1', 922, 'depth', '11,17,18,10', ...
                   'syn_whole', 137, 'syn_whole_depth', '6,6,6,5,5,5,3,7', 'enc_whole', 175, 'enc_whole_depth', 6, ...
                   'dec0', '478,11,17,18', 'dec1', '385,11,17,10', 'cells0', '129,142,207,48');
    case 'hsiao_dictionary'
      st = struct ('ones', 216, 'rows', '30,30,30,26,26,25,25,24', 'largest', 30, 'levels', 5, 'sep', 208, ...
                   'fan', '29,29,29,25,25,24,24,23', 'enc', 200, 'enc_levels', 5, ...
                   'matched', 72, 'nomatch', 183, 'parity', 127, 'only', 56, ...
                   'hand0', 855, 'hand1', 862, 'depth', '9,15,16,9', ...
                   'syn_whole', 164, 'syn_whole_depth', '5,5,5,5,5,5,5,5', 'enc_whole', 156, 'enc_whole_depth', 5, ...
                   'dec0', '549,9,15,16', 'dec1', '457,9,15,10', 'cells0', '139,153,257,39');
    otherwise
      st = struct ('ones', 103, 'rows', '14,14,15,15,15,15,15', 'largest', 15, 'levels', 4, 'sep', 96, ...
                   'fan', '13,13,14,14,14,14,14', 'enc', 89, 'enc_levels', 4, ...
                   'matched', 39, 'nomatch', 88, 'parity', 63, 'only', 25, ...
                   'hand0', 407, 'hand1', 413, 'depth', '8,13,14,8', ...
                   'syn_whole', 80, 'syn_whole_depth', '4,4,4,4,4,4,4', 'enc_whole', 73, 'enc_whole_depth', 4, ...
                   'dec0', '303,8,13,14', 'dec1', '251,8,13,9', 'cells0', '80,91,132,23');
  end
end

function m = entry (all_codes, name)
  % One code's record in the data file, whether jsondecode made a struct array or a cell array.
  for i = 1:numel (all_codes)
    if (iscell (all_codes))
      e = all_codes{i};
    else
      e = all_codes(i);
    end
    if (strcmp (e.code, name))
      m = e;
      return;
    end
  end
  error ('no record for %s in the data file', name);
end

function v = field (s, key)
  v = [s.(key)];
end


% The checks -------------------------------------------------------------------

here = fileparts (mfilename ('fullpath'));
lab = fileparts (here);
package_text = fileread (fullfile (lab, 'hdl', 'secded72', 'secded72_pkg.sv'));
gates = jsondecode (fileread (fullfile (lab, 'data', 'ch30', 'secded72_gates.json')));
C = codes ();
RANDOM_WORDS = 1000;

printf ('Chapter 30, Project: A SEC-DED Encoder and Decoder\n');

section ('The check matrices, rebuilt from chapter 20 and read from the package');
fan_ins = cell (1, numel (C));
for t = 1:numel (C)
  c = C(t);
  st = stated (c.name);
  n = numel (c.cols);
  stored = package_columns (package_text, c.prefix);
  require (sprintf ('%s: package columns equal %s', c.label, c.construction), ...
           isequal (stored(1:n), c.cols) && ~any (stored(n + 1:end)));
  loads = row_loads (c.cols, c.r);
  same ([c.label ': ones in the matrix'], sum (loads), st.ones);
  same ([c.label ': ones in each row'], joined (loads), st.rows);
  same ([c.label ': largest row'], max (loads), st.largest);
  same ([c.label ': XOR levels for the largest row'], levels (max (loads)), st.levels);
  same ([c.label ': syndrome XORs as separate trees'], sum (loads - 1), st.sep);
  fan = check_fan_in (c.cols, c.r, c.checks);
  fan_ins{t} = fan;
  same ([c.label ': data bits in each check equation'], joined (fan), st.fan);
  same ([c.label ': encoder XORs as separate trees'], sum (fan - 1), st.enc);
  same ([c.label ': encoder levels'], levels (max (fan)), st.enc_levels);
end

section ('Two flags, compared on every syndrome and every pattern of weight <= 2');
for t = 1:numel (C)
  c = C(t);
  st = stated (c.name);
  n = numel (c.cols);
  if (strcmp (c.name, 'hamming'))
    require ([c.label ': every column has the parity row set'], all (parity_ok (c.name, c.cols, c.r)));
  else
    require ([c.label ': every column has odd weight'], all (parity_ok (c.name, c.cols, c.r)));
  end
  require ([c.label ': columns distinct and nonzero'], numel (unique (c.cols)) == n && ~any (c.cols == 0));
  s = 0:2 ^ c.r - 1;
  m = ismember (s, c.cols);
  no_match = s ~= 0 & ~m;
  par = s ~= 0 & ~parity_ok (c.name, s, c.r);
  same ([c.label ': syndromes that match a column'], sum (m), st.matched);
  same ([c.label ': syndromes that set the no-match flag'], sum (no_match), st.nomatch);
  same ([c.label ': syndromes that set the parity flag'], sum (par), st.parity);
  same ([c.label ': syndromes that set the no-match flag alone'], sum (no_match & ~par), st.only);
  require ([c.label ': the parity flag never fires alone'], ~any (par & ~no_match));
  pairs = nchoosek (1:n, 2);
  s = bitxor (c.cols(pairs(:, 1)), c.cols(pairs(:, 2)));
  no_match = s ~= 0 & ~ismember (s, c.cols);
  par = s ~= 0 & ~parity_ok (c.name, s, c.r);
  same ([c.label ': single errors whose syndrome is a column'], sum (ismember (c.cols, c.cols)), n);
  same ([c.label ': double errors that set the no-match flag'], sum (no_match), n * (n - 1) / 2);
  same ([c.label ': patterns of weight <= 2 where the flags disagree'], sum (no_match ~= par), 0);
end

section ('What the bench runs');
for t = 1:numel (C)
  c = C(t);
  n = numel (c.cols);
  k = n - c.r;
  patterns = 1 + n + n * (n - 1) / 2;
  words = 2 + k + RANDOM_WORDS;
  if (n == 72)
    stated_patterns = 2629;
    stated_cases = 2802514;
  else
    stated_patterns = 781;
    stated_cases = 807554;
  end
  if (k == 64)
    stated_words = 1066;
  else
    stated_words = 1034;
  end
  same ([c.label ': error patterns of weight <= 2'], patterns, stated_patterns);
  same ([c.label ': data words'], words, stated_words);
  same ([c.label ': cases'], words * patterns, stated_cases);
  if (n == 72)
    same ([c.label ': cases of weight 0, 1 and 2'], ...
          joined ([words, words * n, words * n * (n - 1) / 2]), '1066,76752,2724696');
  end
end

section ('The decoder counted by hand, before any synthesis');
for t = 1:numel (C)
  c = C(t);
  st = stated (c.name);
  n = numel (c.cols);
  r = c.r;
  k = n - r;
  loads = row_loads (c.cols, r);
  syndrome = sum (loads - 1);
  comparators = n * (r - 1);
  base = syndrome + comparators + k + (n - 1) + 1;
  if (strcmp (c.name, 'hamming'))
    flag0 = base + (r - 1);
    flag1 = base + (r - 2);
  else
    flag0 = base + (r - 1);
    flag1 = base + (r - 1) + (r - 1);
  end
  if (n == 72)
    stated_comparators = 504;
  else
    stated_comparators = 234;
  end
  same ([c.label ': comparator AND gates'], comparators, stated_comparators);
  same ([c.label ': gates, no-match flag'], flag0, st.hand0);
  same ([c.label ': gates, parity flag'], flag1, st.hand1);
  ls = levels (max (loads));
  data = ls + levels (r) + 1;
  corrected = ls + levels (r) + levels (n);
  unc0 = corrected + 1;
  if (strcmp (c.name, 'hamming'))
    low = levels (max (loads(1:r - 1)));
    unc1 = max (low + levels (r - 1), levels (loads(r))) + 1;
  else
    unc1 = ls + levels (r) + 1;
  end
  same ([c.label ': depth of data, corrected, both flags'], joined ([data, corrected, unc0, unc1]), st.depth);
end

section ('Measured in the Yosys netlist, data/ch30/secded72_gates.json');
build = regexp (gates.yosys, '\(git sha1 ([0-9a-f]{7})', 'tokens', 'once');
same ('the Yosys build that measured the data file', ['git ' build{1}], 'git 5628c13');
for t = 1:numel (C)
  c = C(t);
  st = stated (c.name);
  m = entry (gates.codes, c.name);
  loads = row_loads (c.cols, c.r);
  fan = fan_ins{t};
  require ([c.label ': the data file''s rows match the matrix'], isequal (m.row_weights(:)', loads));
  require ([c.label ': the data file''s check equations match'], isequal (m.check_fan_in(:)', fan));
  trees = m.syndrome_trees;
  tree_levels = arrayfun (@levels, loads);
  require ([c.label ': each separate tree has W - 1 gates, ceil(log2 W) levels'], ...
           isequal (field (trees, 'gates'), loads - 1) && isequal (field (trees, 'xor'), loads - 1) ...
           && isequal (field (trees, 'depth'), tree_levels));
  same ([c.label ': syndrome XORs, separate trees'], sum (field (trees, 'gates')), st.sep);
  same ([c.label ': syndrome depth, separate trees'], max (field (trees, 'depth')), st.levels);
  enc = m.encoder_trees;
  require ([c.label ': each separate encoder tree has F - 1 gates'], ...
           isequal (field (enc, 'gates'), fan - 1) && isequal (field (enc, 'depth'), arrayfun (@levels, fan)));
  same ([c.label ': encoder XORs, separate trees'], sum (field (enc, 'gates')), st.enc);
  sw = m.syndrome_whole;
  same ([c.label ': syndrome two-input gates, whole block'], sw.two_input, st.syn_whole);
  same ([c.label ': syndrome depth per bit, whole block'], joined (sw.depth), st.syn_whole_depth);
  ew = m.encoder_whole;
  same ([c.label ': encoder two-input gates, whole block'], ew.two_input, st.enc_whole);
  same ([c.label ': encoder depth, whole block'], ew.depth, st.enc_whole_depth);
  keys = {'dec0', 'dec1'};
  for flag = 0:1
    d = m.(sprintf ('decoder_whole_flag%d', flag));
    same (sprintf ('%s: decoder flag %d: gates, data, corrected, flag depth', c.label, flag), ...
          joined ([d.two_input, d.data_depth, d.corrected_depth, d.uncorrectable_depth]), st.(keys{flag + 1}));
  end
  c0 = m.decoder_whole_flag0.cells;
  same ([c.label ': decoder flag 0 cells and, or, xor, not'], ...
        joined ([c0.('and'), c0.('or'), c0.('xor'), c0.('not')]), st.cells0);
end

section ('Hsiao against extended Hamming, computed from the data file');
h = entry (gates.codes, 'hsiao');
e = entry (gates.codes, 'hamming');
hs = sum (field (h.syndrome_trees, 'gates'));
es = sum (field (e.syndrome_trees, 'gates'));
he = sum (field (h.encoder_trees, 'gates'));
ee = sum (field (e.encoder_trees, 'gates'));
hd = max (field (h.syndrome_trees, 'depth'));
ed = max (field (e.syndrome_trees, 'depth'));
same ('syndrome XORs saved by Hsiao, separate trees, percent', pct (es - hs, es), '24.6');
same ('encoder XORs saved by Hsiao, separate trees, percent', pct (ee - he, ee), '13.8');
same ('syndrome levels saved by Hsiao, percent', pct (ed - hd, ed), '28.6');
hw = h.syndrome_whole.two_input;
ew = e.syndrome_whole.two_input;
same ('extra syndrome gates in Hsiao, whole block, percent', pct (hw - ew, ew), '23.4');
hx = h.encoder_whole.two_input;
ex = e.encoder_whole.two_input;
same ('encoder gates saved by Hsiao, whole block, percent', pct (ex - hx, ex), '9.1');
flag_stated = {'15.1', '18.7'};
for flag = 0:1
  key = sprintf ('decoder_whole_flag%d', flag);
  a = h.(key).two_input;
  b = e.(key).two_input;
  same (sprintf ('extra decoder gates in Hsiao, flag %d, percent', flag), pct (a - b, b), flag_stated{flag + 1});
end
a = h.decoder_whole_flag0.data_depth;
b = e.decoder_whole_flag0.data_depth;
same ('data path levels saved by Hsiao, percent', pct (b - a, b), '18.2');
same ('Hsiao: share of the hand count shared away, percent', pct (855 - h.decoder_whole_flag0.two_input, 855), '35.7');
same ('Hamming: share of the hand count shared away, percent', pct (923 - e.decoder_whole_flag0.two_input, 923), '48.2');

section ('Worked examples and exercises');
hc = hsiao_72_64 ();
same ('Hsiao: column at position 0', bits (hc(1), 8), '11100000');
same ('Hsiao: column at position 1', bits (hc(2), 8), '11010000');
same ('Hsiao: syndrome of errors at positions 0 and 1', bits (bitxor (hc(1), hc(2)), 8), '00110000');
p = find (hc == col ([3 5 6])) - 1;
same ('Hsiao: position of the column 00010110', p, 49);
s = bitxor (bitxor (hc(1), hc(p + 1)), hc(71));
same ('Hsiao: syndrome of errors at positions 0, 49 and 70', bits (s, 8), '11110100');
require ('Hsiao: that syndrome has odd weight and is no column', mod (weight (s), 2) == 1 && ~ismember (s, hc));
% With every data bit 1, each check bit is the parity of its equation's size.
same ('Hsiao: check bits of the all-ones data word', char ('0' + mod (fan_ins{1}, 2)), '00000000');
same ('Hamming: check bits of the all-ones data word', char ('0' + mod (fan_ins{2}, 2)), '11111111');
q = 0:71;
index_rows = zeros (1, 7);
for b = 0:6
  index_rows(b + 1) = sum (bitand (bitshift (q, -b), 1));
end
same ('position output: match lines feeding each index bit', joined (index_rows), '36,36,36,32,32,32,8');
same ('position output: OR gates', sum (index_rows - 1), 205);
same ('position output: levels after the match lines', levels (max (index_rows)), 6);

printf ('\nAll chapter 30 numbers reproduced.\n');
