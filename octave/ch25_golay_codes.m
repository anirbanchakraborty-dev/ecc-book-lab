% Reproduce every number in chapter 25 of Error Correction Coding,
% "Golay Codes", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch25_golay_codes.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch25_golay_codes.py and prints the same lines,
% which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in the chapter is a whole number, a binary or
% ternary word or a matrix, so every comparison is exact. Both codes are
% checked exhaustively: all 4,096 codewords of the extended binary Golay code
% and all 729 codewords of the extended ternary Golay code are enumerated,
% every 5-subset of the 24 coordinates is shown to lie in exactly one octad,
% every weight-4 word of length 23 is shown to lie within distance 3 of
% exactly one codeword of the perfect code, and the decoder is run on all
% 2,324 error patterns of weight 1 to 3 and all 10,626 of weight 4.
%
% A binary word of length 24 is held as a whole number whose most
% significant bit is coordinate 0, so dec2bin (x, 24) prints it left to
% right. The left half, coordinates 0 to 11, is the top 12 bits, and the
% right half, coordinates 12 to 23, is the bottom 12 bits. Row i of B is
% held the same way as a 12-bit number, in element i + 1 of brow (). A
% ternary word is a row of 0s, 1s and 2s, and the ternary messages are taken
% in the order Python's itertools.product gives, first coordinate slowest.

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
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', ...
            label, as_repr (computed), as_repr (stated));
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

function s = join (v)
  s = strtrim (sprintf ('%d ', v));
end

function p = pop (x)
  % Number of 1 bits in each element of x, for words of up to 24 bits.
  persistent T;
  if (isempty (T))
    T = zeros (1, 4096);
    for v = 1:4095
      T(v + 1) = T(floor (v / 2) + 1) + mod (v, 2);
    end
  end
  p = T(bitand (x, 4095) + 1) + T(floor (x / 4096) + 1);
end

function s = support (x, n)
  if (nargin < 2)
    n = 24;
  end
  s = find (dec2bin (x, n) == '1') - 1;
end

function m = mask (positions, n)
  if (nargin < 2)
    n = 24;
  end
  m = sum (2 .^ (n - 1 - positions));
end

function s = distribution (w)
  u = unique (w(:))';
  parts = cell (1, numel (u));
  for k = 1:numel (u)
    parts{k} = sprintf ('%d:%d', u(k), sum (w(:) == u(k)));
  end
  s = strjoin (parts, ' ');
end


% The binary code --------------------------------------------------------------

function B = b_matrix ()
  % B as a 12-by-12 matrix of bits: squares pattern, then the border.
  persistent M;
  if (isempty (M))
    S = [0 1 3 4 5 9];
    M = zeros (12, 12);
    for i = 0:10
      for j = 0:10
        M(i + 1, j + 1) = any (mod (i + j, 11) == S);
      end
      M(i + 1, 12) = 1;
    end
    M(12, :) = [ones(1, 11) 0];
  end
  B = M;
end

function r = brow ()
  persistent R;
  if (isempty (R))
    B = b_matrix ();
    R = zeros (1, 12);
    for i = 1:12
      R(i) = sum (B(i, :) .* 2 .^ (11:-1:0));
    end
  end
  r = R;
end

function out = times_b (v)
  % The 12-bit row vector v times B over GF(2), for each element of v.
  persistent TB;
  if (isempty (TB))
    R = brow ();
    TB = zeros (1, 4096);
    for x = 0:4095
      acc = 0;
      for k = 0:11
        if (bitand (floor (x / 2 ^ (11 - k)), 1))
          acc = bitxor (acc, R(k + 1));
        end
      end
      TB(x + 1) = acc;
    end
  end
  out = TB(v + 1);
end

function c = encode (a)
  c = a * 4096 + times_b (a);
end

function s = syndrome (r)
  s = bitxor (floor (r / 4096), times_b (bitand (r, 4095)));
end

function [kase, e] = decode (r)
  % Return the case and the error for a received word r; case 5 means flagged.
  R = brow ();
  s = syndrome (r);
  if (pop (s) <= 3)
    kase = 1;
    e = s * 4096;
    return;
  end
  i = find (pop (bitxor (s, R)) <= 2, 1);
  if (~isempty (i))
    kase = 2;
    e = bitxor (s, R(i)) * 4096 + 2 ^ (12 - i);
    return;
  end
  t = times_b (s);
  if (pop (t) <= 3)
    kase = 3;
    e = t;
    return;
  end
  i = find (pop (bitxor (t, R)) <= 2, 1);
  if (~isempty (i))
    kase = 4;
    e = 2 ^ (12 - i) * 4096 + bitxor (t, R(i));
    return;
  end
  kase = 5;
  e = -1;
end

function w = column_weights ()
  w = sum (b_matrix (), 1);
end

function g = chain_gates (n, m)
  % Two-input AND and OR gates in a thermometer chain that flags weight >= m.
  g = 0;
  for k = 1:n - 1
    g = g + 1;
    for j = 2:min (k + 1, m)
      if (j <= k)
        g = g + 2;
      else
        g = g + 1;
      end
    end
  end
end

function f = chain_flag (x, n, m)
  t = zeros (m + 1, numel (x));
  for k = 0:n - 1
    b = bitand (floor (x / 2 ^ (n - 1 - k)), 1);
    for j = m:-1:2
      t(j + 1, :) = t(j + 1, :) | (t(j, :) & b);
    end
    t(2, :) = t(2, :) | b;
  end
  f = t(m + 1, :);
end

function st = gate_counts ()
  cols = column_weights ();
  st.syndrome_xor = sum (cols);
  st.syndrome_depth = max (arrayfun (@(c) numel (dec2bin (c)), cols));
  st.sb_xor = sum (cols - 1);
  st.sb_depth = max (arrayfun (@(c) numel (dec2bin (c - 1)), cols));
  st.fold_inverters = sum (pop (brow ()));
  st.literal_xor = 12 * 12;
  st.chain3 = chain_gates (12, 3);
  st.chain4 = chain_gates (12, 4);
  R = brow ();
  orr = 0;
  for k = 1:numel (R)
    orr = bitor (orr, R(k));
  end
  st.shared_inverters = pop (orr);
  st.flag_inversions = 2 * (numel (R) + 1);
  st.series_depth = st.syndrome_depth + st.sb_depth;
  st.correction_xor = 2 * numel (R);
end


% The ternary code -------------------------------------------------------------

function B3 = b3_matrix ()
  CHI5 = [0 1 2 2 1];
  B3 = zeros (6, 6);
  B3(1, :) = [0 1 1 1 1 1];
  for i = 1:5
    B3(i + 1, 1) = 1;
    for j = 1:5
      B3(i + 1, j + 1) = CHI5(mod (j - i, 5) + 1);
    end
  end
end


% The checks -------------------------------------------------------------------

function check_b ()
  B = b_matrix ();
  R = brow ();
  stated = {'110111000101', '101110001011', '011100010111', '111000101101', ...
            '110001011011', '100010110111', '000101101111', '001011011101', ...
            '010110111001', '101101110001', '011011100011', '111111111110'};
  meet = {'0 3 4', '1 3 9', '0 1 9', '0 1 5', '0 4 9'};
  section ('The matrix B');
  for i = 0:11
    same (sprintf ('row %d of B', i), dec2bin (R(i + 1), 12), stated{i + 1});
  end
  require ('B is symmetric', isequal (B, B'));
  ok = true;
  for i = 2:11
    ok = ok && isequal (B(i, 1:11), [B(i - 1, 2:11) B(i - 1, 1)]);
  end
  require ('each of rows 1 to 10 is the row above shifted left', ok);
  same ('weights of rows 0 to 10', join (unique (pop (R(1:11)))), '7');
  same ('weight of row 11', pop (R(12)), 11);
  same ('column weights of B', join (column_weights ()), '7 7 7 7 7 7 7 7 7 7 7 11');

  section ('Row overlaps by shift');
  S = [0 1 3 4 5 9];
  ov = @(i, j) sum (B(i + 1, 1:11) & B(j + 1, 1:11));
  for t = 1:5
    m = intersect (S, mod (S - t, 11));
    same (sprintf ('S intersected with S - %d', t), join (m), meet{t});
    same (sprintf ('overlap of rows 0 and %d in columns 0 to 10', t), ov (0, t), 3);
    same (sprintf ('overlap of rows 0 and %d in columns 0 to 10', 11 - t), ov (0, 11 - t), 3);
  end
  P = nchoosek (0:10, 2);
  same ('pairs among rows 0 to 10', size (P, 1), 55);
  d = mod (P(:, 2) - P(:, 1), 11);
  per = zeros (1, 5);
  for t = 1:5
    per(t) = sum (d == t | d == 11 - t);
  end
  same ('pairs whose shift is t or -t, for t = 1 to 5', join (per), '11 11 11 11 11');
  ok = true;
  for k = 1:size (P, 1)
    ok = ok && ov (P(k, 1), P(k, 2)) == 3;
  end
  require ('overlap of any two of rows 0 to 10 is 3', ok);
  same ('overlap of row 11 with each of rows 0 to 10', ...
        join (unique (pop (bitand (R(12), R(1:11))))), '6');
  require ('B times B is the identity over GF(2)', isequal (mod (B * B, 2), eye (12)));
  Q = nchoosek (1:12, 2);
  dists = unique (pop (bitxor (R(Q(:, 1)), R(Q(:, 2)))));
  same ('distances between two different rows of B', join (dists), '6');
end

function check_g24 (CODE24, GROWS)
  B = b_matrix ();
  section ('The extended binary Golay code G24');
  ok = true;
  for i = 1:12
    for j = 1:12
      ok = ok && mod (pop (bitand (GROWS(i), GROWS(j))), 2) == 0;
    end
  end
  require ('G times G transpose is zero', ok);
  same ('weights of the rows of G', join (unique (pop (GROWS))), '8 12');
  same ('codewords enumerated', numel (unique (CODE24)), 4096);
  w = pop (CODE24);
  require ('every codeword weight is divisible by 4', all (mod (w, 4) == 0));
  same ('codewords of weight 4', sum (w == 4), 0);
  same ('minimum distance', min (w(CODE24 > 0)), 8);
  total = 0;
  for i = 1:12
    total = bitxor (total, GROWS(i));
  end
  same ('sum of the 12 rows of G', dec2bin (total, 24), repmat ('1', 1, 24));
  require ('the all-ones word is a codeword', any (CODE24 == 2 ^ 24 - 1));
  same ('weight distribution of G24', distribution (w), '0:1 8:759 12:2576 16:759 24:1');
  same ('A_12 = 4096 - 2 - 2 * 759', 4096 - 2 - 2 * 759, 2576);
  V = dec2bin (0:4095, 12) - '0';
  require ('r_L + r_R B equals r_L + r_R B transpose for every r_R', ...
           isequal (mod (V * B, 2), mod (V * B', 2)) && ...
           isequal ((mod (V * B, 2) * 2 .^ (11:-1:0)')', times_b (0:4095)));
  require ('every codeword has syndrome zero', all (syndrome (CODE24) == 0));
end

function check_octads (CODE24)
  section ('Octads and the Steiner system S(5,8,24)');
  octs = CODE24(pop (CODE24) == 8);
  same ('C(24,5)', nchoosek (24, 5), 42504);
  same ('C(8,5)', nchoosek (8, 5), 56);
  same ('A_8 = C(24,5) / C(8,5)', nchoosek (24, 5) / nchoosek (8, 5), 759);
  same ('octads found', numel (octs), 759);
  sizes = [5 4 3 2 1];
  lam = [1 5 21 77 253];
  for q = 1:5
    sz = sizes(q);
    per = nchoosek (8, sz);
    masks = zeros (numel (octs) * per, 1);
    for k = 1:numel (octs)
      sub = nchoosek (support (octs(k)), sz);
      masks((k - 1) * per + (1:per)) = sum (2 .^ (23 - sub), 2);
    end
    [u, ~, j] = unique (masks);
    cnt = accumarray (j, 1);
    same (sprintf ('%d-subsets reached by the octads', sz), numel (u), nchoosek (24, sz));
    same (sprintf ('octads through each %d-subset', sz), join (unique (cnt)), sprintf ('%d', lam(q)));
  end
  target = mask ([0 1 2 3 4]);
  found = octs(bitand (octs, target) == target);
  same ('octads containing 0 1 2 3 4', numel (found), 1);
  same ('the octad containing 0 1 2 3 4', join (support (found(1))), '0 1 2 3 4 6 14 17');
  R = brow ();
  x = R(1);
  partial = cell (1, 4);
  for i = 2:5
    x = bitxor (x, R(i));
    partial{i - 1} = dec2bin (x, 12);
  end
  same ('partial sums b_0 + b_1, ..., b_0 + ... + b_4', strjoin (partial, ' '), ...
        '011001001110 000101011001 111101110100 001100101111');
  same ('weight of x = b_0 + ... + b_4', pop (x), 7);
  same ('row 6 of B', dec2bin (R(7), 12), '000101101111');
  same ('x + b_6', dec2bin (bitxor (x, R(7)), 12), '001001000000');
  octad = encode (bin2dec ('111110100000'));
  same ('octad of the codeword with message 111110100000', join (support (octad)), '0 1 2 3 4 6 14 17');
end

function check_g23 (CODE24)
  section ('The perfect binary Golay code G23');
  CODE23 = floor (CODE24 / 2);
  w = pop (CODE23);
  same ('distinct words after puncturing', numel (unique (CODE23)), 4096);
  same ('minimum distance of G23', min (w(CODE23 > 0)), 7);
  v = 0;
  for i = 0:3
    v = v + nchoosek (23, i);
  end
  same ('V(23,3) = 1 + 23 + 253 + 1771', v, 2048);
  same ('2^12 * V(23,3)', 4096 * v, 2 ^ 23);
  same ('weight distribution of G23', distribution (w), ...
        '0:1 7:253 8:506 11:1288 12:1288 15:506 16:253 23:1');
  same ('C(23,4)', nchoosek (23, 4), 8855);
  same ('C(7,4)', nchoosek (7, 4), 35);
  same ('A_7 = C(23,4) / C(7,4)', nchoosek (23, 4) / nchoosek (7, 4), 253);
  octs = CODE24(pop (CODE24) == 8);
  same ('octads with a 1 in coordinate 23', sum (bitand (octs, 1) == 1), 253);
  same ('octads with a 0 in coordinate 23', sum (bitand (octs, 1) == 0), 506);
  near = CODE23(w <= 7);
  W = sum (2 .^ (22 - nchoosek (0:22, 4)), 2)';
  nclose = zeros (size (W));
  which = zeros (size (W));
  for k = 1:numel (near)
    hit = pop (bitxor (W, near(k))) <= 3;
    nclose = nclose + hit;
    which(hit) = near(k);
  end
  good = sum (nclose == 1 & pop (which) == 7 & bitand (which, W) == W);
  same ('weight-4 words with one codeword within 3, of weight 7', good, 8855);
  rebuilt = CODE23 * 2 + mod (w, 2);
  require ('extending G23 gives back G24', isequal (sort (rebuilt), sort (CODE24)));
end

function check_ternary ()
  section ('The ternary Golay codes');
  B3 = b3_matrix ();
  stated = {'011111', '101221', '110122', '121012', '122101', '112210'};
  for i = 0:5
    same (sprintf ('row %d of B3', i), sprintf ('%d', B3(i + 1, :)), stated{i + 1});
  end
  require ('B3 is symmetric', isequal (B3, B3'));
  require ('B3 times B3 transpose is -I over GF(3)', isequal (mod (B3 * B3', 3), 2 * eye (6)));
  G = [eye(6) B3];
  require ('G times G transpose is zero over GF(3)', all (all (mod (G * G', 3) == 0)));
  same ('weights of the rows of G', join (unique (sum (G ~= 0, 2))), '6');
  same ('row 1 of B3 with itself, before reducing mod 3', sum (B3(2, :) .^ 2), 11);
  same ('row 1 of B3 with itself, mod 3', mod (sum (B3(2, :) .^ 2), 3), 2);
  same ('row 1 with row 2 of B3, before reducing mod 3', sum (B3(2, :) .* B3(3, :)), 9);
  same ('row 1 with row 2 of B3, mod 3', mod (sum (B3(2, :) .* B3(3, :)), 3), 0);
  chi = [0 1 -1 -1 1];
  tstated = {'0 -1 1 -1 0', '0 -1 -1 0 1'};
  for t = 1:2
    terms = chi(mod ((0:4) + t, 5) + 1) .* chi;
    same (sprintf ('terms chi(x + %d) chi(x), x = 0 to 4', t), join (terms), tstated{t});
    same (sprintf ('their sum for t = %d', t), sum (terms), -1);
  end
  same ('row 0 of G3', sprintf ('%d', G(1, :)), '100000011111');
  same ('row 1 of G3', sprintf ('%d', G(2, :)), '010000101221');
  same ('digit sums of the rows of G3, mod 3', join (mod (sum (G, 2), 3)), '0 2 2 2 2 2');
  same ('digit sum of row 1 of G3', sum (G(2, :)), 8);
  same ('check digit -(c0 + c7 + c8 + c9 + c10) for row 1', mod (-(G(2, 1) + sum (G(2, 8:11))), 3), 1);
  M = dec2base (0:728, 3, 6) - '0';
  C = mod (M * G, 3);
  require ('c11 = -(c0 + c7 + c8 + c9 + c10) on all 729 codewords', ...
           all (C(:, 12) == mod (-(C(:, 1) + sum (C(:, 8:11), 2)), 3)));
  Z = [C(:, 1:11) mod(-sum (C(:, 1:11), 2), 3)];
  require ('appending the zero-sum digit to G11 does not give G12', ...
           ~isequal (sortrows (unique (Z, 'rows')), sortrows (unique (C, 'rows'))));
  wz = sum (Z ~= 0, 2);
  same ('minimum distance of G11 with the zero-sum digit', min (wz(wz > 0)), 5);
  same ('codewords enumerated', size (unique (C, 'rows'), 1), 729);
  w = sum (C ~= 0, 2);
  require ('every codeword weight is divisible by 3', all (mod (w, 3) == 0));
  same ('codewords of weight 3', sum (w == 3), 0);
  same ('minimum distance of G12', min (w(w > 0)), 6);
  same ('weight distribution of G12', distribution (w), '0:1 6:264 9:440 12:24');
  same ('codewords that vanish at a given coordinate', join (unique (sum (C == 0, 1))), '243');
  six = zeros (1, 12);
  for i = 1:12
    six(i) = sum (w == 6 & C(:, i) ~= 0);
  end
  same ('weight-6 codewords nonzero at a given coordinate', join (unique (six)), '132');
  same ('C(11,3) * 2^3', nchoosek (11, 3) * 8, 1320);
  same ('C(11,3) * 2^3 / C(5,3)', nchoosek (11, 3) * 8 / nchoosek (5, 3), 132);
  same ('A_6 = 12 * 132 / 6', 12 * 132 / 6, 264);
  same ('zero entries at a coordinate: 243 - 1 - 132', 243 - 1 - 132, 110);
  same ('A_9 = 12 * 110 / 3', 12 * 110 / 3, 440);
  same ('A_12 = 729 - 1 - 264 - 440', 729 - 1 - 264 - 440, 24);
  v3 = 1 + 11 * 2 + nchoosek (11, 2) * 4;
  same ('V_3(11,2) = 1 + 22 + 220', v3, 243);
  same ('3^6 * V_3(11,2)', 729 * v3, 3 ^ 11);
  C11 = C(:, 1:11);
  w11 = sum (C11 ~= 0, 2);
  same ('distinct words after puncturing', size (unique (C11, 'rows'), 1), 729);
  same ('minimum distance of G11', min (w11(w11 > 0)), 5);
  same ('weight distribution of G11', distribution (w11), '0:1 5:132 6:132 8:330 9:110 11:24');
  same ('weight-9 codewords of G12 with a 0 in coordinate 11', sum (w == 9 & C(:, 12) == 0), 110);
  same ('A_8 of G11 = 440 - 110', 440 - 110, 330);
end

function check_decoder (CODE24)
  section ('Decoding G24');
  same ('patterns of weight 1, 2, 3', join ([nchoosek(24, 1) nchoosek(24, 2) nchoosek(24, 3)]), '24 276 2024');
  same ('patterns of weight 1 to 3 in all', nchoosek (24, 1) + nchoosek (24, 2) + nchoosek (24, 3), 2324);
  same ('patterns of weight 4', nchoosek (24, 4), 10626);
  T = zeros (4, 5);
  k = 0;
  for w = 1:4
    P = nchoosek (0:23, w);
    for row = 1:size (P, 1)
      e = sum (2 .^ (23 - P(row, :)));
      c = CODE24(mod (k, 4096) + 1);
      k = k + 1;
      [kase, found] = decode (bitxor (c, e));
      if (kase < 5 && found ~= e)
        T(w, 5) = T(w, 5) + 1000000;
      end
      T(w, kase) = T(w, kase) + 1;
    end
  end
  same ('weight 1 resolved by cases 1 to 4, then flagged', join (T(1, :)), '12 12 0 0 0');
  same ('weight 2 resolved by cases 1 to 4, then flagged', join (T(2, :)), '66 144 66 0 0');
  same ('weight 3 resolved by cases 1 to 4, then flagged', join (T(3, :)), '220 792 220 792 0');
  same ('weight 4 resolved by cases 1 to 4, then flagged', join (T(4, :)), '0 0 0 0 10626');
  same ('patterns of weight 1 to 3 corrected', sum (sum (T(1:3, 1:4))), 2324);
  same ('patterns of weight 1 to 3 by case', join (sum (T(1:3, 1:4), 1)), '298 948 286 792');
end

function check_examples ()
  section ('Worked decoding examples');
  R = brow ();
  names = {'case 1', 'case 2', 'case 3', 'case 4', 'weight 4'};
  messages = {'101000000000', '010010000001', '111100000000', '000000110001', '100000000000'};
  errors = {[3 7 10], [1 6 16], [12 17 22], [3 14 20], [0 5 13 18]};
  want = { ...
    {'101000000000101011010010', '101100010010', '101011010010', '101000000000', ...
     '000100010010', 3, '8 6 4 10 6 6 6 8 6 6 8 8', '101000010011', 5, ...
     '8 4 4 6 4 4 8 6 8 4 6 8', 1, '3 7 10'}, ...
    {'010010000001100000101110', '000010100001', '100010101110', '100011011010', ...
     '100001111011', 7, '8 6 8 6 2 4 4 6 6 4 6 6', '101001100100', 5, ...
     '6 8 8 4 8 6 6 6 10 4 6 6', 2, '1 6 16'}, ...
    {'111100000000111101110100', '111100000000', '011100110110', '110010010001', ...
     '001110010001', 5, '6 4 4 8 8 6 8 4 4 4 6 8', '100001000010', 3, ...
     '6 6 8 8 4 6 6 8 10 6 6 8', 3, '12 17 22'}, ...
    {'000000110001110001001100', '000100110001', '111001000100', '001010011111', ...
     '001110101110', 7, '8 4 6 6 10 6 4 6 6 8 6 4', '110000100101', 5, ...
     '4 8 6 2 6 4 6 8 6 6 6 8', 4, '3 14 20'}, ...
    {'100000000000110111000101', '000001000000', '100111100101', '001011100100', ...
     '001010100100', 4, '7 7 7 5 11 5 7 5 7 7 5 7', '000101010010', 4, ...
     '7 7 5 11 5 7 5 7 7 5 7 7', 5, 'none'}};
  for q = 1:5
    name = names{q};
    msg = messages{q};
    wq = want{q};
    a = bin2dec (msg);
    c = encode (a);
    r = bitxor (c, mask (errors{q}));
    s = syndrome (r);
    t = times_b (s);
    [kase, found] = decode (r);
    same (sprintf ('%s: codeword for message %s', name, msg), dec2bin (c, 24), wq{1});
    same (sprintf ('%s: received r_L', name), dec2bin (floor (r / 4096), 12), wq{2});
    same (sprintf ('%s: received r_R', name), dec2bin (bitand (r, 4095), 12), wq{3});
    same (sprintf ('%s: r_R B', name), dec2bin (times_b (bitand (r, 4095)), 12), wq{4});
    same (sprintf ('%s: syndrome s', name), dec2bin (s, 12), wq{5});
    same (sprintf ('%s: weight of s', name), pop (s), wq{6});
    same (sprintf ('%s: weights of s + b_i, i = 0 to 11', name), join (pop (bitxor (s, R))), wq{7});
    same (sprintf ('%s: s B', name), dec2bin (t, 12), wq{8});
    same (sprintf ('%s: weight of s B', name), pop (t), wq{9});
    same (sprintf ('%s: weights of s B + b_i, i = 0 to 11', name), join (pop (bitxor (t, R))), wq{10});
    same (sprintf ('%s: case that fires (5 = flagged)', name), kase, wq{11});
    if (found >= 0)
      positions = join (support (found));
    else
      positions = 'none';
    end
    same (sprintf ('%s: error positions found', name), positions, wq{12});
    if (found >= 0)
      same (sprintf ('%s: message recovered', name), dec2bin (floor (bitxor (r, found) / 4096), 12), msg);
    end
  end
end

function check_exercises (CODE24, GROWS)
  section ('Exercises');
  R = brow ();
  same ('row 3 of B', dec2bin (R(4), 12), '111000101101');
  same ('overlap of rows 0 and 3, all 12 columns', pop (bitand (R(1), R(4))), 4);
  g01 = bitxor (GROWS(1), GROWS(2));
  same ('weight of g_0 + g_1', pop (g01), 8);
  same ('support of g_0 + g_1', join (support (g01)), '0 1 13 14 17 20 21 22');
  lam = zeros (1, 4);
  for i = 1:4
    lam(i) = nchoosek (24 - i, 5 - i) / nchoosek (8 - i, 5 - i);
  end
  same ('octads through 1, 2, 3, 4 given points', join (lam), '253 77 21 5');
  r = mask ([5 23]);
  same ('syndrome of an error at coordinates 5 and 23', dec2bin (syndrome (r), 12), '111110111110');
  [kase, found] = decode (r);
  same ('case and positions found for that error', sprintf ('%d %s', kase, join (support (found))), '2 5 23');
  G = [eye(6) b3_matrix()];
  c3 = mod ([1 1 0 0 0 0] * G, 3);
  same ('ternary codeword for message 110000', sprintf ('%d', c3), '110000112002');
  same ('its weight', sum (c3 ~= 0), 6);
  st = gate_counts ();
  same ('XOR gates in the syndrome trees', st.syndrome_xor, 88);
  same ('XOR gates in the s B network', st.sb_xor, 76);
  same ('levels from r to s B, the two networks in series', st.series_depth, 8);
end

function check_hardware ()
  section ('Hardware gate counts');
  st = gate_counts ();
  same ('two-input XOR gates for the syndrome', st.syndrome_xor, 88);
  same ('levels in the deepest syndrome tree', st.syndrome_depth, 4);
  same ('two-input XOR gates for s B', st.sb_xor, 76);
  same ('levels in the deepest s B tree', st.sb_depth, 4);
  same ('XOR gates in one bank of 12 tests, built literally', st.literal_xor, 144);
  same ('inverters in one bank once the constants are folded', st.fold_inverters, 88);
  same ('inverters in one bank with the inverters shared', st.shared_inverters, 12);
  same ('AND and OR gates flagging weight at least 3 of 12', st.chain3, 51);
  same ('AND and OR gates flagging weight at least 4 of 12', st.chain4, 68);
  x = 0:4095;
  require ('both chains agree with a popcount on all 4096 inputs', ...
           isequal (chain_flag (x, 12, 3), pop (x) >= 3) && isequal (chain_flag (x, 12, 4), pop (x) >= 4));
  same ('AND and OR gates in the 24 weight-at-most-2 tests', 24 * st.chain3, 1224);
  same ('AND and OR gates in the 2 weight-at-most-3 tests', 2 * st.chain4, 136);
  same ('AND and OR gates in one bank of 13 tests', 12 * st.chain3 + st.chain4, 680);
  same ('two-input XOR gates in both networks', st.syndrome_xor + st.sb_xor, 164);
  same ('flag inversions in both banks, left to the select', st.flag_inversions, 26);
  same ('two-input XOR gates forming r + e', st.correction_xor, 24);
end


% Main -------------------------------------------------------------------------

printf ('Chapter 25, Golay Codes\n');
CODE24 = encode (0:4095);
GROWS = encode (2 .^ (11:-1:0));
check_b ();
check_g24 (CODE24, GROWS);
check_octads (CODE24);
check_g23 (CODE24);
check_ternary ();
check_decoder (CODE24);
check_examples ();
check_exercises (CODE24, GROWS);
check_hardware ();
printf ('\nAll chapter 25 numbers reproduced.\n');
