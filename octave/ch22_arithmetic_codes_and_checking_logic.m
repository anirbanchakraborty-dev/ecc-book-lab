% Reproduce every number in chapter 22 of Error Correction Coding,
% "Arithmetic Codes and Checking Logic", in GNU Octave.
%
% Run it from the repository root:
%
%     octave --no-window-system --norc --quiet octave/ch22_arithmetic_codes_and_checking_logic.m
%
% It needs only core Octave, with no pkg load and no Forge packages. It is
% the Octave twin of python/ch22_arithmetic_codes_and_checking_logic.py
% and prints the same lines, which tools/compare-routes.sh checks.
%
% Each check prints the value computed here beside the value the chapter
% states, and the program stops with exit status 1 at the first pair that
% disagrees. Every number in the chapter is a whole number except the one
% floating-point sum in the section on checksums for matrix products, which
% is printed to 17 significant digits on both routes.
%
% The exhaustive checks run on whole arrays at once rather than one case at
% a time: an adder is evaluated bit by bit across every operand pair in a
% column, so a stuck carry is injected on every pair together. Residues of
% powers of two are built by doubling modulo A, never from 2^i itself, so no
% value ever exceeds 2^53, where a double still holds every integer exactly.

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

function c = popcount (x)
  % The number of 1 bits in each non-negative whole number of an array.
  c = zeros (size (x));
  while (any (x(:) > 0))
    c = c + mod (x, 2);
    x = floor (x / 2);
  end
end

function p = par (x)
  p = mod (popcount (x), 2);
end

function s = bits_text (x, n)
  s = dec2bin (x, n);
end

function s = low_first (x, n)
  s = char ('0' + mod (floor (x ./ 2 .^ (0:n-1)), 2));
end

function s = join_nums (v)
  parts = arrayfun (@(t) sprintf ('%d', t), v, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function w = arith_weight (e)
  e = abs (e);
  if (e == 0)
    w = 0;
  elseif (e == 1)
    w = 1;
  elseif (mod (e, 2) == 0)
    w = arith_weight (e / 2);
  else
    w = 1 + min (arith_weight ((e - 1) / 2), arith_weight ((e + 1) / 2));
  end
end

function b = bit_length (v)
  b = 0;
  while (v > 0)
    v = floor (v / 2);
    b = b + 1;
  end
end

function [s, c] = ripple_add (x, y, n, c0, fj, fv)
  % Add columns of n-bit words in a ripple-carry adder, one full adder per bit.
  % Column j+1 of c holds the carry c_j. fj = 0 means no fault, otherwise
  % carry line c_fj is stuck at fv.
  x = x(:);
  y = y(:);
  c = zeros (numel (x), n + 1);
  c(:, 1) = c0;
  s = zeros (numel (x), 1);
  for i = 0:n-1
    xi = mod (floor (x / 2 ^ i), 2);
    yi = mod (floor (y / 2 ^ i), 2);
    ci = c(:, i + 1);
    s = s + bitxor (bitxor (xi, yi), ci) * 2 ^ i;
    cn = bitor (bitor (bitand (xi, yi), bitand (xi, ci)), bitand (yi, ci));
    if (fj == i + 1)
      cn(:) = fv;
    end
    c(:, i + 2) = cn;
  end
end

function v = result_value (s, c, n)
  v = s + c(:, n + 1) * 2 ^ n;
end


% A wrong carry in a 16-bit adder ----------------------------------------------

function check_motivation ()
  section ('A wrong carry in a 16-bit adder');
  x = 32760;
  y = 8;
  same ('x = 0x7FF8 in decimal', x, 32760);
  same ('x + y, fault-free', x + y, 32768);
  [s, c] = ripple_add (x, y, 16, 0, 0, 0);
  same ('ripple-carry adder output, fault-free', result_value (s, c, 16), 32768);
  [sf, cf] = ripple_add (x, y, 16, 0, 4, 0);
  bad = result_value (sf, cf, 16);
  same ('output with carry c_4 stuck at 0', bad, 32752);
  same ('output with c_4 stuck at 0, in hexadecimal', sprintf ('%04X', bad), '7FF0');
  same ('sum bits that differ from the right sum', popcount (bitxor (bad, 32768)), 12);
  same ('arithmetic error, faulty minus correct', bad - 32768, -16);
  same ('arithmetic weight of that error', arith_weight (bad - 32768), 1);
  same ('parity of x', par (x), 0);
  same ('parity of y', par (y), 1);
  same ('parity of the correct sum', par (32768), 1);
  same ('parity of the faulty sum', par (bad), 1);
  require ('XOR of operand parities agrees with the faulty sum', bitxor (par (x), par (y)) == par (bad));
  same ('1 + 1: XOR of the operand parities', bitxor (par (1), par (1)), 0);
  same ('1 + 1: parity of the sum 2', par (2), 1);
  same ('carries of 0x7FF8 + 8 that are 1, fault-free', sum (c(2:end)), 12);

  section ('Arithmetic weight');
  same ('arithmetic weight of 16', arith_weight (16), 1);
  same ('arithmetic weight of 7 = 8 - 1', arith_weight (7), 2);
  same ('arithmetic weight of 11', arith_weight (11), 3);
  same ('arithmetic weight of 4080 = 4096 - 16', arith_weight (4080), 2);
  same ('Hamming weight of 4080', popcount (4080), 8);
  same ('arithmetic weight of 255 = 256 - 1', arith_weight (255), 2);

  % One stuck carry adds or removes exactly one power of two.
  n = 6;
  [X, Y] = meshgrid (0:63, 0:63);
  X = X(:);
  Y = Y(:);
  cases = 0;
  wrong = 0;
  single = true;
  for j = 1:n
    for v = 0:1
      [sf, cf] = ripple_add (X, Y, n, 0, j, v);
      e = result_value (sf, cf, n) - (X + Y);
      cases = cases + numel (e);
      ew = e(e ~= 0);
      wrong = wrong + numel (ew);
      if (~all (ew == 2 ^ j | ew == -(2 ^ j)))
        single = false;
      end
    end
  end
  same ('6-bit adder: fault and operand cases tried', cases, 49152);
  same ('6-bit adder: cases with a wrong result', wrong, 24576);
  require ('every wrong result is off by exactly +2^j or -2^j', single);
end


% AN codes ---------------------------------------------------------------------

function k = order_of_two (a)
  k = 1;
  v = mod (2, a);
  while (v ~= 1)
    v = mod (v * 2, a);
    k = k + 1;
  end
end

function ok = correcting (a, m)
  % True when the residues of +2^i and -2^i, 0 <= i < m, are nonzero and distinct mod a.
  seen = [];
  r = mod (1, a);
  ok = true;
  for i = 0:m-1
    for e = [r, mod(a - r, a)]
      if (e == 0 || any (seen == e))
        ok = false;
        return;
      end
      seen(end + 1) = e;
    end
    r = mod (2 * r, a);
  end
end

function m = longest_correctable (a)
  % The largest m with correcting (a, m), found by extending one bit at a time.
  seen = [];
  r = mod (1, a);
  m = 0;
  while (true)
    p = r;
    q = mod (a - r, a);
    if (p == 0 || q == 0 || p == q || any (seen == p) || any (seen == q))
      return;
    end
    seen = [seen, p, q];
    m = m + 1;
    r = mod (2 * r, a);
  end
end

function l = order_limit (a)
  % The longest correctable word the order theorem allows for odd a.
  o = order_of_two (a);
  minus_one = false;
  v = mod (1, a);
  for j = 0:o-1
    if (v == a - 1)
      minus_one = true;
    end
    v = mod (2 * v, a);
  end
  if (minus_one)
    l = o / 2;
  else
    l = o;
  end
end

function check_an_codes ()
  section ('AN codes');
  same ('A = 3, x = 5: stored word 3x', 3 * 5, 15);
  same ('15 in binary', bits_text (15, 5), '01111');
  same ('15 + 4 after one arithmetic error', 15 + 4, 19);
  same ('19 mod 3', mod (19, 3), 1);
  same ('A = 4, x = 5: stored word 4x', 4 * 5, 20);
  same ('20 + 4 mod 4, an undetected error', mod (20 + 4, 4), 0);
  same ('codeword bits for A = 3 and 8-bit x', bit_length (3 * 255), 10);
  [Xs, I] = meshgrid (0:255, 0:9);
  W = [3 * Xs(:) + 2 .^ I(:); 3 * Xs(:) - 2 .^ I(:)];
  same ('A = 3: single errors tried on all 8-bit x', numel (W), 5120);
  same ('A = 3: single errors that go undetected', sum (mod (W, 3) == 0), 0);

  pos = mod (2 .^ (0:8), 19);
  neg = mod (-(2 .^ (0:8)), 19);
  same ('residues of +2^i mod 19, i = 0..8', join_nums (pos), '1 2 4 8 16 13 7 14 9');
  same ('residues of -2^i mod 19, i = 0..8', join_nums (neg), '18 17 15 11 3 6 12 5 10');
  same ('distinct nonzero residues for A = 19, m = 9', numel (unique ([pos, neg])), 18);
  require ('A = 19 corrects every single error in 9 bits', correcting (19, 9));
  require ('A = 19 does not correct every single error in 10 bits', ~correcting (19, 10));
  same ('order of 2 modulo 19', order_of_two (19), 18);
  same ('2^9 mod 19, which is -1', mod (2 ^ 9, 19), 18);
  same ('x = 11: stored word 19x', 19 * 11, 209);
  same ('209 + 32 after an error in bit 5', 209 + 32, 241);
  same ('241 mod 19', mod (241, 19), 13);
  same ('2^5 mod 19', mod (2 ^ 5, 19), 13);
  same ('corrected word (241 - 32) / 19', floor ((241 - 32) / 19), 11);

  agree = true;
  for a = 3:2:99
    if (longest_correctable (a) ~= order_limit (a))
      agree = false;
    end
  end
  require ('order theorem matches a search for every odd A below 100', agree);
  same ('order of 2 modulo 7', order_of_two (7), 3);
  same ('longest correctable word for A = 7', order_limit (7), 3);
  same ('order of 2 modulo 23', order_of_two (23), 11);
  same ('longest correctable word for A = 23', order_limit (23), 11);

  ks = [4, 8, 12, 16, 32];
  expect = [19, 9, 18, 9; 29, 13, 28, 14; 37, 18, 36, 18; 47, 22, 23, 23; 79, 39, 39, 39];
  for t = 1:numel (ks)
    k = ks(t);
    a = 3;
    while (true)
      m = bit_length (a * (2 ^ k - 1));
      if (mod (a, 2) == 1 && correcting (a, m))
        break;
      end
      a = a + 1;
    end
    same (sprintf ('k = %d: smallest single-error-correcting A', k), a, expect(t, 1));
    same (sprintf ('k = %d: codeword length m', k), m, expect(t, 2));
    same (sprintf ('k = %d: check bits m - k', k), m - k, expect(t, 2) - k);
    same (sprintf ('k = %d: order of 2 modulo A', k), order_of_two (a), expect(t, 3));
    same (sprintf ('k = %d: longest correctable length for that A', k), order_limit (a), expect(t, 4));
  end
end


% Residue codes ----------------------------------------------------------------

function x = fold (x, a)
  % Residue modulo 2^a - 1 by summing a-bit chunks with end-around carry, on a whole array.
  mask = 2 ^ a - 1;
  while (any (x(:) > mask))
    total = zeros (size (x));
    while (any (x(:) > 0))
      total = total + mod (x, 2 ^ a);
      x = floor (x / 2 ^ a);
    end
    x = total;
  end
  x(x == mask) = 0;
end

function v = chunks (x, a, n)
  v = mod (floor (x ./ (2 ^ a) .^ (0:(n / a) - 1)), 2 ^ a);
end

function [c1, c0] = canon (p1, p0)
  % A raw bit pair turned into the residue code 00, 01, 10.
  t = bitxor (p1, p0);
  c1 = bitand (p1, t);
  c0 = bitand (p0, t);
end

function [r1, r0] = mod3_add (a1, a0, b1, b0)
  % The ten-gate mod-3 adder on canonical residues.
  za = 1 - bitor (a1, a0);
  zb = 1 - bitor (b1, b0);
  r0 = bitor (bitor (bitand (za, b0), bitand (a0, zb)), bitand (a1, b1));
  r1 = bitor (bitor (bitand (za, b1), bitand (a1, zb)), bitand (a0, b0));
end

function r = residue_tree (x)
  x = x(:);
  hi = cell (1, 8);
  lo = cell (1, 8);
  for j = 0:7
    [hi{j + 1}, lo{j + 1}] = canon (mod (floor (x / 2 ^ (2 * j + 1)), 2), mod (floor (x / 2 ^ (2 * j)), 2));
  end
  while (numel (hi) > 1)
    nh = {};
    nl = {};
    for k = 1:2:numel (hi)
      [h, l] = mod3_add (hi{k}, lo{k}, hi{k + 1}, lo{k + 1});
      nh{end + 1} = h;
      nl{end + 1} = l;
    end
    hi = nh;
    lo = nl;
  end
  r = 2 * hi{1} + lo{1};
end

function [kinds, ins, id] = add_gate (kinds, ins, kind, inputs)
  % Append one gate. Kinds: 1 XOR, 2 AND, 3 NOR, 4 OR. Nets 1 to 16 are x0 to x15.
  kinds(end + 1) = kind;
  ins{end + 1} = inputs;
  id = 16 + numel (kinds);
end

function [kinds, ins, out1, out0] = mod3_tree_netlist ()
  % The 16-bit mod-3 generator as a gate list in build order.
  kinds = [];
  ins = {};
  level = zeros (8, 2);
  for j = 0:7
    p1 = 2 * j + 2;
    p0 = 2 * j + 1;
    [kinds, ins, t] = add_gate (kinds, ins, 1, [p1 p0]);
    [kinds, ins, c1] = add_gate (kinds, ins, 2, [p1 t]);
    [kinds, ins, c0] = add_gate (kinds, ins, 2, [p0 t]);
    level(j + 1, :) = [c1 c0];
  end
  while (rows (level) > 1)
    nxt = zeros (rows (level) / 2, 2);
    for k = 1:2:rows (level)
      a1 = level(k, 1);
      a0 = level(k, 2);
      b1 = level(k + 1, 1);
      b0 = level(k + 1, 2);
      [kinds, ins, za] = add_gate (kinds, ins, 3, [a1 a0]);
      [kinds, ins, zb] = add_gate (kinds, ins, 3, [b1 b0]);
      [kinds, ins, g1] = add_gate (kinds, ins, 2, [za b0]);
      [kinds, ins, g2] = add_gate (kinds, ins, 2, [a0 zb]);
      [kinds, ins, g3] = add_gate (kinds, ins, 2, [a1 b1]);
      [kinds, ins, r0] = add_gate (kinds, ins, 4, [g1 g2 g3]);
      [kinds, ins, g4] = add_gate (kinds, ins, 2, [za b1]);
      [kinds, ins, g5] = add_gate (kinds, ins, 2, [a1 zb]);
      [kinds, ins, g6] = add_gate (kinds, ins, 2, [a0 b0]);
      [kinds, ins, r1] = add_gate (kinds, ins, 4, [g4 g5 g6]);
      nxt((k + 1) / 2, :) = [r1 r0];
    end
    level = nxt;
  end
  out1 = level(1, 1);
  out0 = level(1, 2);
end

function check_residue_codes ()
  section ('Residue codes');
  same ('22 mod 3', mod (22, 3), 1);
  same ('13 mod 3', mod (13, 3), 1);
  same ('22 + 13', 22 + 13, 35);
  same ('35 mod 3', mod (35, 3), 2);
  same ('(1 + 1) mod 3', mod (1 + 1, 3), 2);
  same ('22 times 13', 22 * 13, 286);
  same ('286 mod 3', mod (286, 3), 1);
  same ('(1 times 1) mod 3', mod (1 * 1, 3), 1);
  [X, Y] = meshgrid (0:63, 0:63);
  ok = true;
  for a = 2:16
    ok = ok && all (all (mod (X + Y, a) == mod (mod (X, a) + mod (Y, a), a))) ...
            && all (all (mod (X .* Y, a) == mod (mod (X, a) .* mod (Y, a), a)));
  end
  require ('sum and product residues agree for A = 2 to 16 on all 6-bit pairs', ok);

  same ('181 in binary', bits_text (181, 8), '10110101');
  same ('2-bit chunks of 181, lowest first', join_nums (chunks (181, 2, 8)), '1 1 3 2');
  same ('sum of those chunks', sum (chunks (181, 2, 8)), 7);
  same ('7 folded once, 1 + 3', 1 + 3, 4);
  same ('4 folded once more, 1 + 0', 1 + 0, 1);
  same ('181 mod 3', mod (181, 3), 1);
  same ('fold of 181 for A = 3', fold (181, 2), 1);
  same ('4-bit chunks of 181, lowest first', join_nums (chunks (181, 4, 8)), '5 11');
  same ('5 + 11', 5 + 11, 16);
  same ('181 mod 15', mod (181, 15), 1);
  same ('fold of 181 for A = 15', fold (181, 4), 1);
  xs = 0:4095;
  ok = true;
  for a = 2:4
    if (~all (fold (xs, a) == mod (xs, 2 ^ a - 1)))
      ok = false;
    end
  end
  require ('fold equals x mod 2^a - 1 for a = 2, 3, 4 on 12-bit x', ok);

  [kinds, ins, out1, out0] = mod3_tree_netlist ();
  same ('mod-3 tree netlist: C blocks, one XOR each', sum (kinds == 1), 8);
  same ('mod-3 tree netlist: mod-3 adders, two NOR each', sum (kinds == 3) / 2, 7);
  same ('mod-3 tree netlist: XOR gates', sum (kinds == 1), 8);
  same ('mod-3 tree netlist: AND gates', sum (kinds == 2), 58);
  same ('mod-3 tree netlist: NOR gates', sum (kinds == 3), 14);
  same ('mod-3 tree netlist: OR gates', sum (kinds == 4), 14);
  same ('mod-3 tree netlist: gates in all', numel (kinds), 94);
  depth = zeros (1, 16 + numel (kinds));
  for g = 1:numel (kinds)
    depth(16 + g) = 1 + max (depth(ins{g}));
  end
  same ('mod-3 tree netlist: longest path, in gates', max (depth(out1), depth(out0)), 11);
  xs = (0:65535)';
  V = false (numel (xs), 16 + numel (kinds));
  for i = 0:15
    V(:, i + 1) = mod (floor (xs / 2 ^ i), 2) == 1;
  end
  for g = 1:numel (kinds)
    a = ins{g};
    switch (kinds(g))
      case 1
        V(:, 16 + g) = xor (V(:, a(1)), V(:, a(2)));
      case 2
        V(:, 16 + g) = V(:, a(1)) & V(:, a(2));
      case 3
        V(:, 16 + g) = ~(V(:, a(1)) | V(:, a(2)));
      otherwise
        V(:, 16 + g) = any (V(:, a), 2);
    end
  end
  bad = (V(:, out1) ~= (mod (xs, 3) == 2)) | (V(:, out0) ~= (mod (xs, 3) == 1));
  same ('mod-3 tree netlist: wrong residues over all 65,536 inputs', sum (bad), 0);
  same ('mod-3 tree: wrong residues over all 65,536 inputs', sum (residue_tree (xs) ~= mod (xs, 3)), 0);

  same ('2^16 mod 3', mod (2 ^ 16, 3), 1);
  r = residue_tree (xs);
  ok = true;
  for c16 = 0:1
    [r1, r0] = mod3_add (floor (r / 2), mod (r, 2), 0, c16);
    ok = ok && all (2 * r1 + r0 == mod (xs + c16 * 2 ^ 16, 3));
  end
  require ('tree residue of s, plus c16 in one more mod-3 adder, is r mod 3', ok);
  for n = [6 8]
    top = 2 ^ n;
    [X, Y] = meshgrid (0:top - 1, 0:top - 1);
    alarms = sum (sum (mod (mod (X, 3) + mod (Y, 3), 3) ~= mod (mod (X + Y, top), 3)));
    carries = sum (sum (X + Y >= top));
    same (sprintf ('%d-bit pairs with a carry-out', n), carries, top * (top - 1) / 2);
    same (sprintf ('%d-bit false alarms when the carry-out is ignored', n), alarms, carries);
  end
  pairs = sum (0:65535);
  same ('16-bit pairs with a carry-out, the sum of x over all x', pairs, 2147450880);
  same ('16-bit pairs with a carry-out, 2^16 (2^16 - 1) / 2', 2 ^ 16 * (2 ^ 16 - 1) / 2, 2147450880);
  same ('share of all 2^32 pairs, in ten-millionths', round (pairs * 10 ^ 7 / 2 ^ 32), 4999924);
  same ('mod-3 tree applied to 32768', residue_tree (32768), 2);
  same ('mod-3 tree applied to 32752', residue_tree (32752), 1);
  same ('(32760 + 8) residues: (0 + 2) mod 3', mod (mod (32760, 3) + mod (8, 3), 3), 2);
end


% The 4-bit adder under six carry faults ---------------------------------------

function b = berger (w, n)
  % The Berger check symbol of an n-bit word: its number of 0 bits.
  b = n - popcount (w);
end

function check_adder_faults ()
  section ('Checking a 4-bit adder under six carry faults');
  [s, c] = ripple_add (3, 1, 4, 0, 0, 0);
  same ('3 + 1, fault-free', result_value (s, c, 4), 4);
  same ('carries c_0..c_4 of 3 + 1', join_nums (c), '0 1 1 0 0');
  [sf, cf] = ripple_add (3, 1, 4, 0, 1, 0);
  same ('3 + 1 with c_1 stuck at 0', result_value (sf, cf, 4), 2);
  same ('carries c_0..c_4 with c_1 stuck at 0', join_nums (cf), '0 0 0 0 0');
  same ('residue check: 2 mod 3', mod (2, 3), 2);
  same ('residue check: (3 + 1) mod 3 from the residues', mod (mod (3, 3) + mod (1, 3), 3), 1);
  same ('parity predicted from the fault-free carries', mod (par (3) + par (1) + sum (c(1:4)), 2), 1);
  same ('parity of the faulty sum 0010', par (sf), 1);

  n = 4;
  [X, Y] = meshgrid (0:15, 0:15);
  X = X(:);
  Y = Y(:);
  faults = [1, 0; 1, 1; 2, 0; 2, 1; 3, 0; 3, 1];
  stated = [64, 64, 0, 48, 64, 64; 192, 192, 0, 144, 192, 192; 96, 96, 0, 48, 96, 96; ...
            160, 160, 0, 80, 160, 160; 112, 112, 0, 112, 112, 112; 144, 144, 0, 144, 144, 144];
  [s, c] = ripple_add (X, Y, n, 0, 0, 0);
  good = result_value (s, c, n);
  totals = zeros (1, 6);
  for t = 1:6
    [sf, cf] = ripple_add (X, Y, n, 0, faults(t, 1), faults(t, 2));
    got = result_value (sf, cf, n);
    w = got ~= good;
    mod3 = mod (got, 3) ~= mod (mod (X, 3) + mod (Y, 3), 3);
    ps = par (sf);
    shared = mod (par (X) + par (Y) + sum (cf(:, 1:n), 2), 2);
    dup = mod (par (X) + par (Y) + sum (c(:, 1:n), 2), 2);
    cmp = (ps ~= dup) | any (cf ~= c, 2);
    predicted = berger (X, n) + berger (Y, n) - (n - sum (cf(:, 2:end), 2)) - cf(:, 1) + cf(:, n + 1);
    actual = berger (sf, n);
    row = [sum(w), sum(w & mod3), sum(w & (ps ~= shared)), sum(w & (ps ~= dup)), sum(w & cmp), sum(w & (predicted ~= actual))];
    gap = all (abs (predicted(w) - actual(w)) == 2);
    name = sprintf ('c_%d stuck at %d', faults(t, 1), faults(t, 2));
    same ([name ': wrong results'], row(1), stated(t, 1));
    same ([name ': caught by the mod-3 residue'], row(2), stated(t, 2));
    same ([name ': caught by parity, shared carries'], row(3), stated(t, 3));
    same ([name ': caught by parity, duplicate carries'], row(4), stated(t, 4));
    same ([name ': caught by duplicate carries and compare'], row(5), stated(t, 5));
    same ([name ': caught by Berger check prediction'], row(6), stated(t, 6));
    require ([name ': Berger prediction off by exactly 2'], gap);
    totals = totals + row;
  end
  same ('fault and operand cases in all', 6 * 256, 1536);
  same ('wrong results in all', totals(1), 768);
  same ('caught by the mod-3 residue', totals(2), 768);
  same ('caught by parity prediction, shared carries', totals(3), 0);
  same ('caught by parity prediction, duplicate carries', totals(4), 576);
  same ('escapes with duplicate carries alone', totals(1) - totals(4), 192);
  same ('caught by duplicate carries and comparison', totals(5), 768);
  same ('caught by Berger check prediction', totals(6), 768);
end


% Berger check prediction and parity prediction --------------------------------

function check_predictions ()
  section ('Berger check prediction and parity prediction');
  x = 11;
  y = 6;
  n = 4;
  [s, c] = ripple_add (x, y, n, 0, 0, 0);
  same ('11 + 6 = 17, sum bits', bits_text (s, 4), '0001');
  same ('11 + 6, carry-out c_4', c(5), 1);
  same ('carries c_1..c_4 of 11 + 6', join_nums (c(2:end)), '0 1 1 1');
  same ('ones in x = 1011', popcount (x), 3);
  same ('ones in y = 0110', popcount (y), 2);
  same ('ones among c_1..c_4', sum (c(2:end)), 3);
  same ('ones in the sum, 3 + 2 - 3 + 0 - 1', 3 + 2 - 3 + 0 - 1, 1);
  same ('ones in the sum, counted', popcount (s), 1);
  same ('Berger symbol of x', berger (x, n), 1);
  same ('Berger symbol of y', berger (y, n), 2);
  same ('Berger symbol of the carries', n - sum (c(2:end)), 1);
  same ('predicted Berger symbol, 1 + 2 - 1 - 0 + 1', 1 + 2 - 1 - 0 + 1, 3);
  same ('Berger symbol of the sum, counted', berger (s, n), 3);
  same ('parity of x', par (x), 1);
  same ('parity of y', par (y), 0);
  same ('parity of the carries c_0..c_3', mod (sum (c(1:4)), 2), 0);
  same ('parity of the sum', par (s), 1);

  [X, Y] = meshgrid (0:15, 0:15);
  X = X(:);
  Y = Y(:);
  count = 0;
  ones_ok = true;
  berger_ok = true;
  parity_ok = true;
  for c0 = 0:1
    [s, c] = ripple_add (X, Y, n, c0, 0, 0);
    count = count + numel (X);
    nc = sum (c(:, 2:end), 2);
    if (~all (popcount (s) == popcount (X) + popcount (Y) - nc + c0 - c(:, n + 1)))
      ones_ok = false;
    end
    if (~all (berger (s, n) == berger (X, n) + berger (Y, n) - (n - nc) - c0 + c(:, n + 1)))
      berger_ok = false;
    end
    if (~all (par (s) == mod (par (X) + par (Y) + sum (c(:, 1:n), 2), 2)))
      parity_ok = false;
    end
  end
  same ('4-bit adder inputs, carry-in 0 and 1', count, 512);
  require ('ones in the sum match the count identity', ones_ok);
  require ('Berger symbol of the sum matches the prediction', berger_ok);
  require ('parity of the sum matches the prediction', parity_ok);
end


% Self-checking circuits -------------------------------------------------------

function v = fl (name, value, fname, fval)
  % The value on a line, or the stuck value when the fault sits on it.
  if (strcmp (name, fname))
    v = fval;
  else
    v = value;
  end
end

function [z, zn] = two_rail (a, an, b, bn, fname, fval)
  % The two-rail checker cell: g1 = a b, g2 = a' b', g3 = a b', g4 = a' b,
  % z = g1 + g2, z' = g3 + g4. fname = '' means no fault.
  a = fl ('a', a, fname, fval);
  an = fl ('a''', an, fname, fval);
  b = fl ('b', b, fname, fval);
  bn = fl ('b''', bn, fname, fval);
  g1 = fl ('g1', fl ('a@1', a, fname, fval) & fl ('b@1', b, fname, fval), fname, fval);
  g2 = fl ('g2', fl ('a''@2', an, fname, fval) & fl ('b''@2', bn, fname, fval), fname, fval);
  g3 = fl ('g3', fl ('a@3', a, fname, fval) & fl ('b''@3', bn, fname, fval), fname, fval);
  g4 = fl ('g4', fl ('a''@4', an, fname, fval) & fl ('b@4', b, fname, fval), fname, fval);
  z = double (fl ('z', g1 | g2, fname, fval));
  zn = double (fl ('z''', g3 | g4, fname, fval));
end

function check_two_rail ()
  section ('The two-rail checker cell');
  codes = [0, 0; 0, 1; 1, 0; 1, 1];
  outs = cell (1, 4);
  for k = 1:4
    [z, zn] = two_rail (codes(k, 1), 1 - codes(k, 1), codes(k, 2), 1 - codes(k, 2), '', 0);
    outs{k} = sprintf ('%d%d', z, zn);
  end
  same ('code inputs (a, b)', 4, 4);
  same ('outputs (z, z'') for a b = 00 01 10 11', strjoin (outs, ' '), '10 01 01 10');
  noncode = 0;
  disjoint = true;
  for a = 0:1
    for an = 0:1
      for b = 0:1
        for bn = 0:1
          if (an == 1 - a && bn == 1 - b)
            continue;
          end
          noncode = noncode + 1;
          [z, zn] = two_rail (a, an, b, bn, '', 0);
          if (z == 1 - zn)
            disjoint = false;
          end
        end
      end
    end
  end
  same ('non-code inputs', noncode, 12);
  require ('every non-code input gives a non-code output', disjoint);

  lines = {'a', 'a''', 'b', 'b''', 'a@1', 'b@1', 'a''@2', 'b''@2', 'a@3', 'b''@3', 'a''@4', 'b@4', ...
           'g1', 'g2', 'g3', 'g4', 'z', 'z'''};
  same ('fault sites (lines) in the cell', numel (lines), 18);
  same ('single stuck-at faults', 2 * numel (lines), 36);
  secure = true;
  testable = true;
  per_input = zeros (1, 4);
  sole = zeros (1, 4);
  for li = 1:numel (lines)
    for v = 0:1
      detecting = [];
      for k = 1:4
        a = codes(k, 1);
        b = codes(k, 2);
        [gz, gzn] = two_rail (a, 1 - a, b, 1 - b, '', 0);
        [z, zn] = two_rail (a, 1 - a, b, 1 - b, lines{li}, v);
        if (z == gz && zn == gzn)
          continue;
        end
        if (z == 1 - zn)
          secure = false;
        else
          detecting(end + 1) = k;
        end
      end
      if (isempty (detecting))
        testable = false;
      end
      per_input(detecting) = per_input(detecting) + 1;
      if (numel (detecting) == 1)
        sole(detecting) = sole(detecting) + 1;
      end
    end
  end
  require ('fault-secure: no fault yields a wrong code word', secure);
  require ('self-testing: some code input detects every fault', testable);
  same ('faults detected by a b = 00, 01, 10, 11', join_nums (per_input), '13 13 13 13');
  same ('faults that only a b = 00, 01, 10, 11 detect', join_nums (sole), '5 5 5 5');
  hits = {};
  for k = 1:4
    a = codes(k, 1);
    b = codes(k, 2);
    [gz, gzn] = two_rail (a, 1 - a, b, 1 - b, '', 0);
    [z, zn] = two_rail (a, 1 - a, b, 1 - b, 'g1', 0);
    if (z ~= gz || zn ~= gzn)
      hits{end + 1} = sprintf ('%d%d', a, b);
    end
  end
  same ('code input that detects g1 stuck at 0, as a b', strjoin (hits, ' '), '11');
  [z, zn] = two_rail (1, 0, 0, 1, 'g1', 1);
  same ('a b = 10 with g1 stuck at 1: (z, z'')', sprintf ('%d%d', z, zn), '11');
end

function [f, g] = parity_checker (w, fname, fval)
  % Two parity trees: f = (x0+x1)+(x2+x3), g = NOT((x4+x5)+(x6+x7)), + meaning XOR.
  x = zeros (1, 8);
  for i = 0:7
    x(i + 1) = fl (sprintf ('x%d', i), mod (floor (w / 2 ^ i), 2), fname, fval);
  end
  p1 = fl ('p1', bitxor (x(1), x(2)), fname, fval);
  p2 = fl ('p2', bitxor (x(3), x(4)), fname, fval);
  p3 = fl ('p3', bitxor (x(5), x(6)), fname, fval);
  p4 = fl ('p4', bitxor (x(7), x(8)), fname, fval);
  f = fl ('f', bitxor (p1, p2), fname, fval);
  g = fl ('g', 1 - bitxor (p3, p4), fname, fval);
end

function check_parity_checker ()
  section ('A totally self-checking parity checker');
  words = 0:255;
  code = words(par (words) == 0);
  odd = words(par (words) == 1);
  lines = {'x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'p1', 'p2', 'p3', 'p4', 'f', 'g'};
  same ('even-parity code words of 8 bits', numel (code), 128);
  same ('XOR or XNOR gates in the two trees', 6, 6);
  same ('fault sites (lines)', numel (lines), 14);
  same ('single stuck-at faults', 2 * numel (lines), 28);
  disjoint = true;
  for w = odd
    [f, g] = parity_checker (w, '', 0);
    if (f ~= g)
      disjoint = false;
    end
  end
  require ('every odd-parity word gives a non-code output', disjoint);
  good = zeros (128, 2);
  coded = true;
  for i = 1:128
    [f, g] = parity_checker (code(i), '', 0);
    good(i, :) = [f, g];
    if (f == g)
      coded = false;
    end
  end
  require ('every code word gives a code output', coded);
  secure = true;
  detect = false (28, 128);
  row = 0;
  for li = 1:numel (lines)
    for v = 0:1
      row = row + 1;
      for i = 1:128
        [f, g] = parity_checker (code(i), lines{li}, v);
        changed = f ~= good(i, 1) || g ~= good(i, 2);
        if (changed && f ~= g)
          secure = false;
        end
        detect(row, i) = changed;
      end
    end
  end
  require ('fault-secure: no fault yields a wrong code word', secure);
  require ('self-testing: some code word detects every fault', all (any (detect, 2)));
  pair_covers = false;
  for i = 1:127
    if (any (all (detect(:, i) | detect(:, i + 1:end), 1)))
      pair_covers = true;
    end
  end
  require ('no two code words detect all 28 faults', ~pair_covers);
  tests = [0, 170, 221];
  texts = arrayfun (@(w) low_first (w, 8), tests, 'UniformOutput', false);
  same ('test words, written x_0 x_1 ... x_7', strjoin (texts, ' '), '00000000 01010101 10111011');
  require ('those three words are code words', all (par (tests) == 0));
  idx = arrayfun (@(w) find (code == w), tests);
  require ('those three code words detect all 28 faults', all (any (detect(:, idx), 2)));
end


% Checksums for matrix products ------------------------------------------------

function s = mat_text (m)
  parts = cell (1, size (m, 1));
  for i = 1:size (m, 1)
    parts{i} = join_nums (m(i, :));
  end
  s = strjoin (parts, '; ');
end

function check_abft ()
  section ('Checksums for matrix products');
  A = [2, 1, 0; 1, 3, 1; 0, 2, 4];
  B = [1, 0, 2; 3, 1, 0; 2, 1, 1];
  C = A * B;
  same ('C = AB', mat_text (C), '5 1 4; 12 4 3; 14 6 4');
  col_a = sum (A, 1);
  row_b = sum (B, 2)';
  same ('checksum row of A', join_nums (col_a), '3 6 5');
  same ('checksum column of B', join_nums (row_b), '3 4 4');
  Cf = [A; col_a] * [B, row_b'];
  same ('full checksum matrix, row by row', mat_text (Cf), '5 1 4 10; 12 4 3 19; 14 6 4 24; 31 11 11 53');
  require ('last column holds the row sums of C', all (Cf(1:3, 4) == sum (C, 2)));
  require ('last row holds the column sums of C', all (Cf(4, 1:3) == sum (C, 1)));
  same ('corner entry, (3 6 5) times (3 4 4)', 3 * 3 + 6 * 4 + 5 * 4, 53);

  bad = Cf;
  bad(2, 2) = 9;
  row_gap = sum (bad(1:3, 1:3), 2)' - bad(1:3, 4)';
  col_gap = sum (bad(1:3, 1:3), 1) - bad(4, 1:3);
  same ('row sum minus checksum, rows 1 to 3', join_nums (row_gap), '0 5 0');
  same ('column sum minus checksum, columns 1 to 3', join_nums (col_gap), '0 5 0');
  same ('row 2 sum with the error, 12 + 9 + 3', 12 + 9 + 3, 24);
  same ('column 2 sum with the error, 1 + 9 + 6', 1 + 9 + 6, 16);
  same ('corrected entry 9 - 5', 9 - 5, 4);

  fixed = true;
  for i = 1:3
    for j = 1:3
      for e = [-7, -3, -1, 1, 2, 5, 100]
        m = Cf;
        m(i, j) = m(i, j) + e;
        rg = sum (m(1:3, 1:3), 2)' - m(1:3, 4)';
        cg = sum (m(1:3, 1:3), 1) - m(4, 1:3);
        if (~isequal (find (rg ~= 0), i) || ~isequal (find (cg ~= 0), j) || rg(i) ~= cg(j))
          fixed = false;
          continue;
        end
        m(i, j) = m(i, j) - rg(i);
        if (~isequal (m, Cf))
          fixed = false;
        end
      end
    end
  end
  require ('one wrong entry is located and corrected at all 9 places', fixed);

  left = (0.1 + 0.2) + 0.3;
  right = 0.1 + (0.2 + 0.3);
  same ('(0.1 + 0.2) + 0.3 in double precision', sprintf ('%.17g', left), '0.60000000000000009');
  same ('0.1 + (0.2 + 0.3) in double precision', sprintf ('%.17g', right), '0.59999999999999998');
  same ('their difference', sprintf ('%.17g', left - right), '1.1102230246251565e-16');
  same ('2^-53, the spacing of doubles from 0.5 to 1', sprintf ('%.17g', 2 ^ -53), '1.1102230246251565e-16');
end


% The worked examples and the exercises ---------------------------------------

function v = powers_mod (b, k, a)
  % b^0, b^1, ..., b^(k-1) modulo a, by repeated multiplication.
  v = zeros (1, k);
  t = 1;
  for i = 1:k
    v(i) = t;
    t = mod (t * b, a);
  end
end

function s = carry_names (idx)
  parts = arrayfun (@(j) sprintf ('c_%d', j), idx, 'UniformOutput', false);
  s = strjoin (parts, ' ');
end

function check_examples ()
  section ('The worked examples and the exercises');
  same ('A = 19, N = 11: stored word 19N', 19 * 11, 209);
  same ('209 in binary', bits_text (209, 8), '11010001');
  same ('209 + 32 after the error', 209 + 32, 241);
  same ('19 times 12', 19 * 12, 228);
  same ('241 mod 19', mod (241, 19), 13);
  errs = reshape ([2 .^ (0:8); -(2 .^ (0:8))], 1, []);
  same ('single errors +2^i, -2^i, i < 9, leaving remainder 13', join_nums (errs(mod (errs, 19) == 13)), '32');
  require ('A = 19 corrects every single error in 9-bit words', correcting (19, 9));
  same ('other nonzero remainders modulo 19', 18 - 1, 17);
  same ('(241 - 32) / 19', floor ((241 - 32) / 19), 11);

  same ('powers 2^0, 2^1, 2^2 modulo 7', join_nums (powers_mod (2, 3, 7)), '1 2 4');
  same ('order of 2 modulo 7', order_of_two (7), 3);
  same ('longest correctable length for A = 7', order_limit (7), 3);
  same ('powers 2^0..2^10 modulo 23', join_nums (powers_mod (2, 11, 23)), '1 2 4 8 16 9 18 13 3 6 12');
  same ('2^11', 2 ^ 11, 2048);
  same ('89 times 23 plus 1', 89 * 23 + 1, 2048);
  same ('order of 2 modulo 23', order_of_two (23), 11);
  same ('longest correctable length for A = 23', order_limit (23), 11);
  same ('remainders of 1, -1, 2, -2, 4, -4 modulo 7, sorted', join_nums (sort (mod ([1 -1 2 -2 4 -4], 7))), '1 2 3 4 5 6');

  same ('32752 in hexadecimal', dec2hex (32752), '7FF0');
  ch = chunks (32752, 2, 16);
  same ('2-bit chunks of 32752, lowest first', join_nums (ch), '0 0 3 3 3 3 3 1');
  [h, l] = canon (floor (ch / 2), mod (ch, 2));
  same ('those chunks after the C blocks', join_nums (2 * h + l), '0 0 0 0 0 0 0 1');
  same ('mod-3 tree applied to 32752', residue_tree (32752), 1);
  same ('32768 in hexadecimal', dec2hex (32768), '8000');
  same ('2-bit chunks of 32768, lowest first', join_nums (chunks (32768, 2, 16)), '0 0 0 0 0 0 0 2');
  same ('mod-3 tree applied to 32768', residue_tree (32768), 2);
  same ('3 times 10920', 3 * 10920, 32760);
  same ('32760 mod 3', mod (32760, 3), 0);
  same ('8 mod 3', mod (8, 3), 2);

  [s, c] = ripple_add (3, 1, 4, 0, 1, 0);
  [s_ok, c_ok] = ripple_add (3, 1, 4, 0, 0, 0);
  same ('3 + 1 with c_1 stuck at 0: sum bits', bits_text (s, 4), '0010');
  same ('3 + 1 with c_1 stuck at 0: carries c_0..c_4', join_nums (c), '0 0 0 0 0');
  same ('3 + 1, healthy carries c_0..c_4', join_nums (c_ok), '0 1 1 0 0');
  same ('Berger symbol of x = 0011', berger (3, 4), 2);
  same ('Berger symbol of y = 0001', berger (1, 4), 3);
  same ('Berger symbol of the faulty carries c_1..c_4', 4 - sum (c(2:5)), 4);
  pred = berger (3, 4) + berger (1, 4) - (4 - sum (c(2:5))) - c(1) + c(5);
  same ('predicted Berger symbol, 2 + 3 - 4 - 0 + 0', pred, 1);
  same ('Berger symbol of the sum 0010', berger (s, 4), 3);
  same ('Berger miss, 3 - 1', berger (s, 4) - pred, 2);
  same ('parity prediction, shared faulty carries c_0..c_3', bitxor (bitxor (par (3), par (1)), mod (sum (c(1:4)), 2)), 1);
  same ('parity of the sum 0010', par (s), 1);
  same ('parity of the healthy carries c_0..c_3', mod (sum (c_ok(1:4)), 2), 0);
  same ('parity prediction, duplicated healthy carries', bitxor (bitxor (par (3), par (1)), mod (sum (c_ok(1:4)), 2)), 1);
  same ('carries that differ between the chains, c_0..c_4', carry_names (find (c ~= c_ok) - 1), 'c_1 c_2');
  same ('mod-3 prediction (3 mod 3 + 1 mod 3) mod 3', mod (mod (3, 3) + mod (1, 3), 3), 1);
  same ('faulty result 0010 with c_4 = 0, mod 3', mod (result_value (s, c, 4), 3), 2);

  same ('arithmetic weight of 255', arith_weight (255), 2);
  same ('ones in 255', popcount (255), 8);
  same ('256 - 1', 256 - 1, 255);
  same ('arithmetic weight of 4080', arith_weight (4080), 2);
  same ('ones in 4080', popcount (4080), 8);
  same ('4096 - 16', 4096 - 16, 4080);
  same ('4080 in binary', bits_text (4080, 12), '111111110000');
  same ('arithmetic weight of 4, a run of one 1', arith_weight (4), 1);
  ok = true;
  for run = 2:16
    for shift = 0:16 - run
      ok = ok && arith_weight ((2 ^ run - 1) * 2 ^ shift) == 2;
    end
  end
  require ('a run of 2 to 16 ones, anywhere in 16 bits, has arithmetic weight 2', ok);

  same ('powers 2^0..2^3 modulo 3', join_nums (powers_mod (2, 4, 3)), '1 2 1 2');
  same ('order of 2 modulo 3', order_of_two (3), 2);
  same ('longest correctable length for A = 3', order_limit (3), 1);
  same ('+2 mod 3', mod (2, 3), 2);
  same ('-1 mod 3', mod (-1, 3), 2);
  same ('+1 mod 3', mod (1, 3), 1);
  same ('-2 mod 3', mod (-2, 3), 1);
  ok = true;
  for m = 2:16
    e = [2 .^ (0:m - 1), -(2 .^ (0:m - 1))];
    ok = ok && sum (mod (e, 3) == 1) >= 2 && sum (mod (e, 3) == 2) >= 2;
  end
  require ('A = 3, m = 2 to 16: each nonzero remainder names two or more errors', ok);

  big = 47 * 65535;
  same ('47 times 65535', big, 3080145);
  same ('bit length of 3080145', bit_length (big), 22);
  require ('2^21 <= 3080145 < 2^22', 2 ^ 21 <= big && big < 2 ^ 22);
  same ('order of 2 modulo 47', order_of_two (47), 23);
  same ('longest correctable length for A = 47', order_limit (47), 23);
  require ('A = 47 corrects every single error in 22-bit words', correcting (47, 22));
  same ('remainders of +2^i, -2^i for i < 22', 2 * 22, 44);
  same ('bound 2m + 1 for m = 22', 2 * 22 + 1, 45);

  same ('3-bit chunks of 181, lowest first', join_nums (chunks (181, 3, 9)), '5 6 2');
  same ('sum of those chunks', sum (chunks (181, 3, 9)), 13);
  same ('3-bit chunks of 13, lowest first', join_nums (chunks (13, 3, 6)), '5 1');
  same ('fold of 181 for A = 7', fold (181, 3), 6);
  same ('7 times 25 plus 6', 7 * 25 + 6, 181);

  x = 32760;
  y = 8;
  [s, c] = ripple_add (x, y, 16, 0, 4, 0);
  [s_ok, c_ok] = ripple_add (x, y, 16, 0, 0, 0);
  same ('16-bit adder, c_4 stuck at 0: faulty sum', s, 32752);
  same ('carries that differ between the chains, c_0..c_16', sum (c ~= c_ok), 12);
  require ('healthy carries c_4..c_15 are all 1', all (c_ok(5:16) == 1));
  require ('faulty carries c_4..c_15 are all 0', all (c(5:16) == 0));
  same ('parity prediction from the healthy chain', bitxor (bitxor (par (x), par (y)), mod (sum (c_ok(1:16)), 2)), 1);
  same ('parity of the correct sum 32768', par (s_ok), 1);
  same ('ones in the faulty sum 32752', popcount (s), 11);
  same ('parity of the faulty sum', par (s), 1);

  cm = [5 1 4; 12 4 3; 14 6 4];
  same ('row sums of C', join_nums (sum (cm, 2)'), '10 19 24');
  same ('row 2 sum minus the corrupted checksum 21', sum (cm(2, :)) - 21, -2);
  same ('column sum minus checksum, columns 1 to 3', join_nums (sum (cm, 1) - [31 11 11]), '0 0 0');
end


check_motivation ();
check_an_codes ();
check_residue_codes ();
check_adder_faults ();
check_predictions ();
check_two_rail ();
check_parity_checker ();
check_abft ();
check_examples ();
printf ('\nAll chapter 22 numbers reproduced.\n');
