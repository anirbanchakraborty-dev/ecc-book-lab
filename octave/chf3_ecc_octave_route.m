% Reproduce every number in front-matter chapter f3 of Error Correction Coding, "The Octave Route".
%
% The Octave twin of python/chf3_ecc_octave_route.py. Needs GNU Octave 7 or
% later, for jsonencode and jsondecode, and no package. Run it from the
% repository root:
%
%     octave --no-window-system --norc --quiet octave/chf3_ecc_octave_route.m
%     octave --no-window-system --norc --quiet octave/chf3_ecc_octave_route.m --handoff
%
% The Python script computes the Octave side of each pair in the chapter
% with a small Python model of Octave's rule. This script computes it with
% Octave itself, and the Python side with Octave code that follows Python's
% rule, so the two print the same lines when every model is faithful.
% tools/compare-routes.sh runs both and compares what they print.
%
% Every value computed here is printed beside the value the chapter states.
% The script stops with exit status 1 at the first pair that disagrees by
% more than half a unit in the last digit the chapter prints, and at the
% first count, list or string that differs at all.
%
% With --handoff the script also runs Python as a separate program, the way
% a reader moving between the two routes would. It writes the chapter's
% three data words and 100,000 draws of rand after rand ("state", 1) to
% files, has Python compute the words' parity bits and compare the draws
% with random.Random(1), and reads Python's answer back from JSON. It also
% hands Python 3,146 quotients of Octave's idivide, in its default rounding
% and with "fix", and of fix(a ./ b) on doubles, which Python checks against
% exact truncation toward zero; those counts were measured with Octave 11.3.0
% and are reported as a note on any other version. The
% PYTHON environment variable names the interpreter, python3 by default.
% That mode exits 1 when Python cannot be run.

1;
addpath (fullfile (fileparts (mfilename ('fullpath')), 'common'));

function check (label, computed, stated)
  value = str2double (regexprep (strrep (stated, ',', ''), '%+$', ''));
  computed = double (computed);
  if (stated(end) == '%')
    computed = 100 * computed;
  end
  printf ('  %-64s %-16s %s\n', label, py_g (computed), stated);
  if (abs (computed - value) > half_unit (stated) * (1 + 1e-9))
    printf ('MISMATCH in ''%s'': computed %.17g, the chapter states %s\n', label, computed, stated);
    exit (1);
  end
end

function s = text_of (x)
  % The text Python's str() gives for the values this script compares.
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

function same (label, computed, stated)
  a = text_of (computed);
  b = text_of (stated);
  printf ('  %-64s %-16s %s\n', label, a, b);
  if (! strcmp (a, b))
    printf ('MISMATCH in ''%s'': computed %s, the chapter states %s\n', label, a, b);
    exit (1);
  end
end

function s = spaced (v)
  % A row of whole numbers, one space apart, as Python's " ".join prints it.
  s = strjoin (arrayfun (@(x) sprintf ('%d', x), v(:)', 'UniformOutput', false), ' ');
end

function check_indexing ()
  disp ('Indexing and ranges');
  v = [5 6 7];
  same ('first element of [5 6 7]: v[0] in Python, v(1) in Octave', v(1), 5);
  same ('last element: v[-1] in Python, v(end) in Octave', v(end), 7);
  same ('elements in the Octave range 0:6', numel (0:6), 7);
  same ('elements in Python''s range(7)', numel (0:6), 7);
  same ('elements in the Octave range 1:2:9', numel (1:2:9), 5);
  same ('elements in the Octave range 1:2:10', numel (1:2:10), 5);
  same ('bitget(100, 1:7), the bits of 100 from the lowest', spaced (bitget (100, 1:7)), '0 0 1 0 0 1 1');
  same ('bitget(100, 3), Python (100 >> 2) & 1', double (bitget (100, 3)), 1);
end

function check_doubles ()
  disp ('Doubles');
  check ('eps, the gap from 1 to the next double', eps, '2.22e-16');
  same ('eps equals 2^-52', eps == 2^-52, true);
  same ('flintmax, the integer 2^53', flintmax, 9007199254740992);
  same ('2^53 + 1 stored as a double', 2^53 + 1, 9007199254740992);
  same ('0.1 + 0.2 == 0.3', 0.1 + 0.2 == 0.3, false);
  same ('0.1 + 0.2 written with 17 significant digits', sprintf ('%.17g', 0.1 + 0.2), '0.30000000000000004');
  same ('0.1 + 0.2 written with 16 significant digits', sprintf ('%.16g', 0.1 + 0.2), '0.3');
end

function check_integer_types ()
  disp ('Integer types saturate');
  same ('uint8(250) + 10', uint8 (250) + 10, 255);
  same ('uint8(5) - 10', uint8 (5) - 10, 0);
  same ('int8(-100) - 100', int8 (-100) - 100, -128);
  same ('int8(127) + 1', int8 (127) + 1, 127);
  same ('250 + 10 wrapped modulo 256, as C''s uint8_t would', mod (250 + 10, 256), 4);
  same ('uint8(300)', uint8 (300), 255);
  same ('uint8(2.5), a tie rounds away from zero', uint8 (2.5), 3);
  same ('int8(-2.5)', int8 (-2.5), -3);
  same ('uint8(1) + 1.5, whose result is still a uint8', uint8 (1) + 1.5, 3);
  same ('uint8(200) * 2', uint8 (200) * 2, 255);
  same ('(uint8(250) + 10) - 10, saturating at each step', (uint8 (250) + 10) - 10, 245);
  same ('uint8(250) + (10 - 10)', uint8 (250) + (10 - 10), 250);
end

function check_integer_division ()
  disp ('Integer division');
  same ('uint8(7) / 2, where 3.5 rounds to 4', uint8 (7) / 2, 4);
  same ('uint8(5) / 2', uint8 (5) / 2, 3);
  same ('int32(-7) / int32(2)', int32 (-7) / int32 (2), -4);
  same ('Python 7 // 2, floor(7 / 2) in Octave', floor (7 / 2), 3);
  same ('Python -7 // 2, floor(-7 / 2) in Octave', floor (-7 / 2), -4);
  same ('idivide(int32(-7), int32(2), ''floor'')', idivide (int32 (-7), int32 (2), 'floor'), -4);
  same ('int(-7 / 2), fix(-7 / 2) in Octave', fix (-7 / 2), -3);
end

function check_rounding ()
  disp ('Rounding');
  same ('Octave round(2.5), a tie away from zero', round (2.5), 3);
  same ('Octave round(-2.5)', round (-2.5), -3);
  same ('Octave round(0.5)', round (0.5), 1);
  same ('Python round(2.5), roundb(2.5) in Octave, a tie to even', roundb (2.5), 2);
  same ('Python round(-2.5), roundb(-2.5)', roundb (-2.5), -2);
  same ('Python round(0.5), roundb(0.5)', roundb (0.5), 0);
  same ('Python round(-0.5), roundb(-0.5)', roundb (-0.5), 0);
  same ('fix(-3.5) in Octave, math.trunc(-3.5) in Python', fix (-3.5), -3);
  same ('floor(-3.5) in Octave, math.floor(-3.5) in Python', floor (-3.5), -4);
end

function check_remainders ()
  disp ('Remainders');
  same ('mod(-7, 3) in Octave, -7 % 3 in Python', mod (-7, 3), 2);
  same ('rem(-7, 3) in Octave, math.fmod(-7, 3) in Python', rem (-7, 3), -1);
  same ('mod(7, -3), which takes the divisor''s sign', mod (7, -3), -2);
  same ('rem(7, -3), which takes the dividend''s sign', rem (7, -3), 1);
  [a, b] = ndgrid (-20:20, [-7:-1, 1:7]);
  a = a(:);
  b = b(:);
  same ('pairs (a, b) tried, a from -20 to 20, b from -7 to 7, b != 0', numel (a), 574);
  same ('pairs with a = b*floor(a/b) + mod(a, b)', sum (a == b .* floor (a ./ b) + mod (a, b)), 574);
  same ('pairs with a = b*fix(a/b) + rem(a, b)', sum (a == b .* fix (a ./ b) + rem (a, b)), 574);
  same ('pairs where mod(a, b) and rem(a, b) differ', sum (mod (a, b) != rem (a, b)), 180);
  same ('mod(-6:6, 3), Python [a % 3 for a in range(-6, 7)]', spaced (mod (-6:6, 3)), '0 1 2 0 1 2 0 1 2 0 1 2 0');
  same ('rem(-6:6, 3), math.fmod(a, 3) for the same a', spaced (rem (-6:6, 3)), '0 -2 -1 0 -2 -1 0 1 2 0 1 2 0');
end

function check_bits ()
  disp ('Bit operations');
  same ('bitxor(100, 29), 100 ^ 29 in Python', bitxor (100, 29), 121);
  same ('bitand(100, 29), 100 & 29', bitand (100, 29), 4);
  same ('bitor(100, 29), 100 | 29', bitor (100, 29), 125);
  same ('bitshift(100, -2), 100 >> 2', bitshift (100, -2), 25);
  same ('bitshift(100, 3), 100 << 3', bitshift (100, 3), 800);
  same ('dec2bin(100), format(100, ''b'')', dec2bin (100), '1100100');
  same ('dec2bin(100, 8), format(100, ''08b'')', dec2bin (100, 8), '01100100');
  same ('bin2dec(''1100100''), int(''1100100'', 2)', bin2dec ('1100100'), 100);
  same ('ones in 100, (100).bit_count()', sum (dec2bin (100) == '1'), 3);
end

function check_matrices ()
  disp ('Vectors and matrices');
  words = [100 29 121];
  same ('[1 2 3] .* [4 5 6], element by element', spaced ([1 2 3] .* [4 5 6]), '4 10 18');
  same ('[1 2 3] * [4; 5; 6], a row times a column', [1 2 3] * [4; 5; 6], 32);
  same ('mod([1 1] * [1 0 1; 0 1 1], 2), a codeword', spaced (mod ([1 1] * [1 0 1; 0 1 1], 2)), '1 1 0');
  ones_in = sum (dec2bin (words) == '1', 2);
  same ('ones in each of the words [100 29 121]', spaced (ones_in), '3 4 5');
  same ('even-parity bit of each word', spaced (mod (ones_in, 2)), '1 0 1');
end

function check_random ()
  disp ('Random numbers');
  rand ('state', 1);
  draws = rand (100000, 1);
  check ('first draw of random.Random(1), rand after rand(''state'', 1)', draws(1), '0.13436424411240122');
  check ('draw 100,000 of the same stream', draws(end), '0.6355707617477253');
  seed = 2^33 + 1;
  words = [mod(seed, 2^32), floor(seed / 2^32)];
  same ('32-bit words of the seed 2^33 + 1, lowest first', spaced (words), '1 2');
  rand ('state', words);
  check ('first draw of random.Random(2**33 + 1), rand(''state'', [1 2])', rand (), '0.5845141492495678');
end

function check_json ()
  disp ('JSON between the routes');
  text = jsonencode ([100 29 121]);
  same ('JSON text of [100 29 121], written compactly', text, '[100,29,121]');
  w = jsondecode (text);
  same ('sum of the decoded array', sum (w), 250);
  same ('rows and columns of the array as Octave''s jsondecode returns it', spaced (size (w)), '3 1');
  same ('2^53 + 1 read from JSON into a double', jsondecode ('9007199254740993'), 9007199254740992);
  same ('2^53 + 1 sent as a JSON string, read back as text', jsondecode ('"9007199254740993"'), '9007199254740993');
  same ('JSON true, a logical 1 in Octave', double (jsondecode ('true')), 1);
end

function check_examples_and_exercises ()
  disp ('Worked examples and exercises');
  words = [100 29 121];
  disp ('Example: a byte that saturates');
  same ('uint8(200) + uint8(100)', uint8 (200) + uint8 (100), 255);
  same ('the true sum', 200 + 100, 300);
  same ('the same sum wrapped modulo 256', mod (200 + 100, 256), 44);
  disp ('Example: porting a floor division');
  same ('Python -17 // 5, floor(-17 / 5) in Octave', floor (-17 / 5), -4);
  same ('Python -17 % 5, mod(-17, 5) in Octave', mod (-17, 5), 3);
  same ('rem(-17, 5), math.fmod(-17, 5)', rem (-17, 5), -2);
  same ('int32(-17) / int32(5), where -3.4 rounds to -3', int32 (-17) / int32 (5), -3);
  disp ('Example: the parity of three bytes');
  same ('dec2bin([100 29 121]) rows, ones in each', spaced (sum (dec2bin (words) == '1', 2)), '3 4 5');
  same ('bitxor of the three words', bitxor (bitxor (words(1), words(2)), words(3)), 0);
  disp ('Exercises');
  same ('int8(-128) - 1', int8 (-128) - 1, -128);
  same ('int8(-128) * -1, whose true value 128 does not fit', int8 (-128) * -1, 127);
  same ('uint8(9) / 4, where 2.25 rounds to 2', uint8 (9) / 4, 2);
  same ('uint8(10) / 4, where 2.5 rounds to 3', uint8 (10) / 4, 3);
  same ('mod(-1, 8), Python -1 % 8', mod (-1, 8), 7);
  same ('rem(-1, 8)', rem (-1, 8), -1);
  same ('numel(10:-3:0), the range 10, 7, 4, 1', numel (10:-3:0), 4);
  same ('round(-0.5) in Octave', round (-0.5), -1);
  same ('bitshift(1, 10), 1 << 10', bitshift (1, 10), 1024);
  same ('dec2bin(5, 4), format(5, ''04b'')', dec2bin (5, 4), '0101');
end

function remove_folder (folder)
  % Delete the hand-off folder and everything in it, without asking.
  confirm_recursive_rmdir (false, 'local');
  rmdir (folder, 's');
end

function status = handoff ()
  printf ('\nHanding work to Python\n');
  python = getenv ('PYTHON');
  if (isempty (python))
    python = 'python3';
  end
  tmp = tempname ();
  mkdir (tmp);
  cleanup = onCleanup (@() remove_folder (tmp));
  fid = fopen (fullfile (tmp, 'words.json'), 'w');
  fputs (fid, jsonencode ([100 29 121]));
  fclose (fid);
  rand ('state', 1);
  fid = fopen (fullfile (tmp, 'draws.txt'), 'w');
  fprintf (fid, '%.17g\n', rand (100000, 1));
  fclose (fid);
  [A, B] = ndgrid (-60:60, [-13:-1, 1:13]);
  quotients = struct ('a', A(:)', 'b', B(:)', 'q', double (idivide (int32 (A(:)'), int32 (B(:)'))), ...
                      'qfix', double (idivide (int32 (A(:)'), int32 (B(:)'), 'fix')), 'qtrunc', fix (A(:)' ./ B(:)'));
  fid = fopen (fullfile (tmp, 'idivide.json'), 'w');
  fputs (fid, jsonencode (quotients));
  fclose (fid);
  fid = fopen (fullfile (tmp, 'handoff.py'), 'w');
  fputs (fid, [ ...
    "import json, random, sys\n" ...
    "words = json.load(open(sys.argv[1]))\n" ...
    "ref = random.Random(1)\n" ...
    "draws = [float(t) for t in open(sys.argv[2]).read().split()]\n" ...
    "result = {'version': sys.version.split()[0],\n" ...
    "          'parity': [w.bit_count() % 2 for w in words],\n" ...
    "          'equal': sum(x == ref.random() for x in draws)}\n" ...
    "d = json.load(open(sys.argv[4]))\n" ...
    "trunc = [int(abs(a) // abs(b)) * (1 if (a < 0) == (b < 0) else -1) for a, b in zip(d['a'], d['b'])]\n" ...
    "result['wrong'] = [int(a) for a, q, t in zip(d['a'], d['q'], trunc) if q != t]\n" ...
    "result['wrong_fix'] = sum(q != t for q, t in zip(d['qfix'], trunc))\n" ...
    "result['wrong_trunc'] = sum(q != t for q, t in zip(d['qtrunc'], trunc))\n" ...
    "json.dump(result, open(sys.argv[3], 'w'))\n"]);
  fclose (fid);
  cmd = sprintf ('"%s" "%s" "%s" "%s" "%s" "%s"', python, fullfile (tmp, 'handoff.py'), ...
                 fullfile (tmp, 'words.json'), fullfile (tmp, 'draws.txt'), fullfile (tmp, 'result.json'), ...
                 fullfile (tmp, 'idivide.json'));
  [code, output] = system (cmd);
  if (code != 0)
    printf ('  %s exited with status %d:\n  %s\n', python, code, strtrim (output));
    status = 1;
    return;
  end
  result = jsondecode (fileread (fullfile (tmp, 'result.json')));
  printf ('  Python %s, run as %s\n', result.version, python);
  same ('parity bits Python computed for [100 29 121]', spaced (result.parity), '1 0 1');
  same ('draws of rand after rand(''state'', 1) equal to random.Random(1)', result.equal, 100000);
  wrong = result.wrong(:);
  counts = {'idivide default against truncation, wrong of 3146 pairs', numel(wrong), 1184; ...
            'of those wrong, how many have a nonnegative dividend', sum(wrong >= 0), 0; ...
            'idivide ''fix'' against truncation, wrong of 3146 pairs', result.wrong_fix, 1184; ...
            'fix(a ./ b) on doubles against truncation, wrong of 3146 pairs', result.wrong_trunc, 0; ...
            'idivide(int8(-3), int8(4), ''fix''), documented as 0', ...
            double(idivide (int8 (-3), int8 (4), 'fix')), -1};
  for i = 1:rows (counts)
    if (strcmp (version (), '11.3.0') || counts{i, 2} == counts{i, 3})
      same (counts{i, :});
    else
      printf ('  %-64s %-16d %d\n', counts{i, :});
      printf ('  NOTE: the chapter''s value was measured with Octave 11.3.0, and Octave %s differs here.\n', version ());
    end
  end
  disp ('Hand-off checks pass.');
  status = 0;
end

check_indexing ();
check_doubles ();
check_integer_types ();
check_integer_division ();
check_rounding ();
check_remainders ();
check_bits ();
check_matrices ();
check_random ();
check_json ();
check_examples_and_exercises ();
printf ('\nAll chapter f3 numbers reproduced.\n');
if (any (strcmp (argv (), '--handoff')))
  exit (handoff ());
end
