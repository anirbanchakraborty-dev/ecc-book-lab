classdef PyRandom < handle
  % PYRANDOM  The random-number generator of Python's random.Random, in Octave.
  %
  %   rng = PyRandom (seed) builds the generator that Python builds with
  %   random.Random(seed) for a nonnegative integer seed, so the draws below
  %   come out bit for bit the same as Python's on any machine:
  %
  %     rng.random (n)               n draws of Python's rng.random(), a column
  %     rng.getrandbits (k, n)       n draws of rng.getrandbits(k), for k up to 53
  %     rng.randbelow (m)            one draw of Python's internal _randbelow(m)
  %     rng.choice (seq)             one draw of rng.choice(seq)
  %     rng.choices (pop, k, w)      rng.choices(pop, weights=w, k=k); w may be []
  %     rng.gauss (mu, sigma, n)     n draws of rng.gauss(mu, sigma)
  %
  %   Both languages run the Mersenne Twister MT19937 of Matsumoto and
  %   Nishimura (1998), seeded through its init_by_array routine with the
  %   seed split into 32-bit words, least significant first, as CPython's
  %   Modules/_randommodule.c does. Octave's own rand() is a Mersenne
  %   Twister too, and with Octave 11.3.0 rand ('state', s) was measured to
  %   give the same draws as Python's rng.random() for every seed tried
  %   below 2^32, and for the seed 2^33 + 1 passed as its 32-bit words [1 2].
  %   Octave's documentation does not promise that, and rand has no
  %   counterpart of getrandbits, choice or gauss that follows Python's
  %   rules, so the twins use this class instead. Every 32-bit word is held in a double, which
  %   stores any integer below 2^53 exactly, and every product that could
  %   pass 2^53 is split into 16-bit halves first.

  properties (Access = private)
    mt = [];          % the 624 words of the generator's state
    out = [];         % the 624 tempered outputs of the current state
    pos = 624;        % how many of those outputs have been handed out
    gauss_next = [];  % the second normal of the last Box-Muller pair, if unused
  end

  methods
    function obj = PyRandom (seed)
      if (seed < 0 || seed ~= floor (seed))
        error ('PyRandom: the seed must be a nonnegative integer');
      end
      key = [];
      s = double (seed);
      if (s == 0)
        key = 0;
      end
      while (s > 0)
        key(end + 1) = mod (s, 2^32);
        s = floor (s / 2^32);
      end
      obj.mt = PyRandom.init_by_array (key);
      obj.pos = 624;
    end

    function w = words (obj, n)
      % The next n 32-bit outputs of the generator, as a column of doubles.
      w = zeros (n, 1);
      filled = 0;
      while (filled < n)
        if (obj.pos >= 624)
          obj.mt = PyRandom.twist (obj.mt);
          obj.out = PyRandom.temper (obj.mt);
          obj.pos = 0;
        end
        take = min (n - filled, 624 - obj.pos);
        w(filled + 1:filled + take) = obj.out(obj.pos + 1:obj.pos + take);
        filled = filled + take;
        obj.pos = obj.pos + take;
      end
    end

    function r = random (obj, n)
      % Python's random(): 53 random bits, a from one word and b from the next.
      if (nargin < 2)
        n = 1;
      end
      w = obj.words (2 * n);
      a = floor (w(1:2:end) / 32);
      b = floor (w(2:2:end) / 64);
      r = (a * 67108864 + b) / 9007199254740992;
    end

    function r = getrandbits (obj, k, n)
      % Python's getrandbits(k): whole words, least significant first, with
      % the last word shifted right so that exactly k bits remain.
      if (nargin < 3)
        n = 1;
      end
      if (k > 53)
        error ('PyRandom: getrandbits beyond 53 bits does not fit a double');
      end
      if (k == 0)
        r = zeros (n, 1);
        return;
      end
      nw = ceil (k / 32);
      w = reshape (obj.words (nw * n), nw, n);
      r = zeros (n, 1);
      left = k;
      for i = 1:nw
        part = w(i, :).';
        if (left < 32)
          part = floor (part / 2^(32 - left));
        end
        r = r + part * 2^(32 * (i - 1));
        left = left - 32;
      end
    end

    function r = randbelow (obj, m)
      % Python's _randbelow(m): draw bit_length(m) bits until the draw is below m.
      k = 0;
      while (2^k <= m)
        k = k + 1;
      end
      r = obj.getrandbits (k);
      while (r >= m)
        r = obj.getrandbits (k);
      end
    end

    function x = choice (obj, seq)
      % Python's choice(seq), for a vector, a char row or a cell row.
      i = obj.randbelow (numel (seq)) + 1;
      if (iscell (seq))
        x = seq{i};
      else
        x = seq(i);
      end
    end

    function x = choices (obj, pop, k, weights)
      % Python's choices(pop, weights=weights, k=k). Without weights each draw
      % is floor(random() * n); with them, random() * total is placed among
      % the running totals of the weights by bisect_right.
      n = numel (pop);
      if (nargin < 4 || isempty (weights))
        idx = floor (obj.random (k) * n) + 1;
      else
        cum = cumsum (double (weights(:).'));
        total = cum(end);
        r = obj.random (k) * total;
        idx = sum (r >= cum(1:n - 1), 2) + 1;
      end
      x = pop(idx);
      if (ischar (pop))
        x = reshape (x, 1, []);
      end
    end

    function z = gauss (obj, mu, sigma, n)
      % Python's gauss(mu, sigma): the Box-Muller pair cos and sin, the sin
      % kept for the next call.
      if (nargin < 4)
        n = 1;
      end
      z = zeros (n, 1);
      filled = 0;
      if (n > 0 && ~isempty (obj.gauss_next))
        z(1) = obj.gauss_next;
        obj.gauss_next = [];
        filled = 1;
      end
      need = n - filled;
      if (need > 0)
        pairs = ceil (need / 2);
        r = obj.random (2 * pairs);
        x2pi = r(1:2:end) * (2 * pi);
        g2rad = sqrt (-2.0 * log (1.0 - r(2:2:end)));
        both = [cos(x2pi) .* g2rad, sin(x2pi) .* g2rad].';
        both = both(:);
        z(filled + 1:n) = both(1:need);
        if (mod (need, 2) == 1)
          obj.gauss_next = both(end);
        end
      end
      z = mu + z * sigma;
    end
  end

  methods (Static, Access = private)
    function c = mulmod (a, b)
      % (a * b) mod 2^32 for words a and b, exactly, in doubles.
      ah = floor (a / 65536);
      al = a - ah * 65536;
      c = mod (mod (ah * b, 65536) * 65536 + al * b, 2^32);
    end

    function mt = init_genrand (s)
      mt = zeros (1, 624);
      mt(1) = s;
      for i = 2:624
        prev = mt(i - 1);
        mt(i) = mod (PyRandom.mulmod (bitxor (prev, floor (prev / 2^30)), 1812433253) + (i - 1), 2^32);
      end
    end

    function mt = init_by_array (key)
      mt = PyRandom.init_genrand (19650218);
      i = 2;
      j = 1;
      len = numel (key);
      for count = 1:max (624, len)
        prev = mt(i - 1);
        mixed = PyRandom.mulmod (bitxor (prev, floor (prev / 2^30)), 1664525);
        mt(i) = mod (bitxor (mt(i), mixed) + key(j) + (j - 1), 2^32);
        i = i + 1;
        j = j + 1;
        if (i > 624)
          mt(1) = mt(624);
          i = 2;
        end
        if (j > len)
          j = 1;
        end
      end
      for count = 1:623
        prev = mt(i - 1);
        mixed = PyRandom.mulmod (bitxor (prev, floor (prev / 2^30)), 1566083941);
        mt(i) = mod (bitxor (mt(i), mixed) - (i - 1), 2^32);
        i = i + 1;
        if (i > 624)
          mt(1) = mt(624);
          i = 2;
        end
      end
      mt(1) = 2^31;
    end

    function mt = twist (mt)
      % One pass of MT19937's recurrence over all 624 words. Word k is
      % rebuilt from words k, k+1 and k+397, counted around the ring, and
      % the pass runs in four slices so that each slice reads only words a
      % slice before it has already rebuilt.
      a = 1:227;
      mt(a) = PyRandom.step (mt, a, a + 1, a + 397);
      b = 228:454;
      mt(b) = PyRandom.step (mt, b, b + 1, b - 227);
      c = 455:623;
      mt(c) = PyRandom.step (mt, c, c + 1, c - 227);
      mt(624) = PyRandom.step (mt, 624, 1, 397);
    end

    function v = step (mt, dst, nxt, far)
      % The top bit of word dst joined to the low 31 bits of word nxt,
      % shifted right once, with 0x9908b0df mixed in when the dropped bit
      % was 1, and all of it XORed into word far.
      y = bitand (mt(dst), 2^31) + bitand (mt(nxt), 2^31 - 1);
      v = bitxor (bitxor (mt(far), floor (y / 2)), 2567483615 * mod (y, 2));
    end

    function y = temper (mt)
      y = mt(:);
      y = bitxor (y, floor (y / 2^11));
      y = bitxor (y, bitand (mod (y * 2^7, 2^32), 2636928640));   % 0x9d2c5680
      y = bitxor (y, bitand (mod (y * 2^15, 2^32), 4022730752));  % 0xefc60000
      y = bitxor (y, floor (y / 2^18));
    end
  end
end
