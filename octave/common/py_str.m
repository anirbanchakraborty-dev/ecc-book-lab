function s = py_str (x, quoted)
  % PY_STR  Format a value the way Python's str() prints it.
  %
  %   py_str (x) returns the text that Python's str(x) would print for the
  %   value x, so that an Octave script can print the same check lines as
  %   its Python twin. Octave has one number type where Python has two, so
  %   the caller marks which one is meant:
  %
  %     an integer class (int64 (5), uint8 (3), ...)  prints as a Python int: 5
  %     a double (5, 0.1, Inf)                        prints as a Python float: 5.0, 0.1, inf
  %     a logical (true, false)                       prints as True or False
  %     a char row ('abc')                            prints as abc, or 'abc' inside a list
  %     a vector or a cell row                        prints as a list: [1, 0, 1]
  %     struct ('tuple', {{a, b}})                    prints as a tuple: (a, b)
  %
  %   A Python integer too large for int64 is passed as its decimal text.
  %   The second argument, used for the items of a list or tuple, asks for
  %   repr() rather than str(), which only changes how text is quoted.

  if (nargin < 2)
    quoted = false;
  end
  if (isstruct (x) && isfield (x, 'tuple'))
    items = cellfun (@(v) py_str (v, true), x.tuple, 'UniformOutput', false);
    if (numel (items) == 1)
      s = ['(' items{1} ',)'];
    else
      s = ['(' strjoin(items, ', ') ')'];
    end
  elseif (iscell (x))
    items = cellfun (@(v) py_str (v, true), x, 'UniformOutput', false);
    s = ['[' strjoin(items, ', ') ']'];
  elseif (ischar (x))
    if (quoted)
      s = ['''' x ''''];
    else
      s = x;
    end
  elseif (numel (x) ~= 1)
    s = py_str (num2cell (x(:).'), quoted);
  elseif (islogical (x))
    if (x)
      s = 'True';
    else
      s = 'False';
    end
  elseif (isinteger (x))
    s = sprintf ('%d', x);
  else
    s = float_repr (double (x));
  end
end

function s = float_repr (x)
  % Python's repr() of a float: the shortest text that reads back as x,
  % in fixed notation when the decimal exponent lies in [-4, 16) and in
  % exponent notation otherwise.
  if (isnan (x))
    s = 'nan';
    return;
  elseif (isinf (x))
    if (x > 0)
      s = 'inf';
    else
      s = '-inf';
    end
    return;
  elseif (x == 0)
    if (1 / x < 0)
      s = '-0.0';
    else
      s = '0.0';
    end
    return;
  end
  for p = 1:17
    e_text = sprintf ('%.*e', p - 1, x);
    if (str2double (e_text) == x)
      break;
    end
  end
  [mant, ex] = strtok (e_text, 'e');
  ex = str2double (ex(2:end));
  neg = mant(1) == '-';
  digits = strrep (strrep (mant, '-', ''), '.', '');
  digits = regexprep (digits, '0+$', '');
  if (isempty (digits))
    digits = '0';
  end
  if (ex >= -4 && ex < 16)
    if (ex >= 0)
      if (numel (digits) <= ex + 1)
        int_part = [digits repmat('0', 1, ex + 1 - numel (digits))];
        frac = '0';
      else
        int_part = digits(1:ex + 1);
        frac = digits(ex + 2:end);
      end
    else
      int_part = '0';
      frac = [repmat('0', 1, -ex - 1) digits];
    end
    s = [int_part '.' frac];
  else
    if (numel (digits) > 1)
      body = [digits(1) '.' digits(2:end)];
    else
      body = digits;
    end
    if (ex < 0)
      sign_text = '-';
    else
      sign_text = '+';
    end
    s = sprintf ('%se%s%02d', body, sign_text, abs (ex));
  end
  if (neg)
    s = ['-' s];
  end
end
