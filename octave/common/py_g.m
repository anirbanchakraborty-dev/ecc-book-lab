function s = py_g (x, precision)
  % PY_G  Python's format(x, '.Ng'): N significant digits, trailing zeros dropped.
  %
  %   py_g (x) is py_g (x, 9), which is what every check() line prints.
  %   Octave's sprintf ('%.9g', x) gives the same text for every finite
  %   x. The one difference is the spelling of the non-finite values,
  %   which Python writes inf, -inf and nan and Octave Inf, -Inf and NaN.
  if (nargin < 2)
    precision = 9;
  end
  x = double (x);
  if (isnan (x))
    s = 'nan';
  elseif (isinf (x))
    if (x > 0)
      s = 'inf';
    else
      s = '-inf';
    end
  else
    s = sprintf ('%.*g', precision, x);
  end
end
