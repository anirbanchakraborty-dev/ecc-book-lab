function u = half_unit (stated)
  % HALF_UNIT  Half a unit in the last digit a stated number carries.
  %
  %   half_unit ('0.0280') is 0.00005 and half_unit ('1.5e-07') is 5e-09,
  %   the tolerance every check() line allows, as in the Python scripts.
  %   Commas and a trailing percent sign are ignored, so '12,345' and
  %   '9.6%' count their digits the way the chapter prints them.
  text = lower (strrep (stated, ',', ''));
  if (~isempty (text) && text(end) == '%')
    text = text(1:end - 1);
  end
  [mant, rest] = strtok (text, 'e');
  ex = 0;
  if (~isempty (rest))
    ex = str2double (rest(2:end));
  end
  dot = find (mant == '.', 1);
  places = 0;
  if (~isempty (dot))
    places = numel (mant) - dot;
  end
  u = 0.5 * 10 ^ (ex - places);
end
