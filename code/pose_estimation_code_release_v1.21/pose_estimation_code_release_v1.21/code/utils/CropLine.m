function lo = CropLine(l, wh)

% Crops horizontal or vertical line l = [xmin ymin xmax ymax]
% onto image of width wh(1) and height wh(2)
%
% If l is entirely out of image
% -> lo = []
%

% horizontal or vertical ?
if l(2) == l(4)
  % horizontal
  xov = [max([l(1) 1]) min([l(3) wh(1)])];
  if l(2) < 1 || l(2) > wh(2) || xov(1) > xov(2)
    lo = [];
    return;
  end
  lo = [xov(1) l(2) xov(2) l(4)];
elseif l(1) == l(3)
  % vertical
  yov = [max([l(2) 1]) min([l(4) wh(2)])];
  if l(1) < 1 || l(1) > wh(1) || yov(1) > yov(2)
    lo = [];
    return;
  end
  lo = [l(1) yov(1) l(3) yov(2)];
else
  error([mfilename ' works only for horizontal or vertical lines']);
end
