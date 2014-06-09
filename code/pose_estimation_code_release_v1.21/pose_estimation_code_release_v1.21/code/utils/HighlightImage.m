function im = HighlightImage(im, mask, p, coeff)

% COMMENT TO BE WRITTEN
%
% if p>0 -> adds color plane;
% if p<0 -> removes color plane
%

if nargin < 4
  coeff = 128;
end

s = sign(p);
p = abs(p);

imp = im(:,:,p);
imp = int16(imp);
imp = uint8(max(min(imp + s*coeff*int16(mask), 255),0));
im(:,:,p) = imp;
