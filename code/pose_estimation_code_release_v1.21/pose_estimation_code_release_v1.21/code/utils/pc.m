function Hout = pc(X,Y,F)
if nargin == 1
  F = X;
  if any(size(F) == 1)
    F = squeeze(F);
  end
  [m,n,colours]=size(F);
  X = 1:n;
  Y = 1:m;
end
if any(size(F) == 1)
  F = squeeze(F);
end
if isa(F, 'uint8')
  F = double(F);
end
map = colormap;
OK = finite(F);
if all(OK)
  maxval = max(max(max(F)));
  minval = min(min(min(F)));
else
  I = find(OK);
  n_not_ok = length(find(~OK));
  if n_not_ok > 0
    % fprintf(2, 'pc: See %d nans\n', n_not_ok);
  end
  if isempty(I)
    % all inf
    maxval = 1;
    minval = 0;
    F = zeros(size(F));
  else
    maxval = max(F(I));
    minval = min(F(I));
    I = find(~OK);
    F(I) = minval +0*I;
  end
end
% disp(sprintf('Range = [%g %g]\n',minval,maxval));
if size(F,3) == 3
  % RGB
  scale = 1/(maxval - minval);
  F = (F - minval)*scale;
  H = image(X,Y,F);
else
  % Gray or indexed -- use cmap anyway
  if (maxval - minval) == 0
    fprintf('pc: constant image (value %d)\n', minval);
    F = F - minval;
  else
    scale = length(map)/(maxval - minval);
    F = (F - minval)*scale;
  end
  H = image(X,Y,F);
end

% Added aug02 to make more like imshow
axis image
axis off

if nargout > 0
  Hout = H;
end
