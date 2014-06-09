function sticks = FitSticks(segms, pars, verbose)

% Fit straight line segments sticks(:,p)
% to segmentations segms(:,:,p)
%
% Input:
% - segms(:,:,p) = p-th binary segmentation mask
%
% Output:
% - sticks(:,p) = [x1 y1 x2 y2]' = p-th stick
%   in integer coordinates !
%

% parse arguments
if nargin < 3
  verbose = false;
end

% inits
Np = size(segms,3);
W = size(segms,2);
H = size(segms,1);
len = sqrt(W^2+H^2);                  % start with 'infinite' length
sticks = zeros(4,Np);

% get sticks
for p = 1:Np
  % compute moment-based properties
  [xs ys] = find(segms(:,:,p));
  [ctr orient] = MomentBasedProperties([ys'; xs']);
  %
  % determine stick
  x1 = ctr(1) + len*cos(orient);
  y1 = ctr(2) + len*sin(orient);
  x2 = ctr(1) - len*cos(orient);
  y2 = ctr(2) - len*sin(orient);
  %
  % compute true length
  step = 1/len;              % step is 1 pixel
  for t = 0:step:1
    dx = x2-x1; dy = y2-y1;
    cx = round(x1 + dx*t); cy = round(y1 + dy*t);
    if cx > 0 && cx <= W && cy > 0 && cy <= H
      if segms(cy,cx,p)
        break;
      end
    end 
  end
  x1n = cx; y1n = cy;
  step = 1/len;           
  for t = 1:-step:0
    dx = x2-x1; dy = y2-y1;
    cx = round(x1 + dx*t); cy = round(y1 + dy*t);
    if cx > 0 && cx <= W && cy > 0 && cy <= H
      if segms(cy,cx,p)
        break;
      end
    end 
  end
  x2 = cx; y2 = cy;
  x1 = x1n; y1 = y1n;
  %
  % re-determine stick
  COL = XY2COL([x1 y1 x2 y2]');
  ctr = COL(1:2);
  orient = COL(3);
  len = COL(4);
  len = len * pars.rescale_sticks;                % rescale lengths, as sometimes too long/short
  x1 = ctr(1) + 0.5*len*cos(orient);
  y1 = ctr(2) + 0.5*len*sin(orient);
  x2 = ctr(1) - 0.5*len*cos(orient);
  y2 = ctr(2) - 0.5*len*sin(orient);
  % 
  % store stick
  sticks(:,p) = round([x1 y1 x2 y2]');
  %
  % display stick
  if verbose > 1
    figure; imshow(segms(:,:,p)); axis equal; hold on;
    plot(sticks([1 3],p), sticks([2 4],p), 'y');
    title(['stick for part ' int2str(p)]);
    keyboard;
  end
end
