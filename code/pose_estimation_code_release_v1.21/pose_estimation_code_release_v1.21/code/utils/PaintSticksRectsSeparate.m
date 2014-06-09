function [out bodyMinBBcropped limbMinBBcropped] = PaintSticksRectsSeparate(sticks, siz, aspect_ratios, classname, limbclasses, colmask)
% [out bodyMinBBcropped limbMinBBcropped] = PaintSticksRectsSeparate(sticks, siz, aspect_ratios, classname, limbclasses, colmask)
% modyfication of PaintSticksRects returns image withc channels for each limb separately + last for whole body
% paint sticks(:,p) on an image, using cols(p,:),
% and thickness so as to reach length/width ratios aspect_ratios
% (= automatic scale-adaption).
%
% Input:
% - sticks(:,p) = [x1 y1 x2 y2]' = integers
% - siz = [width height]
% - aspect_ratios(p) = length/width of sticks(:,p)
% - classname - optional if not empty ('ubf' 'ubp' 'full') then occlusion handling
% - limbclasses - optional bool false/true
% - colmask - must be given if limbclasses == true - matrix indicating limb membership to a limb class (limbClasses x limbTypes)
%
% Output:
% - out(:,:,:) (h w nlimbs+1) in [0,1]
% - bodyMinBBcropped - min rect around the stickman [minx miny maxx maxy]
% - limbMinBBcropped - min rects around each limb(or limb class) cropped to the size of the image limbMinBBcropped(:,N) = [minx miny maxx maxy]'

if nargin < 4
  classname = [];
end
if nargin < 5
  limbclasses = false;
end

% initialize minBB with width and hight for xmin ymin and 0s for xmax ymax
bodyMinBBcropped = [siz(1) siz(2) 0 0];

nlimbs = size(sticks,2);
out = false(siz(2),siz(1),nlimbs+1);
% limbareas = zeros(1,nlimbs);
  
limbMinBBcropped = zeros(4,nlimbs);
sticks_col = XY2COL(sticks);                 % convert to sticks_col = [ctr_x ctr_y orient length]
for p = 1:nlimbs
  %
  % prepare data
  s = sticks(:,p);                           % current stick in XY format
  if me_isEmptyStick(s); % if the current stick is not present (body part is occluded then do not draw the mask => everything = 0
    continue;
  end
  scol = sticks_col(:,p);                    % current stick in COL format
  perp = scol(3) + pi/2;                     % perpendicular orientation to length (i.e. oriented on width)
  wid = scol(4) / aspect_ratios(p);     % width
  sidv = [cos(perp)*wid; sin(perp)*wid]/2;   % vector oriented like the shortest side, and half its length
  %lidv = (s(3:4)-s(1:2))*(params.longer(p)-1); % extra length, pointing from (x1,y1) to (x2,y2)
  lidv = 0;
  %
  % build rectangle
  R = zeros(2,4);
  R(:,1) = s(1:2) - sidv - lidv;
  R(:,2) = s(1:2) + sidv - lidv;
  R(:,3) = s(3:4) + sidv + lidv;
  R(:,4) = s(3:4) - sidv + lidv;
  %
%   for i=1:4
%     limbareas(p) = limbareas(p) + R(1,i)*R(2,mod(i,4)+1)-R(1,mod(i,4)+1)*R(2,i);
%   end
%   limbareas(p) = abs(limbareas(p)/2);
  % draw rectangle
  out(:,:,p) = roipoly(out(:,:,p), R(1,:), R(2,:));
  %
  % add to output
  out(:,:,nlimbs+1) = out(:,:,nlimbs+1) | out(:,:,p);
  %calculate minimum bb around whole stickman
  if nargout > 1
    limbMinBBcropped(:,p) = ([min(R(1,:)) min(R(2,:)) max(R(1,:)) max(R(2,:))]');

    % calc minimum bb around whole body
    if(limbMinBBcropped(1,p) < bodyMinBBcropped(1))
      bodyMinBBcropped(1) = limbMinBBcropped(1,p);
    end
    if(limbMinBBcropped(2,p) < bodyMinBBcropped(2))
      bodyMinBBcropped(2) = limbMinBBcropped(2,p);
    end
    if(limbMinBBcropped(3,p) > bodyMinBBcropped(3))
      bodyMinBBcropped(3) = limbMinBBcropped(3,p);
    end
    if(limbMinBBcropped(4,p) > bodyMinBBcropped(4))
      bodyMinBBcropped(4) = limbMinBBcropped(4,p);
    end

    % crop minBB to be in the image boundary (due to rounding errors)
    limbMinBBcropped(:,p) = round([max(limbMinBBcropped(1:2,p),1)' min(limbMinBBcropped(3,p),siz(1)) min(limbMinBBcropped(4,p),siz(2))]'); 
  end
end

%handling occlusions if classname given
%torso may be occluded by any upper or lower arm
if ~isempty(classname)
  classid = class_name2id(classname);
  switch classid
    case {1,2}
      out(:,:,1) = out(:,:,1) - (out(:,:,1) & (out(:,:,2) | out(:,:,3) | out(:,:,4) | out(:,:,5)));
      %upper arms may be ocluded by lower arms;
      out(:,:,2) = out(:,:,2) - (out(:,:,2) & (out(:,:,4) | out(:,:,5)));
      out(:,:,3) = out(:,:,3) - (out(:,:,3) & (out(:,:,4) | out(:,:,5)));
    case {3}
      out(:,:,1) = out(:,:,1) - (out(:,:,1) & (out(:,:,2) | out(:,:,3) | out(:,:,6) | out(:,:,7)));
      %upper arms may be ocluded by lower arms;
      out(:,:,2) = out(:,:,2) - (out(:,:,2) & (out(:,:,6) | out(:,:,7)));
      out(:,:,3) = out(:,:,3) - (out(:,:,3) & (out(:,:,6) | out(:,:,7)));
    otherwise
      error('unsupported class');
  end
end

if limbclasses
%   tempareas = zeros(size(limbareas));
  limbTypes = size(colmask,2);
  limbClasses = size(colmask,1);
  assert(limbTypes == nlimbs); % number of sticks should match ncols in colmask
  %convert limbs masks into limb classes masks
  for c=1:limbClasses
    idx = find(colmask(c,:));
    limbMinBBcropped(:,c) = limbMinBBcropped(:,idx(1));
    out(:,:,c) = out(:,:,idx(1));
    for i=2:length(idx)
      out(:,:,c) = out(:,:,c) | out(:,:,idx(i));
      limbMinBBcropped(:,c) = [min(limbMinBBcropped(1:2,c),limbMinBBcropped(1:2,idx(i)))'...
                               max(limbMinBBcropped(3:4,c),limbMinBBcropped(3:4,idx(i)))']';
%       tempareas(c) = tempareas(c) + limbareas((idx(i)));
    end
  end
  out(:,:,limbClasses+1) = out(:,:,end); % on the last plane all body parts at once
  out(:,:,limbClasses+2:end) = [];
  limbMinBBcropped(:,limbClasses+1:end) = [];
%   limbareas = tempareas(1:limbClasses);
end

bodyMinBBcropped = [floor(bodyMinBBcropped(1:2)) ceil(bodyMinBBcropped(3:4))];
bodyMinBBcropped = [max(bodyMinBBcropped(1:2),1) min(bodyMinBBcropped(3),siz(1)) min(bodyMinBBcropped(4),siz(2))];

end
