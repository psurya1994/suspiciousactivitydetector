function [im_all im_indiv] = PaintSticks(sticks, siz, cols, drawfullskeleton)
% paint sticks(:,p) on an image, using cols(p,:).
%
% Input:
% - sticks(:,p) = [x1 y1 x2 y2]' = integers
% - siz = [width height]
% - cols(p,:) = [r g b] in [0,1]
%
% Output:
% - im_all(:,:,:) in [0,255], all sticks in a color image (thickness = 2)
% - im_indiv(:,:,stix) = binary image for stix-th stick (thickness = 1)
%
sticks = round(sticks);

if nargin < 4
  drawfullskeleton = true; % for example when using sticks for evaluating FGHigh
end

Nsticks = size(sticks,2);
not_occluded_parts = find(~me_isEmptyStick(sticks));

im_all = zeros(siz(2),siz(1),3,'uint8');
im_indiv = false(siz(2),siz(1),Nsticks);


for p = not_occluded_parts
  im_all = DrawLineOnImgCol(im_all, sticks(2,p), sticks(1,p), sticks(4,p), sticks(3,p), uint8(cols(p,:)*255), 4); % last param -> thickness
  im_indiv(:,:,p) = DrawLineOnImg(im_indiv(:,:,p), sticks(2,p), sticks(1,p), sticks(4,p), sticks(3,p), 1); % last param -> color
end

if drawfullskeleton == 2
  orderedsticks = true;
else
  orderedsticks = false;
end

if drawfullskeleton
  links = calcSkeletonLinks(sticks,orderedsticks);
  for i=1:size(links,2)
    im_all = DrawLineOnImgCol(im_all, links(2,i), links(1,i), links(4,i), links(3,i), [255 255 255], false); % last param -> thickness
  end
end
  


end
