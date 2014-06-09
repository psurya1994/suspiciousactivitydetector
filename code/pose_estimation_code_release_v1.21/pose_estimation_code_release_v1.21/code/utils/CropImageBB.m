function [imBB, bb] = CropImageBB(im, bb)

% rectangular subimage inside bb of im.
% if bb extends out of the image,
% it is cropped to fit it.
%
% input bb = [minx miny width height]
%
% size(im) = [height width N] (any dimensionality welcome ;)
%

xmin = max(bb(1),1);
xmax = min(bb(1)+bb(3)-1, size(im,2));
ymin = max(bb(2),1);
ymax = min(bb(2)+bb(4)-1, size(im,1));
bb = [xmin ymin xmax-xmin+1 ymax-ymin+1];

imBB = im(ymin:ymax, xmin:xmax, :);
