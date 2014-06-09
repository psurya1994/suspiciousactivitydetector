function mask = bbs2mask(bbs, mask_size)

% output binary mask of size mask_size = [width height],
% with pixels set to 1 if they fall inside any of the input bbs,
% and 0 otherwise.
%
% bbs(:,bbix) = [xmin ymin xmax ymax]'
%

% create mask
mask = false(mask_size([2 1]));

% fill it in
for bbix = 1:size(bbs,2)
  bb = round(bbs(:,bbix));
  bb = CropBB(bb, mask_size);
  mask(bb(2):bb(4), bb(1):bb(3)) = true;
end
