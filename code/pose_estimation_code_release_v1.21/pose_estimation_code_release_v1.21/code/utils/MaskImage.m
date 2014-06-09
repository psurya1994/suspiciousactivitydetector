function im = MaskImage(im, mask)

% COMMENT TO BE WRITTEN
%

for k = 1:size(im,3)
    imk = im(:,:,k);
    imk(~mask) = 0;
    im(:,:,k) = imk;
end
