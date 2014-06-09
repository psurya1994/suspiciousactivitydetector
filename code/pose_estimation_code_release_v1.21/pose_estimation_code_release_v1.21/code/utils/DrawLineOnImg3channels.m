function im = DrawLineOnImg3channels(im, y1, x1, y2, x2, col)

% COMMENT TO BE WRITTEN
%

im(:,:,1) = DrawLineOnImg(im(:,:,1), y1, x1, y2, x2, col(1));
im(:,:,2) = DrawLineOnImg(im(:,:,2), y1, x1, y2, x2, col(2));
im(:,:,3) = DrawLineOnImg(im(:,:,3), y1, x1, y2, x2, col(3));
