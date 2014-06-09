function im = PaintSegmentation(segm, cols)

% paint segmentation outlines from segm(:,:,p) on an image, using cols(p,:).
%
% Input:
% - segm(:,:,p) = logical
% - cols(p,:) = [r g b] in [0,1]
%
% Output:
% - im(:,:) in [0,255]
%

[height width numTypes] = size(segm);
im = zeros(height,width,3,'uint8');
for p = 1:numTypes
  imp = segm(:,:,p);
  if all(imp(:)==0)                    % no body part
    continue;
  end
  % convert to outlines
  imp = bwoutlines(imp, 4);            % outlines of a mask image, with thickness == 4
  %
  ixs = find(imp);
  rp = im(:,:,1);                      % red plane
  rp(ixs) = cols(p,1);
  im(:,:,1) = rp*255;
  %
  gp = im(:,:,2);                      % green plane
  gp(ixs) = cols(p,2);
  im(:,:,2) = gp*255;
  %
  bp = im(:,:,3);                      % blue plane
  bp(ixs) = cols(p,3);
  im(:,:,3) = bp*255;
end
