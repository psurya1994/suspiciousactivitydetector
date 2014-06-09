function [fgP, imHist, unqCols, Pfg, Pbg, fg_segm] = buildHistExps(im, resp)

% Takes a RGB image im as input and responsibilities for limb pixel types in resp(x,y,part);
%
% if size(resp,3) == 6 -> torso, upper arm, upper leg, lower arm, lower leg, head
% (left/right limbs condensed before invoking this function)
%
% if size(resp,3) == 10 -> full body model
%
% Every part is treated independently !
%
% Output:
% - fgP(bin,p) = posterior probability that a pixel of color bin belongs to part p
%   (fgP is the posterior p(fg|c), including already bg model, which is from all other pixels)
% - Pfg(bin,p), Pbg = likelihoods for foreground and background p(c|fg), p(c|bg)
% - imHist(:) = stacked up pixels, with one index in unqCols each
% - unqCols(k) = index of a 16x16x16 bin in RGB colorspace used in this image
% - fg_segm(:,:,p) = pixel-wise color segmentation for part p
%


imHist = imvq16(im);                        % vector-quantize image pixel colors in to 16x16x16 RGB bins; imHist(x,y) in [1,4096]
%figure; imagesc(imHist); axis equal;       % VF: show color quantization
siz = size(im(:,:,1));                      % image size
N = prod(siz);                              % number of pixels in image; unqCols = indeces in the 16x16x16 RGB cube
[unqCols,dummy,imHist] = unique(imHist(:)); % unique colors in image ...
numBin = length(unqCols);                   % ... = number of bins in the output histos
numTypes = size(resp,3);                    % number of limb types
imHist = imHist(:);                         % stack up all pixels; imHist(pix) = index in unqCols
resp = reshape(resp,N,numTypes);            % stack up all pixels in resp
resp = resp./repmat(sum(resp),N,1);         % normalize
fgP = zeros(numBin,numTypes);               % output fgP(bin,limb_type)

% store likelihood prob of a color given fg and given bg
% p(c|fg) and p(c|bg)
% and posterior p(fg|c)
numCols = length(unqCols);
Pfg = zeros(numCols,numTypes);
Pbg = zeros(numCols,numTypes);
for k = 1:numBin
   inds = find(imHist == k);
   fg = sum(resp(inds,:));
   bg = (length(inds) - fg)/(prod(siz) - 1);
   fgP(k,:) = fg./(fg + bg);                % posterior: p(fg|c)
   % fgP(k,limb_type) = posterior prob that a pixel of color k belongs to limb_type rather than to the background
   Pfg(k,:) = fg;                           % likelihood: p(c|fg)
   Pbg(k,:) = bg;
end


% Pbg sums to 1 
% Pfg nearly, due to roundoff errors
% correct this here:
for t = 1:numTypes
  Pfg(:,t) = Pfg(:,t) / sum(Pfg(:,t));
end


% play around with segmentations
if false
for p = 1:numTypes
  figure; imagesc(reshape(fgP(imHist,p),siz)); axis equal;
  title('posterior p(fg|c) = p(c|fg) / (p(c|fg) + p(c|bg))');
   % for each pixel in image, show prob that it is part of limb p,
   % based only on color models
   %figure; imagesc(reshape(Pfg(imHist,p),siz)); axis equal;  % only p(c|fg),
   %title('likelihood p(c|fg)');
   %without accounting for bg
   %figure; imagesc(reshape(fgP(imHist,p),siz)>.5); axis equal;
   %title('posterior p(fg|c) > .5');
   figure; imagesc(reshape(fgP(imHist,p),siz)>.75); axis equal;
   title('posterior p(fg|c) > .75');
end
end


% compute output color segmentations
fg_segm = zeros([siz numTypes]);
for p = 1:numTypes
  fg_segm(:,:,p) = reshape(fgP(imHist,p),siz)>.75;
end
