function partIms = gen_pp_post(respIm, genmodel, part_ids)
% respIm is expected to be anchored middle-middle
% generate per-pixel posterior prob maps for each body part
% essentially this draws the parts in all positions allowed by respIm,
% weighted by the probs in respIm
%
% if part_ids given, expect size(respIm,4) == length(part_ids),
% and only generte per-pixel posteriors for these body parts
% Output is always partIms(:,:,p), with p pointing to part_ids(p)
% 

if nargin < 3
  part_ids = 1:size(respIm,4);
end

% shortcuts
imy = size(respIm,1);
imx = size(respIm,2);
imo = size(respIm,3);

partIms = zeros(imy,imx,length(part_ids));
for pix = 1:length(part_ids)
  p = part_ids(pix);
  len = genmodel.len(p);
  wid = genmodel.wid(p);
  % Condense into 180 degrees
  tmp = respIm(:,:,1:imo/2,pix) + respIm(:,:,imo/2 + 1:end,pix);
  kernal = ones(len*2,wid*2);
  % Spread by convolution
  for i = 1:imo/2
    partIms(:,:,pix) = partIms(:,:,pix) + filter2(imrotate(kernal,(i-1)*360/imo),tmp(:,:,i));
  end
end
