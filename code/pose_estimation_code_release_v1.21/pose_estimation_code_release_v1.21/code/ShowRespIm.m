function fids = ShowRespIm(respIm)

% COMMENT TO BE WRITTEN
%
% Input:
% respIm(y,x,o,pix)
%

fids = zeros(1,size(respIm,4));
for pix = 1:size(respIm,4)
  fids(pix) = ShowRespImForPart(respIm(:,:,:,pix));
  set(fids(pix), 'name', ['Posterior marginal p(y,x,theta) for part ' int2str(pix)]); 
end
