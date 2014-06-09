function [tot_e, es] = TotalPoseEntropy(respIm)

% total entropy of the pose space respIm
% the distributions respIm are assumed normalized,
% i.e. sum(sum(sum(respIm(:,:,:,p)))) = 1 for all p
%
% The lower the better; best is 0
%
% also returns entropy per part es(pix)
%

% shortcuts
[imy imx imo numTypes] = size(respIm);

% compute entropy
tot_e = 0;
es(numTypes) = 0;
for p = 1:numTypes
  curResp = respIm(:,:,:,p);
  curResp = curResp(:);
  % remove 0 vals (avoid log of 0, and compatible with math def of entropy)
  curResp = curResp(curResp>0); 
  es(p) = -sum(curResp(:).*log(curResp(:)));
  tot_e = tot_e + es(p);
end
