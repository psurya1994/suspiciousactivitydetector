function Y = condenseLRResp(resp, mask)

% Condenses the left right versions of a resp = [imy,imx,numTypes]
%

[imy,imx,numTypes] = size(resp);
numOut = size(mask,1);
Y = zeros(imy,imx,numOut);
for i = 1:numOut,
  Y(:,:,i) = sum(resp(:,:,find(mask(i,:))),3);
end
