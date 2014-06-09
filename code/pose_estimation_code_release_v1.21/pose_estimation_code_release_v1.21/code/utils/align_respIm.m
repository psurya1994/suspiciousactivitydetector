function [respIm respImHidden] = align_respIm(respIm, model, respImHidden )

% Shift the respIm to align with center of segment
%

for p = 1:size(respIm,4)
    currResp = partshiftZ0(respIm(:,:,:,p), -model.len(p), 0);
    if nargin < 3
      norm = sum(currResp(:));
      respImHidden(p) = 0;
    else
      norm = sum(currResp(:))+respImHidden(p);
      respImHidden(p) = respImHidden(p)/norm;
    end
    currResp = currResp/norm;
    respIm(:,:,:,p) = currResp;
end
