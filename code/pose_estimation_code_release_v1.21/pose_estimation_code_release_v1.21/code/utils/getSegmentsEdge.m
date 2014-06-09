function resp = getSegmentsEdge(im, maskOrig, orient_ids)
%resp = getSegmentsColEdge(im,patch)

%Useful debugging
%im = zeros(200,100,3);
%im(100:130,30:40,:) = 255;
%im = uint8(im);

currAngles = 15:15:360;
currAngles = currAngles - 15;

[imy,imx,imz] = size(im);

imo = length(currAngles);

if nargin < 3
  orient_ids = 1:imo;
end

resp = zeros([imy imx imo]);

lensc = .5;
widsc = 1;

[len,wid] = size(maskOrig);
for dir = orient_ids
  %currAngle = dir * dtheta;
  currAngle = currAngles(dir);
  %currIm = imrotate(im,-currAngle);
  mask = imrotate(maskOrig,currAngle);
  [mm,nn] = size(mask);

  thet = currAngle*pi/180;
  offset = round([-sin(thet) cos(thet)]*wid);
  imResp = filter2(mask,im);
  %Zero-out border effects
  %imResp([1:floor(mm/2) imy-ceil(mm/2):end],:) = 0;
  %imResp(:,[1:floor(nn/2) imx-ceil(nn/2):end]) = 0;
  resp(:,:,dir) = imResp;
end

%Debug
%[dummy,I] = max(resp(:)); [y,x,ii] = ind2sub(size(resp),I);
