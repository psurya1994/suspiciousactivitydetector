function Y = partshiftZ0(D,len,wid,padval)

% Pads the matrix with zeros
% Assumes D is 0-degree aligned
%

if nargin < 3,
    wid = 0;
end

[imy, imx, imo] = size(D);
stepo = 360/imo;                                                      % stepo = 15 degrees
[ux,vx] = pol2cart(([stepo:stepo:360]+180-stepo+180) * pi/180,wid);   % Width offset
[uy,vy] = pol2cart(([stepo:stepo:360]+90-stepo+180) * pi/180,len);    % Length offset


vx = -vx;
vy = -vy;
xvec = round(ux+uy);
yvec = round(vx+vy);
maxv = round(max(abs([xvec yvec])));
Y = zeros([imy imx imo]);
if nargin < 4
  D = padarray(D,[maxv maxv 0],'both');
else
  D = padarray(D,[maxv maxv 0],padval,'both');
end

for i = 1:imo,
  Y(:,:,i) = D(maxv+yvec(i)+1:imy+maxv+yvec(i),maxv+xvec(i)+1:imx+maxv+xvec(i),i);
end
