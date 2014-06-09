function [m,theta] = mydetGMc(im,sigma,nbins)
% function [m,theta] = mydetGM(im,sigma,nbins)
%
% Compute image gradient magnitude.
%
% INPUT
%	im	Image.
%	sigma	Scale at which to compute image derivatives.
%
% OUTPUT
%	m	Gradient magnitude.
%	theta	Orientation of gradient + pi/2 
%		(i.e. edge orientation).
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

if nargin < 3,
    nbins = 8;
    nbins = 24;
end

%if isrgb(im), im=rgb2gray(im); end
%idiag = norm(size(im));

%if nargin<2, sigma=2; end
%sigma = max(0.5,sigma);

% compute x and y image derivatives
% use elongated Gaussian derivative filters
%fb = cell(2,1);
fb1 = oeFilterc(sigma*[1 1],3,pi/2,1);
fb2 = fb1';
im = double(im);
[sizy sizx sizz] = size(im);

[dx,dy] = deal(zeros(size(im)));
for i = 1:sizz
  dx(:,:,i) = filter2(fb1,im(:,:,i));
  dy(:,:,i) = filter2(fb2,im(:,:,i));
end
m = sqrt(dx.^2 + dy.^2);

if sizz == 3,
    [m,I] = max(m,[],3);
    inds = sizy*sizx*(I - 1) + reshape([1:(sizy*sizx)],sizy,sizx);
    dx = dx(inds);
    dy = dy(inds);
end
theta = mod(atan2(dy,dx)+pi/2,pi);
m = nonmax(m,theta);
step = 2*pi/nbins;
thets = [0:step:2*pi-step];
[binu,binv] = pol2cart(thets,1);
[dummy, I] = min(dx(:)*binu + dy(:)*binv,[],2);
theta = reshape(I,size(dx));

%Zero-out border effects
r = floor(size(fb1,1)/2);
m([1:r end-r:end],:) = 0;
m(:,[1:r end-r:end]) = 0;
