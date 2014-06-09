function g = vgg_gaussian_1d(u,sig)
% GAUSS      Gaussian function.
%            G = GAUSS(u, sigma);
%
%                  1                -u.^2
%            ---------------  exp(---------)
%            sigma*sqrt(2*pi)     2*sigma^2

% Author: Andrew Fitzgibbon, Edinburgh University AI Dept.
% Email: andrewfg@ed.ac.uk
% Date: 05 Oct 95

inside = 1/(2*sig*sig);
g = (sqrt((1/pi)*inside)) * exp(-u.*u.*inside);
