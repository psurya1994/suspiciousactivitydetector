function [g, kernel_width] = vgg_gsmooth(x,sd, opt)
% VGG_GSMOOTH   Gaussian smoothing.
%           VGG_GSMOOTH(X, SIGMA) applies gaussian smoothing of std. dev.
%           SIGMA to the vector or matrix X.
%           VGG_GSMOOTH(X, SIGMA, OPT) allows specification of the options
%           'full', 'same', 'valid' as with CONV2.
%           'clamp' extends image boundaries with boundary values
%           'wrap' assumes periodic
%           'wrap-' assumes periodic, sign flipped.

%           See also CONV2.

% Author: Andrew Fitzgibbon, Edinburgh University AI Dept.
% Email: andrewfg@ed.ac.uk
% Date: 04 Oct 95

if nargin < 3
  opt = 'full';
end

X = [-500:500];
G = vgg_gaussian_1d(X,sd);
G = max(G) * round(G / max(G) * 256);
G = G(find(G));
G = G./sum(G);
% fprintf(2, 'gsmooth: Gaussian Width = %g',Gc);
kernel_width = size(G,2);

[m, n, channels] = size(x);

% For multiple-channel images, smooth each separately
if channels > 1
  for channel=channels:-1:1
    g(:,:,channel) = vgg_gsmooth_aux(x(:,:,channel), sd, opt, G);
  end
else
  g = vgg_gsmooth_aux(x, sd, opt, G);
end

function g = vgg_gsmooth_aux(x, sd, opt, G)

[m, n] = size(x);
[Gr, Gc] = size(G);

if strcmp(opt, 'clamp')
  l = ones(size(G));
  if m == 1
    g = conv([l*x(1) x l*x(n)], G);
    g = g(Gc + floor(Gc/2) + [1:n]);
  elseif n == 1
    g = conv([l'*x(1) ; x ; l'*x(m)], G);
    g = g([1:m] + Gc + floor(Gc/2) );
  else
    error('no 2d clamp');
  end
elseif strcmp(opt, 'wrap') | strcmp(opt, 'wrap-')
  if strcmp(opt, 'wrap-')
    wrapfac = -1;
  else
    wrapfac = 1;
  end
  if Gc > length(x)
    error('zoiks');
  end
  if m == 1
    new_x = wrapfac*[x(n-Gc:n) wrapfac*x x(1:Gc)];
    g = conv(new_x, G);
    g = g(Gc:Gc+n);
  elseif n == 1
    new_x = wrapfac*[x(m-Gc:m); wrapfac*x; x(1:Gc)];
    g = conv(new_x, G);
    g = g([1:m] + Gc + floor(Gc/2) );
  else
    error('no 2d wrap');
  end
else
  if m == 1
    g = conv2(x, G,opt);
  elseif n == 1
    g = conv2(x, G',opt);
  else
    % Reverse order for speed
    g = conv2(conv2(x,G,opt), G', opt); % / sum(sum(conv2(G',G)));
  end
end

% fprintf(2,' done\n');
