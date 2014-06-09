function [Y cmap] = imvq16(im)

% Vector quantizes an image into RGB bins of size 16x16x16
% Values will lie 1 to 4096
%

dim = 16;
n = 1/(2*dim):1/dim:1-1/(2*dim);
[R,G,B] = ndgrid(n,n,n);
cmap = [R(:) G(:) B(:)];
Y = uint16(double(rgb2ind(uint8(im),cmap,'nodither'))+1);

