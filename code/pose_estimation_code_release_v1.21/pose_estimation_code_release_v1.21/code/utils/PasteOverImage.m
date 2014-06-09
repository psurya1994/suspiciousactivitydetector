function imOut = PasteOverImage(source_im, dest_im, bb, transp_src, transp_dst)

% paste image source_im over dest_im so that it occupies bounding-box bb.
%
% dark pixels in source_im are not pasted.
%
% both source_im and dest_im must be uint8 and have 3 color planes.
%
% if transp
% -> semi-transparent pasting: output = trasp_src*source_im + transp_dest*dst_im
%
% bb = [min_x min_y width height]
%
% warning: no cropping performed ! bb must fit inside dest_im !
%

% process arguments
if nargin < 5
  transp = false;
else
  transp = true;
end

% make sure it all works even if integers passed (0, 1)
if transp
  transp_src = double(transp_src);
  transp_dst = double(transp_dst);
end

% check input format
if size(source_im,3) ~= 3 || size(dest_im,3) ~=3
  error([mfilename ': both source and destination images must have 3 color planes.']);
end

% rescale source image to size of bb (if necessary)
if not(size(source_im,1) == bb(4) && size(source_im,2) == bb(3))
  source_im = imresize(source_im, [bb(4) bb(3)], 'bilinear');
end

% dark pixels of souce image -> don't draw them
T = source_im(:,:,1) + source_im(:,:,2) + source_im(:,:,3) < 100;       %indeces of dark pixels (usual)
%T = source_im(:,:,1) + source_im(:,:,2) + source_im(:,:,3) < 15;         % dagstuhl08 final version
K(:,:,1) = T;  K(:,:,2) = T;  K(:,:,3) = T;                             % dark pixels mask
K = uint8(K);                                                           % enables weighted sum below

% paste non-dark pixels of source image over dest image
O = dest_im( bb(2):(bb(2)+size(source_im,1)-1),  bb(1):(bb(1)+size(source_im,2)-1),  :);   % bb portion of dest_im
if ~transp
  imPaste = (O .* K) + ((1-K).*source_im);
else
  imPaste = (O .* K) + ((1-K).* (min(transp_src*source_im + transp_dst*O, 255)) );
end
imOut = dest_im;
imOut( bb(2):(bb(2)+size(source_im,1)-1),  bb(1):(bb(1)+size(source_im,2)-1),  :) = imPaste;
