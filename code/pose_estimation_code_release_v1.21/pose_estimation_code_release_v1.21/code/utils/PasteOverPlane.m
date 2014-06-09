function imOut = PasteOverPlane(source_im, dest_im, bb, plane, transp_src, transp_dst)

% paste binary image source_im over dest_im so that it occupies bounding-box bb.
%
% both source_im and dest_im must be uint8.
% source_im must be just one plane, whereas dest_im must have 3 color planes
%
% if transp
% -> semi-transparent pasting: output = trasp_src*source_im + transp_dest*dst_im
%
% bb = [min_x min_y width height]
%
% if bb goes out of image
% -> resize it accordingly (crush bb to fit, not crop properly !)
%

% process arguments
if nargin < 6
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
if size(source_im,3) ~= 1 || size(dest_im,3) ~=3
  error([mfilename ': source and destination images must have 1 and 3 color planes respectively.']);
end

% make sure bb fits in dest_im
bb(3) = min(bb(3),size(dest_im,2));
bb(4) = min(bb(4),size(dest_im,1));

% rescale source image to size of bb (if necessary)
if not(size(source_im,1) == bb(4) && size(source_im,2) == bb(3))
  source_im = imresize(source_im, [bb(4) bb(3)], 'bilinear');
end

% pasting
K = uint8((source_im == 0));       % pixels not to copy anyway
dest_im = dest_im * transp_dst;    % must 'darken' all planes of dest_im to avoid artefacts !
O = dest_im( bb(2):(bb(2)+size(source_im,1)-1),  bb(1):(bb(1)+size(source_im,2)-1),  plane);   % bb portion of dest_im(:,:,plane)
if ~transp
  imPaste = (O .* K) + ((1-K).*source_im);
else
  imPaste = (O .* K) + ((1-K).* (min(transp_src*source_im + O, 255)) );
end
imOut = dest_im;
imOut( bb(2):(bb(2)+size(source_im,1)-1),  bb(1):(bb(1)+size(source_im,2)-1),  plane) = imPaste;
