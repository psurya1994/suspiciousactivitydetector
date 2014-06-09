function Aout = TrafoBB2img(A, bb)

% transform areas A from rel-to-bb coordinates to image coordinates.
%
% Input:
% A(:,k) = [xmin ymin xmax ymax]' = a set of rectangles defining the area,
%          in coordinates relative to a canonical bounding-box (0,0)->(1,1)
% bb = [xmin ymin width height] = the actual bb in the image at hand
%
% Output:
% Aout(:,k) = like A but in image coordinates
%

if isempty(A)
  Aout = [];
  return;
end

Aout = A .* repmat(bb([3 4 3 4])',1,size(A,2)) + repmat(bb([1 2 1 2])',1,size(A,2));
