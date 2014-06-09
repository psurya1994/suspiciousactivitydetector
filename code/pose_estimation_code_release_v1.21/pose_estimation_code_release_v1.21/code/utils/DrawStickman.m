function hdl = DrawStickman(sticks, img, colors, thickness, drawidx, drawfullskeleton)
% hdl = DrawStickman(sticks, img)
% Overlays segments in 'sticks' on image 'img'
%
% Input:
%  - sticks: matrix [4, nparts]. sticks(:,i) --> (x1, x2, y1, y2)' 
%  - img: image to show. Can be [].
%  - colors - vector 3x1 then use the specified color for all sticks 
%              or 3xN then color per body part
%  - thickness - thickness of sticks (default 4)
%  - drawidx - draw part index next to the stick (detault true)
%  - drawfullskeleton - if 0 then draw just sticks
%                       if 1 then draw connections (try to figure out the tree structure)
%                       if 2 then draw connections (always sticks(1:2,i) is the limb anchor)
%
% Output:
%  - hdl: figure handler.
%
% See also ReadStickmenAnnotationTxt
%
% MJMJ/2008 modified by Eichner/2009
%


hdl = -1;

if nargin < 2,
   img = [];
end

Nsticks = size(sticks, 2);
not_occluded_parts = find(~me_isEmptyStick(sticks));

if nargin < 3 || isempty(colors)
  colors = sampleColors()';
elseif size(colors,2) == 1
  colors = repmat(colors,1,Nsticks);
elseif size(colors,2) ~= Nsticks
  error('if colors > 1 then ncolors == nsticks');
end

if nargin < 4 || isempty(thickness)
  thickness = 4;
end

if nargin < 5 || isempty(drawidx)
  drawidx = true;
end
if nargin < 6 || isempty(drawfullskeleton)
  drawfullskeleton = false;
end

if ~isempty(img),
   hdl = imshow(img,'Border','Tight');
   hold on
end

for p = not_occluded_parts
   X = sticks([1 3],p);
   Y = sticks([2 4],p);   
   hl = line(X', Y');
   set(hl, 'color', colors(:, p) );
   set(hl, 'LineWidth', thickness);
   if drawidx
     ht = text(double(5+mean(X)), double(5+mean(Y)), num2str(p));
     set(ht, 'color', colors(:, p) );
     set(ht , 'FontWeight', 'bold');
   end
end

if drawfullskeleton
    % drawfullskeleton == 1 figure out tree structure yourself
    % drawfullskeleton == 2 follow the tree as anchor is always sticks(1:2,:)
    links = calcSkeletonLinks(sticks,drawfullskeleton-1); 
    for i=1:size(links,2);
      X = links([1 3],i);
      Y = links([2 4],i);
      hl = line(X', Y');
      %set(hl, 'color', [.99 .99 .99] );
      set(hl, 'color', [.7 .7 .7] );
      set(hl, 'LineWidth', thickness/2);
    end
end


 
