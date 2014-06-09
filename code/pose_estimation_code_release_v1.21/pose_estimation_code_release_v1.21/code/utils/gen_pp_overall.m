function colIm = gen_pp_overall(partIms, pixs, rgbs)

% generate nice colorful per-pixel 'posterior' 
% by stacking up all per-part ones.
%
% can choose which parts to draw in pixs
% if not given -> draw them all (supports an occlusion model)
%
% pixs = list of body parts to draw
% rgbs = list of color planes to use (part pixs(p) uses color plane rgbs(p))
%

if nargin < 2
  pixs = 1:size(partIms,3);
end

% turn into binary vector
temp = zeros(1,size(partIms,3));
temp(pixs) = 1;
pixs = temp;

% default rgb planes used:
%
% for full body model:
% red:   torso (part 1)
% blue:  upper arms, upper legs
% green: lower arms, lower legs, head
% colors get superposed, so white means everything, yellow means torso
% and upper limbs, cyan means upper + lower arms or legs, and so on
%
% for upper-body model:
% red:   torso
% blue:  upper arms/legs
% green: lower arms/legs, head

imy = size(partIms,1);
imx = size(partIms,2);
numTypes = size(partIms,3);

if nargin < 3
  switch numTypes
  case 10   % full body
    rgbs = [1 3 3 3 3 2 2 2 2 2];
  case 6    % upper body
    rgbs = [1 3 3 2 2 2];
  otherwise
    error(['unknown number of body parts ' num2str(numTypes)]);      
  end
else
  % turn rgbs into rgbs(pix) = color_plane, with pix absolute
  temp = zeros(1,size(partIms,3));
  temp(logical(pixs)) = rgbs;
  rgbs = temp;
end

colIm = zeros(imy,imx,3);
for i = 1:3
  colIm(:,:,i) = sum(partIms(:,:,rgbs==i & pixs),3);
end
