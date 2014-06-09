function imMask = HighlightImageMultiAreas(im, mask, classname, coeff)

% comment to be written
%

if nargin < 4
  coeff = 128;
end

% find out color planes
switch lower(classname)
    case  'full'
        % limb ids in full body model:
        % 1 = torso; 2-3 = upper arms; 4-5 = upper legs; 6-7 = lower arms; 8-9: lower legs; 10 = head
        mask_color_planes = [1 2 2 2 2 3 3 3 3 3];      % 1 -> red, 2-> green, 3 -> blue
    case 'ubf'
        % limb ids in upper-body model:
        % 1 = torso; 2-3 = upper arms; 4-5 = lower arms; 6 = head
        mask_color_planes = [1 2 2 3 3 3];     
    case '-full'
        mask_color_planes = -[1 2 2 2 2 3 3 3 3 3]; 
    case '-ubf'
        mask_color_planes = -[1 2 2 3 3 3];    
    otherwise
        error([mfilename ': unknown class ' classname]);
end

% highlight
imMask = im;
for lix = 1:length(mask)
    if not(isempty(mask{lix}))
        imMask = HighlightImage(imMask, mask{lix}, mask_color_planes(lix), coeff);
    end
end
