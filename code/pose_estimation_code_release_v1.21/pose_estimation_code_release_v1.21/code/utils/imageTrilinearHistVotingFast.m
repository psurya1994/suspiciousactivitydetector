function [hist, binC] = imageTrilinearHistVotingFast(im,quan,weights)
% labImageTrilinearHistVotingFast(im,quan,weight)
% im - input image in Lab color space
% dims - quantization in L a b dimensions
% weight - for each pixel in the corresponding image <0-1>

% !!!image's channels must be scaled accross interval <0,255> !!!
  if ~isfloat(im)
    error('im - must be floating point array use double(im) ');
  end
  siz = size(im);
  if numel(siz) ~= 3 || siz(3) ~= 3
    error('input image must be 3 channel image');
  end
  if isscalar(quan)
    quan = [quan quan quan];
  end
  if numel(quan) ~= 3
    error('bins - incorrect size');
  end

  N = siz(1)*siz(2); % n of pixels
  
  im = reshape(im,N,siz(3));
  
  if nargin == 3
    weights = reshape(weights,N,size(weights,3));
    [hist, binC] = triLinearVoting(im,quan,weights);
  else
    [hist binC] = triLinearVoting(im,quan);
  end
  
  
  
  
  
  
  
  
  
  
  
  
end



% lower = floor((im(:,:)-repmat(binC(1,:),N,1))./repmat(sizbin,N,1)) + 1;  
   
% closest = lower(:,1)*quan(1) + lower(:,2)*quan(2)+lower(3);


%   index = 1;
%   for b = 1:quan(3)
%       for g = 1:quan(2)
%         for r = 1:quan(1)
%         temp(index) = n{1}(mod(r,quan(1)+1));
%         temp(dim+index) = n{2}(mod(g,quan(2)+1));
%         temp(2*dim+index) = n{3}(mod(b,quan(3)+1));
%         index = index+1;
%         end
%       end
%   end
