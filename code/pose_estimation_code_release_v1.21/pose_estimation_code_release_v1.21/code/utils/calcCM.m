function [cm] = calcCM(img,weights,cm_params,verbose)
% [cm backgr] = calcCM(img,weights,cm_params,what)
% computes color models given image img, corresponding per pixel weights, and color space parameters (in cm_params)
% weights may be a 3 channel image, when the last channel corresponds to number of color models to be generated.
% output is a bins x classes matrix where bins is the color model dimensionality and classes = size(weighs,3)
  if nargin < 4 || isempty(verbose)
    verbose = 0;
  end

%   if verbose > 2 
%     w = sum(weights,3);
%     w = w/max(w(:));
%     siz = size(img);
%     newImg = PasteOverPlane(uint8(w*255), uint8(img), [1 1 siz([2 1])], 2, 1.0, .4);
%     imshow(newImg);
%   end

  classes = size(weights,3);
  
%   if strcmp(cm_params.colorSpace,'rgb') || (~cm_params.normalizeLightness && strcmp(cm_params.colorSpace,'lab'))

    [cm binc] = imageTrilinearHistVotingFast(img ...
                                        ,cm_params.quantiz ...
                                        ,weights); ...     
%   elseif cm_params.normalizeLightness && strcmp(cm_params.colorSpace,'lab')
%     cm = normalizeLightnessAndVote(img ...
%                                      ,cm_params.quantiz ...
%                                      ,weights ...
%                                      ,cm_params.lightnessNormCroppFactor);
%   else
%     error('unsupported color space set in cm_params');
%   end

  for i=1:classes
    if norm(cm(:,i),cm_params.norm) ~= 0
      cm(:,i) = cm(:,i)./norm(cm(:,i),cm_params.norm);
    else
      disp(['WARNING: empty color model of limb class ' num2str(i)]);
    end
  end
    
%%%  nicely visualize color histogram - use for small number of bins
% nbins = size(binc,1);
% rgb_binc = reshape(lab2rgb_vitto(reshape(labImgScale255(binc,'backward'),[1,nbins,3])),nbins,3);
% figure,me_bar_custcolors(rgb_binc,cm(:,1));

  
end

