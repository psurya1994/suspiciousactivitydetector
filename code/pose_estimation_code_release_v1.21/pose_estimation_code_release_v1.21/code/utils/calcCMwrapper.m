function [forgr backgr img] = calcCMwrapper(img,bbt,bbtc,mask,lp_params,cm_params,what,verbose)
% usefull routine computing color models in many ways depending on options
% img - input image
% bbt - enlared bounding box in wh format
% bbtc - cropped enlared bb in wh format
% mask - (optional) if [] lp_params.LP.map is used
% lp_params - location prior and appearance transfer params
% cm_params - color space params
% what - 'fg' - foreground, 'bg' - background, 'fgbg' - both
% cm_params.bgmodel - create background model from whole image 'img', or just from enlared bounding box area 'bb'

  if nargin < 8 || isempty(verbose)
    verbose = 0;
  end
  
  
  imsize = size(img);
  bbtcMM = wh2minmax(bbtc);
 
  if isempty(mask)
    offset = [ bbtcMM(1)-bbt(1)+1 bbtcMM(2)-bbt(2)+1 bbtcMM(3)-bbt(1)+1 bbtcMM(4)-bbt(2)+1]; 
    mask = imresize_vitto(lp_params.LP.map,bbt([4 3]));
    mask = mask(offset(2):offset(4),offset(1):offset(3),:);
  else
    if max(mask(:))-1e-10 > 1 || min(mask(:))+1e-10 < 0
      error('invalid mask range, only values in the interval <0,1> are allowed');
    end
    classes = size(mask,3);
  end
  
  forgr = [];
  backgr = [];
  
  img = img(bbtcMM(2):bbtcMM(4),bbtcMM(1):bbtcMM(3),:);

  % overlay 4 channel LP on the image
%   aaa = PasteOverPlane(uint8(mask(:,:,1)*255),img,[1 1 size(img,2) size(img,1)],1,1,1);
%   aaa = PasteOverPlane(uint8(mask(:,:,2)*255),aaa,[1 1 size(img,2) size(img,1)],2,1,1);
%   aaa = PasteOverPlane(uint8(sum(mask(:,:,[3 4]),3)*255),aaa,[1 1 size(img,2) size(img,1)],3,1,1);
  
  switch cm_params.colorSpace
    case 'lab'
      img = labImgScale255(rgb2lab_vitto(img));
    case 'rgb'
      img = double(img);
  end

  if strcmp(what,'fg') || strcmp(what,'fgbg')
    forgr = calcCM(img,mask,cm_params,verbose);
  end
  if strcmp(what,'bg') || strcmp(what,'fgbg');
    switch cm_params.bgmodel
      case 'bb'
        backgr = calcCM(img,1-mask,cm_params,verbose);
      case 'img'
        error('not supported anymore')
        mask_img = zeros(imsize(1),imsize(2),classes);
        % assign to weights part of the location prior corresponding to croped enlarged bb present in the image
        mask_img(bbtcMM(2):bbtcMM(4),bbtcMM(1):bbtcMM(3),:) = mask;
        backgr = calcCM(img,1-mask_img,cm_params,verbose);
      otherwise
        error('only ''img'', ''bb'' (bgmodel) strings supported');
    end
  end
  if isempty(forgr) && isempty(backgr)
    error('only ''f'', ''b'', ''a'' (what) strings supported');
  end
  if nargout < 2
    if ~isempty(forgr) && ~isempty(backgr)
      warning('calculated both forground and background color models, but only one expected');
    end
    if isempty(forgr) 
      forgr=backgr; % the case when what = 'b'
    end
  end

end