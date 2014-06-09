function im = labImgScale255(im,way)
%im = labScale01(img,way)
% scales all channels of Lab color space image im to <0,1> or rescale them back
% depending upon way value ( 'forward' or 'backward' allowed)
% forward is default operation (if way not given)
% rescaling color map done if im is two dimensional [nOfIndexes,3]

if nargin < 2
  way = 'forward';
end

switch numel(size(im))
  case 3 % normal 3 dimensional image
    switch way
      case 'forward'
        im(:,:,1) = 255/100*im(:,:,1);
        im(:,:,2) = 255/220*(im(:,:,2)+110);    
        im(:,:,3) = 255/220*(im(:,:,3)+110);
      case 'backward'
        im(:,:,1) = 100/255*im(:,:,1);
        im(:,:,2) = 220/255*im(:,:,2)-110;    
        im(:,:,3) = 220/255*im(:,:,3)-110;
      otherwise
        error('incorrect second param only ''forward'' and ''backward'' allowed');
    end
  case 2 % color map
    switch way
      case 'forward'
        im(:,1) = 255/100*im(:,1);
        im(:,2) = 255/220*(im(:,2)+110);    
        im(:,3) = 255/220*(im(:,3)+110);
      case 'backward'
        im(:,1) = 100/255*im(:,1);
        im(:,2) = 220/255*im(:,2)-110;    
        im(:,3) = 220/255*im(:,3)-110;
      otherwise
        error('incorrect second param only ''forward'' and ''backward'' allowed');
    end
end

end