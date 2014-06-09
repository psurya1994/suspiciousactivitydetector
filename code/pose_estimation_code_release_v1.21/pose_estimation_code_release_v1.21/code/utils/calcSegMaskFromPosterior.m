function out = calcSegMaskFromPosterior(im,posterior,quantiz,verbose)

  if nargin < 4
    verbose = 0;
  end
% if verbose > 1
%   fig = figure;
%   imgpanel = uipanel('Title','Image','FontSize',12,...
%               'Position',[.0 0.7 1.0 .3]);
%   axes('Parent',imgpanel);
%   image(im);
%   goodaxis;
%   axis off;
%   segmaskpanel = uipanel('Title','Segmentation masks','FontSize',12,...
%              'Position',[.0 .0 1.0 0.7]);
%   axes('Parent',segmaskpanel);
% end
  
  if isscalar(quantiz)
    quantiz = [quantiz quantiz quantiz];
  end
  if numel(quantiz) ~= 3
    error('bins - incorrect size');
  end
  sizbin = 255./quantiz';
  n=cell(3,1);
  for i=1:3
    n{i} = 0 + sizbin(i)/2 : sizbin(i) : (255-sizbin(i)/2);
  end
  [l,a,b] = ndgrid(n{1},n{2},n{3});
  imsize = size(im);
  classes = size(posterior,2);
  %im = labImgScale255(rgb2lab_vitto(im),'forward');
  im = reshape(im,imsize(1)*imsize(2),imsize(3));
  out = triLinearInterpolation(im,quantiz,posterior);
  out = reshape(out,imsize(1),imsize(2),classes);
   
%   out = zeros(imsize(1)*imsize(2),classes);
%   for clas = 1:classes
%     temp = reshape(posterior(:,clas),quantiz(1),quantiz(2),quantiz(3));
%     out(:,clas) = interp3(a,l,b,temp,im(:,2),im(:,1),im(:,3));
%     out(isnan(out(:,clas)),clas) = 0;
%   end
%   out = reshape(out,imsize(1),imsize(2),classes);
  
  if verbose > 2
    fig = showImgCh(out,[],[],'Jet','off',{'Torso','Upper Arms','Lower Arms', 'Head'});
    set(fig,'Name','Image segmentation using posterior color model');
  end
    
end
