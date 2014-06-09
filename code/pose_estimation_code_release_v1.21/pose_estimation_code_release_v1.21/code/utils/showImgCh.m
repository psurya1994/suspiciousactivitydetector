function fig = showImgCh(img,fig,clim,cmap,ax,titles)
% fig = showImgCh(img,fig,fun,cmap)
%
% displays each channel (last dimension of 3 dimensional matrix)
% fun is a function handler used for displaying if not provided imagesc is assumed 
%
% fig - if given showImgCh returns fig else returns handle to new figure
% clim - [low high] - limits on image values;
% cmap - if given specifies what color map should be used for the figure
% ax - string passed to axis option for every subplot eg. 'on' or 'off'
% titles - cell array with titles of subplots
% default value assumed if [] passed as parameter
  if nargin < 6
    titles = [];
  end
  if nargin < 5 || isempty(ax)
    ax = 'on';
  end
  if nargin < 4 || isempty(cmap)
    cmap = 'default';
  end
  if nargin < 3 || isempty(clim)
    clim = [];
  end
  if nargin < 2 || isempty(fig)
    fig = figure;
  end
  
  if isempty(clim)
    clim = [min(img(:)) max(img(:))];
  end
  
  figure(fig);
  N = size(img,3);
  for i=1:N
    subplot(ceil(sqrt(N)),ceil(sqrt(N)),i)
    if isempty(clim)
      imagesc(img(:,:,i));
    else
      imagesc(img(:,:,i),clim);
    end
    goodaxis;
    if ~isempty(titles)
      title(titles{i});
    end
    axis(ax);
  end
  colormap(cmap);

end