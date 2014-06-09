function [pm, cm] = ParseImage(im, CM, fg, mask, class_id, pars, verbose)

% parse image im using Deva's limb pose estimation from NIPS 2006
%
% returns posterior map pm and color model cm
%
% if color model cm given, skip parse1 phase and directly
% go to parse2 phase using it
%
% if pars.use_fg_high
% -> use foreground highlighting image fg to
%    discard all pixels of im not in fg;
%    facilitates parsing as remove lots of bg clutter
%
% if limb location constraints mask{lix} given
% -> prepare a separate image for each limb,
%    discarding all pixels not in the mask.
%    This option combines with fg: only the pixels in the intersection
%    of fg and mask are considered for a limb.
%
% if pars.cm_from_org_im
% -> learn color models from the original image,
%    without fghigh regardless of pars.use_fg_high.
%
% if pars.img_lik_only
% -> don't do full parsing, stop at extracting the
%    image likelihood and return it inside pm;
%    only works when cm given
%
% CM.cm(:,end)        = indeces of colors into RGB cube binned into 16x16x16 bins
% CM.cm(:,1:(end-1))  = color models for each limb
%                    (a column per limb type, each entry is prob that that color occurs on limb)
%

% process arguments
if nargin < 7
  verbose = false;
end

if isstruct(im) %motion has been passed but nothing will by done with it anyway (maybe a code update)
  im = im.im;
end

classname = class_id2name(class_id);
iterModel = load_iterModel(classname,pars);


% iterModel.www is the weight for the edge filter when used in combo with the color filter
iterModel.www = ones(size(iterModel.len));

% show original image
if verbose > 1
  figh = figure;
  subplot(2,2,1);
  imagesc(im);
  axis equal; axis tight; axis off;
  title('image');
  drawnow;
end

% get edge map
if verbose
  disp('computing edge map');
end
% even in case of fg high, det edges from org img,
% to avoid spurious edges on fg boundaries
[m, trash] = mydetGMc(im,2); m = double(imdilate(m > iterModel.edgeThresh, ones(3)));
if verbose > 1
  figure(figh); subplot(2,2,2);
  imshow(m);
  axis equal; axis tight; axis off;
  title('edges');
  drawnow;
end

% foreground highlighting
imOrg = im;                                   % useful for learning color models if pars.cm_from_org_im
if pars.use_fg_high
  fg = imdilate(fg,ones(9));                  % dilate to make sure all edges are in fg
  m(~fg) = 0;                                 % remove non-fg pixels from the edgemap m
  im = MaskImage(im, fg);                     % remove non-fg pixels from the actual image im
  if verbose
    disp('keep only pixels and edgels on foreground highlit area');
  end
  %
  if verbose > 1
    figure(figh); subplot(2,2,1);
    imagesc(im);
    axis equal; axis tight; axis off;
    subplot(2,2,2);
    imshow(m);
    axis equal; axis tight; axis off;
    drawnow;
  end
else
  fg = true(size(im,1),size(im,2));
end

% limb location constraints
if not(islogical(mask))
  % dilate and intersect with fg
  for lix = 1:length(mask)
    if not(isempty(mask{lix}))
      %moutl{lix} = imdilate(mask{lix},ones(9));  % CVPR08 figures
      %moutl{lix} = bwoutlines(moutl{lix}, 4); 
      %moutl{lix} = imdilate(moutl{lix},ones(3));      % add some thickness
      mask{lix} = imdilate(mask{lix},ones(9)) & fg;
    end
  end
  %
  if verbose > 1
    imMask = HighlightImageMultiAreas(im, mask, classname);
    %imMask = SuperimposeImageMultiAreas(im, moutl, classname); % CVPR08 figs; 1 -> not transparent
    figure(figh); subplot(2,2,1);
    imagesc(imMask);
    axis equal; axis tight; axis off;
    drawnow;
  end % verbose
end % given limb masks ?



% color model not given ?
if islogical(CM)

% parse1 phase
% Do first parse
% get first posterior map for each limb
if verbose
  disp('color model not given -> computing first parse using edges only');
end
%keyboard;
res1 = expected_genmodel_FHedgesN(m, mask, iterModel, pars.part_ids{class_id}, pars.orient_ids{class_id}, pars.relor_inhibit);
if verbose
  disp(['total pose entropy ' num2str(res1.e)]);
  disp(['total pixel confidence ' num2str(res1.p)]);
end
% Show it
if verbose > 1
  figure(figh); subplot(2,2,3); imagesc(uint8(res1.a*2e3));
  axis equal; axis tight; axis off;
  title('parse on edges');
  drawnow;
end

% Build region appearance model (given posterior maps from first parse)
% condenseLRResp -> fuses probs for left/right limbs into just 6 classes
if ~pars.cm_from_org_im
  if verbose && pars.use_fg_high
    disp('learning color models from foreground highlit image');
  end        
  [fgP imHist unqCols Pfg Pbg] = buildHistExps(im, condenseLRResp(res1.b, iterModel.colmask));
else
  % learn limb color models from full image -> better lrs vs bg ?
  [fgP imHist unqCols Pfg Pbg] = buildHistExps(imOrg, condenseLRResp(res1.b, iterModel.colmask));
  if pars.use_fg_high
    if verbose
      disp('learning color models from original image');
    end
    % transfer color model to the image that will actually be used for parsing
    % important: otherwise indeces of colors don't match !
    cm = [fgP double(unqCols)];
    [fgP imHist unqCols] = TransferColorModel(cm, im);
  end
end

else  % end of first parse and building color model
    
  % color model given
  if verbose
    disp('color model given -> transferring it to new image and skipping first parse');
  end
  [fgP imHist unqCols] = TransferColorModel(CM.cm, im);
  
end

if isfield(pars,'init_stage_only') && pars.init_stage_only
  pm = res1;
  cm.cm = [fgP double(unqCols)];
  cm.Pfg = Pfg;
  cm.Pbg = Pbg;
  return;
end


% parse2 phase
% Run the region-based body model (one iteration)
% only the region appearance model is given (in fgP),
% not the posterior map from first parse
if verbose
  disp('computing second parse using edges and color models');
  if pars.img_lik_only
    disp('only compute image likelihood');
  end
end
% select colors which are more likely to be foreground than background for a limb group
iterModel.fgP = fgP > .5;
res2 = expected_genmodel_FHedgecolsW(im, m, mask, iterModel, pars.img_lik_only, pars.part_ids{class_id}, pars.orient_ids{class_id}, pars.relor_inhibit);
% a second iter of edge+col parsing
if pars.img_lik_only
  pm = res2;
  cm = [];
  return;
end
if verbose
  disp(['total pose entropy ' num2str(res2.e)]);
  disp(['total pixel confidence ' num2str(res2.p)]);
end
% Show final parse
if verbose > 1
  figure(figh); subplot(2,2,4); imagesc(uint8(res2.a*2e3));
  axis equal; axis tight; axis off;
  title('parse on edges and color');
  drawnow;
  keyboard;
  if isfield(res2,'MAP') && ~isempty(res2.MAP)
    DrawStickman(res2.MAP.sticks);
  end
end

% finalize color models based on new posterior maps from second parse
if ~pars.cm_from_org_im
  [fgP imHist unqCols Pfg Pbg] = buildHistExps(im, condenseLRResp(res2.b, iterModel.colmask));
else
  [fgP imHist unqCols Pfg Pbg] = buildHistExps(imOrg, condenseLRResp(res2.b, iterModel.colmask));
end


% assign outputs
pm = res2;
clear cm;
cm.cm = [fgP double(unqCols)];
cm.Pfg = Pfg;
cm.Pbg = Pbg;
