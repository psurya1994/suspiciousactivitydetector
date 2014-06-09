function [imProc, mask, bb_tot, bb_tot_stdscale, bb_tot_stdscale_nocrop, fid] = PrepareFGHigh(im, bb, pars, class_id, verbose, fr)

% prepare image imProc and labels mask for foreground highlighting for bounding-box bb in image im.
%
% Input:
% - bb = [minx miny width height]
% - pars.abg_area .bg_area .fg_area .afg_area -> areas of init for grabcut, corresp to codes 64, 96, 160, 192 (128 == neutral)
% - pars.border = list of area labels to give to each side of a 2% border around the actual area to be processed
%   (useful to have some label at the border even if the enlarged bb goes out of the image
%   cell array of pairs of strings, e.g. pars.border = {'left' 'bg' 'top' 'bg' 'right' 'bg'}
%   all mask pixels on the border gets reassigned (not just neutral ones)
%   the border overwrites any other previous assignment
% - verbose > 1 -> open nice figure with all preparations; must provide fr too !
%
% Output:
% imProc = standard-sized image, corresp to the area bb_tot in the original image im
% mask = labels to be input to grabcut, for image imProc
% bb_tot = area of im where fg highlighting it to be applied = [xmin ymin width height]
% bb_tot_stdsize = bb in the standard-sized image (directly contains imProc)
% bb_tot_stdsize_nocrop = before cropping
% fid = figure id that it opened if verbose > 1
%

% process arguments
if nargin < 5
  verbose = false;
end
fid = -1;

% complete pars if necessary
pars = setfields_empty(pars, {'abg_area', 'bg_area', 'afg_area', 'fg_area', 'border'},cell(1,3));

% transform fb/bg areas from rel-to-bb coordinates to image coordinates
abg_area = TrafoBB2img(pars.abg_area{class_id}, bb);
bg_area = TrafoBB2img(pars.bg_area{class_id}, bb);
fg_area = TrafoBB2img(pars.fg_area{class_id}, bb);
afg_area = TrafoBB2img(pars.afg_area{class_id}, bb);

% resize the whole image so that det bb heigh = bb_rescale
sizOrg = [size(im,2) size(im,1)];
imscale = pars.bb_rescale(class_id) / bb(4);
im = imresize(im,imscale,'bilinear');

% trafo bb, bg_area, and fg_area to resized image
bb = bb * imscale;
abg_area = abg_area * imscale;
bg_area = bg_area * imscale;
fg_area = fg_area * imscale;
afg_area = afg_area * imscale;

% produce fg/bg masks (includes cropping if necessary)
siz = [size(im,2) size(im,1)];                                    % [width height]
abg_mask = bbs2mask(abg_area, siz);
bg_mask = bbs2mask(bg_area, siz);
fg_mask = bbs2mask(fg_area, siz);
afg_mask = bbs2mask(afg_area, siz);

% determine total area where to apply fg highlightning
enlarge = pars.bb_enlarge(:,class_id)';
bb_tot = enlarge .* bb([3 4 3 4]) + bb([1 2 1 2]);
bb_tot = round(minmax2wh(bb_tot));                                % keep same format as bb

% show rescaled image, input afg/bg masks, and bb_tot
% everything not in given areas is fed as unknown =
% don't learn color models and not clamped
if verbose > 1
    fid = figure;
    subplot(3,2,1);
    imagesc(im);
    axis equal; axis tight; axis off; drawnow;
    title('rescaled image');
    %
    subplot(3,2,2);
    imV = im;
    imV = HighlightImage(imV, abg_mask, 1, 255);                  % show abg in deep red
    imV = HighlightImage(imV, bg_mask, 1, 128);                   % show bg in light red
    imV = HighlightImage(imV, fg_mask, 2, 128);                   % show fg in light green
    imV = HighlightImage(imV, afg_mask, 2, 255);                  % show afg in deep green
    % unknown area is simply not highlit
    imagesc(imV);
    %
    DrawBB(wh2mxmxmymy(bb_tot), [0 0 1], 2);                      % show expanded bb
    DrawBB(wh2mxmxmymy(bb), [0 0 1], 2);                          % show detection bb
    axis equal; axis tight; axis off; drawnow;
    title(['frame ' num2str(fr)]);
end

% crop rescaled image to bb_tot
% if bg/afg masks don't fit bb_tot, they are cropped
bb_tot_stdscale_nocrop = bb_tot;
[imProc bb_tot] = CropImageBB(im, bb_tot);                        % crop also bb_tot
bb_tot_stdscale = bb_tot;
abg_maskProc = CropImageBB(abg_mask, bb_tot);
bg_maskProc = CropImageBB(bg_mask, bb_tot);
fg_maskProc = CropImageBB(fg_mask, bb_tot);
afg_maskProc = CropImageBB(afg_mask, bb_tot);

% show cropped image and afg/bg masks to be processed
% anything not in afg/bg is set to unknown (= erodable fg)
if verbose > 1
    subplot(3,2,3);
    imV = imProc;
    imV = HighlightImage(imV, abg_maskProc, 1, 255);              % show abg in deep red
    imV = HighlightImage(imV, bg_maskProc, 1, 128);               % show bg in light red
    imV = HighlightImage(imV, fg_maskProc, 2, 128);               % show fg in light green
    imV = HighlightImage(imV, afg_maskProc, 2, 255);              % show afg in deep green
    imagesc(imV);
    axis equal; axis tight; axis off; drawnow;
    title('overall areas');
end

% prepare input labels for grabcut
mask = zeros(size(bg_maskProc)) + 128;                            % everything not in given masks is unknown: learn no color models from it, and not clamped either
mask(abg_maskProc) = 64;                                          % abs background: learn col distr, and clamped
mask(bg_maskProc) = 96;                                           % background:     learn col distr, but not clamped = grabcut can change its label to fg
mask(fg_maskProc) = 160;                                          % foreground:     learn col distr, but not clamped = grabcut can change its label to bg
mask(afg_maskProc) = 192;                                         % absolute fg:    learn col distr, and clampled

% add border to final mask
sizMask = [size(mask,2) size(mask,1)];
for bix = 1:length(pars.border{class_id})/2
  bord_area = side2area(pars.border{class_id}{bix*2-1});
  bord_area = bord_area .* sizMask([1 2 1 2])';
  bord_mask = bbs2mask(bord_area, sizMask);
  area_type = pars.border{class_id}{bix*2};
  switch lower(area_type)
      case 'abg',  mask(bord_mask) = 64;
      case 'bg',   mask(bord_mask) = 96;
      case 'fg',   mask(bord_mask) = 160;
      case 'afg',  mask(bord_mask) = 192;
  end
end

% trafo bb_tot to original image (before rescaling);
% cropping is necessary, as large imscale + rounding can lead to out-of-image
% bb_tot = [xmin ymin width height]
bb_tot = minmax2wh(CropBB(wh2minmax(round(bb_tot * (1/imscale))), sizOrg));
