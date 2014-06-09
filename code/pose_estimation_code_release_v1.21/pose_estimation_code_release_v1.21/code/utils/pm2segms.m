function out = pm2segms(pm, class_id, pars, verbose)

% Find segmentation for each part pm.b(:,:,p).
%
% Respects optional non-overlap preferences pars.non_overlap:
% pars.non_overlap(:,k) = [pi pj]' -> prefers part pj to not overlap with part pi
% these soft constraints are limited to a tree-fashion: each pj can only
% non-overlap with at most one pi. The order of part selection must be
% specified in pars.order (although it's in theory derivable automatically).
%
% Output:
% out(:,:,p) = binary segmentation mask for part p
%

% robustly find candidate segmentations for each limb (= regions)
% based on watershedding the PMs for each limb, after very low thresh.
% + find all possible limb alternatives
% + no risk of mergers
% + no under-estimation of limb size
%
% ++ get more than rectangles: full limb outlines
% ++ same-side upper/lower arm exclusion constraints
%    by preferring a non-overlapping alternative, if one is available,
%    while allowing for overlapping solutions when no alternative available.
% ++ same-side upper/lower arm connectedness constraints
%    by preferring an alterative for the lower arm close to its upper arm
%

global colors;


% initialization
pm.b = pm.b .* (pm.b > pars.abs_min);      % extremely low threshold: make sure we find every possible hint of a limb ;)
b = pm.b;
[height width numTypes] = size(b);
segm(numTypes).all = [];
segm(numTypes).regs = [];
% show pose-map after thresholding
if verbose > 1
  ShowParseResult(pm, false, true);
  set(gcf,'name','pose maps (per-pixel probs for each part)');
end

% water-shed based segmentation:
% find limb alternatives and rank them by their maximum pixel
for p = 1:numTypes
  im = b(:,:,p);                           % per-pixel prob map (each entry in [0,1]) for pth body part
  if isnan(im(1,1)) || all(im(:)==0)       % no body part
    segm(p).regs = [];
    segm(p).all = zeros(height,width);
    segm(p).scores = [];
    continue;
  end
  if pars.smooth_sigma > 0                 % smooth out minor ripples that cause watershed to break contiguous regions
    im = vgg_gsmooth(im, pars.smooth_sigma, 'same');
  end
  ws = watershed(1-im);                    % watershed
  numRegs = max(max(ws));                  % number of regions
  segm(p).regs = false(height,width);      % candidate regions for each limb
  for rix = 1:numRegs
    imr = im .* (ws == rix);               % region only
    mx = max(max(imr));                    % maximum
    imrt = (imr > min(pars.rel_to_max*mx,pars.segm_min)); 
    segm(p).regs(:,:,rix) = imrt;          % store results
    segm(p).scores(rix) = mx;
  end
  % sort regions by their score
  [trash order] = sort(-segm(p).scores);
  segm(p).regs = segm(p).regs(:,:,order);
  segm(p).scores = segm(p).scores(order);
  % compose results on .all
  segm(p).all = zeros(height,width);
  for rix = 1:numRegs
    segm(p).all = segm(p).all + segm(p).regs(:,:,rix)*rix;
  end
end % loop over parts

% show candidate regions
if verbose > 1
    figure;
    nsp = ceil(sqrt(numTypes));
    for p = 1:numTypes
        subplot(nsp,nsp,p);
        image(segm(p).all+1);                            % +1 -> direct entries into 'colors', otherwise duplicates 0 and 1 to 1
        colormap([0 0 0; colors]);
        axis equal; axis tight; axis off;
    end
    set(gcf,'name','candidate regions');
end

% compose final output
%%%%%%%%%%%%%%%%%%%%%%
% determine part order from non-overlappiness preferences
out = false(height,width,numTypes);                      % final output segmentation matrix
for p = pars.order{class_id}
  nopix = find(pars.non_overlap{class_id}(2,:)==p);                % part preferred not to overlap with
  nop = pars.non_overlap{class_id}(1,nopix);
  %
  clpix = find(pars.close{class_id}(2,:)==p); 
  clp = pars.close{class_id}(1,clpix);                             % part preferred to be close to
  scores = zeros(1,length(segm(p).scores));
  for rix = 1:length(segm(p).scores)
    if ~isempty(nop)
      nisect = segm(p).regs(:,:,rix) .* ~out(:,:,nop);   % non-intersection with nop 
    else
      nisect = segm(p).regs(:,:,rix);
    end
    %
    if ~isempty(clp)
      [ys1 xs1 trash] = find(out(:,:,clp)); P1 = [xs1 ys1]';
      [ys2 xs2 trash] = find(segm(p).regs(:,:,rix)); P2 = [xs2 ys2]';
      dist = sqrt(min(min(vgg_nearest_neighbour_dist(P1, P2))));
      dist_err = max(dist-pars.close_thresh,0);
    else
      dist_err = 0;
    end
    %
    if ~isempty(dist_err)
      scores(rix) = sum(sum(nisect .* b(:,:,p))) * exp(-(dist_err)/pars.close_sigma); % score for this region
    else
      scores(rix) = 0;                                   % in the rare case the segmentation for this limb is entirely empty
    end
  end
  %
  [trash wix] = max(scores);                             % winner region
  if isempty(wix)                                        % no body part
    continue;
  end
  if pars.chop && ~isempty(nop)
    % chops away chunk overlapping with forbidden area (non-overlap part)
    out(:,:,p) = segm(p).regs(:,:,wix) .* ~out(:,:,nop);
  else
    out(:,:,p) = segm(p).regs(:,:,wix);
  end
  %
  % keep only the largest conn comp
  % (sometimes happens to have two in nasty cases of chopping)
  out(:,:,p) = LargestCC(out(:,:,p));
  %
  % remove cusps and smooth out in general
  out(:,:,p) = imclose(imopen(out(:,:,p),ones(3,3)),ones(3,3));
end

% show final result
if verbose > 1
    figure;
    for p = 1:numTypes
        subplot(nsp,nsp,p);
        imshow(out(:,:,p));
        axis equal; axis tight; axis off;
    end
    set(gcf,'name','output segmentation');
    keyboard;
    close all;
end
