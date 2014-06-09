function expWWAll = expected_genmodel_FHedgesN(im, limb_mask, genmodel, part_ids, orient_ids, relor_inhibit)

% Takes edge image im.
% Assumes genmodel has a allfgP and colmask field.
%
% Every mention of 'limbs', refers to all of head, torso, upper/lower arms, legs.
% More generally, each of the 'types' in genmodel
%
% if limb_mask{p} given
% -> only search for limb lix in the area limb_mask{p}
% give limb_mask == false if no such constraints,
% or give limb_mask{p} = logical image when constr on p, or limb_mask{p} =
% [] when no constr on p
%
% limb ids in full body model:
% 1 = torso; 2-3 = upper arms; 4-5 = upper legs; 6-7 = lower arms; 8-9: lower legs; 10 = head
%
% limb ids in upper-body model:
% 1 = torso; 2-3 = upper arms; 4-5 = lower arms; 6 = head
%
% all fields of genmodel are documented in genmodel_fields.m
%


% shortcuts
[imy imx imz] = size(im);                             % must keep third argout, otherwise second is corrupted
numClass = max(genmodel.equiv_class);
numTypes = size(genmodel.dag,1);                      % number of body parts
if nargin < 4
  part_ids = 1:numTypes;
end
imo = size(genmodel.orient,2);                        % number of possible limb orientations
if nargin < 5
  orient_ids = ones(imo,max(part_ids));               % which orientations to consider orient_ids(o, pix) == 1 -> consider orient o for part pix
end
if nargin < 6
  relor_inhibit = false;
end
respIm = zeros([imy imx imo numTypes]);               % pixel-to-limb responsibility maps
assert(numClass == numTypes);
%lens = genmodel.len*2; wids = genmodel.wid*2;        % unused ?
LOG_MAX = 500;
LOG_MIN = -500;

% OLD LLC - WRONG ONE !
% limb location constraints
% generate numTypes images im{p},
% with only the masked pixels in them
% -> each limb has its own edge image where to be found !
%if not(islogical(limb_mask))
%  orgim = im;  clear im;
%  for p = 1:numTypes
%    im{p} = MaskImage(orgim, limb_mask{p});
%  end
%else
%  orgim = im; clear im;
%  for p = 1:numTypes
%    im{p} = orgim;
%  end
%end

% new llc - right one
if islogical(limb_mask)
  clear limb_mask;
  for p = 1:numTypes
    limb_mask{p} = [];           % all locations are valid
  end
end

% loop over limbs (num_Types)
% estimate unary prob that a limb appears at a certain loc and orient
% based only on the edge model (pre-learnt)
% there is no search over scales
% output is respIm(x,y,orient,limb) and pNN (normalization constant)
pNN = zeros(numTypes,1);
for p = part_ids
    ww = genmodel.ww(:,:,p);
    % resize mask to 2*len x  2*2*wid (in width mask was trained to capture context around the limb)
    mask = imresize(ww,2*[genmodel.len(p) 2*genmodel.wid(p)]); % originally 'nearest'    
    if any(mask(:)),
      mask = mask*(sum(abs(ww(:)))/sum(abs(mask(:))));
    end  

    resp = getSegmentsEdge(im, mask, find(orient_ids(:,p))');
    resp = partshiftZ0(resp,genmodel.len(p),0);  % Shift so rectangles are anchored at the top
    resp(resp > LOG_MAX) = LOG_MAX; resp(resp < LOG_MIN) = LOG_MIN;
    resp = exp(resp);
    % location constraint; set to exactly 0 prob out of llc mask (and not to some not-so-small value as other places in image with no evidence)
    % careful with partshift !! must set it exactly as follows, don't call
    % with just one plane !
    if not(isempty(limb_mask{p}))
      temp = partshiftZ0(repmat(limb_mask{p},[1 1 imo]),genmodel.len(p),0);
      resp(~temp) = 1;    % not set-hard-0 (when limb_mask == all ones -> changes nothing)
    end
    NN = sum(resp(:));
    pNN(p) = log(NN);
    respIm(:,:,:,p) = resp/NN;
end
%figure; imagesc(respIm(:,:,1,1)); axis equal; % VF: show prob for torso at first orientation

% Convert orient to exp scale (and normalize)
genmodel.orient = exp(genmodel.orient);
genmodel.orient = genmodel.orient./repmat(sum(genmodel.orient,2),1,imo);

% play around with rel-orient prior
if relor_inhibit
  inhibit = ones(1,imo);
  inhibit(1:3) = [0.01 0.1 0.5];
  inhibit((imo-1):imo) = [0.5 0.1];
  %genmodel.orient([4 5],[1:12 14:24]) = 0; % must keep same orient (i.e. won't superimpose)
  %genmodel.orient([4 5], 2:24) = 0;        % forces to flip orient (i.e. must superimpose)
  if numTypes == 6 % upper body model
    relor_pidx = [4 5];
  elseif numTypes == 10 % full body model
    relor_pidx = [6 7 8 9];
  else
    error('relor inhibit undefined for loaded PS model');
  end
  for p=relor_pidx
    genmodel.orient(p,:) = genmodel.orient(p,:) .* inhibit;  % smoothly recommends not to flip orient -> rarely will superimpose
    genmodel.orient(p,:) = genmodel.orient(p,:) / sum(genmodel.orient(p,:));
  end
end

% run message-passing optimization
respIm = quick_msgpass_tree_corr(respIm, genmodel, true);

% % exclude overlap between upper and lower arms of same side
% if false
%     respIm(:,:,:,5) = respIm(:,:,:,5) .* (respIm(:,:,:,3)<10^-7);
%     respIm(:,:,:,5) = respIm(:,:,:,5) .* (respIm(:,:,:,5) > 10^-7); % remove its own near-0 -> avoid unstable normalizations
%     respIm(:,:,:,5) = respIm(:,:,:,5) / (sum(sum(sum(respIm(:,:,:,5)))+eps));
%     %
%     respIm(:,:,:,4) = respIm(:,:,:,4) .* (respIm(:,:,:,2)<10^-7);
%     respIm(:,:,:,4) = respIm(:,:,:,4) / (sum(sum(sum(respIm(:,:,:,4))))+eps);
%     respIm(:,:,:,4) = respIm(:,:,:,4) .* (respIm(:,:,:,4) > 10^-7);
% end

% generate per-pixel posterior prob maps for each body part
partIms = gen_pp_post(respIm, genmodel);

% fuse part-specific maps into a single colorful one
colIm = gen_pp_overall(partIms);

% save output
clear expWWAll;
expWWAll.a = colIm;
expWWAll.b = partIms;
expWWAll.NN = NN;

% VF: total pose entropy helps evaluating the confidence in the
% estimated poses
expWWAll.e = TotalPoseEntropy(respIm);
% VF: total pixel confidence, evaluating spatial separateness of the
% limbs and their confidences
expWWAll.p = TotalPixelConfidence(partIms);

