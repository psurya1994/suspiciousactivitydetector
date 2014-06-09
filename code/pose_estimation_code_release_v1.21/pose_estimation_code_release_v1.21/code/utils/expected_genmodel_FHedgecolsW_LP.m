function [expWWW samples] = expected_genmodel_FHedgecolsW_LP(seg, m, limb_mask, genmodel, img_lik_only, part_ids, orient_ids, relor_inhibit,motion,lastiter)
% 26/01/09 this file is changed comparing to the oryginal expected_genmodel_FHedgecolsW,
% takes segmentation image (seg) of the same size as edgemap image (m)
% the number of channels (size(seg,3)) must correspond to number of parts or number of part classes
% doesn't use genmodel.fgP field anymore as seg image is passed.
% Assumes genmodel has colmask field
%
% ALSO assumes genmodel has a .www field, which is a scalar that weights the
% .www masks; this is the only thing that is tweaked!
%
% Cleaned up and extended by Vittorio Ferrari (VF)
%
% if limb_mask{p} given
% -> only search for limb lix in the area limb_mask{p}
% give limb_mask == false if no such constraints,
% or give limb_mask{p} = logical image when constr on p, or limb_mask{p} =
% [] when no constr on p
%
% output in expWWW:
% .a = pose map mixed over all limb types
% .b = pose maps for each limb type
% .e, .p = total pose entropy, total pixel confidence (measures of the quality of the result)
% .respIm = compressed version of respIm, with subfields:
%           .vals = entries of original respIm >10^-5, in 'single' format.
%           .ixs = corresponding indeces (in uint32 format)
%           .siz = original size of respIm
%
% Every mention of 'limbs', refers to all of head, torso, upper/lower arms, legs.
% More generally, each of the 'types' in genmodel
%
% part ids in full-body model:
% 1 = torso; 2-3 = upper arms; 4-5 = upper legs; 6-7 = lower arms; 8-9: lower legs; 10 = head
% (for symmetric parts: the one that's left in the image comes first)
%
% part ids in upper-body model:
% 1 = torso; 2-3 = upper arms; 4-5 = lower arms; 6 = head
%
% all fields of genmodel are documented in genmodel_fields.m
%
% if img_lik_only == 1
% -> don't do full parsing, stop at extracting the
%    image likelihood and return it in expWWW.respIm (uncompressed)
%    in this case can choose which parts in part_ids
%    returned likelihoods are shifted so that expWWW.respIm aligns with middle-top of segment
% if img_lik_only == 2
%    same as above, but likelihoods aligned with center of segment
%

% process arguments


% shortcuts
[imy imx imc] = size(seg);                             % must keep third argout, otherwise second is corrupted
numClass = max(genmodel.equiv_class);
numTypes = size(genmodel.dag,1);                      % number of body parts

if nargin < 5
  img_lik_only = false;
end
if nargin < 6
  part_ids = 1:numTypes;
end
imo = size(genmodel.orient,2);                        % number of possible limb orientations
if nargin < 7
  orient_ids = ones(imo,max(part_ids));               % which orientations to consider orient_ids(o, pix) == 1 -> consider orient o for part pix
end
if nargin < 8
  relor_inhibit = false;
end
if nargin < 9
  motion = [];
end
if nargin < 10
  lastiter = false;
end
respIm = zeros([imy imx imo numTypes]);               % pixel-to-limb responsibility maps
assert(numClass == numTypes);
%lens = genmodel.len*2; wids = genmodel.wid*2;        % unused ?
LOG_MAX = 500;
LOG_MIN = -500;

samples = [];

% Get out color histogram
% this is how it can interpret genmodel.fgP
% (the appearance model for each limb)
% wasteful: imHist might have been computed already before invoking this function
% (but negligible cost compared to estimations that follow)
%imHist = imvq16(im);
% unique colors in imHist must match the order in genmodel.fgP
%[dummy1,dummy2,imHist] = unique(imHist);

% OLD LLC - WRONG ONE !
% ONLY WORKING ON EDGEMAP, NOT ON COLORS,
% AND ONLY BLACKENING EDGEMAP, NOT FORBIDDING LOCATIONS IN STATE SPACE !
% REALLY BAD !
% limb location constraints
% generate numTypes images im{p},
% with only the masked pixels in them
% -> each limb has its own edge image where to be found !
%if not(islogical(limb_mask))
%  orgm = m; clear m;
%  for p = part_ids
%     m{p} = MaskImage(orgm, limb_mask{p});
%  end
%else
%  orgm = m; clear m;
%  for p = part_ids
%     m{p} = orgm;
%  end
%end

% new llc - right one
if islogical(limb_mask)
  clear limb_mask;
  for p = 1:numTypes
    limb_mask{p} = [];           % all locations are valid
  end
end

% all the stages below aim at progressively improving respIm

% loop over limb types (num_Types)
% estimate unary prob that a limb appears at a certain loc and orient
% based on the edge model (pre-learnt) and limb color information
% (from first parse, in genmodel.fgP)
% there is no search over scales
% output is respIm(x,y,orient,limb_type) and pNN (normalization constant)
pNN = zeros(numTypes,1);
%tic

[edgeResp colorResp] = deal(zeros([imy imx imo numTypes]));

   
for p = part_ids
    if numTypes ~= imc
      %working with limb classes location priors
      currLp = find(genmodel.colmask(:,p));
    else
      %working with limb lication priors
      currLp = p;
    end
      
    % take the edge masks
    
    ww = genmodel.ww(:,:,p);
    % resize mask to 2*len x  2*2*wid (in width mask was trained to capture context around the limb)
    mask = imresize(ww,2*[genmodel.len(p) 2*genmodel.wid(p)]);
    if any(mask(:)),
      mask = mask*(sum(abs(ww(:)))/sum(abs(mask(:))));
    end
    currOrient = find(orient_ids(:,p))';
    if isfield(genmodel,'unaryAndriluka')
      for i=currOrient
        edgeResp(:,:,i,p) = genmodel.unaryAndriluka(:,:,i,p);
      end
    else
      edgeResp(:,:,:,p) = getSegmentsEdge(m, mask, currOrient);
    end
    resp = edgeResp(:,:,:,p) * genmodel.cues_weight.edge(p);
    % Treat this as a feature that is weighted by genmodel.www
    resp = resp*genmodel.www(p);
    
    
    % take the color masks
    ww = genmodel.wc(:,:,p);
    mask = imresize(ww,2*[genmodel.len(p) 2*genmodel.wid(p)]);
    if any(mask(:)),
      mask = mask*(sum(abs(ww(:)))/sum(abs(mask(:))));
    end
    %fgP = reshape(genmodel.fgP(imHist,find(genmodel.colmask(:,p))),imy,imx);
    colorResp(:,:,:,p) = getSegmentsEdge(seg(:,:,currLp), mask, currOrient);    
    resp = resp + colorResp(:,:,:,p) * genmodel.cues_weight.color(p);
    
    resp = partshiftZ0(resp,genmodel.len(p),0);        % Shift so rectangles are anchored at the middle-top
    resp(resp > LOG_MAX) = LOG_MAX; resp(resp < LOG_MIN) = LOG_MIN;
    resp = exp(resp); 
    % location constraint; set to exactly 0 prob out of llc mask (and not to some not-so-small value as other places in image with no evidence)
    % careful with partshift !! must set it exactly as follows, don't call
    % with just one plane !
    if not(isempty(limb_mask{p}))
      temp = partshiftZ0(repmat(limb_mask{p},[1 1 imo]),genmodel.len(p),0);
      resp(~temp) = 1e-10;    % not set-hard-0 (when limb_mask == all ones -> changes nothing) % originally 1 but when andriluka used it must be smaller
    end
    NN = sum(resp(:));
    pNN(p) = log(NN);  % not used !
    respIm(:,:,:,p) = resp/NN;
    maxrespIm(p) = max(resp(:)/NN);
end

% if user wants only the image likelihood, return it and quit
if img_lik_only && lastiter
  if img_lik_only == 2
    % shift respIm to align with center of segment
    expWWW.respIm = align_respIm(respIm, genmodel);
  else
    expWWW.respIm = respIm;
  end
  
  return;
end

% Convert orient to normalized probability distribution
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

% 
if lastiter
%     tic
    [MAP samples trash] = sampleFromRespIm(respIm,genmodel); % get best pose from posterior marginals
%     toc
end

[respIm] = quick_msgpass_tree_corr(respIm, genmodel, true);
% now respIm is back with middle-middle anchor
 
% compress and output respIm
% it's a large matrix, but allows to sample puppets later on
temp = CompressRespIm(respIm, 10^-5);        % only keep elements larger than 10^-5
expWWW.respIm = temp;                        % total compression is very large !

% generate per-pixel posterior prob maps for each body part
%keyboard; % fake setting a single best on max (so can see body part lengths)
if false
siz = size(respIm);
for p = part_ids
  temp = respIm(:,:,:,p);
  temp = temp(:);
  [trash m] = max(temp);
  temp = zeros(size(temp));
  temp(m) = 1;
  respIm(:,:,:,p) = reshape(temp, siz(1:3));
end
end % if false



partIms = gen_pp_post(respIm, genmodel);

% fuse part-specific maps into a single colorful one
colIm = gen_pp_overall(partIms);

expWWW.a = colIm;
expWWW.b = partIms;
expWWW.pNN = pNN;
expWWW.e = TotalPoseEntropy(respIm);       % confidence in the estimated poses
expWWW.p = TotalPixelConfidence(partIms);  % spatial separateness of the limbs and their confidences
% save out all fields
if lastiter
  expWWW.MAP = MAP;
end


end