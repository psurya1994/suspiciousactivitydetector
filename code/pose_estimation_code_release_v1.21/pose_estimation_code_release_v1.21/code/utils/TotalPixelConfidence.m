function tot_c = TotalPixelConfidence(partIms)

% total confidence over all pixels and all limbs partIms
% the per-pixel probs partIms(x,y,limb) are assumed normalized
% (i,e, range in [0,1])
%

maxperp = max(partIms, [], 3);
tot_c = sum(sum(maxperp));
