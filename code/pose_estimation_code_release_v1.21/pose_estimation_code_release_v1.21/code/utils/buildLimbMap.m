function map =  buildLimbMap(PM,model,lp_params)
  % converts the part specific segmentations (obtained from parsing routine PM.b, 
  % into limb specyfic segmentations if necessary
  
  types = size(model.colmask,2);
  classes = size(model.colmask,1);
  
  if size(lp_params.LP.map,3) == types % working with limbs
    map = respb;
  elseif size(lp_params.LP.map,3) == classes % working with limb classes
    map = condenseLRResp(PM.b,model.colmask);
    normfact = sum(model.colmask,2);
    for p = 1:classes, % additional normalization is required
      map(:,:,p) = map(:,:,p)/normfact(p);
    end
   else
      error('location prior do not correspond to number of limbs nor to number of limb classes');
   end