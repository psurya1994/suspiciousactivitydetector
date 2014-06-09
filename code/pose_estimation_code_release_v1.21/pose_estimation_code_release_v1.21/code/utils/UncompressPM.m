function PM = UncompressPM(compPM)

% uncompress pose map compPM
%
% if compPM is already in uncompressed form,
% do nothing.
%

% is it already uncompressed ?
if isfield(compPM, 'a') && strcmp(class(compPM.a), 'double')
  PM = compPM;
  return;
elseif isfield(compPM, 'b') && strcmp(class(compPM.b), 'double')
  PM = compPM;
  return;
end

% uncompress
PM = compPM;
if isfield(PM, 'a')
  PM.a = double(PM.a)/255;
end
if isfield(PM, 'b')
  PM.b = double(PM.b)/255;
end
if isfield(PM,'mb')
  PM.mb = double(PM.mb)/255;
end