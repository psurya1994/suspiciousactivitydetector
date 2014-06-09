function compPM = CompressPM(PM)

% compress pose map PM to save storage space.
% 
% Minor loss in precision.
%
% if compPM is already in uncompressed form
% -> do nothing.
%


% is it already compressed ?
if strcmp(class(PM.a), 'uint8')
  compPM = PM;
  return;
end

% compress PM
compPM = PM;
compPM.a = uint8(PM.a*255);
compPM.b = uint8(PM.b*255);
if isfield(PM,'mb')
  compPM.mb = uint8(PM.mb*255);
end
