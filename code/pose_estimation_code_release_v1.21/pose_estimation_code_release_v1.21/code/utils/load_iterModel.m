function iterModel = load_iterModel(classname,pars,verbose)
% load iterModel data structure (encodes part models and part structures)
% if pars passed to the routing contain field customModelName then
% class specyfic customModel is loaded
% these mat-files are found in Matlab/ThirdParties/parse_matlab/util

if nargin < 3
  verbose = false;
end

if nargin < 2
  pars = [];
end

varName = 'model';
varNameAlt = 'iterModel';
class_id = class_name2id(classname);

if isfield(pars,'customModelName') && ~isempty(pars.customModelName{class_id})
  fileName = pars.customModelName{class_id};
else
  switch lower(classname)
  case 'full'
     fileName = 'peopleAllRunPostStoch';
  case 'ubf'
     %fileName = 'peopleUpperBody';   % new experiental model with longer, thinner arms, and wider torso
     %fileName = 'peopleUpperBodyCVPR08';    % model valid until Dec 2008, then changed because of .cont_pts y-flipped head 
     fileName = 'peopleUpperBodyPostCVPR09';
  otherwise
     error([filename ': unknown class ' classname]);
  end
end

try
  iterModel = load(fileName,varName);
  iterModel = iterModel.(varName);
catch
  iterModel = load(fileName,varNameAlt);
  iterModel = iterModel.(varNameAlt);
end


if verbose > 0
  disp([fileName ' model loaded...']);
end