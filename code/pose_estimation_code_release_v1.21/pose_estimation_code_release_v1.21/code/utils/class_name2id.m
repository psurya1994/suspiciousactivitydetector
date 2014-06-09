function id = class_name2id(name)
% id = class_name2id(name)
%
% object class name-to-id convertor
%
% names are case-insensitive
%
% See also class_id2name
%

switch lower(name)
    case 'ubp'
        id = 1;
    case 'ubf'
        id = 2;
    case 'full'
        id = 3;
    case 'ubfreg' %if it is a ubf regressed from face detection (treat the same as normal ubf)
        id = 2;
    case 'calvinubf' %if it is obtained using Felzenszwalb model (pff) retrainded with the upper body frontal deta from Marin
                  %size is exactly the same as in original marin's detector so we can treat it the same
        id = 2;
    case 'ubfpff' %if it is obtained using Felzenszwalb model (pff) retrainded with the upper body frontal deta from Marin
                  %size is exactly the same as in original marin's detector so we can treat it the same
        id = 2;
    case 'ubpffvoc2008' %if it is obtained using Felzenszwalb model (pff) trained on the voc2008 data (component 1- upper body)
                        % different size -> we need a new class
        id = 4;
    otherwise
        error([mfilename ': unknown class name ' name]);
end
