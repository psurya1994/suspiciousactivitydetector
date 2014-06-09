function name = class_id2name(id)
% name = class_id2name(id)
% object class id2name convertor
%
% See also class_name2id
%

switch id
    case 1
        name = 'ubp';
    case 2
        name ='ubf';
    case 3
        name ='full';
    otherwise
        error([mfilename ': unknown class id ' num2str(id)]);
end
