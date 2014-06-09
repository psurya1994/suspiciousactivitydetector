function D = rmfield_safe(D, field_name)

% safe rmfield: does not crash if
% field_name is not a field of D
% (it just returns D unchanged in that case)
%

if isfield(D, field_name)
  D = rmfield(D, field_name);
end
