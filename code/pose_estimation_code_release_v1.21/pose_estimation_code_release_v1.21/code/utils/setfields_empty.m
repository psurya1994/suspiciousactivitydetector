function s = setfields_empty(s, fs, value)

% set unexisting fields of struct s to empty.
%
% Input:
% - s  : any struct
% - fs : cell arrary, with fs{fix}: name of a field (char)
%
% Output:
% s.(fs{fix}): empty if field fs{fix} didn't exist before, otherwise unchanged
%

if nargin < 3
  value = [];
end
  


for f = fs
  t = char(f);
  if ~isfield(s, t)
    s.(t) = value;
  end
end
