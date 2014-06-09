function trg = copy_fields(src, trg, fs, n)

% copy data in src.(fs{fix}){n} to trg.(fs{fix}),
% for all fs{fix}
%
% Input:
% - src,trg  : two structures
% - fs   : cell arrary, with fs{fix}: name of a field (char);
%          src.(fs{fix}) must be a cell array !
% - n    : index into src.(fs) (integer)
%
% Output:
% trg.(fs{fix}) = src.(fs{fix}){n}
%                 (of course this notation is for me, not valid in matlab)
%

for f = fs
  t = char(f);
  data = src.(t);
  trg.(t) = data{n};
end
