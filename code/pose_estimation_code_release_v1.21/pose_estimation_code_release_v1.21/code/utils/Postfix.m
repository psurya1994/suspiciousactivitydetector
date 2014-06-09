function [so temp] = Postfix(s, t)

% substrig of s succeeding last
% occourrence of t in s
%
% t can be any string
%

% indexes of occurrences of t in s
ocix = strfind(s, t);
temp = [];
if isempty(ocix)
  so = [];
  return;
elseif ocix(end)+length(t) == length(s)+1
  so = [];
  return;
else
  so = s( (ocix(end)+length(t)) : end);
  temp = s(1:ocix(end)-1);
end
