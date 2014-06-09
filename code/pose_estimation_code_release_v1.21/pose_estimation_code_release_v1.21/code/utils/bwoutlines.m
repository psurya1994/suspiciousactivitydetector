function out = bwoutlines(bwi, thick)

% easy ;)
%

outi = imerode(bwi,ones(thick));
out = bwi - outi;           
