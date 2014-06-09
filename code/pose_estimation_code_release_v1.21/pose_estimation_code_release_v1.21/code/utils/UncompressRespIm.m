function respIm = UncompressRespIm(comp)

% COMMENT ME !
%

% if already uncompressed
if not(isstruct(comp))
  respIm = comp;
  return;
end

respIm = zeros(comp.siz);
respIm(comp.ixs) = comp.vals;
