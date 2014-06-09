function comp = CompressRespIm(respIm, min_val)

% COMMENT TO BE WRITTEN
%

comp.ixs = uint32(find(respIm>min_val));
comp.vals = single(respIm(comp.ixs));
comp.siz = size(respIm);
