function [cx, cy] = ClipPoint(x, y, siz)

% COMMENT TO BE WRITTEN
%

cx = min(max(x, siz(1)),siz(3));
cy = min(max(y, siz(2)),siz(4));
