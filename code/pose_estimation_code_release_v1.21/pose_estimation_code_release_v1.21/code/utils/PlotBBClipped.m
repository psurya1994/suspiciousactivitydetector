function PlotBBClipped(BB, col, thick, limits)

% BB in form:
% [minx maxx miny maxy]
%
% limits = [width height]
%

% draw left line
if BB(1) > 0 && BB(1) <= limits(1) && BB(3) <= limits(2) && BB(4) > 0
  miny = max(BB(3),1); maxy = min(BB(4),limits(2));
  plot([BB(1) BB(1)], [miny maxy], 'Color', col, 'LineWidth', thick);
end

% draw right line
if BB(2) > 0 && BB(2) <= limits(1) && BB(3) <= limits(2) && BB(4) > 0
  miny = max(BB(3),1); maxy = min(BB(4),limits(2));
  plot([BB(2) BB(2)], [miny maxy], 'Color', col, 'LineWidth', thick);
end

% draw top line
if BB(3) > 0 && BB(3) <= limits(2) && BB(1) <= limits(1) && BB(2) > 0
  minx = max(BB(1),1); maxx = min(BB(2),limits(1));
  plot([minx maxx], [BB(3) BB(3)], 'Color', col, 'LineWidth', thick);
end

% draw bottom line
if BB(4) > 0 && BB(4) <= limits(2) && BB(1) <= limits(1) && BB(2) > 0
  minx = max(BB(1),1); maxx = min(BB(2),limits(1));
  plot([minx maxx], [BB(4) BB(4)], 'Color', col, 'LineWidth', thick);
end
