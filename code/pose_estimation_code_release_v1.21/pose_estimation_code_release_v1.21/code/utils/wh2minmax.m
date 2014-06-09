function bb = wh2minmax(bb)

% convert bounding-box from [xmin ymin width height] 
% to [xmin ymin xmax ymax] format
%
% assumes it's pixels, so '-1'
%

%bbmm = [bbwh(1:2) bbwh(1:2)+bbwh(3:4)-1];
bb(3:4) = bb(1:2)+bb(3:4)-1;
