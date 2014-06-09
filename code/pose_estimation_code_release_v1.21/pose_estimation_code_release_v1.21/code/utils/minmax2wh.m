function bb = minmax2wh(bb)

% convert bounding-box from [xmin ymin xmax ymax]
% to [xmin ymin width height] format
%
% assumes it's pixels, so '+1'
%

%bbwh = [bbmm(1:2) bbmm(3:4)-bbmm(1:2)+1];
bb(3:4) = bb(3:4)-bb(1:2)+1;
