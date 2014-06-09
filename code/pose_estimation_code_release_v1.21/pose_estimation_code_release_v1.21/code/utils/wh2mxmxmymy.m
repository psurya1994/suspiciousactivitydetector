function bb = wh2mxmxmymy(bbwh)

% convert bounding-box from [xmin ymin width height]
% to [xmin xmax ymin ymax] format
%
% assumes it's pixels, so '-1'
%

%bbmm = [bbwh(1) bbwh(1)+bbwh(3)-1  bbwh(2) bbwh(2)+bbwh(4)-1];
bb(4) = bbwh(2)+bbwh(4)-1;
bb(3) = bbwh(2);
bb(2) = bbwh(1)+bbwh(3)-1;
bb(1) = bbwh(1);