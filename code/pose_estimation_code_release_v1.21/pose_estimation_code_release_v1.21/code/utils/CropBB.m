function bbc = CropBB(bb, siz)

% crop BB inside image with size siz = [width height]
%
% input bb = [xmin ymin xmax ymax]
%

xmin = max(bb(1),1);
xmax = min(bb(3), siz(1));
ymin = max(bb(2),1);
ymax = min(bb(4), siz(2));

bbc = [xmin ymin xmax ymax];
