function bbc = CropBBwh(bb, siz)

% crop BB inside image with size siz = [width height]
%
% input bb = [xmin ymin width height]
%

xmin = max(bb(1),1);
xmax = min(bb(1)+bb(3)-1, siz(1));
ymin = max(bb(2),1);
ymax = min(bb(2)+bb(4)-1, siz(2));

bbc = [xmin ymin xmax-xmin+1 ymax-ymin+1];
