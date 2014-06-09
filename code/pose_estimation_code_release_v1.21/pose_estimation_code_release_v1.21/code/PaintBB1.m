function im = PaintBB1(im, bb, col)

% Paint bounding-box bb = [xmin ymin width height]
% on image im, using color col.
%
% bb must be integer (pixel values) and its lines are cropped to fit im
%
% col should be a 3x1 vector, each element in [0,1]
%
%
%

% from [0,1] to [1 255]

if max(im(:))>1
    col = round(col*255);
end

% image width and height
iw = size(im, 2);
ih = size(im, 1);

% top-left -> top-right
l = CropLine([bb(1) bb(2) bb(1)+bb(3)-1 bb(2)], [iw ih]);  % l = [xmin ymin xmax ymax]
if not(isempty(l))
  im(l(2):l(4), l(1):l(3), 1) = col(1);
  im(l(2):l(4), l(1):l(3), 2) = col(2);
  im(l(2):l(4), l(1):l(3), 3) = col(3);
end

% bottom-left -> bottom-right
l = CropLine([bb(1) bb(2)+bb(4)-1 bb(1)+bb(3)-1 bb(2)+bb(4)-1], [iw ih]);  % l = [xmin ymin xmax ymax]
if not(isempty(l))
  im(l(2):l(4), l(1):l(3), 1) = col(1);
  im(l(2):l(4), l(1):l(3), 2) = col(2);
  im(l(2):l(4), l(1):l(3), 3) = col(3);
end

% top-left -> bottom-left
l = CropLine([bb(1) bb(2) bb(1) bb(2)+bb(4)-1], [iw ih]);  % l = [xmin ymin xmax ymax]
if not(isempty(l))
  im(l(2):l(4), l(1):l(3), 1) = col(1);
  im(l(2):l(4), l(1):l(3), 2) = col(2);
  im(l(2):l(4), l(1):l(3), 3) = col(3);
end

% top-right -> bottom-right
l = CropLine([bb(1)+bb(3)-1 bb(2) bb(1)+bb(3)-1 bb(2)+bb(4)-1], [iw ih]);  % l = [xmin ymin xmax ymax]
if not(isempty(l))
  im(l(2):l(4), l(1):l(3), 1) = col(1);
  im(l(2):l(4), l(1):l(3), 2) = col(2);
  im(l(2):l(4), l(1):l(3), 3) = col(3);
end
