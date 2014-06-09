function im = DrawLineOnImgCol(im, y1, x1, y2, x2, col, thick)

% col = [r g b]
%

% draw line
im = DrawLineOnImg3channels(im, y1, x1, y2, x2, col);

% draw second thickness
if thick
  nx1 = x1+1; nx2=x2+1;
  im = DrawLineOnImg3channels(im, y1, nx1, y2, nx2, col);
  %
  nx1 = x1-1; nx2=x2-1;
  im = DrawLineOnImg3channels(im, y1, nx1, y2, nx2, col);
  %
  ny1 = y1+1; ny2=y2+1;
  im = DrawLineOnImg3channels(im, ny1, x1, ny2, x2, col);
  %
  ny1 = y1-1; ny2=y2-1;
  im = DrawLineOnImg3channels(im, ny1, x1, ny2, x2, col);
end
