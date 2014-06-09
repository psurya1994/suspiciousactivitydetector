function im = PaintBB(im, bb, col, thick, nr, score, pos)

% Paint bounding-box bb = [xmin ymin width height]
% on image im, using color col.
%
% bb must be integer (pixel values) and its lines are cropped to fit im
%
% col should be a 3x1 vector, each element in [0,1]
%
% if thick(i) = offset given
% -> offsets BB to the list, generating different thicknesses
% >0 offsets -> outward, <0 offsets -> inward
% thick is a LIST of offsets ! if you want a 3-thick bb use:
% thick = [-1 0 1]
%

% process arguments
if nargin < 7
  pos = [3 3];
end
if nargin < 6
  score = [];
end
if nargin < 5
  nr = [];
end
if nargin < 4
  thick = 0;
end

text = '';
if ~isempty(nr) && ~ischar(nr)
  text = num2str(nr);
end
if ~isempty(score)
  text = [text ':' sprintf('%.2f',score)];
end

% draw all layers of bb
% !! BUG INTRODUCED BY MARCIN ! BROKE MY POSE SEARCH DRAWING ROUTINES !
% NOW FIXED IT, BUT MIGHT BREAK HIS ROUTINES
for t = thick
  bbt = [bb(1)-t bb(2)-t bb(3)+2*t bb(4)+2*t];
  im = PaintBB1(im, bbt, col);
end
if ~isempty(text)
  im = double(im)/255;
  im = PaintText(im,text,bb(1)+pos(1),bb(2)+pos(2)-1,[0 0 0],'letterssmall'); %text shadow
  im = PaintText(im,text,bb(1)+pos(1),bb(2)+pos(2)+1,[0 0 0],'letterssmall'); %text shadow
  im = PaintText(im,text,bb(1)+pos(1),bb(2)+pos(2),col,'letterssmall');
  im = uint8(im*255);
end

