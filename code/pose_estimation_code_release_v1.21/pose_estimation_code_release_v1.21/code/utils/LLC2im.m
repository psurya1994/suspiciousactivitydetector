function mask = LLC2im(im, LLC, bb)

% transform LLC areas from rel-to-bb coordinates to image coordinates of bb
% and produce LCC masks (includes cropping if necessary)
%

siz = [size(im,2) size(im,1)];                                    % [width height]
Nlimbs = length(LLC);
mask{Nlimbs} = [];
for lix = 1:Nlimbs
  if not(isempty(LLC{lix}))
    temp = TrafoBB2img(LLC{lix}, bb);
    mask{lix} = bbs2mask(temp, siz);
  end
end
