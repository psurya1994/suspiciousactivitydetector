function coor = convertSticksToImgCoor(sticks,frameSize,enl_bb)
% sticks - [4x6] array with sticks endpoints in some frame coordinates
% frameSize - [width height]
% enl_bb - position of the frame wrt to the image [x y width height]
  sticks = double(sticks);
  x1s = ((sticks(1,:)-1)/frameSize(1))*enl_bb(3) + enl_bb(1);
  x2s = ((sticks(3,:)-1)/frameSize(1))*enl_bb(3) + enl_bb(1);
  y1s = ((sticks(2,:)-1)/frameSize(2))*enl_bb(4) + enl_bb(2);
  y2s = ((sticks(4,:)-1)/frameSize(2))*enl_bb(4) + enl_bb(2);
  coor = [x1s; y1s; x2s; y2s];
end