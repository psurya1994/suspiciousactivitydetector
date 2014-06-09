function fgP = applyCMTransfer(fg,bg,transfer_weights)

if nargin < 3
  pred = fg;
else
  pred = fg*transfer_weights;
end
fgP = pred./(pred+bg);
fgP(isnan(fgP(:))) = 0;