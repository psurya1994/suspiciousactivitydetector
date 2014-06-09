function [fgP, imHist, unqCols] = TransferColorModel(cm, im)

% Transfers a previously learnt color model (cm) to a novel image (im).
%
% Input:
% - cm(:,end) = indeces of colors into RGB cube binned into 16x16x16 bins
% - cm(:,1:(end-1)) = color models for each limb
%                     (a column per limb type, each entry is prob that that color occurs on limb)
% - im = a novel image
% 

imHist = imvq16(im);                               % vector-quantize im pixel colors
[unqCols,dummy,imHist] = unique(imHist(:));        % unique colors in im
numTypes = size(cm,2)-1;                           % number of limb types

fgPu = zeros(16^3, numTypes);                      % unfold color model to direct indeces into quantized RGB cube
fgPu(uint16(cm(:,end)),:) = cm(:,1:(end-1));

fgP = fgPu;                                        % transfer to unique colors of im
fgP = fgP(unqCols,:);
