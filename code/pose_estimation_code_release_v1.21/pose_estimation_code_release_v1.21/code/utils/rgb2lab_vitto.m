function Ilab = rgb2lab_vitto(I)

% convenient wrapper around RGB2Lab
%

I = double(I);
[Ilab(:,:,1) Ilab(:,:,2) Ilab(:,:,3)] = RGB2Lab(I(:,:,1), I(:,:,2), I(:,:,3));
