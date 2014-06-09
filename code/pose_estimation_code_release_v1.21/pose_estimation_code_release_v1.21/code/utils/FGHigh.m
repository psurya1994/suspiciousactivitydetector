function Aout = FGHigh(STD, STM, BBTSC1, BBTC, pars, fids, verbose)

% invoke grabcut for spatio-temporal chunk STD, and initial labels mask STM
% then crop and rescale A to the original BBs BBA and BBT
%
% Input:
% - STD(y,x,col_chan,frm) = 3D image data, with frm going from 1 to dix_end-dix_begin+1
% - STM(y,x,frm) = 3D label map for initializing grabcut
% - BBTSC1(1:4,frm) = coordinates of cropped bb in the standard-scaled coordinates of STD/STM
% - BBTC(1:4,frm) = coordinates of cropped bb in the original image; all bbs in [min_x max_x width height] format
% - fids(frm) = figure with preparations for frm (if verbose > 1)
% - verbose>1 -> expects there is a figure for each frame in fids(frm), and adds a segmentation each in subplot(3,2,4)
%
% Output:
% - Aout{frm} = binary segmentation matrix
%   if no STM is foreground -> set all output to bg
%

% process arguments
if nargin < 7
  verbose = false;
end

% shortcut
F = size(STD,4);                                          % number of frames

% segment STD with grab-cut
opts.fgModelTimeWindow = F-(1-mod(F,2));                  % foreground color models are global 
opts.bgModelTimeWindow = min(pars.bg_model_temp_win,F);   % background model adapts over time
opts.bgModelTimeWindow = opts.bgModelTimeWindow - (1-mod(opts.bgModelTimeWindow,2));
disp('Running grabcut ...');
if any(any(STM == 160 | STM == 192)) && any(any(STM == 64 | STM == 96)) % fg and bg regions to build cm present
  tic;
  A = vg_grabCut3D(STD, STM, 'single', opts);
  toc;
  A = logical(A);                                           % don't want it in double
elseif any(any(STM == 160 | STM == 192)) % no background pixels 
  disp('fg pixels present but no bg pixels -> setting everything to fg');
  A = false(size(STM)); % put all fg or undefinded regions to fg
  A(STM > 0) = true;
else % no foreground pixels
  disp('no fg pixels -> setting everything to bg');
  A = false(size(STM));                                     % if no STM is fg -> set all output to bg
end

% show output segmentations
if verbose > 1
    for frm = 1:F
      figure(fids(frm));
      subplot(3,2,4);
      imShow = HighlightImage(STD(:,:,:,frm), A(:,:,frm), 2);
      imagesc(imShow);
      axis equal; axis tight; axis off; drawnow;
      title('overall foreground');
    end
    keyboard;
    close(fids);
end

% recrop and rescale A back to original images
Aout{F} = [];
for frm = 1:F
  %
  % fetch data
  bbtc = BBTC(:,frm);
  bbtsc1 = BBTSC1(:,frm);
  Acur = A(:,:,frm);
  %
  % recrop
  Acur = Acur(bbtsc1(2):(bbtsc1(2)+bbtsc1(4)-1), bbtsc1(1):(bbtsc1(1)+bbtsc1(3)-1), :);
  %
  % rescale
  % beware: imresize in matlab 2007b is bugged !
  % not only slightly different results from 2006,
  % but also crashes on some image/sizes combinations !
  Acur = imresize_vitto(Acur, bbtc([4 3]));
  %
  % store
  Aout{frm} = Acur;
end
