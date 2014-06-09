function [STD, STM, BBTSC1, BBTC, fids] = PrepareFGHighRun(T, dix_begin, dix_end, image1, pars, verbose)

% prepare 3D image data STD and 3D label map STM
% for spatio-temporal grabcut
%
% Output:
% - STD(y,x,col_chan,frm) = 3D image data, with frm going from 1 to dix_end-dix_begin+1
% - STM(y,x,frm) = 3D label map for initializing grabcut
% - BBTSC1(1:4,frm) = coordinates of cropped bb in the standard-scaled coordinates of STD/STM
% - BBTC(1:4,frm) = coordinates of cropped bb in the original image; all bbs in [min_x max_x width height] format
% - fids(frm) = figure with preparations for frm (if verbose > 1)
%

% class info
class_id = T.D(9,dix_begin);
classname = class_id2name(class_id);
display(['Class: ' classname]);

% find out ideal width W and height H of enlarged and rescaled window
ar = T.D(4,dix_begin)' / T.D(5,dix_begin)';         % aspect-ratio (width/height)
H = round((pars.bb_enlarge(4,class_id) - pars.bb_enlarge(2,class_id)) * pars.bb_rescale(class_id));
W = round((pars.bb_enlarge(3,class_id) - pars.bb_enlarge(1,class_id)) * pars.bb_rescale(class_id) * ar);
STD = zeros(H,W,3,(dix_end-dix_begin+1),'uint8');   % whole spatio-temporal data
STM = zeros(H,W,(dix_end-dix_begin+1),'double');    % whole spatio-temporal masks
BBTC = zeros(4,dix_end-dix_begin+1);                % all bbtc
BBTSC1 = zeros(4,dix_end-dix_begin+1);              % all bbtsc1

fids = zeros(dix_end-dix_begin+1);                  % figure indeces (in verbose mode)
for cur_dix = dix_begin:dix_end
    % load image
    fr = T.D(1,cur_dix);                            % absolute frame index
    im = gray2rgb(image1);
    %
    % prepare det_bb = [xmin ymin width height]
    det_bb = T.D(2:5, cur_dix)';
    %
    % prepare foreground highlight
    [imProc mask bbtc bbtsc bbts fids(cur_dix-dix_begin+1)] = PrepareFGHigh(im, det_bb, pars, class_id, verbose, fr);
    % all bbt* are in [min_x max_x width height] format
    % bbtsc corresponds exactly to imProc
    %
    % set bbtsc1 to right place within bbts
    bbtsc1 = bbtsc;
    bbtsc1(1:2) = bbtsc(1:2) - bbts(1:2) + 1;
    %
    % pad imProc and mask to make them the size of bbt_scal
    imProcPadded = zeros(bbts(4), bbts(3), 3, class(imProc));
    imProcPadded(bbtsc1(2):(bbtsc1(2)+bbtsc1(4)-1), bbtsc1(1):(bbtsc1(1)+bbtsc1(3)-1), :) = imProc;
    maskPadded = zeros(bbts(4), bbts(3), 1, class(mask));                                 % padded pixels are 0-labeled -> clamped, background, don't learn color distr (essentially discard them ;)
    maskPadded(bbtsc1(2):(bbtsc1(2)+bbtsc1(4)-1), bbtsc1(1):(bbtsc1(1)+bbtsc1(3)-1), :) = mask;
    %
    % add to total spatio-temporal data
    % (resizing just to set 1-2 pixels roundoff errors)
    STD(:,:,:,cur_dix-dix_begin+1) = imresize_vitto(imProcPadded,[H W],'nearest');        % nearest: avoid blurring, as it's just 1-2 pixels diff !
    STM(:,:,cur_dix-dix_begin+1) = imresize_vitto(maskPadded,[H W],'nearest');
    % adjust bbtsc1 too, otherwise it might go out of STD/STM ! (due to the above resizing)
    s = [H W] ./ [size(imProcPadded,1) size(imProcPadded,2)];                             % scale factors for height and width
    bbtsc1(3:4) = floor(bbtsc1(3:4) .* s([2 1]));                                         % floor -> make sure stay within limits when inverse trafo after grabcut
    bbtsc1(1:2) = floor(1 + (bbtsc1(1:2)-1).*s([2 1]));
    %
    % store bb info for later re-cropping
    BBTC(:,cur_dix-dix_begin+1) = bbtc';                                                  % cropped bb in coords of org image
    BBTSC1(:,cur_dix-dix_begin+1) = bbtsc1';                                              % cropped bb in abs coords of STD/STM
end % loop over frames in the temporal window

% show output
if verbose > 1
    % show 'videos' of the aligned image data and label maps
    if false
    [t t t t t fg_segm] = buildHistExps(STD, STM==192);     % check what a direct pixel-labeling based on colors would give
    STMnew = STM;
    STMnew(fg_segm) = 192;                                     % help grabcut, by showing more samples of fg colors
    STMnew(STM==192) = 192;
    STM = STMnew; clear STMnew;
    end
    video_fid = ShowFGHighRun(STD, STM, 0.5);
    keyboard;
    close(video_fid);
end
