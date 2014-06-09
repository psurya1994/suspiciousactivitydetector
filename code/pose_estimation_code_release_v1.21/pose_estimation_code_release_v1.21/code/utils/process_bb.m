function [bbt, bbtc, bbts, bbtsc, bbts1, bbtsc1] = process_bb(det_bb, siz, class_id, pars, verbose)

% Computes all important transformations
% of the detection bounding-box det_bb
%
% bbt     =  det_bb after enlarging
% bbtc    =  bbt after crop
% bbts    =  bbt after rescaling so that det_bb's height == pars.bb_rescale(class_id)
% bbtsc   =  bbtc after rescaling by the same factor
% bbts1   =  bbts starting at (1,1)
% bbtsc1  =  bbtsc in coordinate system starting at (1,1) for bbts1
%
% siz(1:2) = [W H] = width and height of image
%
% all input/output bb* = [minx miny width height]
%

% process arguments
if nargin < 5
  verbose = false;
end

% find enlarged bounding-box before cropping to fit image (bbt)
enlarge = pars.bb_enlarge(:, class_id);
bbt = enlarge .* det_bb([3 4 3 4])' + det_bb([1 2 1 2])';   % [minx miny maxx maxy]
bbt = minmax2wh(bbt');                                      % [minx miny width height]
bbt = round(bbt);
% crop bbt
bbtc = CropBBwh(bbt, siz);

% scale bbt
s = pars.bb_rescale(class_id) / det_bb(4);                  % scale factor
bbts = bbt;
bbts(3:4) = bbts(3:4)*s;
bbts(1:2) = (bbts(1:2)-1)*s + 1;                            % rescale min,max so that origin always (1,1)

% crop bbts
bbtsc = CropBBwh(bbts, siz*s);

% set origin of rescaled bbs to (1,1)
bbts1 = bbts;
bbts1(1:2) = 1;
bbtsc1 = bbtsc;
bbtsc1(1:2) = bbtsc(1:2) - bbts(1:2) + 1;

% round everything
bbtc = round(bbtc);
bbts = round(bbts);
bbtsc = round(bbtsc);
bbts1 = round(bbts1);
bbtsc1 = round(bbtsc1);

% enforce constraints to counter roundoff errors
bbtc = CropBBwh(bbtc, siz);

% display all produced BBs, to verify correctness
if verbose > 2
    figure; hold;
    DrawBB([1 siz(1) 1 siz(2)], [0 0 0], 1);                % image border (black)
    DrawBB(wh2mxmxmymy(det_bb), [0 1 0], 3);                % detection bb (green thick)
    DrawBB(wh2mxmxmymy(bbt), [0 1 0], 3);                   % enlarged bb (green thick)
    DrawBB(wh2mxmxmymy(bbtc), [0 0 1], 1);                  % cropped bb (blue thin)
    axis equal; axis tight; axis ij;
    %
    figure; hold;
    DrawBB(wh2mxmxmymy(bbts1), [0 1 0], 3);                 % rescaled enlarged bb at (1,1) (green thick)
    DrawBB(wh2mxmxmymy(bbtsc1), [0 0 1], 1);                % rescaled enlarged cropped bb in coord sys starting at (1,1) before cropping (blue thin)
    axis equal; axis tight; axis ij;
end
