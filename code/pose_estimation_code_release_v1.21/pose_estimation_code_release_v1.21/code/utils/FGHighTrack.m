function FGH = FGHighTrack(frames_dir, format, T, pars, fghigh_dir, verbose)

% foreground highlight all detections in the track, based on a conservative grab-cut.
%
% Idea is to remove as much of the bg as possible, without losing any of the fg.
%
% T is a track and T.D contains all the info (frame number, BB coords, class, etc.)
%
% if fghigh_dir given
% -> output a .jpg for each frame with the estimated fg high superimposed;
% if image already exists, add to it (allows for multiple tracks on same shot)
%
% if pars.temp_win > 1 -> apply spatio-temporal grabcut
%
% Output:
% FGH(dix).bb = [xmin ymin xmax ymax] = BB where FG high applied (in img coordinates)
%         .fg = fg segm of bb (directly firs .bb, no rescaling needed)
%
% parameters:
% pars.abg_area(:,rix) = set of rects; absolute bg: clamped and learn color models
% pars.bg_area(:,rix)   = set of rects defined relative to det BB, as if it were [0,0] to [1,1]
%                         bg: learn col models but not clamped
% pars.fg_area(:,rix)   = set of rects, or empty if not used (coord sys = [xmin ymin xmax ymax])
% pars.afg_area(:,rix)  = afg and fg are analo to abg and bg, but for foreground
% pars.bb_enlarge = factor to enlarge BB to obtain region to apply fg highlight to
% pars.bb_rescale = height to which to normalize BB to a constant height
%                   + approach parsing conditions, + keep processing time constant
% pars.temp_win         = number of frames to process at the same time
%
% supports multi-grabcut:
% .abg, .bg, .fg, .afg are cell arrays, .abg{ix} etc. contain the rectangles for the ix-th grabcut
% final output is union of all individual grabcuts; if want just good old
% single grabcut, just use a length 1 cell
%

% parse arguments

if nargin < 6
  verbose = false;
end
if nargin < 5 
  fghigh_dir = false;
end
if not(islogical(fghigh_dir))
  mkdir(fghigh_dir);
end

% init outputs
Ndets = size(T.D,2);
if Ndets == 0
  FGH = [];
  return;
end
FGH(Ndets).fg = [];
FGH(Ndets).bb = [];

% estabilish temporal schedule
% run R{rix} = list of frames in it
R = ContRuns(T.D(1,:));
display([num2str(length(R)) ' runs in this track.']);

% process runs
for rix = 1:length(R)
  % spatio-temporal window
  % continuous, at most length pars.temp_win
  dixs = R{rix};
  display(['Foreground highlighting run ' num2str(T.D(1,dixs(2,[1 end])))]);
  
  for frm_begin = 1:pars.temp_win:size(dixs,2)
    dix_begin = dixs(2,frm_begin);         % dix point to T.D; whereas frm point to the current chunk
    dix_end = min(dixs(2,end),dix_begin+pars.temp_win-1);
    
    newline;
    display(['Foreground highlighting chunk ' num2str(T.D(1,[dix_begin dix_end]))]);
    
    % loop over the grabcuts to carry out (multi-grabcut)
    Ngcs = length(pars.afg_area);          % number of grabcuts requested
    for gcix = 1:Ngcs
      if Ngcs > 1
        display(['Grabcut #' int2str(gcix)]);
      end
      %
      % prepare data for foreground highlighting of this time chunk
      cur_pars = pars;                     % parameters of current grabcut
      cur_pars = copy_fields(pars, cur_pars, {'abg_area', 'bg_area', 'afg_area', 'fg_area', 'border'}, gcix);
      [STD STM BBTSC1 BBTC fids] = PrepareFGHighRun(T, dix_begin, dix_end, frames_dir, format, cur_pars, verbose);
      %
      % run grabcut
      Atemp = FGHigh(STD, STM, BBTSC1, BBTC, pars, fids, verbose);
      clear STD STM;
      %
      % union to previous grabcuts if necessary
      if gcix > 1
        for fix = 1:(dix_end-dix_begin+1)
          A{fix} = (A{fix} | Atemp{fix});
        end
      else
         A = Atemp;
         clear Atemp;
      end
    end
    clear Atemp;
    
    % set outputs
    % now A{fix} == binary segmentation for frame fix, resulting from union
    % of all grabcuts at that frame
    for frm = 1:length(A)
      FGH(dix_begin+frm-1).bb = wh2minmax(BBTC(:,frm));
      FGH(dix_begin+frm-1).fg = A{frm};
    end
    
    % paste fghighed area A over original frames
    % if an fghighed image file already exists (from another track) -> add to it
    if not(islogical(fghigh_dir))
        for frm = 1:length(A)
            fr = T.D(1,dix_begin+frm-1);
            fghigh_fname = fullfile(fghigh_dir, sprintf(format, fr));
            imFG = false;
            if exist(fghigh_fname, 'file')
                imFG = try_imread(fghigh_fname);
            end
            if islogical(imFG)
                im_fname = fullfile(frames_dir, sprintf(format, fr));
                imFG = gray2rgb(imread(im_fname));
            end
            % data
            %keyboard;
            bbtc = BBTC(:,frm);
            Acur = A{frm};
            % semi-transparent pasting
            imFG = PasteOverPlane(uint8(Acur)*255, imFG, bbtc, 2, 0.5, 1.0);   % 2 = green plane; 0.5 = visibility of Acur; 1.0 = visibility of org im
            %
            % paint bbs
            det_bb = T.D(2:5,dix_begin+frm-1);
            imFG = PaintBB(imFG, det_bb, [0 1 0], -1:1);                % paint detection bb
            imFG = PaintBB(imFG, bbtc, [0 1 0], -1:1);                  % paint enlarged and cropped bb where fghigh took place
            % writeout
            imwrite(imFG, fghigh_fname);
        end % loop over frames in this chunk
    end % export fghigh to images ?
    
  end % loop over chunks
  
end % loop over runs
