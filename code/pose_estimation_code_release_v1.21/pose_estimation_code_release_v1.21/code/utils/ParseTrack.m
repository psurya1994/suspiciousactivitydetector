function [PM, CM, angles, joints] = ParseTrack(frames_dir, format, motion, T, reparse, frixs, pars, poses_dir, verbose)

% estimate limb pose using Ramanan's image parser.
%
% T is a track and T.D contains all the info (frame number, BB coords, class, etc.)
%
% if poses_dir given
% -> output a .jpg for each frame with the estimated human body poses
% superimposed; if image already exists, add to it (allows for multiple tracks on same shot)
%
% if reparse
% -> use color model T.SCM.cm to parse all frms; this implies skipping the parse1 phase
%
% if frixs
% -> only process frames in the intersection between what the track covers and frixs
%
% if pars.use_fg_high
% -> use T.FGH(dix).bb, .fg has the area where to perform parsing;
%    this foreground highlit area substantially facilitates parsing
%
% if pars.only_if_improved
% -> PM(dix) = T.PM(dix) if PM(dix).p > .p from this parsing
%    this is useful during reparsing, where the color model cm might or might
%    not lead to a better pose map. This option needs T.PM to exist.
%
% parameters:
% pars.histeq     = contrast normalization through adaptive histogram equalization
% pars.bb_enlarge = factor to enlarge BB; make sure the whole person is inside BB
% pars.bb_rescale = height to which to normalize BB to a constant height -> easier life for parser (as know limb lengths)
%                   the BB is first enlarged and then rescaled and cropping is taken into account,
%                   so that height of BB before enlarging and cropping = bb_rescale
%

% process arguments
if nargin < 9
  verbose = false;
end
if nargin < 8
  poses_dir = false;
end

if (islogical(poses_dir) || ~isempty(poses_dir))
  if ~exist(poses_dir,'dir')
    mkdir(poses_dir);
  end
end

if reparse
  cm = T.SCM;
else
  cm = false;
end
if islogical(frixs)
  dixs = 1:size(T.D,2);   % just process all frames in the track
else
  [trash dixs trash] = intersect(T.D(1,:), frixs);
end

% want to do full parsing
pars.img_lik_only = false;

%assert(not(pars.use_fg_high && pars.naked), 'Contraddictory parameters !');

% process dets
if isempty(dixs)
  PM = []; CM = [];      % in the case dixs = []
  return;
end

for dix = dixs

    
  % info
  class_id = T.D(9,dix); classname = class_id2name(class_id);
  if verbose
    newline;
    display(['Parsing frame ' num2str(T.D(1,dix)')]);
    display(['Class: ' classname]);
  end
  
  % load image
  fr = T.D(1,dix);
  im = gray2rgb(imread(fullfile(frames_dir, sprintf(format, fr))));
  
  %imwrite(im,'DataSet/o5.jpg','jpg');
  
  % select motion for current frame
  curr_motion = [];
  
  % prepare det bb = [xmin ymin width height]
  bb = T.D(2:5, dix)';     
  
  % parsing
  if ~isfield(pars,'lp_use')
    pars.lp_use = 0;
  end
  if      pars.use_fg_high && ~pars.naked
    [curPM curCM] = ParseImageBB(im, bb, T.FGH(dix), cm, class_id, pars, verbose);
  elseif ~pars.use_fg_high && ~pars.naked
    [curPM curCM] = ParseImageBB(im, bb, false,      cm, class_id, pars, verbose);
  elseif ~pars.use_fg_high &&  pars.naked && ~pars.lp_use     % parse with no aid whatsoever
    [curPM curCM] = ParseImage(im, false, false, false, class_id2name(class_id), pars, verbose);
    curPM.bb = [1 1 size(im,2) size(im,1)];
  else
    error('Contraddictory parameters !');
  end
  % 
  
  % keep previous PM ?
  if isfield(pars,'only_if_improved') && pars.only_if_improved
    sf = pars.pm_qual_fct;
    if (strcmp(sf,'p') && curPM.p < T.PM(dix).p) || (strcmp(sf,'e') && curPM.e > T.PM(dix).e)
      if verbose
        display(['previous PM score better (prev = ' num2str(T.PM(dix).(sf)) ', new = ' num2str(curPM.(sf)) ') -> keep previous PM.']);
      end
      curPM = UncompressPM(T.PM(dix));
      curCM = T.CM(dix);
    else
      if verbose
        display(['new PM score better (prev = ' num2str(T.PM(dix).(sf)) ', new = ' num2str(curPM.(sf)) ') -> keep new PM.']);
      end
    end
  end
  %
  % store PM and CM
  PM(dix) = CompressPM(curPM);
  CM(dix) = curCM;
  
  if isstruct(im)
    im = im.im;
  end

  % paste pose over entire image
  % if image file already exists -> add to it
  if not(islogical(poses_dir))
    

    curPose = uint8(curPM.a*2e3);
    pose_fname = fullfile(poses_dir, sprintf(format, fr));
    imPose = try_imread(pose_fname);
    if islogical(imPose)
    	imPose = im;
    end
    %
    if reparse && ismember(dix, T.SCM.best_dix)                     % paint dets contrib to SCM in red
      col = [1 0 0];
    else
      col = [0 1 0];
    end
    imPose = PasteOverImage(curPose, imPose, curPM.bb, 0.6, 0.5);  % last params -> visibility of curPose and imPose respectively
    if isfield(curPM,'MAP') && isstruct(curPM.MAP)
      curSticks = PaintSticks(round(curPM.MAP.sticks), [size(curPM.a,2) size(curPM.a,1)], class_id2cols(class_id),2);
      imPose = PasteOverImage(curSticks, imPose, curPM.bb);
    end

    safe_imwrite(imPose, pose_fname);
  end
  
  % wait before going to next det
  if verbose > 1
    keyboard;
  end
 
  
  imtool(curSticks);
  
  %% FInding exact angles 
  jointsX = [];
  jointsY = [];
  % body
  body = (curSticks(:,:,1)==255).*(curSticks(:,:,2)==0).*(curSticks(:,:,3)==0);
  [k num] =  my_hough(double(body),2);
  theta = k(1).theta+90;
  
  if(k(1).point1(2)>k(1).point2(2))
      jointsX = [jointsX k(1).point2(1)];
      jointsY = [jointsY k(1).point2(2)];
  else
      jointsX = [jointsX k(1).point1(1)];
      jointsY = [jointsY k(1).point1(2)];
  end
  

  %upper limbs
  ulimb = (curSticks(:,:,1)==0).*(curSticks(:,:,2)==255).*(curSticks(:,:,3)==0);
  [k1 num1] = my_hough(ulimb,8);
  pointX=[];
  for i=1:num1
    pointX = [pointX k1(i).point1(1)];
  end
  
  avgX=mean(pointX);

    for i=1:num1
        if(pointX(i)<avgX)
            lultheta = k1(i).theta+90;
            if(distance(jointsX(1),jointsY(1),k1(i).point1(1),k1(i).point1(2))>distance(jointsX(1),jointsY(1),k1(i).point2(1),k1(i).point2(2)))
                jointsX = [jointsX k1(i).point1(1)];
                jointsY = [jointsY k1(i).point1(2)];
            else
                jointsX = [jointsX k1(i).point2(1)];
                jointsY = [jointsY k1(i).point2(2)];
            end
            break;
        end;
    end;
  
  
    for i=1:num1
        if(pointX(i)>avgX)
            rultheta = k1(i).theta+90;
            if(distance(jointsX(1),jointsY(1),k1(i).point1(1),k1(i).point1(2))>distance(jointsX(1),jointsY(1),k1(i).point2(1),k1(i).point2(2)))
                jointsX = [jointsX k1(i).point1(1)];
                jointsY = [jointsY k1(i).point1(2)];
            else
                jointsX = [jointsX k1(i).point2(1)];
                jointsY = [jointsY k1(i).point2(2)];
            end
            break;
        end;
    end;

    
    
    % lower limbs
    
    llimb = (curSticks(:,:,1)==255).*(curSticks(:,:,2)==255).*(curSticks(:,:,3)==0);
    
    [k2 num2]= my_hough(llimb,10);
    pointY=[];
    for i=1:num2
        pointY = [pointY k2(i).point1(1)];
    end
  
    avgY=mean(pointY);

    
    
    for i=1:num2
        if(pointY(i)<avgY)
            llltheta = k2(i).theta+90;
            if(distance(jointsX(2),jointsY(2),k2(i).point1(1),k2(i).point1(2))>distance(jointsX(2),jointsY(2),k2(i).point2(1),k2(i).point2(2)))
                jointsX = [jointsX k2(i).point1(1)];
                jointsY = [jointsY k2(i).point1(2)];
            else
                jointsX = [jointsX k2(i).point2(1)];
                jointsY = [jointsY k2(i).point2(2)];
            end
            break;
        end;
    end;
    
    
    for i=1:num2
        if(pointY(i)>avgY)
            rlltheta = k2(i).theta+90;
            if(distance(jointsX(3),jointsY(3),k2(i).point1(1),k2(i).point1(2))>distance(jointsX(3),jointsY(3),k2(i).point2(1),k2(i).point2(2)))
                jointsX = [jointsX k2(i).point1(1)];
                jointsY = [jointsY k2(i).point1(2)];
            else
                jointsX = [jointsX k2(i).point2(1)];
                jointsY = [jointsY k2(i).point2(2)];
            end
            break;
        end;
    end;
    joints = [jointsX;jointsY];
    
    
    if(jointsY(2)<jointsY(1))
        lultheta = lultheta - 180;
    end
    if(jointsY(3)<jointsY(1))
        rultheta = rultheta - 180;
    end
    if(jointsY(4)<jointsY(2))
        llltheta = llltheta - 180;
    end
    if(jointsY(5)<jointsY(3))
        rlltheta = rlltheta - 180;
    end
    
    angles = [theta lultheta rultheta llltheta rlltheta];
    %angles_n = angles - theta;
end % loop over dets in the track

