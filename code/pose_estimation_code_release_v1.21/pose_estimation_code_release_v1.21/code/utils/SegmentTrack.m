function PM = SegmentTrack(frames_dir, format, T, pars, segms_dir, verbose)

% Compute hard segmentations for each part in each frame of the pose-maps PM(dix).
% Also fits sticks to the segments.
%
% Only operates on frames for which T.PM(dix) isn't empty.
%
% Output:
% - PM(dix).segm(:,:,p) = binary segmentation for part p
% - PM(dix).sticks(:,p) = [x1 y1 x2 y2]' = stick of part p
% - org image with segmentation and sticks overlaid
%

% process arguments
if nargin < 6
  verbose = false;
end
if nargin < 5
  segms_dir = false;
else
  mkdir(segms_dir);
end

% process dets
for dix = 1:size(T.D,2)
%for dix = find(ismember(T.D(1,:), 32383))     % Glory swinging
%for dix = 230  % DEVELOP; very cool frame on Tara track 239-1
%for dix = 27  % DEVELOP: double-head on Joyce track 39-1; and overlap ua/la
%for dix = 37  % DEVELOP: re-parsing on Joyce track 39-1
%for dix = 14 % DEVELOP: double-head frame on Buffy gym track 46-3
%for dix = 44 % DEVELOP: re-parsing on Buffy track 46-3
%for dix = 77 % DEVELOP: double-head on Riley track 373-1; and also overlap ua/la 
%for dix = 18 % DEVELOP: lower/upper arm overlap Joyce track 328-1
%for dix = 5 % DEVELOP: re-parsing on Willow track 359-1
%for dix = find(ismember(T.D(1,:), 60190))   % DEBUG

  % skip previously unprocessed frames
  if dix > length(T.PM) || isempty(T.PM(dix).b)
    continue;
  end

  % info
  newline;
  display(['Segmenting detection ' num2str(T.D(:,dix)')]);
  class_id = T.D(9,dix); classname = class_id2name(class_id);
  display(['Class: ' classname]);
  
  % load image
  fr = T.D(1,dix);
  im = gray2rgb(imread(fullfile(frames_dir, sprintf(format, fr))));

  % compute segmentation
  T.PM(dix).segm = pm2segms(UncompressPM(T.PM(dix)), class_id, pars, verbose);
  
  % fit sticks to segmentation
  T.PM(dix).sticks = FitSticks(T.PM(dix).segm, pars, verbose);
  
  % paste segmentation over entire image
  % if image file already exists -> add to it
  if not(islogical(segms_dir))
    segm_fname = fullfile(segms_dir, sprintf(format, fr));
    if isfield(pars,'outimg_type')
      [type name] = Postfix(segm_fname,'.');
      segm_fname = [name '.' pars.outimg_type];
    end
    %segm_fname = [segm_fname(1:end-4) '.tiff'];  % better image quality
    imSegm = try_imread(segm_fname);
    if islogical(imSegm)
      imSegm = im;
    end
    if ~isfield(pars,'draw_segms') || pars.draw_segms
      curSegm = PaintSegmentation(T.PM(dix).segm, class_id2cols(class_id));
      imSegm = PasteOverImage(curSegm, imSegm, T.PM(dix).bb);
    end
    if isfield(pars,'draw_sticks') && pars.draw_sticks
      curSticks = PaintSticks(T.PM(dix).sticks, [size(T.PM(dix).segm,2) size(T.PM(dix).segm,1)], class_id2cols(class_id));
      imSegm = PasteOverImage(curSticks, imSegm, T.PM(dix).bb);
    end
    if isfield(pars,'draw_bb') && pars.draw_bb
      imSegm = PaintBB(imSegm, T.D(2:5,dix), pars.bb_col, -1:1);      % paint detection bb
    end
    safe_imwrite(imSegm, segm_fname);
  end
  % wait before going to next det
  if verbose > 1
    keyboard;
  end

end % loop over dets in the track

% set output
PM = T.PM;
