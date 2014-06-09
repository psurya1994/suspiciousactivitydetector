function [T sticks_imgcoor angles joints] = PoseEstimStillImage(base_dir, img_dir, img_fname_format, img_number, classname, bb, fghigh_pars, parse_pars, addinf_pars, segm_params, verbose)

% Run a complete pose estimation cycle,
% starting from a detected bounding-box bb = [x y width height]'
%
% base_dir - directory where the results should be stored
% img_dir - name of the directory with images relative to the base_dir
% img_fname_format - format of the image files in the img_dir, e.g. %06d.jpg means that images are in the format xxxxxx.jpg (where x is a digit)
% img_number - image number, used together with img_fname_format to determine the filename of the image to process
% classname - either 'ubf', or 'ubs', or 'full' - used as a selector for proper parameters.
%             'ubf' = upper body front, 'ubs' = upper body side, 'full' = full body
%             in this release, only 'ubf' is supported cleanly.
% bb - window around object of interest [x y width height]';
%      when classname == 'ubf' -> please set bb to be around the head and shoulders (see [2-4] for examples); e.g. use our upper-body detector [5]
% fghight_pars - parameters of the foreground highlighting algorithm (used if it is switched on)
% parse_params - parameters of the parsing algorithm
% addinf_params - parameters of the repulsive model (not used from version 1.2)
% segm_params - parameters of the routine that derives hard segmentations and stickmen from posterior marginals (called 'pose maps' in our code)
% verbose - 0 = no output, 1 = text output, 2 = displaying intermediate figures
%
% some pars.<field>(class_id) = parameters can vary per class
% example bb for the joyce image = [343 27 242 217]'


% process arguments
if nargin < 11
  verbose = true;           % set to verbose=2 to have figures
end

T.D = [img_number; bb];
T.D(9) = class_name2id(classname);
  
% foreground highlighting
if parse_pars.use_fg_high
  T.FGH = FGHighTrack(fullfile(base_dir, img_dir), img_fname_format, T, fghigh_pars, [base_dir '/fghigh_' classname], verbose);
else
  T.FGH = false;
end

% parsing
[T.PM T.CM angles joints] = ParseTrack(fullfile(base_dir, img_dir), img_fname_format, [], T, false, false, parse_pars, [base_dir '/poses_' classname], verbose);
sticks_imgcoor=0;
%T.PM = SegmentTrack(fullfile(base_dir, img_dir),img_fname_format, T, segm_params, [base_dir '/segms_' classname], verbose);
%sticks_imgcoor = convertSticksToImgCoor(T.PM.sticks,[size(T.PM.a,2) size(T.PM.a,1)], T.PM.bb);
end