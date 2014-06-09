function [PM, CM, imProc] = ParseImageBB(im, det_bb, FGH, cm, class_id, pars, verbose)

if isstruct(im)
  motion = im.motion;
  im = im.im;
else
  motion = [];
end

% parse image starting from detection bb, according to params pars.
%
% det bb = [xmin ymin width height]
%
% FGH.bb + FGH.fg:
% foreground highlighting
%
% Returns uncompressed PM;
% also includes PM.bb = [xmin ymin width height] = image area where parsing applied
%

% limb location constraints (LLC)
% transform LLC areas from rel-to-bb coordinates to image coordinates of det_bb
% and produce LCC masks (includes cropping if necessary)

%not anymore
%assert(not(pars.use_fg_high && pars.lp_use), 'Contradictory parameters !');

Nlimbs = length(pars.part_ids{class_id});
limb_mask = LLC2im(im, pars.limb_areas{class_id}, det_bb);


% enlarge det BB
org_height = det_bb(4);                                           % keep org_height to make resizing properly
enlarge = pars.bb_enlarge(:, class_id);
bb = enlarge .* det_bb([3 4 3 4])' + det_bb([1 2 1 2])';
bb = round(minmax2wh(bb'));                                        % enlarged bb = [xmin ymin width height]

% crop enlarged BB on image and constrain search within it only
if (pars.use_fg_high)     % always use bb from FGH, regardless of what pars.bb_enlarge says !
      bb = FGH.bb;                                                  % use fg high bb = [xmin ymin xmax ymax]
    bb = minmax2wh(bb);
end
[imProc bb] = CropImageBB(im, bb);


% crop limb location constraint masks
% at this point bb = [xmin ymin width height] is already the correct one
% indep of pars.use_fg_high
for lix = 1:Nlimbs
   if not(isempty(limb_mask{lix}))
     limb_mask{lix} = CropImageBB(limb_mask{lix}, bb);
   end
end

% contrast normalization
if pars.histeq(class_id)
    imProc = vitto_clahe(imProc);
end

% resize
imscale = pars.bb_rescale(class_id) / org_height;                 % rescale so that height of det_BB(before enlarging and cropping) = bb_rescale
imProc = imresize(imProc, imscale, 'bilinear');

%handle motion
if ~isempty(motion)
  temp = imProc;
  imProc = [];
  imProc.im = temp;
  switch pars.motioncue_type
    case 'klt'
      [imProc.motion.mag imProc.motion.theta imProc.motion.x imProc.motion.y] = vgg_klt_tracks_to_image(motion,false,wh2minmax(bb),imscale);
    case 'ofTVL1'
      motionCrop = CropImageBB(motion, bb);
      motionProc = imresize(motionCrop, imscale, 'bilinear') * imscale;%rescale motion map the motion itselfs too
      imProc.motion.x = motionProc(:,:,1);
      imProc.motion.y = motionProc(:,:,2);
      imProc.motion.mag = sqrt(imProc.motion.x.^2+imProc.motion.y.^2);
      imProc.motion.theta = atan2(imProc.motion.y,imProc.motion.x);
     otherwise
      error('unsupported motioncue_type string');
  end
  % move motion thetas into range <0,2pi>
  imProc.motion.theta = mod(imProc.motion.theta + 2*pi,2*pi);
  % warning motion.x y and mag are not normalized to perframe values yet.
end


if pars.use_fg_high
    imFG = FGH.fg;
    siz = size(imProc);
    imFG = imresize(imFG, siz(1:2), 'bilinear');
else
    imFG = false;
end
for lix = 1:Nlimbs
   if not(isempty(limb_mask{lix}))
     limb_mask{lix} = imresize(limb_mask{lix}, imscale, 'bilinear');
   end
end

% parse image (find limb poses)
if pars.lp_use
  [PM CM] = ParseImageUsingLP(imProc, imFG, cm, im, det_bb, limb_mask, class_id, pars, verbose);
else
  [PM CM] = ParseImage(imProc, cm, imFG, limb_mask, class_id, pars, verbose);
end
PM.bb = bb;       % image area on which parsing applied
if isstruct(imProc)
  imProc = imProc.im;
end