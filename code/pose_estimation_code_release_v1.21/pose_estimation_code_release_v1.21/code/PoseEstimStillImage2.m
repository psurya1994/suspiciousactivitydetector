function [T sticks_imgcoor curPose curSticks curPM angles joints] = PoseEstimStillImage2(image1, classname, bb, fghigh_pars, parse_pars, addinf_pars, segm_params, verbose)

if nargin < 11
  verbose = true;           % set to verbose=2 to have figures
end

T.D = [0;bb];
T.D(9) = class_name2id(classname);

 if parse_pars.use_fg_high
   T.FGH = FGHighTrack2(image1, T, fghigh_pars, verbose);
 else
   T.FGH = false;
 end

[T.PM T.CM curPose curSticks curPM angles joints] = ParseTrack2(image1, [], T, false, false, parse_pars, verbose);
sticks_imgcoor=0;

end