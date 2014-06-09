% script compiling platform-specific mex files
% run this before executing any routines
utilsdir = './utils';
try
  eval(['mex -outdir ' utilsdir ' ' fullfile(utilsdir,'mexDGC.cpp')]);
  eval(['mex -outdir ' utilsdir ' ' fullfile(utilsdir,'nema_lognorm_fast.cxx')]);
catch
  disp('WARNING: Unable to compile mex files required by foregrou nd highlighting, please switch it off by setting the use_fg_high field in parsing parameters to false, or try again with a different compiler');
end
try
  eval(['mex -outdir ' utilsdir ' ' fullfile(utilsdir,'triLinearInterpolation.cpp')]);
  eval(['mex -outdir ' utilsdir ' ' fullfile(utilsdir,'triLinearVoting.cpp')]);
  eval(['mex -outdir ' utilsdir ' ' fullfile(utilsdir,'vgg_nearest_neighbour_dist.cxx')]);
catch
  error('Unable to compile mex files required by the application, please contact the authors');
end

