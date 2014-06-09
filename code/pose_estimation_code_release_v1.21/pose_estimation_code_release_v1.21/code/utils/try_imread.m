function im = try_imread(fname)

% tries to load fname
%
% if fname doesn't exist
% -> return im = false
%
% if fname exists
% but corrupt image file (or not an image file at all)
% -> delete fname and return im = false
%

if exist(fname, 'file')
  try 
    im = imread(fname);
  catch
    disp(['WARNING: ' fname ' is not a healthy image file:']);
    trash = lasterror;
    disp(trash.message);
    disp(['deleting ' fname]);
    system(['rm ' fname]);
    im = false;
  end
else
  im = false;
end
