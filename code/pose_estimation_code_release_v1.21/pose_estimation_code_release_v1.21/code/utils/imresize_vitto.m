function Iout = imresize_vitto(I, varargin)

% imresize compatible with matlab 2006 upwards
%

temp = version('-release');
if str2double(temp(1:4)) > 2006
  Iout = imresize_old(I, varargin{:});
else
  Iout = imresize(I, varargin{:});
end
