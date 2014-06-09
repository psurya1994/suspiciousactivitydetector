function A = side2area(s)

% converts string s to a border area of 2% within a standard (0,0)->(1,1) square
%

switch lower(s)
    case 'left'
        A = [0     0    0.02  1]';
    case 'right'
        A = [0.98  0    1.00  1]';
    case 'top'
        A = [0     0    1  0.02]';
    case 'bottom'
        A = [0   0.98   1     1]';
    otherwise
        error([mfilename ': unknown area named ' s]);
end
