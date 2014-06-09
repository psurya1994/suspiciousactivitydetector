function Scol = XY2COL(S)

% converts segments S(:,ix) = [x1 y1 x2 y2]'
% to Scol(:,ix) = [ctr_x ctr_y orient length]'
%
% with orient in [0,pi]
%

ctrs = [mean(S([1 3], :)); mean(S([2 4], :))];
dx = S(3,:)-S(1,:);
dy = S(4,:)-S(2,:);
lens = sqrt(dx.^2+dy.^2);
orients = atan(dy./dx);  % in [-pi/2,pi/2]
neg = orients < 0;
orients(neg) = orients(neg) + pi;

Scol = [ctrs; orients; lens];
