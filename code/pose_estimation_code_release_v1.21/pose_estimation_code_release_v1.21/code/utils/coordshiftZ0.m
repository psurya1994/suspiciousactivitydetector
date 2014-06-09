function newcoord = coordshiftZ0(coord,siz,len,wid)

if nargin < 4,
    wid = 0;
end

imo = siz(3);
stepo = 360/imo;                                                      % stepo = 15 degrees
[ux,vx] = pol2cart(([stepo:stepo:360]+180-stepo+180) * pi/180,wid);   % Width offset
[uy,vy] = pol2cart(([stepo:stepo:360]+90-stepo+180) * pi/180,len);    % Length offset

y=coord(1);
x=coord(2);
o=coord(3);

vx = -vx;
vy = -vy;
xvec = round(ux+uy);
yvec = round(vx+vy);

newcoord = [-yvec(o)+y,-xvec(o)+x,o];

