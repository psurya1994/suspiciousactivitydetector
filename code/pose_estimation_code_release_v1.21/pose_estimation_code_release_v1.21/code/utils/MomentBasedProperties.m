function [center, orient, elongation, max_theta] = MomentBasedProperties(C)

% computes the center, principal axes of inertia (orient),
% and elongation of region C.
%
% Input:
% C(:,p) = [x y]'
%

% just in case ! If C is integer, the whole operations below are screwed up !
% (especially if it's unsigned !)
C = double(C);    

% region center
center = mean_matrix(C')';

% moments
m11 = sum( (C(1,:)-center(1)).*(C(2,:)-center(2)) );
m20 = sum( (C(1,:)-center(1)).^2 );
m02 = sum( (C(2,:)-center(2)).^2 );

% orientation
theta = 0.5*atan(2*m11/(m20-m02));
thetaP = theta+pi/2;

% is it the theta of the max or min energy ?
Mtheta = m20*sin(theta)^2 + m02*cos(theta)^2 -2*m11*sin(theta)*cos(theta);
MthetaP = m20*sin(thetaP)^2 + m02*cos(thetaP)^2 -2*m11*sin(thetaP)*cos(thetaP);

% want orientation of minimal energy
% (as opposed to greylevel images, where want orient of max energy)
if MthetaP < Mtheta 
  orient = thetaP;
  max_theta = MthetaP;
else
  orient = theta;
  max_theta = Mtheta;
end

% force orient to be in [0,pi]
if orient < 0
  orient = orient + pi;
end
if orient > pi
  orient = orient - pi;
end

% elongation
if MthetaP==0 || Mtheta==0
  elongation = 0;
else
  elongation = min(Mtheta/MthetaP, MthetaP/Mtheta);
end
