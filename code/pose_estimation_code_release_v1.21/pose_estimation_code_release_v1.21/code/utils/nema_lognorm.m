function p = nema_lognorm( x, mu, sigma )

% NEMA_LOGNORM A function that calculates the log of the normal
% function directly with evaluating the exponential or the logarithm.
%
%    LOGP = NEMA_LOGNORM( X, MU, SIGMA) calculates the log of the
%    normal function LOGP of the N-D dimensional columnwise points
%    in X according to
%
%      LOGP(X) = log( 1 / sqrt( (2*pi)^N * det(SIGMA) ) ) -
%                0.5 * (X-MU)' * SIGMA^(-1) * (X-MU)
%
%    MU is a Dx1 mean column vector and SIGMA is the DxD covariance
%    matrix.



% Author: Nicholas Apostoloff <nema@robots.ox.ac.uk>
% Date: 03 Jun 05

d = length( mu );

if size( x, 1 ) ~= d | size( sigma, 1 ) ~= size( sigma, 2 ) ...
      | size( sigma, 1 ) ~= d
  error( 'Inconsistent input dimensions' );
end

n = size(mu,1);
npoints = size(x,2);

t1 = ( (2*pi)^n * abs(det(sigma)) ).^(-0.5);
t2 = inv(sigma);
%p = zeros( 1, npoints );
p=nema_lognorm_fast( x, mu, t2, t1 );
