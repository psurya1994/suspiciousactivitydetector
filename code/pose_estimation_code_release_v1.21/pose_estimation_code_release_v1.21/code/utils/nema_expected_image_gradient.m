function OUT = nema_expected_image_gradient( C, N )

% NEMA_EXPECTED_IMAGE_GRADIENT A function that calculates the
% expected squared gradient magnitude over an image.
%
%    E = NEMA_EXPECTED_IMAGE_GRADIENT( C, N ) calculates the expected
%    image gradient over the entire image C.  For grayscale images
%    it computes
%      E = <(z_i - z_j)^2>
%    where i and j are neighbouring (N-connected) pixels.  In the case
%    of colour images, we estimate the squared norm of the colour
%    gradient
%      E = <||z_i - z_j||^2>
%
%    N can be either 4 or 8 connected.
%
% See also NEMA_ALPHA_MATTE_GRABCUT for usage.

% Author: Nicholas Apostoloff <nema@robots.ox.ac.uk>
% Date: 13 Apr 06

[ nrows, ncols, ndims ] = size( C );
t1 = sqrt(2);
if nargin < 2 | N == 4
  OUT = ...
      nema_sum( ( C( 2:end, 1:end, : ) - C( 1:end-1, 1:end, : ) ).^2 ) + ...
      nema_sum( ( C( 1:end, 2:end, : ) - C( 1:end, 1:end-1, : ) ).^2 );
  OUT = OUT / ( ( nrows - 1 ) * ncols + nrows * ( ncols - 1 ) );
elseif N == 8
  OUT = ...
      nema_sum( ( C( 2:end, 1:end, : ) - C( 1:end-1, 1:end, : ) ).^2 ) + ...
      nema_sum( ( C( 1:end, 2:end, : ) - C( 1:end, 1:end-1, : ) ).^2 ) + ...
      nema_sum( ( C( 2:end, 2:end, : ) - C( 1:end-1, 1:end-1, : )).^2)/t1 + ...
      nema_sum( ( C( 2:end, 1:end-1, : ) - C( 1:end-1, 2:end, : )).^2)/t1;
  
  OUT = OUT / ( (nrows-1)*ncols + nrows*(ncols-1) + 2*(nrows-1)*(ncols-1) );
  
else
  error( 'N can only be 4 or 8.' );
end
