function OUT = nema_sum( IN, DIM, CLASS )

% NEMA_SUM A function that overrides the native matlab SUM function
% such that if presented with no dimension to sum over, all
% elements of the input array are summed and the output is scalar.
%
% See also SUM for the interface to NEMA_SUM.

% Author: Nicholas Apostoloff <nema@robots.ox.ac.uk>
% Date: 13 Apr 06

if nargin == 2
  if ~isnumeric( DIM )
    OUT = sum( IN(:), DIM );
  else
    OUT = sum( IN, DIM );
  end
elseif nargin == 3
  OUT = sum( IN, DIM, CLASS );
else
  OUT = sum( IN(:) );
end
