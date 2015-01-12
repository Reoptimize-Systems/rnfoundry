function len = chordlen (r, theta)
% returns the length of the chord with given radius and segment angle
%
% Input
%
%   r - matrix of radii at which the chord length is to be calculated
%
%   theta - matrix of span angles
%
%
% Output
%
%   len - the chord length
%

    if any (theta(:) <= 0 | theta(:) > 180)
        error ('theta must be between 0 and 180')
    end
    
    len = 2 .* r .* sin ( theta ./ 2 );

end