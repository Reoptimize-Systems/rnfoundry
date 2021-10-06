function magn = centrifugalbodyforce(xyz, omega, density)
% calculates the centrifugal force per unit volume throughout a mass
% rotating about the z axis
%
% Syntax
% 
% magn = centrifugalbodyforce(xyz, v, density)
%
% Input
%
%   xyz - coordinate at which the force is to be calculated
%
%   omega - velocity in radians/s
%
%   density - density of the mrotating material
%
% Output
%
%   magn - a vector of forces per unit mass acting at the point xyz in the
%          mass
%

    % get the radial distance from the center
    [theta,R] = cart2pol(xyz(1), xyz(2));
    
    % calculate the acceleration at this distance
    a = omega^2 * R;
    
    % the force per unit volume is density times the acceleration (the
    % total force would be F = ma = density * V * a, so force per unit
    % volume is just F = density * a)
    Fcentrifugal = density * a;
    
    magn = radialforcevec(xyz, Fcentrifugal);

end