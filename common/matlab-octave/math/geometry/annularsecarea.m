function area = annularsecarea(Ri, Ro, theta)
% calculates the area of an annular section, or an angular subportion of
% the annulus
%
% Syntax
%
% area = annularsecarea(Ri, Ro, theta)
%
% Input
%
%   Ri - Inner radius of the annulus
%
%   Ro - Outer radius of the annulus
%
%   theta - (optional) an angle representing a small section of the
%     annulus, the area returned will be the area of this portion of the
%     annulus. If not supplied the full annulus area is returned
%     (equivalent to setting theta = 2*pi). Angle in radians.
%
% Output
%
%   area - area of the annulus, or portion of the annulus
%

    if nargin < 3
        theta = 2*pi;
    end
    
    area = (theta / 2) * (Ro.^2 - Ri.^2);
    
end