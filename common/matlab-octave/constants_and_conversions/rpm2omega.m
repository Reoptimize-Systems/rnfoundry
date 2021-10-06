function omega = rpm2omega(rpm)
% converts revolutions per minute to angular velocity (radians/s)
%
% Syntax
%
% omega = rpm2omega(rpm)
%

% Created by Richard Crozier 

    omega = 2 * pi * rpm2freq(rpm);
    
end