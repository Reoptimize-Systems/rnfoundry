function T = power2torque(omega, P)
% converts a power developed by a system at a given angular velocity to a torque
%
% Syntax
%
% T = power2torque(omega, P)
%
% Input
%
%   omega - angular velocity in rad/2
%
%   P - developed power
%
% Output
%
%   T - torque
%

% Created by Richard Crozier 

    % get the torque due to the developed power, defined as acting in the
    % opposite direction to the motion of the machine
    T = -P ./ omega;
    
    % set the forces at zero velocity to be zero
    T(omega == 0) = 0;

end