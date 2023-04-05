function P = torque2power(omega, T)
% converts a torque at a given angular velocity to a power
%
% Syntax
%
% P = torque2power(omega, T)
%
% Input
%
%   omega - angular velocity in rad/2
%
%   T - torque
%
% Output
%
%   P - developed power
%

% Created by Richard Crozier 

    % get the power developing the torque
    P = T .* omega;

end