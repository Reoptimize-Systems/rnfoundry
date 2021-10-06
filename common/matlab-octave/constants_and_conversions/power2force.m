function F = power2force(v, P)
% converts a power developed by a system at a given velocity to a force
%
% Syntax
%
% F = power2force(v, P)
%
% Input
%
%   v - velocity
%
%   P - developed power
%
% Output
%
%   F - force
%

% Created by Richard Crozier 

    % get the force due to the developen power, defined as acting in the
    % opposite direction to the motion of the machine
    F = -P ./ v;
    
    % set the forces at zero velocity to be zero
    F(v == 0) = 0;

end