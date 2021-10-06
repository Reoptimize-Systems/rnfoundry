function vel = rpm2vel(rpm, r)
% converts a revolutions per minute to a tangential velocity at a given
% radius 
%
% Syntax
%
% vel = rpm2vel(rpm, r)
%

% Created by Richard Crozier 

    % get the number of revolutions per second
    rps = rpm / 60;
    
    % multiply revolutions per second by the distance of that revolution to
    % get the tangential velocity
    vel = rps .* (2 * pi * r);

end