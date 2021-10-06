function omega = vel2omega(v, r)
% converts a tangential velocity at a given radius to radians per second
%
% Syntax
%
% omega = vel2omega(v, r)
%

% Created by Richard Crozier 

    % get the number of revolutions per second by dividing the velocity by
    % the length of one revolution
    rps = v ./ (2 * pi * r);
    
    % multiply by 2 * pi to get radians per second
    omega = rps * 2 * pi;

end