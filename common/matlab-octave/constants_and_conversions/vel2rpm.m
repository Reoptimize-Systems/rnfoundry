function rpm = vel2rpm(v, r)
% converts a velocity at a given radius to revolutions per minute
%
% Syntax
%
% rpm = vel2rpm(v, r)
%

% Created by Richard Crozier 

    % get the number of revolutions per second by dividing the velocity by
    % the length of one revolution
    rps = v ./ (2 * pi * r);
    
    % multiply by 60 to get revolutions per minute
    rpm = rps * 60;

end