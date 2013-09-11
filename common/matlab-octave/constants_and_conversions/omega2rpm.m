function rpm = omega2rpm(omega)
% converts an angular velocity in radians per second to rpm

% Created by Richard Crozier 

    rpm = 60 * (omega / (2*pi));


end