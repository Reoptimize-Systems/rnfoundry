function T = omega2period (omega)
% convert frequency in radians/s to period in seconds

    T = 1./(omega ./ (2 .* pi));

end