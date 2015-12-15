function hz = omega2hz (omega)
% convert frequency in radians/s to Hz

    hz = omega ./ (2 .* pi);

end