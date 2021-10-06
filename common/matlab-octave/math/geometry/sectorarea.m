function area = sectorarea(R, theta)
% calculates the area of a circular sector

    area = (R.^2  ./ 2) .* (theta - sin(theta));

end