function theta = chord2angle (c, R)
% calculate angle (in radians) from chord length and radius

    theta = 2 .* asin (c ./ (2 .* R));

end