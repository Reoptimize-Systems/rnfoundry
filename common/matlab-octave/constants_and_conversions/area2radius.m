function r = area2radius(area)
% calculates the radius of a circle from its cross-sectional area

    % Area = pi r^2, so 
    %
    % Area / pi = r^2
    % 
    % sqrt(Area / pi) = r
    r = sqrt(area / pi);

end