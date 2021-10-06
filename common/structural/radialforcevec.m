function forcevec = radialforcevec(xyz, forcemag)
% calculates a vector of forces ponting in the radial direction on a plane
% given the x-y coordinate of the point of application in the plane

    % calculate the force vector, this will be the force times a unit
    % vector pointing in the radial direction
    %
    % this unit vector can be found from the cross product of the vector in
    % pointing from the origin to the x and y location of the point in the
    % z = 0 plane, and a vector pointing along the z axis
    forcevec = forcemag * unit([xyz(1:2), 0]);
    

end