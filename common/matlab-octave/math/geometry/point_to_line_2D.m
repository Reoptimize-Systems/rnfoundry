function d = point_to_line_2D(pt, v1, v2)
% find the shortest distance of a point from a line in 2D coordinates

    a = v1 - v2;
    b = pt - v2;
    d = norm(det([a,b])) / norm(a);

end