function isin = incylinder(centre, r, h, x, y, z)
% determines if points are in a cylinder oriented with its axis along the z
% direction
%
% Syntax
%
% isin = incylinder(centre, r, h, x, y, z)
%
% Input
%
%  centre - 3 element vector containing the x, y and z coodinate of the
%    center of the cylinder
%
%  r - radius of the cylinder
% 
%  h - height, or length of the cylinder along its axis
% 
%  x,y,z - vectors or matrices of the same size containing the x,y and z
%    positions of the coordinates to test
%
% Output
%
%  isin - matrix of the same size as the input x,y and z matrices
%   containing boolean flags, true if the corresponding point lies in the
%   cylinder (or on its surface), or false otherwise.
%


    isbelowz = (z <= (centre(3) + h/2));
    isabovez = (z >= (centre(3) - h/2));
    
    radialdistsquared = (x - centre(1)).^2 + (y - centre(2)).^2;
    
    isinrad = (radialdistsquared <= r^2);
    
    isin = (isbelowz & isabovez & isinrad);

end
