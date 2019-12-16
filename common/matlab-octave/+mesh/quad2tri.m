function trifaces = quad2tri (quadfaces)
% creates a triangle mesh from a quad mesh
%
% Syntax
%
% trifaces = mesh.quad2tri (quadfaces)
%
% Description
%
% mesh.quad2tri creates a triangle mesh from a quad mesh by splitting each
% quad into two triangles. It only takes face as input as no new vertices
% are created.
%
% Input
%
%  quadfaces - (4 x n) matrix of indices into a vertices matrix where each
%   column represents one quadrilateral element.
%
% Output
%
%  trifaces - (3 x n) matrix of indices into the same vertices matrix as 
%   referred to by quadfaces, where each column represents one triangular
%   element.
%
%

    trifaces = [ quadfaces(1:3,:), ...
                 [ quadfaces(3,:);
                   quadfaces(4,:);
                   quadfaces(1,:) ] ...
                ];

end