function [v, f, n, name] = read (fileName)
% stl.read reads an STL file (binary or ASCII)
%
% Syntax
%
% [v, f, n, name] = stl.read (fileName)
%
% Input
%
%  fileName - path of stl file to be read
% 
% Output
% 
%  v - (n x 3) matrix with the x, y and z coordinates of the mesh vertices
%
%  f - (m x 3) matrix of connectivities representing the faces of the mesh
%
%  n - (m x 3) matrix of face normals
%
%  name - the name of the STL object stored in the file (NOT the name of
%   the STL file)
%
%

    format = stl.getFormat(fileName);
    
    if strcmp(format,'ascii')
        
        [v,f,n,name] = stl.readAscii(fileName);
        
    elseif strcmp(format,'binary')
        
        [v,f,n,name] = stl.readBinary(fileName);
        
    end

end