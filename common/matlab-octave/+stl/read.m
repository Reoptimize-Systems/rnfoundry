function [v, f, n, name] = read(fileName)
%STLREAD reads any STL file not depending on its format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

format = stl.getFormat(fileName);
if strcmp(format,'ascii')
  [v,f,n,name] = stl.readAscii(fileName);
elseif strcmp(format,'binary')
  [v,f,n,name] = stl.readBinary(fileName);
end