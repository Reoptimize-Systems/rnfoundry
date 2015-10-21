function [xq,yq,zq] = mesh2dtable(x, y, data, varargin)
% creates a mesh plot of the data in a 2D table
%
% Syntax
%
% mesh2dtable(x, y, data, ...)
%
% Inputs
% 
%  x - vector of values representing the independent variable variables of
%   the rows of the table, must be linearly spaced
%
%  y - vector of values representing the independent variable variables of
%   the columns of the table, must be linearly spaced
%
%  data - (2 x 2) matrix of data containing the dependent variables
%    corresponding to each combination of the independent values provided
%    in the x and y vectors
%
% Example
%
% x = 1:10
% y = 0.5:0.5:5
% data = bsxfun(@times, x, y')
% mesh2dtable(x, y, data)
%
% See also, mesh
%

% Copyright Richard Crozier 2012-2015

    [newx, newy, newz] = makevars2dtable23dplot(x, y, data);
    
    % now mesh the data
    [xq,yq] = meshgrid(x,y);
    
    zq = griddata(newx, newy, newz, xq, yq);
    
    % and create the mesh plot
    mesh(xq,yq,zq, varargin{:});

end