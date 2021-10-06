function [newx,newy] = scatterplot2dtable(x, y, data, varargin)
% creates a scatter plot of the data in a 2D table
%
% Syntax
%
% scatterplot2dtable(x, y, data, ...)
%
% Inputs
% 
%  x - vector of values representing the independent variable variables of
%   the rows of the table
%
%  y - vector of values representing the independent variable variables of
%   the columns of the table
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
% scatterplot2dtable(x, y, data)
%
% See also, scatter3
%

% Copyright Richard Crozier 2012

    [newx, newy, newz] = makevars2dtable23dplot(x, y, data);
    
    clear data;
    
    % now plot the data using scatter3
    scatter3(newx, ...
             newy, ...
             newz, ...
             varargin{:});

end