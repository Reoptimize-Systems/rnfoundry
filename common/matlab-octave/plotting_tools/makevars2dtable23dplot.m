function [newx, newy, newz] = makevars2dtable23dplot(x, y, data)
% convert 2d table data to a format suitible for plotting functions
%
% Syntax
%
% [newx, newy, newz] = makevars2dtable23dplot(x, y, data)
%
% Description
%
% converts data formatted as a 2d table (row heading, column headings and
% data) to a format suitible for the standard 3d plotting functions.
%
% Input
%
%  x - vector of m values the same length as the number of rows in 'data'.
%
%  y - vector of n values the same length as the number of columns in
%    'data'.
%
%  data - (m x n) matrix of numeric values. 
%
% Output
%
%  newx - column vector of x coordinates at which to plot the corresponding
%    value in 'newz' (see below)
%
%  newy - column vector of y coordinates at which to plot the corresponding
%    value in 'newz' (see below)
%
%  newz - column vector of data points (the data input reshaped into a
%    column vector).
%
% See Also:  contour2dtable.m  contourf2dtable.m
%

    % do some input checking
    if ~isvector(x) || numel(x) ~= size(data, 1)
        error('x must be a vector of the same number of elements as the number of rows in data')
    elseif ~isvector(y) || numel(y) ~= size(data, 2)
        error('y must be a vector of the same number of elements as the number of columns in data')
    end
    
    newx = repmat (x(:), numel (y), 1);
    newy = reshape (repmat (y(:)', numel (x), 1), [], 1);
    newz = data(:);

end