function [x,y,z] = table2vectors(colh,rowh,values,doprunenans)
% converts a 2d table of values to a set of x, y and z coordinates in the
% function space
%
% Syntax
%
% [x,y,z] = table2vectors(colh,rowh,values)
% [x,y,z] = table2vectors(colh,rowh,values,doprunenans)
%
% Input
%
%  colh - a vector of n column headers from the table
% 
%  rowh - a vector of m row headers from the table
% 
%  values - an (m x n) matrix of values from the table, or a 
% 
%  doprunenans - optional boolean determining if members of the table whose
%    values are nan are removed from the output. Default is false, keeping
%    the nans
%
% Output
%
%  x - a vector of (m x n) values containing the column headings repeated
%    for each column. This may be less than (m x n) if nan values are
%    removed
%
%  y - a vector of (m x n) values containing the corresponding row
%    heading values. This may be less than (m x n) if nan values are
%    removed
%
%  z - a vector of (m x n) vavalueslues containing the table data values
%    for each combination of x and y from the headings. This may be less
%    than (m x n) if nan values are removed
%

% Created by Richard Crozier 2013

    if nargin < 4
        doprunenans = false;
    end
    
    y = repmat(rowh(:), numel(colh), 1);
    x = reshape(repmat(colh(:)', numel(rowh), 1), [], 1);
    z = values(:);

    if doprunenans
        
        nonnaninds = find(~isnan(z));

        x = x(nonnaninds);
        y = y(nonnaninds);
        z = z(nonnaninds);
    
    end

end