function [newx, newy, newz] = makevars2dtable23dplot(x, y, data)
% converts data formatted as a 2d table (row heading, column headings and
% data) to a format suitible for the standard 3d plotting functions
%

    % do some input checking
    if ~isvector(x) || numel(x) ~= size(data, 1)
        error('x must be a vector of the same number of elements as the number of rows in data')
    elseif ~isvector(y) || numel(y) ~= size(data, 2)
        error('y must be a vector of the same number of elements as the number of columns in data')
    end
    
    newx = repmat(x(:), numel(y), 1);
    newy = reshape(repmat(y(:)', numel(x), 1), [], 1);
    newz = data(:);

end