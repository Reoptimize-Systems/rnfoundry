function C = mat2cellstr(A)
% converts a matrix into a cell array of strings of the same dimensions
%
% Syntax
%
% C = mat2cellstr(A)
%
% Input
%
% A - numeric matrix to be converted to a cell array of strings
%
% Output
%
% C - cell array of strings of the same dimensions as A where each string
%   is the corresponding number in A converted to a string (with num2str)
%
% See also: num2str

    C = reshape(cellstr(num2str(A(:))), size(A));

end