function C = str2cellstr (S)
% split string or character vector on newlines and return as a cell aray
%
% Syntax
%
% S = str2cellstr (C)
%
% Description
%
% str2cellstr splits a string or character vector on newlines and return
% the resulting lines as a string array, or a cell array of character
% vectors. The newline characters are removed.
%
% Input
%
%  S - string object or character vector
%
% Output
%
%  C - string array or cell array of character vectors containing each line
%   of text from the input.
%
%
%
% See Also: 
%

    C = regexp (S, '.*', 'match', 'dotexceptnewline')';

end