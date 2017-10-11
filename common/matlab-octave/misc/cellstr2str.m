function s = cellstr2str(C)
% converts a cell array of strings to a single string with newlines
%
% Syntax
%
%
% Description
%
% converts a cell array of strings to a single string with newline
% characters between the string in each cell.
%
% Input
%
%  C - cell string array
%
% Output
%
%  s - string composed of the strings in each cell of the input C cell
%   array concatenated with a newline character inserted between each
%   substring.
%
% See Also: cellstr2txtfile, printcellstr
%

    if ~iscellstr (C)
        error ('C must a cell array of strings');
    end
    
    s = sprintf ('%s\n', C{:});

end