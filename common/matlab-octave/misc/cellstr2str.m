function s = cellstr2str(C, splitchar)
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
%  splitchar - optional character or string to insert between each
%   character vector in the cell array. Devault is '\n' if not supplied. 
%
% Output
%
%  s - string composed of the strings in each cell of the input C cell
%   array concatenated with a newline character inserted between each
%   substring.
%
% See Also: cellstr2txtfile, printcellstr
%
    if nargin < 2
        splitchar = '\n';
        numsplitchar = 1;
    end
    
    numsplitchar = numel(sprintf (splitchar));

    if ~iscellstr (C)
        error ('C must a cell array of strings');
    end
    
    s = sprintf (['%s', splitchar], C{:});
    
    s((end-numsplitchar+1):end) = [];

end