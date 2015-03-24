function cellstrs = indentstrs (cellstrs, indent)
% prepend an indent (or other string) to all strings in a cell array of
% strings
%
% Syntax
%
% cellstrs = indentstrs (cellstrs)
% cellstrs = indentstrs (cellstrs, indent)
%
% Input
%
%   cellstrs - cell array of strings. Each string will have an indent
%     prepended. By defaut the indent is four spaces ('    ').
%
%   indent - options alternative indent for the strings, this can in fact
%     be any string
%
% Output
%
%   cellstrs - the input cell array of strings but with the desired indent
%     prepended to every string
%
%

    if nargin < 2
        indent = '    ';
    end
    
    if ~ischar (indent)
        error ('indent must be a char array (string) to prepend to the supplied strings in cellstrs');
    end
    
    if ~iscellstr (cellstrs)
        error ('cellstrs must be a cell array of strings');
    end
    
    cellstrs = cellfun (@(x) sprintf('%s%s', indent, x) , ...
                    cellstrs, 'UniformOutput', false);

end