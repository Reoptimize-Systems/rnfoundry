function printcellstr(C)
% print all the strings in a cell string array to the command line, each on
% a new line
%
% Syntax
%
% printcellstr(C)
%
% Input
%
%   C - a cell array of strings to be printed
%

    if ~iscellstr(C)
        error('input must be a cell array of strings');
    end

    for ind = 1:numel(C)
        
       fprintf(1, '%s\n', C{ind});
       
    end

end