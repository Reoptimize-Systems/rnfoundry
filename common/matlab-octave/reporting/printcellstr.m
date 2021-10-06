function printcellstr(C, fid)
% print all the strings in a cell string array to the command line, or a
% file, each on a new line
%
% Syntax
%
% printcellstr(C)
% printcellstr(C,fid)
%
% Input
%
%   C - a cell array of strings to be printed
%
%   fid - optional file id to print to, will print to comamnd line if not
%     supplied.
%

    if ~iscellstr(C)
        error('input must be a cell array of strings');
    end
    
    if nargin < 2
        fid = 1;
    end

    for ind = 1:numel(C)
        
       fprintf(fid, '%s\n', C{ind});
       
    end

end