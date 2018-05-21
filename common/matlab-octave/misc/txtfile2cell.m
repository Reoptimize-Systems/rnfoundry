function [strs, nlines] = txtfile2cell(fname)
% read an entire text file into a cell array of strings
%
% Syntax
%
% [strs, nlines] = txtfile2cell(fname)
%
% Description
%
% txtfile2cell ead an entire text file into a cell array of strings. The
% text file is split on the newline character, so each line in the file is
% in a cell.
%
% Input
%
%  fname - character vector containing th name of the file to read in
%
% Output
%
%  strs - cell array of character vectors containing each line of text from
%   the file
%
%  nlines - the number of lines read in 
%
%
%
% See Also: cellstr2txtfile, cellstr2str
%

    assert (ischar (fname), 'fname must be a character vector');
    
    % try to open the file
    [fid, message] = fopen(fname, 'r');
    
    if fid == -1
        error(message);
    end
    
    nlines = 0;
    firstrun = true;
    tline = [];
    
    % read it in line by line
    while firstrun || ischar(tline)
        
        firstrun = false;
        nlines = nlines + 1;
        tline = fgetl(fid);
        
        if ischar(tline)
            strs{nlines,1} = tline;
        end
        
    end

    fclose(fid);

end