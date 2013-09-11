function [strs, nlines] = txtfile2cell(fname)
% read an entire text file into a cell array of strings
%
% Syntax
%
% [strs, nlines] = txtfile2cell(fname)
%
% 

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