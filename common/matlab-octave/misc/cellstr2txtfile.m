function cellstr2txtfile (file, strs)
% Write one or more strings or character vectors to a text file
%
% Description
%
% write all strings in an array, or all character vectors in a cell array
% to a text file, each string or character vector will be separated by a
% newline character
%
% Syntax
%
% cellstr2txtfile (file, strs)
%

    [fid, message] = fopen (file, 'w');
    
    if fid == -1
        error (message);
    end
    
    CC = onCleanup (@() fclose (fid));
    
    if isstring (strs)
        
        for ind = 1:numel (strs)

            fprintf (fid, '%s\n', strs(ind));

        end
        
    elseif iscellstr (strs)
        
        for ind = 1:numel (strs)

            fprintf (fid, '%s\n', strs{ind});

        end
        
    else
        error ('strs must be a string array or a cell array of character vectors');
    end
    
end