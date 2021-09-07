function usages = searchfileforexp(filename, exp)
% search a file line by line for the regular expression exp 

    usages = {};
    
    [fid, msg] = fopen(filename);
    
    [pathstr, name, ext] = fileparts(filename);
        
    linenum = 0;
    
    while 1
    
        linenum = linenum + 1;
        
        tline = fgets(fid);

        if tline ~= -1

            [matchstart,...
                matchend,...
                tokenindices,...
                matchstring,...
                tokenstring,...
                tokenname ] = regexp(tline, exp);

            if ~isempty(matchstart)
                % string newlines 
                tline = strrep(tline, sprintf('\r\n'), ''); % windows line ending
                tline = strrep(tline, sprintf('\n'), ''); % unix line ending
                usages = [ usages; {[name, ext], linenum, tline, matchstart, matchend} ];
            end
        else
            break;
        end
    end
    
    fclose(fid);

end