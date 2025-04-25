function usages = searchfileforexp(filename, exp)
% search a file line by line for the regular expression exp 
%
% Syntax
%
% usages = searchfileforexp (filename, exp)
%
% Description
%
% Searches a file line-by-line for the regular expression provided in exp. 
%
% Input
%
%  filename - path to file to be searched
%
%  exp - regular expression
%
% Output
%
%  usages - (n x 5) cell array where each row contains information about an
%   expression match found in the file, i.e:
%   { name of the file, 
%     line number match was found,
%     copy of the full line from the file where the match was found,
%     start indices of matches as returned by the regexp function,
%     end indices of vectors as returned by the regexp function }
%
% See Also: regexp, regexprepfile
%

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