function changecount = regexprepfile(filename, exp, repstr, display, dryrun)
% Replace string in one or more files using regular expressions.
%
% Syntax
%
% regexprepfile(filename, exp, repstr)
% regexprepfile(..., display)
% regexprepfile(..., dryrun)
%
% Description
% 
% regexprepfile(filename,exp,repstr) performs a regular expression
% replacement using the pattern exp in the file filename with the string
% repstr. exp, repstr may also be cell arrays of strings of the same size,
% in which case the replacement is performed for each pair by performing a
% regexprep using corresponding elements of the inputs. Alternatively
% repstr may be a string and exp a cell array, in this case the single
% string repstr will replace all the patterns in exp
% 
% filename may also be a cell array of strings, each a file name, in which
% case the replacement is perfromed for each file.
%
% Input
%
%  filename - string or cell array of strings being file paths to be
%   checked by regexprepfile for the expression exp
%
%  exp - regular expression replacement using the pattern exp, in the same
%   format as required by regexprep.m
%
%  repstr - string with which to replace any occurances of the pattern in
%   exp, in all files in 'filename'
%
%  display - optional T/F flag determining whether to display output about
%   the replacement process. Default is false
% 
%  display - optional T/F flag determining whether to do a dry run test
%   where the replacement is not actually performed. Default is false
%
% Output
% 
%  changecount - number 
%
%
% See also regexprep, strrepfile
%

    if nargin < 4
        display = false;
    end
    
    if nargin < 5
        dryrun = false;
    end
    
    if ~iscellstr(filename) && ischar(filename)
        filename = {filename};
    elseif ~iscellstr(filename)
        error('filename must be a string or cell array of strings containing file locations.')
    end
    
    changecount = zeros (size (filename));
    
    for i = 1:numel(filename)
        
        if exist(filename{i}, 'file') ~= 2
            error('file %s not found', filename{i});
        end
       
        if iscellstr(exp) && iscellstr(repstr) 
            
            if ~samesize(exp, repstr)
                error('If two cell arrays of strings are supplied, all dimensions must be the same size')
            end
            
        elseif iscellstr(exp) && ischar(repstr)
            
            % replicate string in exp to make cell array same size as repstr
            repstr = repmat({repstr}, size(exp));
            
        elseif ischar(exp) && iscellstr(repstr)
            
            % throw an error
            error('exp cannot be a string if repstr is a cell array of strings.')
            
        elseif ischar(exp) && ischar(repstr)
            % convert strings to scalar cell arrays
            exp = {exp};
            repstr = {repstr};
            
        else
            error('exp and repstr must be two cell arrays of strings of the same size, or a combination of cell string array and strings, or two strings')
        end
        
        [pathstr, name, ext] = fileparts(filename{i});
        
        [origfid, msg] = fopen(fullfile(pathstr, [name, ext]), 'r');
        
        % create a temporary file to write the replaced string version of
        % the file to
        tempfilename = [tempname, ext];
        
        [tempfid, msg] = fopen(tempfilename, 'w');

        % check the files were opened properly
        if tempfid == -1 && origfid == -1
            error('error opening both search and temporary file');
        elseif tempfid == -1
            origstatus = fclose(origfid);
            if origstatus == -1
                error('Error opening temporary file, and also now unable to close original file.');
            else
                error('Error opening temporary file.');
            end
        elseif origfid == -1
            error('Error opening search file.');
        end
        
        linenum = 1;
        while 1
            tline = fgets(origfid);
            
            if tline ~= -1
                
                for ii = 1:numel(exp)
                    newtline = regexprep(tline, exp{ii}, repstr{ii});
                end
                
                % print the new line to the file (may be unchanged)
                fprintf(tempfid, '%s', newtline);
                
                if ~strcmp (tline, newtline)
                    changecount(i) = changecount(i) + 1;
                    
                    if display
                        
                        fprintf (1, 'match found on line %d of file %s\n', linenum, filename{i});
                        fprintf (1, 'orig: %s\n', strtrim (tline));
                        fprintf (1, 'new:  %s\n', strtrim (newtline));
                        
                    end
                end
                
            else
                tempstatus = fclose(tempfid);
 
                origstatus = fclose(origfid);
                
                if tempstatus == -1 && origstatus == -1
                    error('error closing both search and temporary file');
                elseif tempstatus == -1
                    error('Error closing temporary file.');
                elseif origstatus == -1
                    error('Error closing search file.');
                end
                
                break;
            end
            
            linenum = linenum + 1;
            
        end
        
        if changecount(i) > 0 && ~dryrun
            % copy over the new file to the original location
            [status,message,messageid] = copyfile(tempfilename, fullfile(pathstr, [name, ext]));
            
            % if there was a problem rethrow the error
            if status == 0
                rethrow(struct('message', message, 'identifier', messageid));
            end
            
            if display
                fprintf (1, '^^^ changes above found in file %s\n', filename{i});
            end
            
        end
        
        try
            % attempt to delete the temporary file
            delete (tempfilename);
        end
        
    end

end