function usages = find_in_files(exp, topdir, dodisplay)
% search for useages of an expression in mfiles in a directory and all its
% subfolders recursively, or the entire path
%
% Syntax
%
% find_in_files(exp)
% find_in_files(exp, topdir, dodisplay)
%
% Description
%
% find_in_files(exp) searches for the regular expression exp in all mfiles
% in the entire path. Infomation on the usages is returned on the command
% line, and in the cell array usages, described below.
%
% find_in_files(exp, topdir) searches for the regular expression exp in all
% mfiles in the folder topdir and all its subfolders recursively.
% Infomation on the usages is returned on the command line, and in the cell
% array usages, described below. If topdir is empty, the entire path is
% searched.
%
% find_in_files(exp, topdir, dodisplay) performs as
% find_in_files(exp,topdir) but allows control over whether results are
% displayed on the command line. If dodisplay evaluates to true, results
% are displayed, otherwise they are not.
%
% Output
%
% usages - cell array containing information on the usages found, will be
%   empty if no usages are found. Otherwise usages will be an (n x 4) cell
%   array where n is the number of usages found. Each row of the cell array
%   will contain the file name in which the usage was found, a copy of the
%   line from the file in which it was found, the start position of the
%   regular expression match in the line, and the end position of the match
%   in the line in that order.
%
% 
% See also: regexprepfile
%

% Created by Richard Crozier 2013

    if ~ischar(exp)
        error('fcnname must be a string');
    end
    
    if nargin < 2 || isempty(topdir)
        thepath = path2cell(path);
    else
        if ~ischar(topdir)
            error('topdir must be a string');
        end
        if exist(topdir, 'file') ~= 7
            error('supplied directory name does not exist.')
        end
        thepath = path2cell(genpath(topdir));
    end
    
    if nargin < 3
        dodisplay = true;
    end
    
    usages = {};

    for indi = 1:numel(thepath)

        mfiles = dir(fullfile(thepath{indi}, '*.m' ));

        for indii = 1:numel(mfiles)
            % find expression in file using regexp
            usages = [ usages; searchfile(fullfile(thepath{indi}, mfiles(indii).name), exp) ];
        end

    end
    
    if dodisplay
        
        nusages = size(usages,1);

        fprintf(1, 'Found %d occurances of pattern %s in files:\n', nusages, exp); 
        if nusages > 0
            for i = 1:nusages
                fprintf(1, '%s  | line: %d  |  %s\n', usages{i,1}, usages{i,2}, usages{i,3});
            end
        end
        
    end

end


function usages = searchfile(filename, exp)
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