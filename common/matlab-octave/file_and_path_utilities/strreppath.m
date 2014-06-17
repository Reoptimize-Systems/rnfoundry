function strreppath(S1, S2, topdir)
% replace strings in all files in the matlab path or a sub-path recursively
%
% Syntax
%
% strreppath(S1, S2)
% strreppath(S1, S2, topdir)
%
% Description
% 
% strreppath(S1,S2) replaces all occurrences of the string S1 in all files
% on the matlab path with the string S2. S1, S2 may also be cell arrays of
% strings of the same size, in which case the replacement is performed for
% each pair by performing a STRREP using corresponding elements of the
% inputs. Alternatively S2 may be a string and S1 a cell array, in this
% case the single string S2 will replace all the strings in S1.
%
% strreppath(S1,S2,topdir) replaces the strings the directory supplied in
% 'topdir', and all subdirectories recursively.
%
%

    if nargin < 3 || isempty(topdir)
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
    
    for indi = 1:numel(thepath)

        mfiles = dir([ thepath{indi}, '\*.m' ]);

        for indii = 1:numel(mfiles)
            % replace string in file using regexprepfile, with a regex
            % looking for word boundaries to avoid functions with similar
            % substrings
            strrepfile(fullfile(thepath{indi}, mfiles(indii).name), S1, S2);
        end

    end

end