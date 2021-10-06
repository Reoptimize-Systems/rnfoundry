function strreppath(S1, S2, varargin)
% replace strings in all files in the matlab path or a sub-path recursively
%
% Syntax
%
% strreppath(S1, S2)
% strreppath(S1, S2, 'Parameter', Value)
%
% Description
% 
% strreppath(S1,S2) replaces all occurrences of the string S1 in all mfiles
% (with a '.m' extension) on the matlab path with the string S2. S1, S2 may
% also be cell arrays of strings of the same size, in which case the
% replacement is performed for each pair by performing a STRREP using
% corresponding elements of the inputs. Alternatively S2 may be a string
% and S1 a cell array, in this case the single string S2 will replace all
% the strings in S1.
%
% Additional arguments can also be supplied as parameter-value pairs.
%
% strreppath(S1,S2,'topdir',path) replaces the strings the directory
% supplied in path, and all its subdirectories recursively.
%
% strreppath(S1,S2,'ext',ext) replaces the strings the directory supplied
% in 'topdir', and all subdirectories recursively, but in files with the
% extension supplied in ext, instead of files with a '.m' extension. The
% extension should be supplied without the '.', e.g. just 'm', or 'txt', or
% 'cpp'.
%

    Inputs.ext = 'm';
    Inputs.topdir = '';
    
    Inputs = parseoptions (Inputs, varargin);
    
    if isempty(Inputs.topdir)
        thepath = path2cell(path);
    else
        if ~ischar(Inputs.topdir)
            error('topdir must be a string');
        end
        if exist(Inputs.topdir, 'file') ~= 7
            error('supplied directory name does not exist.')
        end
        thepath = path2cell(genpath(Inputs.topdir));
    end
    
    for indi = 1:numel(thepath)

        mfiles = dir(fullfile (thepath{indi}, ['*.', Inputs.ext ]));

        for indii = 1:numel(mfiles)
            % replace string(s) in file using strrepfile
            strrepfile(fullfile(thepath{indi}, mfiles(indii).name), S1, S2);
        end

    end

end
