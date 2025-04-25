function refactor_fcn_name(fcnname, newfcnname, topdir, domove, varargin)
% refactor the name of a function, changing all references to the function
% name in the path and moving the function file to the new named file
%
% Syntax
%
% refactor_fcn_name(fcnname, newfcnname)
% refactor_fcn_name(fcnname, newfcnname, topdir)
% refactor_fcn_name(fcnname, newfcnname, topdir, domove)
% refactor_fcn_name(..., 'Parameter', Value)
%
%
% Description
%
% refactor_fcn_name(fcnname, newfcnname) finds all uses of the function
% name in fcnname and replaces it with the string in newfcnname in the
% entire matlab path, except the matlabroot paths (see discussion of topdir
% option below). The function must be function in an m-file on the Matlab
% path. By default the mfile containing the fuction is also then renamed to
% the new name.
%
% refactor_fcn_name(fcnname, newfcnname, topdir) performs the same action
% but searching the folder provided in 'topdir' and all it's
% subdirectories. If topdir is the string  'all', or is empty, the entire
% path will be searched, including matlab's own toolbox directories
% contained in the directory pointed to by matlabroot(). If topdir is the
% string 'allbutroot' all direcotries on the path will be searched, which
% are not sibdirectories of the matlab root path. In practice this means
% built-in toolbox directories will be ignored.
%
% refactor_fcn_name(fcnname, newfcnname, topdir, domove) performs the same
% action but the 'domove' flag determines whether the function file is
% actually moved or not. If this evaluates to true the file is moved to the
% mfile with the new name, if false nothing is moved.
%
%
% Input
%
%  fcnname - name of the function to be replaced, can be part of a package,
%   in which case use the fully qualified name (mypackage.myfunction)
%
%  newfcnname - the new name of the function to be replaced
%
%  topdir - character vector containing either the string 'allbutroot'
%
%  domove - true/false flag indicating whether to move the origninal
%   function file name to the new file name. It will be kept in the same
%   directory.
%
% Addtional arguments may be supplied as parameter-value pairs. The available options are:
%
%  'StripVersionControl' - true/false Flag indicating whether to remove
%    paths containing '.svn', '.hg', '.git', '.cvs'.
%
%  'AllowDotBeforeName' - true/false Flag indicating whether names preceded
%    by a '.' character will be modified. This is really only useful if
%    coercing this code to refactor somehting like a structure variable
%    rather than a function name.
%
%  'WarnIfFcnNotOnPath' - true/false Flag indicating whether to issue a
%    warning if the function to be modified does not appear to be on the
%    Matlab path.
%
%
% See also: regexprepfile
%

    options.StripVersionControl = true;
    options.AllowDotBeforeName = true;
    options.WarnIfFcnNotOnPath = true;
    options.UseGenpathForSubdirs = true;
    
    options = parse_pv_pairs (options, varargin);

    if ~ischar (fcnname)
        error ('fcnname must be a string');
    end
    
    if ~ischar (newfcnname)
        error ('newfcnname must be a string');
    end
    
    % refactor_fcn_name
    fullpath = which (fcnname);
    
    if isempty (fullpath) && options.WarnIfFcnNotOnPath
        warning('refactor_fcn_name:notonpath', 'fcnname is not a function or class on the path')
    end

    [pathstr, ~, ext] = fileparts (fullpath);

    fcnname_dots = strfind (fcnname, '.');

    if isempty (fcnname_dots)
        fcnname_no_package = fcnname;
    else
        fcnname_no_package = fcnname (fcnname_dots(end)+1:end);
    end
    
    if isempty(strcmpi(ext, '.m'))
       warning('function is not an m-file on the path.') 
    end
    
    if nargin < 3 
        topdir = 'allbutroot';
    elseif isempty (topdir)
        topdir = 'all';
    end

    if ~ischar (topdir)
        error ('topdir must be a character vector');
    end

    if strcmpi (topdir, 'all')
        thepath = path2cell (path);
    elseif strcmpi (topdir, 'allbutroot')
        thepath = path2cell (path);
        thepath(strncmpi (thepath, matlabroot, numel (matlabroot))) = [];
    else
        if exist (topdir, 'file') ~= 7
            error ('supplied directory name does not exist.')
        end
        if options.UseGenpathForSubdirs
            thepath = path2cell (genpath (topdir));
        else
            dirlist = dir(fullfile(topdir, '**'));  % get list of files and folders reursively
            dirlist = dirlist(...
                [dirlist.isdir] ...
                    & ~arrayfun(@(x) string(x.name).matches('.'), dirlist)' ...
                    & ~arrayfun(@(x) string(x.name).matches('..'), dirlist)' ...
            );  %remove non directories from list

            thepath = cell(numel(dirlist),1);
            for ind = 1:numel(dirlist)
                thepath{ind} = fullfile(dirlist(ind).folder, dirlist(ind).name);
            end
        end
    end

    if nargin < 4
        domove = true;
    end
    
    if options.StripVersionControl
        rminds = [];
        vclist = {'.svn', '.hg', '.git', '.cvs'};
        for indi = 1:numel(thepath)
            for indii = 1:numel(vclist)
                if ~isempty ( strfind (thepath{indi}, vclist{indii}) )
                    rminds = [rminds, indi];
                    break;
                end
            end
        end
        thepath(rminds) = [];
    end
    
    for indi = 1:numel(thepath)

        mfiles = dir(fullfile (thepath{indi}, '*.m'));

        for indii = 1:numel(mfiles)
            % replace string in file using regexprepfile, with a regex
            % looking for word boundaries to avoid functions with similar
            % substrings
            if options.AllowDotBeforeName
                regexprepfile(...
                    fullfile(thepath{indi}, mfiles(indii).name), ...
                    ['(?<=^|\W)', fcnname, '(?=((\s*\()|(\s*;)|(\s*$)))'], ...
                    newfcnname ...
                );
            else
                regexprepfile(...
                    fullfile(thepath{indi}, mfiles(indii).name), ...
                    ['(?<=^|\W)(?<![.])', fcnname, '(?=((\s*\()|(\s*;)|(\s*$)))'], ...
                    newfcnname ...
                );
            end
        end

    end

    if domove
        % finally move the function file to the new m-file with the correct
        % name, but fix the function file if it's in a package

        new_fcnname_dots = strfind (fcnname, '.');

        if isempty (new_fcnname_dots)
            new_fcnname_no_package = newfcnname;
        else
            new_fcnname_no_package = newfcnname (new_fcnname_dots(end)+1:end);
        end

        if ~isempty (new_fcnname_dots)
            % function was in a package so fix the function line
            strrepfile (fullpath, fcnname_no_package, new_fcnname_no_package);
        end

        if ispc
            % On windows we can't use the built-in movefile function as the
            % file system is not case sensitive, and will refuse to move
            % the file to what it beleives is the same location, so we take
            % the roundabout route
            
            % move to a temporary location
            tmpdest = [tempname, '_', newfcnname, ext];
            movefile(fullfile(pathstr, [fcnname, ext]), tmpdest);
            
            % move it again to where we really want it
            movefile(tmpdest, fullfile(pathstr, [newfcnname, ext]));
        
        else
            % on other systems, just move it!
            movefile(fullfile(pathstr, [fcnname_no_package, ext]), fullfile(pathstr, [new_fcnname_no_package, ext]));  
        end

    end

end
