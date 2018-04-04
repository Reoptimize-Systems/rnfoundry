function refactor_code_string (expression, replace, varargin)
% refactor a code string, changing all occurances of the string name in the
% path to the new string
%
% Syntax
%
% replace_fcnname (expression, replace)
% replace_fcnname (..., 'Parameter', Value)
%
% Description
%
% refactor_code_string(expression, replace) finds all uses of expression
% and replaces it with the string in newstr in all .m files in the entire
% matlab path. 'expression' is a regular expression as understood by the
% regexp family of functions. The expression is applied to each line of
% all mfiles searched, and not to the entire file as sngle string. An
% alternative list of directories can be searched by using the TopDir
% option. By default version control directoies 
%
% Input
%
%  expression - character vector containing regular expression to search
%   for.
%
%  replace - character vector containing the string with which to replace
%   any instances of expression which are found
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'IgnoreVCDirs' - true/false flag determining whether common source
%    control system directories will be ignored. This currently includes
%    directories starting with '.svn', '.hg', '.git', '.cvs'. Default is
%    true, so these directories will be ignored.
%
%  'TopDir' - character vector used to control what paths are searched,
%    instead of searching the entire Matlab path. If TopDir is specified,
%    the folder provided in TopDir' and all it's subdirectories are
%    searched rather than the entire Matlab path. If topdir is the string
%    'allbutroot' all direcotries on the path will be searched, which are
%    not subdirectories of the matlab root path (as returned by
%    matlabroot.m). In practice this means built-in toolbox directories
%    will be ignored.
%
%  'DryRun' - true/false flag indicating whether to perform a dry run only,
%    where no changes are really made to the files. The output of what
%    would be changed is still displayed on the command line.
%
% See Also: refactor_fcn_name.m, regexprepfile.m
%

% Copyright Richard Crozier 2013-2018

    options.IgnoreVCDirs = true;
    options.TopDir = [];
    options.DryRun = false;
    
    options = parse_pv_pairs (options, varargin);
    
    check.isLogicalScalar (options.IgnoreVCDirs, true, 'IgnoreVCDirs');
    check.isLogicalScalar (options.DryRun, true, 'DryRun');

    if ~ischar (expression)
        error ('fcnname must be a string');
    end
    
    if ~ischar (replace)
        error ('newfcnname must be a string');
    end
    
    if isempty (options.TopDir)
        thepath = path2cell (path);
    else
        if ~ischar (options.TopDir)
            error ('topdir must be a string');
        end
        if strcmpi (options.TopDir, 'allbutroot')
            thepath = path2cell (path);
            thepath(strncmpi (thepath, matlabroot, numel (matlabroot))) = [];
        else
            if exist (options.TopDir, 'file') ~= 7
                error ('supplied directory name does not exist.')
            end
            thepath = path2cell (genpath (options.TopDir));
        end
    end
    
    % +package directories are not on the path, but we want to search them
    % too, so we add them here
    pkgpaths = {};
    for ind = 1:numel(thepath)
        
        pkgpaths = [ pkgpaths; addpkgpaths(thepath{ind}) ];
        
    end
    
    thepath = [thepath; pkgpaths];

    if nargin < 4
        stripvc = true;
    end
    
    if options.IgnoreVCDirs
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
            % replace string in file using regexprepfile
            regexprepfile(fullfile(thepath{indi}, mfiles(indii).name), expression, replace, true, options.DryRun);
        end

    end

end

function newpaths = addpkgpaths (testpath)

    pkgpaths = dir (fullfile (testpath, '+*'));
    
    newpaths = {};

    for ind = 1:numel (pkgpaths)
        
        if pkgpaths(ind).isdir
            
            newpaths = [ newpaths; fullfile(pkgpaths(ind).folder, pkgpaths(ind).name) ];
            
            newpaths = [ newpaths; addpkgpaths(newpaths{end}) ];
        
        end
        
    end

end


function pathlist = path2cell(pathstr)
%PATH2CELL Convert search path to cell array.
%
%   PATH2CELL returns a cell array where each element is a directory
%   in the search path.
%
%   PATH2CELL(MYPATH) converts MYPATH rather than the search path.
%
%   Empty directories are removed, so under UNIX, PATH2CELL('foo::bar')
%   returns {'foo', 'bar'} rather than {'foo', '', 'bar'}.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-10-18 21:23:19 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
%    error(nargchk(0, 1, nargin));

   % use MATLAB's search path if no input path is given
   if ~nargin
      pathstr = path;
   end

   k = find(pathstr == pathsep);            % find path separators
   k = [0 k length(pathstr)+1];             % find directory boundaries
   ndirs = length(k)-1;                     % number of directories
   pathlist = cell(0);                      % initialize output argument
   for i = 1 : ndirs
      dir = pathstr(k(i)+1 : k(i+1)-1);     % get i'th directory
      if ~isempty(dir)                      % if non-empty...
         pathlist{end+1,1} = dir;           % ...append to list
      end
   end
end

