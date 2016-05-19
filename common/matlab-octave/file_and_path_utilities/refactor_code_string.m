function refactor_code_string (expression, replace, topdir, stripvc)
% refactor a code string, changing all occurances of the string name in the
% path to the new string
%
% Syntax
%
% refactor_fcn_name(expression, replace)
% refactor_fcn_name(..., topdir)
% refactor_fcn_name(..., strinvc)
%
% Description
%
% refactor_code_string(expression, replace) finds all uses of expression
% and replaces it with the string in newstr in all .m files in the entire
% matlab path. 'expression' is a regular expression as understood by the
% regexp family of functions. The expression is applied to each line of
% all mfiles searched, and not to the entire file as sngle string.
%
% refactor_fcn_name(expression, replace, topdir) performs the same action
% but searching the folder provided in 'topdir' and all it's subdirectories
% rather than the entire Matlab path. If topdir is the string 'allbutroot'
% all direcotries on the path will be searched, which are not
% sibdirectories of the matlab root path (as returned by matlabroot.m). In
% practice this means built-in toolbox directories will be ignored.
%
% refactor_fcn_name(expression, replace, topdir, strinvc) performs the same
% action but the 'strinvc' flag determines whether the directories used by
% common version control systems will be ignored. This currently includes
% directories starting with '.svn', '.hg', '.git', '.cvs'. Default is true,
% so these directories will be ignored.
%
%
% See also: regexprepfile
%

% Copyright Richard Crozier 2013-2015


    if ~ischar (expression)
        error ('fcnname must be a string');
    end
    
    if ~ischar (replace)
        error ('newfcnname must be a string');
    end
    
    if nargin < 3 || isempty (topdir)
        thepath = path2cell (path);
    else
        if ~ischar (topdir)
            error ('topdir must be a string');
        end
        if strcmpi (topdir, 'allbutroot')
            thepath = path2cell (path);
            thepath(strncmpi (thepath, matlabroot, numel (matlabroot))) = [];
        else
            if exist (topdir, 'file') ~= 7
                error ('supplied directory name does not exist.')
            end
            thepath = path2cell (genpath (topdir));
        end
    end
    
    if nargin < 4
        stripvc = true;
    end
    
    if stripvc
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
            regexprepfile(fullfile(thepath{indi}, mfiles(indii).name), expression, replace);
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

