function refactor_fcn_name(fcnname, newfcnname, topdir, domove, varargin)
% refactor the name of a function, changing all references to the function
% name in the path and moving the function file to the new named file
%
% Syntax
%
% refactor_fcn_name(fcnname, newfcnname)
% refactor_fcn_name(fcnname, newfcnname, topdir)
% refactor_fcn_name(fcnname, newfcnname, topdir, domove)
%
% Description
%
% refactor_fcn_name(fcnname, newfcnname) finds all uses of the function
% name in fcnname and replaces it with the string in newfcnname in the
% entire matlab path. The function must be function in an m-file on the
% Matlab path. By default the mfile containing the fuction is also then
% renamed to the new name.
%
% refactor_fcn_name(fcnname, newfcnname, topdir) performs the same action
% but searching the folder provided in 'topdir' and all it's subdirectories
% rather than the entire Matlab path. If topdir is the string 'allbutroot'
% all direcotries on the path will be searched, which are not
% sibdirectories of the matlab root path (as returned by matlabroot.m). In
% practice this means built-in toolbox directories will be ignored.
%
% refactor_fcn_name(fcnname, newfcnname, topdir, domove) performs the same
% action but the 'domove' flag determines whether the function file is
% actually moved or not. If this evaluates to true the file is moved to the
% mfile with the new name, if false nothing is moved.
%
%
% See also: regexprepfile
%

% Copyright Richard Crozier 2013-2015

    options.StripVersionControl = true;
    options.AllowDotBeforeName = true;
    options.WarnIfFcnNotOnPath = true;
    
    options = parse_pv_pairs (options, varargin);

    if ~ischar (fcnname)
        error ('fcnname must be a string');
    end
    
    if ~ischar (newfcnname)
        error ('newfcnname must be a string');
    end
    
    % refactor_fcn_name
    loc = which (fcnname);
    
    if isempty (loc) && options.WarnIfFcnNotOnPath
        warning('refactor_fcn_name:notonpath', 'fcnname is not a function or class on the path')
    end
    
    [pathstr, ~, ext] = fileparts (which (fcnname));
    
    if isempty(strcmpi(ext, '.m'))
       warning('function is not an m-file on the path.') 
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
                regexprepfile(fullfile(thepath{indi}, mfiles(indii).name), ['(?<=^|\W)', fcnname, '(?=((\s*\()|(\s*;)|(\s*$)))'], newfcnname);
            else
                regexprepfile(fullfile(thepath{indi}, mfiles(indii).name), ['(?<=^|\W)(?<![.])', fcnname, '(?=((\s*\()|(\s*;)|(\s*$)))'], newfcnname);
            end
        end

    end

    if domove
        % finally move the function file to the new m-file with the correct
        % name, 
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
            movefile(fullfile(pathstr, [fcnname, ext]), fullfile(pathstr, [newfcnname, ext]));  
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

