function nslavefiles = countslaveIDfiles (dirpath)
% count the number of slave ID fles in a directory
%
% Syntax
%
% nslavefiles = mcore.countslaveIDfiles ()
% nslavefiles = mcore.countslaveIDfiles (dirpath)
%
% Input
%
%  dirpath - directory in which to search for slave ID files. If not
%   specified the directory returned by mcore.defaultmulticoredir will be
%   used.
%
% Output
%
%  nslavefiles - the number of slave ID files found in the specified
%  directory
%
%
%
% See Also: 
%

    if nargin < 1
        dirpath = mcore.defaultmulticoredir ();
    end

    nslavefiles = mcore.countmatfilesstartingwith ('slaveID_', dirpath);

end


