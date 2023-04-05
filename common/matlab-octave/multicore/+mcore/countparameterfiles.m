function nparamfiles = countparameterfiles (dirpath)
% count the number of multicore job parameter files
%
% Syntax
%
% nparamfiles = mcore.countparameterfiles (dirpath)
%
%
% Input
%
%  dirpath - directory in which to search for parameter files. If not
%   specified the directory returned by mcore.defaultmulticoredir will be
%   used.
%
% Output
%
%  nparamfiles - number of parameter files
%
%
%
% See Also: 
%

    if nargin < 1
        dirpath = mcore.defaultmulticoredir ();
    end

    nparamfiles = mcore.countmatfilesstartingwith ('parameters_', dirpath);
    
end


