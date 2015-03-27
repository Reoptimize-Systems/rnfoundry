function top = topdir (path)
% return the top level directory in a path
%
% Syntax
%
% top = topdir (path)
% 
% 

    C = strsplit (path, filesep ());
    
    top = C{end};

end