function filepath = getmfilepath(mfile)
% getmfilepath: gets the directory containing an mfile function
%
% Syntax
%
% filepath = getmfilepath(mfile)
%
% Description
%
% mfile is a string containing the name of the mfile for which the location
% is to be found, the .m extension is optional. A package member may also
% be specified, e.g.
%
%   getmfilepath('package.function')
%
% In fact, getmfilepath will work with any syntax supported by the 'which'
% function.
% 
% Output
%
% filepath is the directory path to the file containing the specified
% function
%
%
% See also: which.m
%
 
    loc = which(mfile);

    if isempty(loc)
        error('UTILS:nofile', 'm-file or class does not appear to exist')
    else
        filepath = fileparts(loc);
    end
    
end