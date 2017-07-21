function [status,result] = cleansystem (cmd, restoreuser)
% equivalent of 'system' but removes Matlab's modifications to the library
% path
%  and restores the users environment if necessary
%
% Syntax
%
% [status,result] = cleansystem (cmd)
% [status,result] = cleansystem (cmd, restoreuser)
%
% Description
%
% cleansystem is identical to 'system' but ensures the command runs in a
% clean environment. On linux, before running the command, LD_LIBRARY_PATH
% is removed (using unset) the user's bashrc is rerun to restore the
% user environment.
%
% Input
%
%   cmd - command to be executed on the host computer
%
%   restoreuser - optional logical flag determining whether to restore the
%     user's environment. Only relevant for non-windows. This is achieved
%     by running source ~/.bashrc before running the command (after
%     unsetting LD_LIBRARY_PATH). Default is true if not supplied.
%
% Output
%
%  see help for system.m
% 
% See also: system.m

    if nargin < 2
        restoreuser = true;
    end

    if ispc
        precmd = '';
    else
        precmd = 'unset LD_LIBRARY_PATH ; ';
        
        if restoreuser
            precmd = [precmd, 'source ~/.bashrc ;'];
        end
        
    end
    
    cmdline = sprintf ('%s %s', precmd, cmd);

	[status,result] = system ( cmdline );
            
end