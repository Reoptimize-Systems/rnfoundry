function pattern = mcorewinreboot(directory, timeout)

    if nargin < 2
        timeout = 30;
    end
    
    % get the computer name it should be used in the slave id of any slaves
    % associated with this machine
    [status, pattern] = system('echo %COMPUTERNAME%');
    
    fprintf(1, 'pattern is %s\n', pattern);
    
    % delete all slave files for this machine based on a supplied
    % pattern, e.g. *3224*.mat
    delete(fullfile(directory, ['*', pattern, '*']));
    
    fprintf(1, 'deleting %s\n', ['*', pattern, '*']);
    
    % reboot the computer in timeout seconds, forcing all programs to
    % close and not waiting for the result of the command
    [status] = system(sprintf('shutdown -r -t %d -f &', timeout));
    
    fprintf('issued shutdown command\n')

end