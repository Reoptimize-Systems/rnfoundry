function sendwinrebootsignal(multicoredir)
% creates a parameter file with instruction to reboot a windows machine
%

    dateStr = sprintf('%04d%02d%02d%02d%02d%02d', round(clock));
    
    parameterfile = fullfile(multicoredir, sprintf('parameters_%s_0001.mat', dateStr));
    
    functionHandles = 'mcore.winreboot';
    parameters      = {{multicoredir}};
    nresultargs     = 1;
    
    % Generate a file semaphore to prevent workers attempting to access
    % parameter file as we generate it
    sem = mcore.setfilesemaphore(parameterfile, true, true, 0, 0);

    pause(0.01);
    
    % save the parameter file
    save(parameterfile, 'functionHandles', 'parameters', 'nresultargs'); 
    
    % remove the semaphore so workers can access the file
    mcore.removefilesemaphore2(sem);

end