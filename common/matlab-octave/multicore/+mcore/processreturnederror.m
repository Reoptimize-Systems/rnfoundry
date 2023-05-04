function settings = processreturnederror(settings, errorstruct, i, parameterCell, evalfcn)

    if strcmp(errorstruct.identifier, 'MATLAB:nomem')
        
        % quit all the condor slaves as they have less memory at their
        % disposal than manual slaves, and change the settings so no new
        % slaves are launched
        
        mcore.quitallslaves (settings.MulticoreSharedDir);
        
        % remove the monitor function so we don't spawn any new slaves
        settings.monitorFunction = {};
        
        % pause for some time to make sure new slaves have not spawned
        % after quitting the existing ones
        pause (settings.DeletePauseTime);
        
        % call quitallslaves again to quit any new ones that have spawned
        % while we waited
        mcore.quitallslaves (settings.MulticoreSharedDir);
        
        % give them some time to quit
        pause(5);
        
    end
    
    % send an email warning about the error
    mcore.errormail (settings, errorstruct, i, parameterCell, evalfcn);

end