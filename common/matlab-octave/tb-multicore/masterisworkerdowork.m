function [mcstate, resultCell] = masterisworkerdowork(mcstate, settings, evalfcns, fcnparams, resultCell)
% looks for results and performs jobs if they have timed out in a multicore
% process
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % work down the file list from top %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if settings.debugMode
        fprintf(1,'********** 1. Working from top to bottom (file nr %d)\n', mcstate.lastFileNrMaster);
    end

    % curFileNr = mcstate.lastFileNrMaster; % for simpler copy&paste
    % construct the parameter, result and working file names from
    % the templates and the current master file number
    paramfname = strrep(mcstate.parameterFileNameTemplate, 'XX', sprintf('%04d', mcstate.lastFileNrMaster));
    resultfname    = strrep(paramfname, 'parameters', 'result' );
    workfname   = strrep(paramfname, 'parameters', 'working');

    if mcstate.multicoreCancelled, return; end

    if settings.debugMode, t1 = mbtime(); end

    % get the indices of the set of parameters in the parameter
    % file we are about to examine
    mcstate.parIndex = ((mcstate.lastFileNrMaster-1)*settings.nrOfEvalsAtOnce+1) : min(mcstate.lastFileNrMaster*settings.nrOfEvalsAtOnce, mcstate.nrOfEvals);

    % create a semaphore for the file so the slaves will not try to
    % access it while the master works on it
    sem = setfilesemaphore(paramfname, settings.showWarnings, settings.debugMode, 0, mcstate.dirTimeDiff);

    if settings.debugMode, mcstate.setTime = mcstate.setTime + mbtime() - t1; end

    % check if the parameter file is still in the multicore
    % directory
    parameterFileExisting = existfile(paramfname);

    if parameterFileExisting
        % If the parameter file is still in the directory, no other
        % process has started working on that job --> Remove
        % parameter file, so that no slave process can load it. The
        % master will do the current job.
        mbdelete2(paramfname, settings.showWarnings);
        if settings.debugMode
            fprintf(1,'Parameter file nr %d deleted by master as it is going to do the job.\n', ...
                mcstate.lastFileNrMaster);
        end
    end

    % check if the current parameter set was evaluated before by a
    % slave process
    resultLoaded = false;

    if parameterFileExisting
        % If the master has taken the parameter file, there is no
        % need to check for a result. Semaphore will be removed
        % below.
        if settings.debugMode
            fprintf(1,'Not checking for result because parameter file nr %d was existing.\n', mcstate.lastFileNrMaster);
        end

    else
        % The master has not taken the parameter file, and it also
        % is not present in the multicore directory, this means a
        % slave must have taken it. Now we will look for a result
        % file, or if a result file can not be found, wait for the
        % job to time out.
        [mcstate, resultCell, sem, resultLoaded] = ...
            mcoremasterwaitforresult(mcstate, settings, resultCell, ...
                sem, paramfname, resultfname, workfname);

    end % if parameterFileExisting

    if settings.debugMode, t1 = mbtime(); end

    % remove file semaphore
    removefilesemaphore2(sem);

    if settings.debugMode, mcstate.removeTime = mcstate.removeTime + mbtime() - t1; end

    if mcstate.multicoreCancelled, return; end

    if ~resultLoaded
        % have the master evaluate function if the result could not
        % be loaded for some reason. The most common reason will be
        % that the parameter file was still present, as it has not
        % been evaluated yet by a slave
        [mcstate, resultCell] = mcoremastereval(mcstate, settings, evalfcns, fcnparams, resultCell);
    end

    % move to next file
    mcstate.lastFileNrMaster = mcstate.lastFileNrMaster + 1;
    if settings.debugMode
        fprintf(1,'Moving to next file (%d -> %d).\n', mcstate.lastFileNrMaster, mcstate.lastFileNrMaster + 1);
    end
            
end