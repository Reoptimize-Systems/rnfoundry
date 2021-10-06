function [mcstate, resultCell, sem, resultLoaded] = mcoremasterwaitforresult(mcstate, settings, resultCell, sem, paramfname, resultfname, workfname)
% looks for results files in a multicore process and waits for jobs to time
% out if they are not present
%
%

    mcstate.curPauseTime = settings.initCheckPauseTime;
    firstRun = true;
    
    while 1
        % this while-loop will be left if result was loaded or
        % job timed out

        if firstRun
            % use the semaphore passed in
            firstRun = false;
        else
            % set semaphore
            if settings.debugMode, t1 = mbtime(); end

            sem = setfilesemaphore(paramfname, settings.showWarnings, settings.debugMode, 0, mcstate.dirTimeDiff);

            if settings.debugMode, mcstate.setTime = mcstate.setTime + mbtime() - t1; end
        end

        % Check if the result is available. The semaphore file of the
        % parameter file is used for the following file accesses of the
        % result file.
        if existfile(resultfname)
            
            % attempt to load the result from result file
            [result, resultLoaded] = loadResultFile2(resultfname, settings.showWarnings);

            if resultLoaded && isempty(result) && ~settings.expectEmptyResults
                % if result is empty and we're not expecting that,
                % something has gone wrong and we should mark the result
                % load as failed and do the work
                resultLoaded = false;
                unexpectedEmptyResult = true;

                if settings.debugMode
                    fprintf(1,'Result file nr %d has empty result matrix.\n', mcstate.lastFileNrMaster);
                end
%             else
%                 unexpectedEmptyResult = false;
            end

            if resultLoaded && settings.debugMode
                fprintf(1,'Result file nr %d loaded.\n', mcstate.lastFileNrMaster);
            end

        else
            % the result file is not present, so flag it as not loaded
            resultLoaded = false;
            unexpectedEmptyResult = false;
            if settings.debugMode
                fprintf(1,'Result file nr %d was not found.\n', mcstate.lastFileNrMaster);
            end
        end

        if resultLoaded

            % Save result in the result cell array
            resultCell(mcstate.parIndex) = result;

            % increment the count of the number of slaves that have been
            % evaluated
            mcstate.nrOfFilesSlaves = mcstate.nrOfFilesSlaves + 1;

            % Leave while-loop immediately after result was loaded. The
            % calling function must handle removal of the file semaphore.
            break;

        else % if resultLoaded
            % The result was not loaded, so check if the processing time
            % (current time minus time stamp of working file) exceeds the
            % maximum wait time. Still using the semaphore of the parameter
            % file from above.
            if existfile(workfname)

                if settings.debugMode
                    fprintf(1,'Master found working file nr %d.\n', mcstate.lastFileNrMaster);
                end

                % Check if the job timed out by getting the time when the
                % slave started working on that file. If the job has timed
                % out, the master will do the job, .
                jobTimedOut = mbtime() - getfiledate(workfname) * 86400 > mcstate.maxMasterWaitTime;
            else
                % No working file has been found. The loop is immediately
                % left and the master will do the job.
                if settings.showWarnings
                    fprintf(1,'Warning: Working file %s not found.\n', workfname);
                end
                jobTimedOut = true;
            end

            if jobTimedOut || unexpectedEmptyResult
                if settings.debugMode && jobTimedOut
                    fprintf(1,'Job nr %d has timed out.\n', mcstate.lastFileNrMaster);
                elseif settings.debugMode && unexpectedEmptyResult
                    fprintf(1,'Job nr %d returned unexpected empty result.\n', mcstate.lastFileNrMaster);
                end

                % Reset empty result flag
%                 unexpectedEmptyResult = false;

                % As the slave process seems to be dead or too slow, or the
                % result was unexpectedly empty, the master will do the job
                % itself (semaphore will be removed by the calling
                % function).
                break;
            else
                if settings.debugMode
                    fprintf(1,'Job nr %d has NOT timed out.\n', mcstate.lastFileNrMaster);
                end
            end

            % If the job did not time out, remove semaphore and wait a
            % moment before checking again
            if settings.debugMode, t1 = mbtime(); end

            removefilesemaphore2(sem);

            if settings.debugMode
                mcstate.removeTime = mcstate.removeTime + mbtime() - t1; 
                fprintf(1,'Waiting for result (file nr %d).\n', mcstate.lastFileNrMaster);
            end
            
            pause(mcstate.curPauseTime);

            % make the wait time bigger up to a given maximum value on
            % each loop
            mcstate.curPauseTime = min(settings.maxCheckPauseTime, mcstate.curPauseTime + settings.initCheckPauseTime);

            % run the multicore monitor function before beginning
            % evaluation if it has been supplied
            if ~isempty(settings.monitorFunction)
                mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData, mcstate.monitorstate);
            end

        end % if resultLoaded
        
    end % while 1

end