function lastFileNrMaster = domasterwork(parameterFileNameTemplate, lastFileNrMaster, lastFileNrSlave, resultCell, functionHandleCell, showWarnings, debugMode, dirtdiff)

            if debugMode
                fprintf(1,'********** 1. Working from top to bottom (file nr %d)\n', lastFileNrMaster);
            end

            curFileNr = lastFileNrMaster; % for simpler copy&paste
            parameterFileName = strrep(parameterFileNameTemplate, 'XX', sprintf('%04d', curFileNr));
            resultFileName    = strrep(parameterFileName, 'parameters', 'result' );
            workingFileName   = strrep(parameterFileName, 'parameters', 'working');
            parIndex = ((curFileNr-1)*nrOfEvalsAtOnce+1) : min(curFileNr*nrOfEvalsAtOnce, nrOfEvals);
% 
%             if multicoreCancelled
%                 multicoreCancel2(lastFileNrMaster, lastFileNrSlave);
%                 return
%             end

            if debugMode, t1 = mbtime; end

            sem = setfilesemaphore(parameterFileName, showWarnings, debugMode, 0, dirtdiff);

            if debugMode, setTime = setTime + mbtime - t1; end

            parameterFileExisting = existfile(parameterFileName);

            if parameterFileExisting
                % If the parameter file is existing, no other process has started
                % working on that job --> Remove parameter file, so that no slave
                % process can load it. The master will do the current job.
                mbdelete2(parameterFileName, showWarnings);
                if debugMode
                    fprintf(1,'Parameter file nr %d deleted by master.\n', curFileNr);
                end
            end

            % check if the current parameter set was evaluated before by a slave process
            resultLoaded = false;
            if parameterFileExisting
                % If the master has taken the parameter file, there is no need to check
                % for a result. Semaphore will be removed below.
                if debugMode
                    fprintf(1,'Not checking for result because parameter file nr %d was existing.\n', curFileNr);
                end

            else
                % Another process has taken the parameter file. This branch is
                % entered if master and slave "meet in the middle", i.e. if a slave
                % has taken the parameter file of the job the master would have done
                % next. In this case, the master will wait until the job was finished
                % by the slave process or until the job has timed out.
                curPauseTime = startPauseTime;
                firstRun = true;
                while 1 
                    % this while-loop will be left if result was loaded or
                    % job timed out
                    
                    if firstRun
                        % use the semaphore generated above
                        firstRun = false;
                    else
                        % set semaphore
                        if debugMode, t1 = mbtime; end
                        
                        sem = setfilesemaphore(parameterFileName, showWarnings, debugMode, 0, dirtdiff);
                        
                        if debugMode, setTime = setTime + mbtime - t1; end
                    end

                    % Check if the result is available. The semaphore file of the
                    % parameter file is used for the following file accesses of the
                    % result file.
                    if existfile(resultFileName)
                        % attempt to load the result from result file
                        [result, resultLoaded] = loadResultFile2(resultFileName, showWarnings);
                        
                        if resultLoaded && isempty(result) && ~settings.expectEmptyResults
                            % if result is empty and we're not expecting
                            % that, something has gone wrong and we should
                            % mark the result load as failed and do the
                            % work
                            resultLoaded = false;
                            unexpectedEmptyResult = true;
                            
                            if debugMode, fprintf(1,'Result file nr %d has empty result matrix.\n', curFileNr); end
                        end
                        
                        if resultLoaded && debugMode, fprintf(1,'Result file nr %d loaded.\n', curFileNr); end
                        
                    else
                        resultLoaded = false;
                        if debugMode, fprintf(1,'Result file nr %d was not found.\n', curFileNr); end
                    end

                    if resultLoaded
                        % Save result
                        resultCell(parIndex) = result;
                        nrOfFilesSlaves = nrOfFilesSlaves + 1;

                        % Update waitbar
%                         multicoreWaitbar('update2', nrOfFiles, nrOfFilesMaster, nrOfFilesSlaves);

                        % Leave while-loop immediately after result was loaded. Semaphore
                        % will be removed below.
                        break
                    end

                    % Check if the processing time (current time minus time stamp of
                    % working file) exceeds the maximum wait time. Still using the
                    % semaphore of the parameter file from above.
                    if existfile(workingFileName)
                        if debugMode, fprintf(1,'Master found working file nr %d.\n', curFileNr); end

                        % Check if the job timed out by getting the time when the slave
                        % started working on that file. If the job has timed out, the
                        % master will do the job, .
                        jobTimedOut = mbtime - getfiledate(workingFileName) * 86400 > maxMasterWaitTime;
                    else
                        % No working file has been found. The loop is immediately left
                        % and the master will do the job.
                        if showWarnings, fprintf(1,'Warning: Working file %s not found.\n', workingFileName); end
                        jobTimedOut = true;
                    end

                    if jobTimedOut || unexpectedEmptyResult
                        if debugMode && jobTimedOut
                            fprintf(1,'Job nr %d has timed out.\n', curFileNr);
                        elseif debugMode && unexpectedEmptyResult
                            fprintf(1,'Job nr %d returned unexpected empty result.\n', curFileNr);
                        end
                        
                        unexpectedEmptyResult = false;
                        
                        % As the slave process seems to be dead or too
                        % slow, or the result was unexpectedly empty, the
                        % master will do the job itself (semaphore will be
                        % removed below).
                        break;
                    else
                        if debugMode, fprintf(1,'Job nr %d has NOT timed out.\n', curFileNr); end
                    end

                    % If the job did not time out, remove semaphore and wait a moment
                    % before checking again
                    if debugMode, t1 = mbtime; end
                    
                    removefilesemaphore2(sem);
                    
                    if debugMode, removeTime = removeTime + mbtime - t1; end

                    if debugMode, fprintf(1,'Waiting for result (file nr %d).\n', curFileNr); end

                    pause(curPauseTime);
                    curPauseTime = min(maxPauseTime, curPauseTime + startPauseTime);
                end % while 1
            end % if parameterFileExisting

            % remove semaphore
            if debugMode, t1 = mbtime; end
            
            removefilesemaphore2(sem);
            
            if debugMode, removeTime = removeTime + mbtime - t1; end

%             if multicoreCancelled
%                 multicoreCancel2(lastFileNrMaster, lastFileNrSlave);
%                 return
%             end

            % evaluate function if the result could not be loaded
            if ~resultLoaded
                
                if debugMode
                    fprintf(1,'Master evaluates job nr %d.\n', curFileNr);
                    t0 = mbtime;
                end
                
                for k = parIndex
                    if debugMode
                        %fprintf(' %d,', k);
                    end
                    
                    if iscell(parameterCell{k})
                        resultCell{k} = feval(functionHandleCell{k}, parameterCell{k}{:});
                    else
                        resultCell{k} = feval(functionHandleCell{k}, parameterCell{k});
                    end

                    if multicoreCancelled
                        multicoreCancel2(lastFileNrMaster, lastFileNrSlave);
                        return
                    end
                end
                
                nrOfFilesMaster = nrOfFilesMaster + 1;

                % Run postprocessing function
                if ~isempty(settings.postProcessHandle)
                    postProcStruct.state               = 'after master evaluation'; % no copy & paste here!!
                    postProcStruct.lastFileNrReady     = lastFileNrMaster;          % no copy & paste here!!
                    postProcStruct.lastFileNrMaster    = lastFileNrMaster;
                    postProcStruct.lastFileNrSlave     = lastFileNrSlave;
                    postProcStruct.nrOfFilesMaster     = nrOfFilesMaster;
                    postProcStruct.nrOfFilesSlaves     = nrOfFilesSlaves;
                    postProcStruct.resultCell          = resultCell;
                    postProcStruct.parIndex            = parIndex;
                    feval(settings.postProcessHandle, postProcStruct);
                end

                % Update waitbar
%                 multicoreWaitbar('update2', nrOfFiles, nrOfFilesMaster, nrOfFilesSlaves);

                if debugMode
                    fprintf(1,'Master finished job nr %d in %.2f seconds.\n', curFileNr, mbtime - t0);
                end
            end

            % move to next file
            lastFileNrMaster = lastFileNrMaster + 1;
            
            if debugMode, fprintf(1,'Moving to next file (%d -> %d).\n', curFileNr, curFileNr + 1); end

end