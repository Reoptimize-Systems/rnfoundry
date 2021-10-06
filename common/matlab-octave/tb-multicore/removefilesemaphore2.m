function removefilesemaphore2(semaphore)
%REMOVEFILESEMAPHORE  Remove semaphore after file access.
%   REMOVEFILESEMAPHORE(SEMAPHORE) removes the semaphore(s) set by function
%   SETFILESEMAPHORE to allow file access for other Matlab processes.
%
%		Example:
%		sem = setfilesemaphore('test.mat');
%		% access file test.mat here
%		dir test.mat.semaphore.*
%		removefilesemaphore(sem);
%
%		Markus Buehren
%		Last modified 07.04.2008
%
%   See also SETFILESEMAPHORE.

    checkWaitTime = 0.1;

    showWarnings  = 1;

    % remove semaphore files
    for fileNr = 1:length(semaphore)
        
        if existfile(semaphore{fileNr})

            % do not use function deletewithsemaphores.m here!

            % sometimes deletion permission is not given, so try several times to
            % delete the file
            if ~isoctave
                warnID = 'MATLAB:DELETE:Permission';
                warnState = warning('query', warnID);
                warning('off', warnID);
            end

            fileDeleted = false;
            for attemptNr = 1:10
                lastwarn('');
                try
                    delete(semaphore{fileNr}); %% file access %%
                    if isempty(lastwarn)
                        fileDeleted = true;
                        break;
                    elseif ~existfile(semaphore{fileNr})
                        if showWarnings
                            disp('Delete was not successful, but semaphore file is not found anyway? Marking as deleted.')
                        end
                        fileDeleted = true;
                        break;
                    end
                catch
                    % wait before checking again
                    pause(checkWaitTime);
                end
            end
            
            if ~isoctave
                warning(warnState.state, warnState.identifier);
            end

            if showWarnings && ~fileDeleted
                % try one last time with display of warning message
                delete(semaphore{fileNr});
            end
            
        end % if existfile(semaphore{fileNr})
        
    end % for fileNr = 1:length(semaphore)

end
