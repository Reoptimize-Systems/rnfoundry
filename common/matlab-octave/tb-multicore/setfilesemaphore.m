function semaphore = setfilesemaphore(fileList, showWarnings, debugMode, slaveID, dirtdiff)
%SETFILESEMAPHORE  Set semaphore for file access.
%   SEMAPHORE = SETFILESEMAPHORE(FILENAME) sets a semaphore to get
%   exclusive access on file FILE. The semaphore is realized by generating
%   a simple Matlab data file after checking that no other semaphores are
%   existing. The function exits if the semaphore is set. Exclusive file
%   access is of course only guaranteed if all other Matlab processes use
%   semaphores to access the same file.
%
%   The output variable SEMAPHORE is needed to correctly remove the file
%   semaphore after file access. It is an error to call function
%   SETFILESEMAPHORE without output arguments.
%
%   SEMAPHORE = SETFILESEMAPHORE(FILELIST) sets semaphores for all files
%   given in cell array FILELIST. Note that function SETFILESEMAPHORE waits
%   for exclusive file access on ALL files in the list before exiting.
%
%		Note: A semaphore older than 20 seconds is considered as invalid and
%		will immediately be deleted.
%
%		Example:
%		sem = setfilesemaphore('test.mat');
%		% access file test.mat here
%		dir test.mat.semaphore.*
%		removefilesemaphore(sem);
%
%		Markus Buehren
%		Last modified 05.04.2009
%
%   See also REMOVEFILESEMAPHORE.

% Todo: What about system time differences between the local machine and
% the machine where the temporary directory lies?

  
    if nargin < 2
        showWarnings = 0;
    end
    
    if nargin < 3
        debugMode = false;
    end
    
    if nargin < 4
        slaveID = 0;
    end
    
    if nargin < 5
        dirtdiff = 0;
    end

    persistent filesToIgnore

    if nargout ~= 1
        error('Function %s must be called with one output argument!', mfilename);
    end

    if ischar(fileList)
        % single file given
        fileList = {fileList};
        if debugMode, disp('Single file name received in setfilesemaphore'); end
    end
    
    if debugMode, disp(['Local to remote time diff is: ', num2str(dirtdiff)]); end
    
     % set times (all in seconds)
    semaphoreOldTime = 60 + dirtdiff;
    fixedWaitTime    = 0.05; % wait fixedWaitTime after generating semaphore
    checkWaitTime    = 0.02;
    waitInfoPeriod   = 5;
    maxRandomTime    = 0.3;
    
    nOfFiles = length(fileList);
    semaphore = cell(nOfFiles, 1);
    hostname = gethostname;
    for fileNr = 1:nOfFiles
        fileName = fileList{fileNr};

        % check if given file is itself a semaphore file
        if ~isempty(regexp(fileName, '\.semaphore\.\w+\.\d+\.mat$', 'once'))
            %disp('Warning: Trying to generate a semaphore file for a semaphore file!! Will be ignored.');
            semaphore{fileNr, 1} = '';
            continue
        end

        % generate semaphore file pattern of current file
        semaphorePattern     = [fileName, '.semaphore.*.mat'];
        semaphorePatternPath = fileparts(semaphorePattern);

        startWaitTime   = now;
        displayWaitInfo = true;
        while 1
            % get list of existing semaphores
            dirStruct = dir(semaphorePattern);
            
            if debugMode, disp('dir returns:'); dir(semaphorePattern); end

            semaphoreExisting = false;

            if ~isempty(dirStruct)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % other semaphore file existing %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % check if any existing semaphore file has to be respected
                for k=1:length(dirStruct)

                    % get file date
                    if isempty(dirStruct(k).date)
                        % it happens that the info in dirStruct is damaged, file is
                        % checked again later in that case
                        semaphoreExisting = true;
                        if debugMode, disp('file dirStruct Info Damaged, empty date'); end
                        continue
                    elseif isfield(dirStruct, 'datenum')
                        fileDatenum = dirStruct(k).datenum;
                    else
                        % in older Matlab versions, the field datenum seems not to exist
                        fileDatenum = datenum2(dirStruct(k).date);
                        if debugMode, disp(dirStruct(k).date); end
                    end

                    % check if current semaphore is very old and can be ignored
                    if  (now - fileDatenum) * 86400 <= semaphoreOldTime && ...
                            (now - startWaitTime) * 86400 <= semaphoreOldTime
                        % semaphore file is not old and must be respected
                        if debugMode, disp('Existing semaphore NOT old'); end
                        semaphoreExisting = true;
                    else
                        if debugMode, disp('Existing semaphore IS old'); end
                        
                        % Get the old semaphore file name
                        oldSemaphoreFileName = concatpath(semaphorePatternPath, dirStruct(k).name);

                        % avoid to issue more than one warning for each file
                        if ~isempty(filesToIgnore) && ismember(oldSemaphoreFileName, filesToIgnore)
                            % ignore file
                            continue
                        end

                        % add file to ignore list
                        filesToIgnore{end+1} = oldSemaphoreFileName; %#ok
                        disp(textwrap2(sprintf('Ignoring old semaphore of file %s.', fileName)));
                        % limit the number of saved files
                        if length(filesToIgnore) > 200
                            filesToIgnore = filesToIgnore(end-100:end);
                        end

                        % try to remove old semaphore file
                        deleteSuccessful = mbdelete2(oldSemaphoreFileName, showWarnings, 0); %% file access %%
                        
                        if debugMode
                            if deleteSuccessful
                               disp(textwrap2(sprintf('Successfully deleted old semaphore file: %s.', oldSemaphoreFileName)));
                            else
                               disp(textwrap2(sprintf('Failed to delete old semaphore file: %s.', oldSemaphoreFileName)));
                            end
                        end

                    end
                end % k=1:length(dirStruct)

                % display info
                if semaphoreExisting && displayWaitInfo && (now - startWaitTime) * 86400 >= waitInfoPeriod
                    disp(textwrap2(sprintf('Waiting for semaphore of file %s to disappear.', fileName)));
                    if debugMode
                        disp(textwrap2(sprintf('now - startWaitTime: %f', (now - startWaitTime) * 86400)));
                    end
                    displayWaitInfo = false;
                end

                % wait before checking again
                pause(checkWaitTime + maxRandomTime * generaterandomnumber);

            end % if ~isempty(dirStruct)

            if ~semaphoreExisting
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % set own semaphore file %
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                for attemptNr = 1:10
                    % build semaphore file name
                    [randomNr, randomStr] = generaterandomnumber; 
                    semaphoreFileName = [fileName, '.semaphore.', hostname, '.', randomStr, '.mat'];
                    
                    if ~isoctave
                        lasterror('reset'); % WHY?
                    end
                    
                    try
                        generateemptyfile(semaphoreFileName); %% file access %%
                        break;
                    catch
                        % in very very very unlikely cases two processes might have
                        % generated the same semaphore file name
                        if showWarnings
                            disp(textwrap2(sprintf('Warning: An error occured while generating semaphore file %s:', semaphoreFileName)));
                            displayerrorstruct;
                        end

                        % wait random time and try again
                        pause(checkWaitTime + maxRandomTime * randomNr);
                    end
                end

                % Two semaphore files might have been created at the same time -->
                % wait fixed time, then check if any other semaphore file is
                % existing. This pause time unfortunately causes overhead but seems
                % to be necessary.
                pause(fixedWaitTime + fixedWaitTime * generaterandomnumber);

                removeOwnSemaphore = false;
                dirStruct = dir(semaphorePattern);
                if length(dirStruct) > 1
                    for k=1:length(dirStruct)
                        currFileName = dirStruct(k).name;
                        if ~isempty(currFileName) && ...
                                ~strcmp(currFileName, semaphoreFileName) && ...
                                (isempty(filesToIgnore) || ~ismember(currFileName, filesToIgnore))

                            % at least one of the semaphores found may not be ignored
                            removeOwnSemaphore = true;
                            break
                        end
                    end
                end

                if ~removeOwnSemaphore
                    % exclusive file access is guaranteed
                    % save semaphore file name in output cell and leave while loop
                    semaphore{fileNr, 1} = semaphoreFileName;
                    break
                else
                    % remove own semaphore file
                    mbdelete2(semaphoreFileName, showWarnings, 0); %% file access %%

                    % wait random time before checking everything again in while-loop
                    pause(maxRandomTime * generaterandomnumber);
                end
            end % if ~semaphoreExisting
        end % while 1
    end % for fileNr = 1:nOfFiles
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [randomNr, randomStr] = generaterandomnumber
% %GENERATERANDOMNUMBER
% %   in very unlikely cases, it might happen that the random states of rand
% %   and randn are equal in two Matlab processes calling function
% %   SETFILESEMAPHORE. For this reason, the system and cpu time are used to
% %   create different random numbers even in this unlikely case.
% %
% %		This all were not necessary if it were possible to get some sort of a
% %		Matlab process ID.
% 
%     nOfDigits = 8; % length of random string will be 4*nOfDigits
% 
%     randNr    = rand;
%     randnNr   = mod(randn+0.5, 1);
%     cputimeNr = mod(cputime, 100)/100;
%     nowNr     = mod(rem(now,1)*3600*24, 100)/100;
% 
%     % random number is used for random pause after conflict
%     randomNr = 0.25 * (randNr + randnNr + cputimeNr + nowNr);
% 
%     % random string is used for the semaphore file name
%     if nargout > 1
%         ee = 10^nOfDigits;
%         randomStr = sprintf('%.0f%.0f%.0f%.0f', ...
%             ee * randNr,    ...
%             ee * randnNr,   ...
%             ee * cputimeNr, ...
%             ee * nowNr      ...
%             );
%     end
% 
% end
