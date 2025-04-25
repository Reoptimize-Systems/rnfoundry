function startslave(varargin)
% mcore.startslave: Start multi-core processing slave process.
%
% Syntax
%
% mcore.startslave()
% mcore.startslave('Param', value)
%
% Description
%
% mcore.startslave starts a slave process which looks in a given directory
% for files created by a companion function, mcore.startmaster, running in
% another matlab client. mcore.startslave runs the jobs specified in the
% communication files created by more.startmaster and saves the output
% to results files in the same directoyr for collection by
% mcore.startmaster. If called with no arguments mcore.startslave
% uses the default directory 'Temp_Files' expected to be in the same
% directory as the mcore.startslave function file. The slave also
% produces an individual file which it moitors between processing of job
% files. This file is given the default name 'slaveID_XXX' where XXX is any
% number of random digits, specific to the slave. This number is also
% appended to some other files produced by the slave to identify their
% origin. If this slaveID file is not found in the multicore directory the
% slave quits.
%
% mcore.startslave(MulticoreSharedDir) starts the slave in the directory
% specified by MulticoreSharedDir.
%
% If supplied, quitmode is a scalar value determining how the slave behaves
% on exit. If zero, the function simply exits as normal. If 1, the matlab
% process is also ended by calling the 'quit' function. If quitmode is
% empty a default value of zero is used.
%
% quitdatenum is an optional matlab date vector specifying a date and time
% when the slave should quit, regardless if further jobs are available.
% This date is only checked between each job evalation, so long jobs may
% cause the slave to overrun this date. If an empty matrix is supplied, no
% end date is assumed.
%
% IDfilestart is a string containing a replacement for the starting part of
% the slave ID file name. i.e. 'slaveID_' in the file name is replaced with
% the contents of this string. The random number identifying the slave will
% still be appended to this string. 
%
%
% See also: mcore.startmaster, mcore.startnslaves, mcore.quitallslaves
%

% Orignial Author: Markus Buehren
% Last modified: 10.04.2009
%
% Subsequently:
%
% Modified by         Date         Description
%
% Richard Crozier     2011-12-05   Added ability to specify a time when
%                                  slave should quit, ability to specify
%                                  whether matlab is quitted on slave
%                                  quitting. Added ability to quit client
%                                  by deleting a monitored file and the
%                                  ability to specify the starting name of
%                                  the file.
%
% Richard Crozier     2011-12-07   Added more detailed error checking for
%                                  multicore directoy 
%
% Further changes recorded in hg and git logs etc


% some day we will get the slave to clean up after itself when interrupted
%     cleanupinfo = struct('SlaveIDFile', [], ...
%                          'WorkingFile', [], ...
%                          'WorkingIDFile', [], ...
%                          'Sem', []);
    
    options.DebugMode = false;
    options.ShowWarnings = false;
    options.MulticoreSharedDir = mcore.defaultmulticoredir ();
    options.QuitMatlabOnExit = true;
    options.QuitDatenum = [];
    options.IDfilestart = 'slaveID_';
    options.ThrowErrors = false;
    options.SlaveID = [];
    options.LoadIntoBaseWorkspace = '';
    
    options = parse_pv_pairs (options, varargin);

    debugMode = options.DebugMode;
    showWarnings = options.ShowWarnings;
    throwerrors = options.ThrowErrors;
    IDfilestart = options.IDfilestart;
    MulticoreSharedDir = options.MulticoreSharedDir;
    quitdatenum = options.QuitDatenum;
    quitmode = options.QuitMatlabOnExit;

    maxID = 99999;

    if options.DebugMode
        fprintf(1, 'Slave options are: \n'); 
        dispstruct (options);
    end

    if isempty (options.SlaveID)
        try
            % try and get a really random number from online
            options.SlaveID = randi_org(maxID, 1);
        catch
            % if that doesn't work, try generating one locally, using a new
            % seed based on the system clock
            if showWarnings, fprintf(1, 'Could not access online random number generator, using randi()\n'); end
            options.SlaveID = randi(99999);
        end
    else
        check.isScalarInteger (options.SlaveID, true, 'SlaveID');
        if options.SlaveID <= 1 || options.SlaveID >= maxID
            error ('SlaveID must be <= 1 and >= %d, but supplied value was %d', maxID, options.SlaveID);
        end
        slaveID = options.SlaveID;
    end
    
    if ~ischar(MulticoreSharedDir)
        
        error('MulticoreSharedDir must be a character array specifying a directory location');
        
    elseif isempty(MulticoreSharedDir)
        
        % default directory for multicore
        MulticoreSharedDir = mcore.defaultMulticoreSharedDir ();
        disp(['No dir specified, starting in: ', MulticoreSharedDir])
        
    elseif ~exist(MulticoreSharedDir, 'dir')
        
        error('slave file directory %s does not exist.', MulticoreSharedDir);
        
    end

    if debugMode
        % activate all messages
        showWarnings = true;
    end

    % parameters
    firstWarnTime = 10;
    startWarnTime = 10*60;
    maxWarnTime   = 24*3600;
    startWaitTime = 1.5;
    maxWaitTime   = 13;

    if debugMode
        firstWarnTime = 10;
        startWarnTime = 10;
        maxWarnTime   = 60;
        maxWaitTime   = 1;
    end

    persistent lastSessionDateStr

    
    startWaitTime = startWaitTime + (slaveID/100000);
    maxWaitTime = maxWaitTime + (3*slaveID/100000);

    % initialize variables
    lastEvalEndClock = clock;
    lastWarnClock    = clock;
    firstRun         = true;
    curWarnTime      = firstWarnTime;
    curWaitTime      = startWaitTime;
    
%     C = onCleanup(@() onfinish(slaveID));
    fprintf (1, 'This is a multicore slave process with ID: %d\n', slaveID)
    
    slaveIDfile = fullfile(MulticoreSharedDir, [IDfilestart, num2str(slaveID), '.mat']);
    fprintf (1, 'Slave control file is:\n%s\nDelete this file to quit slave process.\n', slaveIDfile);
    
    % try to delete the slave ID file when we quit the slave for any reason
    % using an onCleanup object
    slaveIDcleanup = onCleanup(@() delete(slaveIDfile));
    
    currentTime = clock;
    
    % create the the slave ID file
    save(slaveIDfile, 'slaveID', 'currentTime');
    
%     cleanupinfo.SlaveIDFile = slaveIDfile;

    % Determine if there is a difference in local computer time and the
    % time of the multicore directory
    dirtdiff = mcore.localvfiledirtimediff(MulticoreSharedDir, slaveID);
    
    fprintf(1,'Slave ID number: %d\n', slaveID)

    if isempty(quitdatenum)
        willexitbydate = false;
    else
        willexitbydate = true;
    end
    
    settoquitasnodiraccess = false;

    if ~isempty(options.LoadIntoBaseWorkspace)
        evalin("base", "load('" + options.LoadIntoBaseWorkspace + ");");
    end
    
    while 1

        % Check if a job end datetime has been specified
        if willexitbydate
            if datenum(clock) > quitdatenum
                disp(['Clock states: ', datestr(clock)])
                disp(['quitdatenum is: ', datestr(quitdatenum)])
                disp('Quit time and date has passed, exiting multicoreslave')
                % Delete the slave ID file
                delete(slaveIDfile);
                break;
            end
        end
        
        % check the directory can be accessed
        if exist(MulticoreSharedDir, 'dir')
            
            % check for the existence of an active file for this slave
            if exist(slaveIDfile,'file') ~= 2
                % Exit the loop if it is not found
                disp('File not found for slave, exiting mcore.startslave')
                break;
            else
                % overwrite the file so we can look at the time stamp
                currentTime = clock;
                try
                    save(slaveIDfile, 'slaveID', 'currentTime');
                end
            end
            % if we have previously been unable to access the directory for
            % some reason the slave will have been set to quit, reset the
            % conditions to their original values
            if settoquitasnodiraccess
                settoquitasnodiraccess = false;
                if nargin < 3
                    willexitbydate = false;
                end
                curWaitTime = startWaitTime;
            end
        else
            
            willexitbydate = true;
            
            % keep trying for half a day before quitting if no quit date
            % is already specified
            if nargin < 3
                quitdatenum = datenum(clock) + 0.5;
            end
            
            settoquitasnodiraccess = true;
            
            if showWarnings
                disp(['Multicore directory not currently accessible, will keep trying until: ', datestr(quitdatenum)])
            end
            
            pause(curWaitTime);
            
            curWaitTime = min(curWaitTime + 0.5, maxWaitTime);
            
            continue;
        end
            
        parameterFileList = mcore.findfiles(MulticoreSharedDir, 'parameters_*.mat', 'nonrecursive');

        % get last file that is not a semaphore file
        parameterFileName = '';
        for fileNr = length(parameterFileList):-1:1
            % Look for the word 'semaphore in the file name'
            if isempty(strfind(parameterFileList{fileNr}, 'semaphore'))
                % if 'semaphore' is not found, the file is a paramaters file, so
                % choose this file name to store in parameterFileName
                parameterFileName = parameterFileList{fileNr};
                break % leave the for-loop
            end
        end

        % Check if a parameter file was found and stored in 'parameterFileName'
        if ~isempty(parameterFileName)

            % Attempt to load the parameters and function handles
            [loadSuccessful, sem, functionHandles, parameters, nresultargs, workingFile] = ...
                mcore.loadslaveparams(parameterFileName, debugMode, showWarnings, slaveID);

            % remove semaphore and continue if loading was not successful
            if ~loadSuccessful
                mcore.removefilesemaphore(sem);
                continue;
            else
%                 cleanupinfo.Sem = sem;
            end

            % Generate a temporary file which shows when the slave started working.
            % Using this file, the master can decide if the job timed out.
            % Still using the semaphore of the parameter file above.
            workingIDFileName = strrep(parameterFileName, 'parameters', sprintf('%d_', slaveID));

%             cleanupinfo.WorkingFile = workingFile;
            
            if debugMode
                fprintf(1,'Working file nr %d generated.\n', fileNr);
            end

            % remove semaphore file
            mcore.removefilesemaphore(sem);
            
%             cleanupinfo.Sem = [];

            % Now produce a file stating what this slave is working on
            mcore.generateemptyfile(workingIDFileName);
            
%             cleanupinfo.WorkingIDFile = workingIDFileName;
            
            % show progress info
            if firstRun
                fprintf(1,'First function evaluation (%s)\n', datestr(clock, 'mmm dd, HH:MM'));
                firstRun = false;
            elseif etime(clock, lastEvalEndClock) > 60
                fprintf(1,'First function evaluation after %s (%s)\n', ...
                    mcore.formattime(etime(clock, lastEvalEndClock)), datestr(clock, 'mmm dd, HH:MM'));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % *****************    evaluate function    *******************   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if debugMode
                fprintf(1,'Slave evaluates job nr %d.\n', fileNr);
                t0 = mcore.mbtime();
            end

            % Check if date string in parameter file name has changed. If yes, call
            % "clear functions" to ensure that the latest file versions are used,
            % and there are no older versions in Matlab's memory.
            sessionDateStr = mcore.regexptokens(parameterFileName, 'parameters_(\d+)_\d+\.mat');
            if ~strcmp(sessionDateStr, lastSessionDateStr) && ~mcore.isoctave
                clear functions

                if debugMode
                    disp('New multicore session detected, "clear functions" called.');
                end
            end
            lastSessionDateStr = sessionDateStr;

            % Now evaluate each of the functions described by the
            % functionHandles and parameters value pairs

            % first create a cell array to hold each of the return values from
            % the function evaluations
            result = cell(size(parameters)); %#ok

            % Now evaluate the functions
            for k = 1:numel(parameters)
                
                % preallocate a cell array to hold the function
                % evaluation return arguments
                thisresult = cell(1, nresultargs(k));

                if iscell(parameters{k})
                    
                    % if the parameter cell contains another cell array, we use the
                    % colon to pass all the input arguments
                    
                    if debugMode
                        disp('Starting evaluation of functions (cell array input)');
                    end
                    
                    if throwerrors
                        
                        [thisresult{1:nresultargs(k)}] = feval(functionHandles{k}, parameters{k}{:}); %#ok
                        
                    else
                        
                        try
                            [thisresult{1:nresultargs(k)}] = feval(functionHandles{k}, parameters{k}{:}); %#ok
                            % result{k} = functionHandles{k}(parameters{k}{:});
                        catch
                            
                            % get the error thrown
                            ME = lasterror;
                            
                            % display the error so we can find which slaves
                            % are finding the error
                            mcore.displayerrorstruct(ME);
                            
                            % return the error structure instead
                            thisresult = {ME}; %#ok<LERR>
                            
                        end
                        
                    end
                    
                else
                    
                    % If the parameter cell contains a matrix, pass this value
                    % to the function
                    
                    if debugMode
                        disp('Starting evaluation of functions (cell array input)');
                    end
                    
                    if throwerrors
                        [thisresult{1:nresultargs(k)}]  = feval(functionHandles{k}, parameters{k}); %#ok
                    else
                        try
                            [thisresult{1:nresultargs(k)}]  = feval(functionHandles{k}, parameters{k}); %#ok
                            % result{k} = functionHandles{k}(parameters{k});
                        catch
                            
                            % get the error thrown
                            ME = lasterror;
                            
                            % display the error so we can find which slaves
                            % are finding the error
                            mcore.displayerrorstruct(ME);
                            
                            % return the error structure instead
                            thisresult = {ME}; %#ok<LERR>
                            
                        end
                    end
                    
                end             
                
                % now put the result in the results cell array. If
                % we requested only a single argument from the
                % evaluation function, we extract the contents of
                % the cell array, otherwise, the whole cell array
                % goes in the results cell
                if numel(thisresult) > 1
                    % pass back the whole cell array containing the return
                    % arguments
                    result{k} = thisresult;
                else
                    % only one result, so just put the contents of the the
                    % first (and only) cell of the cell array of result
                    % arguments in the results cell array
                    result{k} = thisresult{1};
                end
                
            end

            if debugMode
                fprintf(1,'Slave finished job nr %d in %.2f seconds.\n', fileNr, mcore.mbtime() - t0); 
            end

            % Save result. Use file semaphore of the parameter file to reduce the
            % overhead.
            sem = mcore.setfilesemaphore(parameterFileName, showWarnings, debugMode, slaveID, dirtdiff);
            resultFileName = strrep(parameterFileName, 'parameters', 'result');
            try
                save(resultFileName, 'result'); %% file access %%
                if debugMode
                    fprintf(1,'Result file nr %d generated.\n', fileNr);
                end
            catch
                if showWarnings
                    fprintf(1,'Warning: Unable to save file %s.\n', resultFileName);
                    mcore.displayerrorstruct;
                end
            end

            % remove working file
            mcore.mbdelete2(workingFile, showWarnings); %% file access %%
            if debugMode
                fprintf(1,'Working file nr %d deleted.\n', fileNr);
            end
            
%             cleanupinfo.WorkingFile = [];
            
            % remove parameter file (might have been re-generated again by master)
            mcore.mbdelete2(parameterFileName, showWarnings); %% file access %%
            if debugMode
                fprintf(1,'Parameter file nr %d deleted.\n', fileNr);
            end

            % remove semaphore
            mcore.removefilesemaphore(sem);

            % save time
            lastEvalEndClock = clock;
            curWarnTime = startWarnTime;
            curWaitTime = startWaitTime;

            % remove info file as we're no longer working on it
            mcore.mbdelete2(workingIDFileName, showWarnings);
%             cleanupinfo.WorkingIDFile = [];
            
            % remove variables before next run
            clear result functionHandle parameters

        else %~isempty(parameterFileName)
            
            % There are no parameter files, so continue to idle, some
            % messages are displayed if idle for long time
            timeSinceLastEvaluation = etime(clock, lastEvalEndClock);
            
            if min(timeSinceLastEvaluation, etime(clock, lastWarnClock)) > curWarnTime
                
                if timeSinceLastEvaluation >= 10*60
                    % round to minutes
                    timeSinceLastEvaluation = 60 * round(timeSinceLastEvaluation / 60);
                end
                
                fprintf(1,'Warning: No slave files found during last %s (%s).\n', ...
                    mcore.formattime(timeSinceLastEvaluation), datestr(clock, 'mmm dd, HH:MM'));
                
                lastWarnClock = clock;
                
                if firstRun
                    curWarnTime = startWarnTime;
                else
                    curWarnTime = min(curWarnTime * 2, maxWarnTime);
                end
                
                curWaitTime = min(curWaitTime + 0.5, maxWaitTime);
                
            end
            % wait before next check
            pause(curWaitTime);

        end

    end

    % We have exited the loop, so either the slave id file was deleted or
    % inaccessible for too long, or the quit time was specified and has 
    % passed
    if quitmode == true
        disp('QuitMatlabOnExit is set to true, so quitting matlab client')
        quit;
    else
        disp('QuitMatlabOnExit is set to false, not quitting client')
    end
    
end

% function onfinish(cleanupinfo)
% % deletes the temporary files created by the slave on exit
%                      
%     if ~isempty(cleanupinfo.SlaveIDFile) && exist(cleanupinfo.SlaveIDFile, 'file')
%         delete(cleanupinfo.SlaveIDFile);
%     end
%     
%     if ~isempty(cleanupinfo.WorkingFile) && exist(cleanupinfo.WorkingFile, 'file')
%         delete(cleanupinfo.WorkingFile);
%     end
%     
%     if ~isempty(cleanupinfo.WorkingIDFile) && exist(cleanupinfo.WorkingIDFile, 'file')
%         delete(cleanupinfo.WorkingIDFile);
%     end
%     
%     if ~isempty(cleanupinfo.Sem) && exist(cleanupinfo.Sem, 'file')
%         delete(cleanupinfo.Sem);
%     end
%     
% end
