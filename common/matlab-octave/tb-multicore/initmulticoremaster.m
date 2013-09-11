function [mcstate, settings, evalfcns, postProcStruct, parameters] = initmulticoremaster(evalfcns, fcnparams, settings)
% initialisation function for startmulticoremaster2

    % internal parameters, will be stored in a structure which contains the
    % internal state of the algorithm
    multicorerootdir = fileparts(which('startmulticoremaster2'));

    % default settings
    settingsDefault.multicoreDir        = fullfile(multicorerootdir,'Temp_Files');
    settingsDefault.nrOfEvalsAtOnce     = 1;
    settingsDefault.maxEvalTimeSingle   = 60;
    settingsDefault.initCheckPauseTime  = 0.1;
    settingsDefault.maxCheckPauseTime   = 2;
    settingsDefault.masterIsWorker      = 1;
    settingsDefault.postProcessHandle   = '';
    settingsDefault.postProcessUserData = {};
    settingsDefault.expectEmptyResults  = 0;
    settingsDefault.maxMasterEvals      = 3;
    settingsDefault.preProcessHandle    = [];
    settingsDefault.preProcessUserData  = {};
    settingsDefault.debugMode           = false;
    settingsDefault.showWarnings        = false;
    settingsDefault.nResults            = 1;
    settingsDefault.monitorFunction     = [];
    settingsDefault.monitorUserData     = [];
    settingsDefault.clearExistingFiles  = true;
    
    %%%%%%%%%%%%%%%%
    % check inputs %
    %%%%%%%%%%%%%%%%
    error(nargchk(2, 3, nargin, 'struct'))

    % check function handle cell
    if isa(evalfcns, 'function_handle') || ischar(evalfcns)
        % expand to cell array if a single function handle, or string is
        % supplied       
        evalfcns = repmat({evalfcns}, size(fcnparams));
    else
        if ~iscell(evalfcns)
            error('First input argument must be a function handle or a cell array of function handles.');
        elseif any(size(evalfcns) ~= size(fcnparams))
            error('Input cell arrays evalfcns and fcnparams must be of the same size.');
        end
    end

    % check parameter cell
    if ~iscell(fcnparams)
        error('MULTICORE:badfcnparams', ...
            'fcnparams must be a cell array each containing arguments for the evaluation functions.');
    end

    % get settings
    if ~isstruct(settings)
        error('settings must be a structure.');
    else
        % set default values where fields are missing
        fieldNames = fieldnames(settingsDefault);
        for k=1:length(fieldNames)
            if ~isfield(settings, fieldNames{k})
                settings.(fieldNames{k}) = settingsDefault.(fieldNames{k});
            end
        end
    end
    
    if settings.debugMode
        fprintf(1,'*********** Start of function %s **********\n', mfilename);
        mcstate.startTime    = mbtime;
        settings.showWarnings = true;
        mcstate.setTime      = 0;
        mcstate.removeTime   = 0;
    end
    
    % check the number of arguments desired to be returned for each function
    % in evalfcns
    if isscalar(settings.nResults)
        settings.nResults = repmat(settings.nResults, size(evalfcns));
    elseif ~samesize(evalfcns, settings.nResults)
        error('MULTICORE:nresults', ...
            'If not scalar, settings.nResults must be a matrix the same size as evalfcns');
    end
    
    % now look for inf values which indicates all results should be
    % extracted
    for ind = 1:numel(settings.nResults)
        if isinf(settings.nResults(ind))
            settings.nResults(ind) = nargout(evalfcns{ind});
            if settings.nResults(ind) == -1
                error('MULTICORE:nresults', ...
                    'Number of results for function handle %d has been set to inf, but nargout returns -1', ind);
            end
        end
    end

    % check number of evaluations at once
    mcstate.nrOfEvals = numel(fcnparams);
%     nrOfEvalsAtOnce = settings.nrOfEvalsAtOnce;
    if settings.nrOfEvalsAtOnce > mcstate.nrOfEvals
        settings.nrOfEvalsAtOnce = mcstate.nrOfEvals;
    elseif settings.nrOfEvalsAtOnce < 1
        error('Parameter nrOfEvalsAtOnce must be greater or equal one.');
    end
    settings.nrOfEvalsAtOnce = round(settings.nrOfEvalsAtOnce);


    if ~exist(settings.multicoreDir, 'dir')
        error('Slave file directory %s does not exist.', settings.multicoreDir);
    end
    
    % Determine if there is a difference in local computer time and the
    % time of the directory for determining correctly if file semaphores
    % are out of date or not
    mcstate.dirTimeDiff = localvfiledirtimediff(settings.multicoreDir, 0);

    % check maxEvalTimeSingle
    if settings.maxEvalTimeSingle < 0
        error('Parameter maxEvalTimeSingle must be greater than or equal to zero.');
    end

    % compute the maximum waiting time for a complete job
    mcstate.maxMasterWaitTime = settings.maxEvalTimeSingle * settings.nrOfEvalsAtOnce;

    % compute number of files/jobs
    mcstate.nrOfFiles = ceil(mcstate.nrOfEvals / settings.nrOfEvalsAtOnce);
    if settings.debugMode
        fprintf(1,'nrOfFiles = %d\n', mcstate.nrOfFiles);
    end

    % Initialize structure for postprocessing function
    if ~isempty(settings.postProcessHandle)
        postProcStruct.state     = 'initialization';
        postProcStruct.nrOfFiles = mcstate.nrOfFiles;
        postProcStruct.evalfcns  = evalfcns;
        postProcStruct.fcnparams = fcnparams;
        postProcStruct.userData  = settings.postProcessUserData;
        feval(settings.postProcessHandle, postProcStruct);
    else
        postProcStruct = struct;
    end

    % remove all existing temporary multicore files
    existingMulticoreFiles = [...
        findfiles(settings.multicoreDir, 'parameters_*.mat', 'nonrecursive'), ...
        findfiles(settings.multicoreDir, 'working_*.mat',    'nonrecursive'), ...
        findfiles(settings.multicoreDir, 'result_*.mat',     'nonrecursive')];
    
    if settings.clearExistingFiles
        if settings.debugMode, fprintf(1, 'Master is deleting previous mcore files\n'); end
        deletewithsemaphores2(existingMulticoreFiles);
    end
    
    try
        % try and get a really random number from online
        mcstate.ID = randi_org(99999, 1);
    catch
        % if that doesn't work, try generating one locally, using a new
        % seed based on the system clock
        fprintf(1, 'Could not access online random number generator, using rand()\n')
        mcstate.ID = round(randMat(1,99999,0,1));
    end
    
    fprintf(1, 'Master process ID is %d\n', mcstate.ID);

    % build parameter file name (including the date is important because slave
    % processes might still be working with old parameters)
    dateStr = sprintf('%04d%02d%02d%02d%02d%02d', round(clock));
    mcstate.parameterFileNameTemplate = fullfile(settings.multicoreDir, sprintf('parameters_%s_XX_mID_%d.mat', dateStr, mcstate.ID));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate parameter files %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % generate and save parameter files for all the parameter sets
    mcstate.multicoreCancelled = false;
    
    for curFileNr = mcstate.nrOfFiles:-1:1
        
%         curFileNr = mcstate.lastFileNrMaster; % for simpler copy&paste

        % Generate the parameter file name
        parameterFileName = strrep(mcstate.parameterFileNameTemplate, 'XX', sprintf('%04d', curFileNr));
        % Get the correct parameters and function handles for this set of
        % parameters from the cell arrays
        mcstate.parIndex = ((curFileNr-1)*settings.nrOfEvalsAtOnce+1) : min(curFileNr*settings.nrOfEvalsAtOnce, mcstate.nrOfEvals);
        functionHandles = evalfcns(mcstate.parIndex); %#ok
        parameters      = fcnparams(mcstate.parIndex); 
        nresultargs     = settings.nResults(mcstate.parIndex);

        if settings.debugMode, t1 = mbtime; end
        
        % Generate a file semaphore to prevent workers attempting to access
        % parameter file as we generate it
        sem = setfilesemaphore(parameterFileName, settings.showWarnings, settings.debugMode, 0, mcstate.dirTimeDiff);
        
        pause(0.01);
        
        if settings.debugMode, mcstate.setTime = mcstate.setTime + mbtime - t1; end

        % Now attempt to save the parameter file
        try
            save(parameterFileName, 'functionHandles', 'parameters', 'nresultargs'); %% file access %%
            if settings.debugMode
                fprintf(1,'Parameter file nr %d generated.\n', curFileNr);
            end
        catch
            ME = lasterror;
            if settings.showWarnings
                disp(sprintf('Warning: Unable to save parameter file No. %d at\n%s.', curFileNr, parameterFileName));
                displayerrorstruct;
            end
        end

        if settings.debugMode, t1 = mbtime; end
        
        % remove the semaphore so workers can access the file
        removefilesemaphore2(sem);
        
        if settings.debugMode, mcstate.removeTime = mcstate.removeTime + mbtime - t1; end

        if mcstate.multicoreCancelled
            return;
        end

    end
    
    mcstate.preProcState = [];
    
    % run the multicore monitor function if supplied
    if ~isempty(settings.monitorFunction)
        mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData);
    else
        mcstate.monitorstate = [];
    end

    mcstate.lastFileNrMaster = 1;         % start working down the list from top to bottom
    mcstate.lastFileNrSlave = mcstate.nrOfFiles; % check for results from bottom to top
    mcstate.parameterFileFoundTime  = NaN;
    mcstate.parameterFileRegCounter = 0;
    mcstate.nrOfFilesMaster = 0;
    mcstate.nrOfFilesSlaves = 0;
    

end