function resultCell = startmulticoremaster2(evalfcns, fcnparams, settings)
%STARTMULTICOREMASTER  Start multi-core processing master process.
%   RESULTCELL = STARTMULTICOREMASTER(FHANDLE, fcnparams)
%   starts a multi-core processing master process. The function specified
%   by the given function handle is evaluated with the parameters saved in
%   each cell of fcnparams. Each cell may include parameters in any
%   form or another cell array which is expanded to an argument list using
%   the {:} notation to pass multiple input arguments. The outputs of the
%   function are returned in cell array RESULTCELL of the same size as
%   fcnparams. Only the first output argument of the function is
%   returned. If you need to get multiple outputs, write a small adapter
%   that puts the outputs of your function into a single cell array.
%
%   To make use of multiple cores/machines, function STARTMULTICOREMASTER
%   saves files with the function handle and the parameters to a temporary
%   directory (default: <TEMPDIR2>/multicorefiles, where <TEMPDIR2> is the
%   directory returned by function TEMPDIR2). These files are loaded by
%   function STARTMULTICORESLAVE running in other Matlab processes which
%   have access to the temporary directory. The slave processes evaluate
%   the given function with the saved parameters and save the result in
%   another file. The results are later collected by the master process.
%
%   Note that you can make use of multiple cores on a single machine or on
%   different machines with a commonly accessible directory/network share
%   or a combination of both.
%
%   RESULTCELL = STARTMULTICOREMASTER(FHANDLE, fcnparams, SETTINGS)
%   The additional input structure SETTINGS may contain any of the
%   following fields:
%
%   settings.multicoreDir:
%     Directory for temporary files (standard directory is used if empty)
%   settings.nrOfEvalsAtOnce:
%     Number of function evaluations gathered to a single job.
%   settings.maxEvalTimeSingle:
%     Timeout for a single function evaluation. Choose this parameter
%     appropriately to get optimum performance.
%   settings.masterIsWorker:
%     If true, master process acts as worker and coordinator, if false the
%     master acts only as coordinator.
%   settings.preProcessHandle:
%     A handle to a function to be called before evaluation of the jobs,
%     default is an empty matrix meaning no prepreocessing is performed
%   settings.preProcessUserData:
%     A cell array of arguments to be passed to the function in
%     settings.preProcessHandle, default is empty cell array
%   settings.debugMode:
%     Display detailed info for debugging if true, or not if false
%   settingsDefault.showWarnings:
%     Display more info about progress, a subset of the info supplied if
%     settings.debugMode is true
%   settingsDefault.nResults:
%     Integer value determining the number of return arguments from the
%     functions in the evalfcns to be captued and returned in the
%     results cell array.
%
%   Please refer to the heavily commented demo function MULTICOREDEMO for
%   details and explanations of the settings.
%
%   RESULTCELL = STARTMULTICOREMASTER(FHANDLECELL, fcnparams, ...),
%   with a cell array FHANDLECELL including function handles, allows to
%   evaluate different functions.
%
%   Example: If you have your parameters saved in parameter cell
%   fcnparams, the for-loop
%
%   	for k=1:numel(fcnparams)
%   		RESULTCELL{k} = FHANDLE(fcnparams{k});
%   	end
%
%   which you would run in a single process can be run in parallel on
%   different cores/machines using STARTMULTICOREMASTER and
%   STARTMULTICORESLAVE. Run
%
%   	RESULTCELL = STARTMULTICOREMASTER(FHANDLE, fcnparams, DIRNAME)
%
%   in one Matlab process and
%
%   	STARTMULTICORESLAVE(DIRNAME)
%
%   in one or more other Matlab processes.
%
%		Markus Buehren
%		Last modified 19.06.2009
%
%   See also STARTMULTICORESLAVE, FUNCTION_HANDLE.
    
    if nargin < 3
        settings = struct;
    end
    
    [mcstate, settings, evalfcns, postProcStruct, parameters] = ...
        initmulticoremaster(evalfcns, fcnparams, settings);
    
    % preallocate a cell array to hold the results
    resultCell = cell(size(fcnparams));

    % Call "clear functions" to ensure that the latest file versions are used,
    % no older versions in Matlab's memory.
    if ~isoctave, clear functions; end
    
    while 1 % this while-loop will be left if all work is done
        
        % check the directory can be accessed
        if ~exist(settings.multicoreDir, 'dir')

            if settings.showWarnings
                disp('Multicore directory not currently accessible')
            end
            
            if mcstate.multicoreCancelled
                return;
            else
                pause(mcstate.curPauseTime);
                mcstate.curPauseTime = min(settings.maxCheckPauseTime, mcstate.curPauseTime + settings.initCheckPauseTime);
            end

            continue;
        end
        
        % run a supplied preprocessing function if present
        if ~isempty(settings.preProcessHandle)
            mcstate.preProcState = feval(settings.preProcessHandle, settings.preProcessUserData{:}, mcstate.preProcState);
        end
        
        % run the multicore monitor function before beginning evaluation if
        % it has been supplied
        if ~isempty(settings.monitorFunction)
            mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData, mcstate.monitorstate);
        end
        
        if settings.masterIsWorker 
            
            % if the master is also a worker, it looks for results, and
            % also does jobs which have timed out here
            [mcstate, resultCell] = masterisworkerdowork(mcstate, settings, evalfcns, fcnparams, resultCell);
            
        else % if settings.masterIsWorker
            
            % pause a little time to let the slaves catch up (this is
            % different from the main loop pause time stored in
            % curPauseTime
            pause(min(settings.maxCheckPauseTime, max(settings.initCheckPauseTime, 0.01 * settings.maxEvalTimeSingle)));
            
            % run the multicore monitor function if supplied
            if ~isempty(settings.monitorFunction)
                mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData, mcstate.monitorstate);
            end
            
        end % if settings.masterIsWorker

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if all work is done %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % (mcstate.lastFileNrMaster - 1) is the number of the file/job that
        % was last computed/loaded when working down the list from top to
        % bottom.
        %
        % (mcstate.lastFileNrSlave + 1) is the number of the file/job that
        % was last computed/loaded when checking for results from bottom to
        % top.
        if (mcstate.lastFileNrMaster - 1) + 1 == (mcstate.lastFileNrSlave + 1)
            % all results have been collected, leave big while-loop
            if settings.debugMode
                disp('********************************');
                fprintf(1,'All work is done (mcstate.lastFileNrMaster = %d, mcstate.lastFileNrSlave = %d).\n', mcstate.lastFileNrMaster, mcstate.lastFileNrSlave);
            end
            break;
        end

        mcstate.curPauseTime = settings.initCheckPauseTime;

        % check for results from the slaves
        [mcstate, resultCell] = mcorecheckforresults(mcstate, settings, evalfcns, fcnparams, resultCell);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check again if all work is done %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (see comment between the two while-loops)
        if (mcstate.lastFileNrMaster - 1) + 1 == (mcstate.lastFileNrSlave + 1)
            % all results have been collected, leave big while-loop
            if settings.debugMode
                disp('********************************');
                fprintf(1,'All work is done (mcstate.lastFileNrMaster = %d, mcstate.lastFileNrSlave = %d).\n', ...
                        mcstate.lastFileNrMaster, mcstate.lastFileNrSlave);
            end
            break
        end

%         firstRun = false;
    end % while 1

    if settings.debugMode
        fprintf(1,'\nSummary:\n--------\n');
        fprintf(1,'%2d jobs in total\n',         mcstate.nrOfFiles);
        fprintf(1,'%2d jobs done by master\n', mcstate.nrOfFilesMaster);
        fprintf(1,'%2d jobs done by slaves\n', mcstate.nrOfFilesSlaves);
        %disp('No jobs done by slaves. (Note: You need to run function startmulticoreslave.m in another Matlab session?)');

        overallTime = mbtime() - mcstate.startTime;
        fprintf(1,'Processing took %.1f seconds.\n', overallTime);
        fprintf(1,'Overhead caused by setting  semaphores: %.1f seconds (%.1f%%).\n', ...
            mcstate.setTime,    100*mcstate.setTime    / overallTime);
        fprintf(1,'Overhead caused by removing semaphores: %.1f seconds (%.1f%%).\n', ...
            mcstate.removeTime, 100*mcstate.removeTime / overallTime);
        fprintf(1,'\n*********** End of function %s **********\n', mfilename);
    end

end % function startmulticoremaster2




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancelCallback(hObject, eventdata, multicoreCancelHandle) %#ok
    %CANCELCALLBACK  Cancel multicore computation on user request.
    %   This callback function is used to avoid having to use "hObject" and
    %   "eventdata" in function multicoreCancel above.

    % remove cancel button
    delete(hObject);

    % call cancel function
    multicoreCancelHandle();

end % function






