function [mcstate, resultCell] = mcoremastereval(mcstate, settings, evalfcns, fcnparams, resultCell)
% evaluates funtion in a multicore process for the master process
%
% 

    if settings.debugMode
        fprintf(1,'Master evaluates job nr %d.\n', mcstate.nrOfFilesMaster);
        t0 = mbtime;
    end

    for k = mcstate.parIndex

        if settings.debugMode
            fprintf('Master evaluates parameter index %d from parameter file %d\n', k, mcstate.nrOfFilesMaster);
        end
        
        % run the multicore monitor function before beginning evaluation if
        % it has been supplied
        if ~isempty(settings.monitorFunction)
            mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData, mcstate.monitorstate);
        end

        firstmasterfuneval = true;
        evalcount = 0;
        
        while firstmasterfuneval || ((~settings.expectEmptyResults && isempty(resultCell{k})) && evalcount < settings.maxMasterEvals)

            firstmasterfuneval = false;

            % preallocate a cell array to hold the function
            % evaluation return arguments
            thisresult = cell(1, settings.nResults(k));

            if iscell(fcnparams{k})
                [thisresult{1:settings.nResults(k)}] = feval(evalfcns{k}, fcnparams{k}{:});
            else
                [thisresult{1:settings.nResults(k)}] = feval(evalfcns{k}, fcnparams{k});
            end

            % now put the result in the results cell array. If
            % we requested only a single argument from the
            % evaluation function, we extract the contents of
            % the cell array, otherwise, the whole cell array
            % goes in the results cell
            if numel(thisresult) > 1
                resultCell{k} = thisresult;
            else
                resultCell{k} = thisresult{1};
            end

            evalcount = evalcount + 1;
        end

        if evalcount >= settings.maxMasterEvals
            % the function could not be evaluated
            % despite trying pretty hard, return empty
            % matrix in the cell so the user can deal with
            % the result
            resultCell{k} = [];
        end

        if mcstate.multicoreCancelled, return; end
    end
    
    mcstate.nrOfFilesMaster = mcstate.nrOfFilesMaster + 1;

    % run the multicore monitor function before beginning evaluation of the
    % post-processing function, if it has been supplied
    if ~isempty(settings.monitorFunction)
        mcstate.monitorstate = feval(settings.monitorFunction, settings.monitorUserData, mcstate.monitorstate);
    end
        
    % Run postprocessing function
    if ~isempty(settings.postProcessHandle)
        postProcStruct.state               = 'after master evaluation'; % no copy & paste here!!
        postProcStruct.lastFileNrReady     = mcstate.lastFileNrMaster;  % no copy & paste here!!
        postProcStruct.lastFileNrMaster    = mcstate.lastFileNrMaster;
        postProcStruct.lastFileNrSlave     = mcstate.lastFileNrSlave;
        postProcStruct.nrOfFilesMaster     = mcstate.nrOfFilesMaster;
        postProcStruct.nrOfFilesSlaves     = mcstate.nrOfFilesSlaves;
        postProcStruct.resultCell          = resultCell;
        postProcStruct.parIndex            = mcstate.parIndex;
        feval(settings.postProcessHandle, postProcStruct);
    end

    if settings.debugMode
        fprintf(1,'Master finished job nr %d in %.2f seconds.\n', mcstate.nrOfFilesMaster - 1, mbtime - t0);
    end

end