function [istate] = mcoreresubmitparams(istate, settings, functionHandleCell, parameterCell, parameterFileName)
% resubmit multicore parameter files
%

    functionHandles = functionHandleCell(istate.parIndex); %#ok
    parameters      = parameterCell(istate.parIndex);
    nresultargs     = settings.nResults(istate.parIndex);
    
    try
    
        save(parameterFileName, 'functionHandles', 'parameters', 'nresultargs'); %% file access %%

        if settings.debugMode
            
            fprintf(1,'Parameter file nr %d was generated again (%d. time).\n', ...
                curFileNr, istate.parameterFileRegCounter);
            
        end
        
    catch
        
        if settings.showWarnings
            disp(textwrap2(sprintf('Warning: Unable to save file %s.', parameterFileName)));
            displayerrorstruct;
        end
        
    end
    
    istate.parameterFileRegCounter = istate.parameterFileRegCounter + 1;
                        
end