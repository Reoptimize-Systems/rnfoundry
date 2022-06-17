function [result, resultLoaded] = loadResultFile2(resultFileName, showWarnings)
% loadResultFile2: loads a result file produced by multicore slave process

    % reset warnings and errors
    lastwarn('');
    lasterror('reset');

    % try to load file
    try
        result = []; % (only for M-Lint)
        load(resultFileName, 'result'); %% file access %%
        resultLoaded = true;
    catch %#ok
        resultLoaded = false;
        if showWarnings
            fprintf('Warning: Unable to load file %s.\n', resultFileName);
            mcore.displayerrorstruct;
        end
    end

    % display warning (if any)
    if showWarnings
        lastMsg = lastwarn;
        if ~isempty(lastMsg)
            fprintf('Warning issued when trying to load file %s:\n%s\n', ...
                resultFileName, lastMsg);
        end
    end

    % check if variable 'result' is existing
    if resultLoaded && ~exist('result', 'var')
        if showWarnings
            fprintf('Warning: Variable ''%s'' not existing after loading file %s.\n', ...
                'result', resultFileName);
        end
        resultLoaded = false;
    end

    if resultLoaded
        % it seems that loading was successful
        % try to remove result file
        mcore.mbdelete2(resultFileName, showWarnings); %% file access %%
    end

end % function