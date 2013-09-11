function [loadSuccessful, sem, functionHandles, parameters, nresultargs] = ...
            loadslaveparams(parameterFileName, debugMode, showWarnings, slaveID)
% loadslaveparams: Loads the function handles and input variables created
% by startmulticoremaster for use by startmulticoreslave2
%
% Input:
%
%   parameterFileName
%
%   debugMode
%
%   showWarnings
%
%
% Output:
%
%   loadSuccessful - boolean scalar reporting whether the variables were
%                    successfully loaded
%
%   sem - semaphore file
%
%   functionHandles - 
%
%   parameters - 
%
%

    functionHandles = [];
    parameters = [];
    nresultargs = [];

    if debugMode
        % get parameter file number for debug messages
        fileNr = str2double(regexptokens(parameterFileName,'parameters_\d+_(\d+)\.mat'));
        disp(sprintf('****** Slave %d is checking file nr %d *******', slaveID, fileNr));
    end

    % Set a semaphore to get exclusive access on parameter file
    sem = setfilesemaphore(parameterFileName, showWarnings, debugMode, slaveID);

    % load and delete last parameter file
    loadSuccessful = true;

    if existfile(parameterFileName)

        if ~isoctave
            lastwarn(''); % Set the last matlab warning to an empty string, WHY??
            lasterror('reset'); %#ok<LERR> % reset the last matlab error message, WHY??
        end
        
        % Now attempt to load the parameter file, displaying some
        % information if this cannot be achieved
        try
            load(parameterFileName, 'functionHandles', 'parameters', 'nresultargs'); %% file access %%
        catch
            loadSuccessful = false;
            if showWarnings
                disp(sprintf('Warning: Unable to load parameter file %s.', parameterFileName));
                lastMsg = lastwarn;
                if ~isempty(lastMsg)
                    disp(sprintf('Warning message issued when trying to load:\n%s', lastMsg));
                end
                displayerrorstruct;
            end
        end

        % Check that the variables and function handles have in fact
        % been loaded correctly from the parameter file
        if loadSuccessful && (~exist('functionHandles', 'var') || ~exist('parameters', 'var') || ~exist('nresultargs', 'var'))
            % if they haven't set loadSuccessful to false and display
            % some info if showWarnings is set to true
            loadSuccessful = false;
            if showWarnings
                disp(textwrap2(sprintf(['Warning: Either variable ''%s'', ''%s''', ...
                    'or ''%s'' not existing after loading file %s.'], ...
                    'functionHandles', 'parameters', 'nresultargs', parameterFileName)));
            end
        end

        % If in debugMode, display some further information about what
        % has occured up to this point
        if debugMode
            if loadSuccessful
                disp(sprintf('Successfully loaded parameter file nr %d.', fileNr));
            else
                disp(sprintf('Problems loading parameter file nr %d.', fileNr));
            end
        end

        % Remove the parameter file from the multicore directory
        deleteSuccessful = mbdelete2(parameterFileName, showWarnings); %% file access %%
        if ~deleteSuccessful
            % If deletion is not successful it can happen that other slaves or
            % the master also use these parameters. To avoid this, ignore the
            % loaded parameters
            loadSuccessful = false;
            if debugMode
                disp(sprintf('Problems deleting parameter file nr %d. It will be ignored', fileNr));
            end
        end

    else % existfile(parameterFileName)

        % the parameter file could not be found and the variables cannot be
        % loaded
        loadSuccessful = false;
        if debugMode
            disp('No parameter files found.');
        end

    end
    
    if debugMode
        disp('Finished in loadslaveparams')
    end
        

end