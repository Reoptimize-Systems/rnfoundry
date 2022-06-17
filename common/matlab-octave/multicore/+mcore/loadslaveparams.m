function [loadSuccessful, sem, functionHandles, parameters, nresultargs, workingFile] = ...
            loadslaveparams (parameterFileName, debugMode, showWarnings, slaveID)
% loadslaveparams: Loads the function handles and input variables created
% by startmulticoremaster2 for use by mcore.startmulticoreslave2
%
% Syntax
%
% [loadSuccessful, sem, functionHandles, parameters, nresultargs, workingFile] = ...
%            mcore.loadslaveparams(parameterFileName, debugMode, showWarnings, slaveID)
%
% Description
%
% loadslaveparams Loads the function handles and input variables created
% by startmulticoremaster2 for use by startmulticoreslave2. It creates a 
% semaphore to prevent access by multiple slave processes and attempts to load 
% the variables from the file. If successful a working file name is created for 
% use to indicate that this slave is working on the file. If the load is 
% successful the parameter file is moved to this working file name, preserving 
% its contents. This aids investigation in the event of a failure of any kind.
%
% Input
%
%  parameterFileName - string containing the path of the parameter file to be 
%    loaded. If successfully loaded, the parameter file will be removed and a
%    a corresponding working file will be created.
%
%  debugMode - flag determining whether to operate n a debugging mode with 
%    display of warnings and detailed information about the progress of the 
%    process
%
%  showWarnings - flag determining whether to show warning messages
%
%  slaveID - ID number of the slave process calling the function
%
% Output
%
%  loadSuccessful - boolean scalar reporting whether the variables were
%    successfully loaded from the parameter file
%
%  sem - name of semaphore file created to prevent other slaves accessing the 
%    parameter file while this slave process works with it
%
%  functionHandles - function handles to be evaluated which were loaded from the 
%    parameter file (empty if this was not successful)
%
%  parameters - input arguements for the functionHandles to be evaluated (empty
%    if loading was not successful)
%
%  nresultargs - the number of output arguments to be stored from the output of 
%    the functions to be evaluated
%
%  workingFile - the name of the working file created 
%
%
% See also: mcore.startmulticoreslave2
%

    functionHandles = [];
    parameters = [];
    nresultargs = [];

    if debugMode
        % get parameter file number for debug messages
        fileNr = str2double (mcore.regexptokens (parameterFileName,'parameters_\d+_(\d+)\.mat'));
        disp (sprintf ('****** Slave %d is checking file nr %d *******', slaveID, fileNr));
    end

    % Set a semaphore to get exclusive access on parameter file
    sem = mcore.setfilesemaphore (parameterFileName, showWarnings, debugMode, slaveID);

    % load and delete last parameter file
    loadSuccessful = false;
    
    % generate the working file name, the parameter file name will be
    % renamed to this on successful load of the parmeters file
    workingFile = strrep (parameterFileName, 'parameters', 'working');

    if mcore.existfile (parameterFileName)

        lastwarn (''); % Set the last matlab warning to an empty string, WHY??
        lasterror ('reset'); %#ok<LERR> % reset the last matlab error message, WHY??
        
        % Now attempt to load the parameter file, displaying some
        % information if this cannot be achieved
        try
            load (parameterFileName, 'functionHandles', 'parameters', 'nresultargs'); %% file access %%
            loadSuccessful = true;
        catch
            loadSuccessful = false;
            if showWarnings
                disp (sprintf ('Warning: Unable to load parameter file %s.', parameterFileName));
                lastMsg = lastwarn;
                if ~isempty (lastMsg)
                    disp (sprintf ('Warning message issued when trying to load:\n%s', lastMsg));
                end
                mcore.displayerrorstruct;
            end
        end

        % Check that the variables and function handles have in fact
        % been loaded correctly from the parameter file, and if so create a 
        % working file
        if loadSuccessful 
        
            if (~exist ('functionHandles', 'var') ...
                  || ~exist ('parameters', 'var') ...
                  || ~exist ('nresultargs', 'var'))
                            
                % if they haven't, set loadSuccessful to false and display
                % some info if showWarnings is set to true
                loadSuccessful = false;
                if showWarnings
                    disp (mcore.textwrap2 (sprintf (['Warning: Either variable ''%s'', ''%s''', ...
                        'or ''%s'' not existing after loading file %s.'], ...
                        'functionHandles', 'parameters', 'nresultargs', parameterFileName)));
                end
            
            else
            
                if debugMode, disp(sprintf('Successfully loaded parameter file nr %d.', fileNr)); end
                
                % Move the parameter file to the working file location
                moveSuccessful = movefile (parameterFileName, workingFile); %% file access %%

                if ~moveSuccessful
                    % If move is not successful it can happen that other slaves or
                    % the master also use these parameters. To avoid this, ignore the
                    % loaded parameters
                    loadSuccessful = false;
                    if debugMode
                        disp(sprintf('Problems deleting parameter file nr %d. It will be ignored', fileNr));
                    end
                end
            
            end
        else % if loadSuccessful 
            if debugMode, disp (sprintf ('Problems loading parameter file nr %d.', fileNr)); end
        end

    else % mcore.existfile(parameterFileName)

        % the parameter file could not be found and the variables cannot be
        % loaded
        loadSuccessful = false;
        if debugMode
            disp ('No parameter files found.');
        end

    end
    
    if debugMode
        disp ('Finished in loadslaveparams')
    end

end