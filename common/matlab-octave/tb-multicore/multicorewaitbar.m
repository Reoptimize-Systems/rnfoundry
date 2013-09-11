function multicorewaitbar(command, varargin)
%MULTICOREWAITBAR  Handle multicore waitbar.

    persistent useWaitbar initialized waitbarHandle waitbarMessage fractionReady
    persistent clockStart1 clockStart2 clockInit0 clockUpdate multicoreCancelHandle
    if isempty(initialized)
        initialized   = 0;
        useWaitbar    = 0;
        fractionReady = 0;
    end

    switch command
        case 'init0'
            % called at the very beginning

            % save persistent variables
            multicoreCancelHandle = varargin{1};
            useWaitbar = 1;

            % remember time
            clockInit0  = clock;
            clockUpdate = clock;

            % nothing else to do now
            return

        otherwise
            % nothing to do
    end

    if ~useWaitbar
        return
    end

    tag = 'Multicore waitbar';
    waitbarExisting = ~isempty(waitbarHandle) && ishandle(waitbarHandle);
    updateNow = etime(clock, clockUpdate) > 1.0;

    switch command
        case 'init1'
            % called before parameter file generation

            % change waitbar message
            istate.nrOfFiles = varargin{1};
            waitbarMessage = sprintf('Generating parameter files.\n0/%d done.\n\n', istate.nrOfFiles);
            fractionReady = 0;

            % remember time
            clockStart1 = clock;
            clockUpdate = clock;

        case 'update1'
            % update during parameter file generation
            if updateNow
                istate.nrOfFiles = varargin{1};
                istate.lastFileNrMaster    = varargin{2};
                nrOfFilesReady = istate.lastFileNrMaster - 1;
                fractionReady = 1 - (nrOfFilesReady / istate.nrOfFiles); % istate.lastFileNrMaster decreases from nrOfFiles to 1
                if fractionReady > 0
                    timeLeft = etime(clock, clockStart1) * (1 - fractionReady) / fractionReady;
                    waitbarMessage = sprintf('Generating parameter files.\n%d/%d done.\nest. time left: %s\n', ...
                        istate.nrOfFiles - nrOfFilesReady, istate.nrOfFiles, formatTime(round(timeLeft), 'short'));
                else
                    waitbarMessage = '';
                end
            end

        case 'init2'
            % called before start of big while-loop

            % change waitbar message
            waitbarMessage = sprintf('0.0%% done by master\n0.0%% done by slaves\n0.0%% done overall\n');
            fractionReady  = 0;

            % force update, as first function evaluation may take a long time
            updateNow = 1;

            % remember time
            clockStart2 = clock;
            clockUpdate = clock;

        case 'update2'
            % update after function evaluation
            if updateNow
                istate.nrOfFiles       = varargin{1};
                istate.nrOfFilesMaster = varargin{2};
                istate.nrOfFilesSlaves = varargin{3};
                nrOfFilesReady  = istate.nrOfFilesMaster + istate.nrOfFilesSlaves;

                if 1%nrOfFilesReady > 1 % master should have checked for results at least once

                    fractionReady = nrOfFilesReady / istate.nrOfFiles;

                    if istate.nrOfFilesSlaves > 0
                        nrOfFilesSlavesTmp = istate.nrOfFilesSlaves + 0.5; % tweak for better estimation
                    else
                        nrOfFilesSlavesTmp = istate.nrOfFilesSlaves;
                    end
                    filesPerSecond = ...
                        istate.nrOfFilesMaster    / etime(clock, clockStart2) + ...
                        nrOfFilesSlavesTmp / etime(clock, clockStart1);
                    timeLeft = (istate.nrOfFiles - nrOfFilesReady) / filesPerSecond;

                    waitbarMessage = sprintf('%.1f%% done by master\n%.1f%% done by slaves\n%.1f%% done overall\nest. time left: %s', ...
                        100 * istate.nrOfFilesMaster / istate.nrOfFiles, 100 * istate.nrOfFilesSlaves / istate.nrOfFiles, ...
                        100 * nrOfFilesReady  / istate.nrOfFiles, formatTime(round(timeLeft), 'short'));
                else
                    waitbarMessage = '';
                end
            end

        case 'init3'
            % called after user cancellation

            % change waitbar message
            waitbarMessage = sprintf('Removing parameter files.\n0.0%% done.\n\n');
            fractionReady = 0;

            % remember time
            clockStart1 = clock;
            clockUpdate = clock;

        case 'update3'
            % update during deleting parameter files (if user cancelled)
            if updateNow
                minFileNr = varargin{1};
                maxFileNr = varargin{2};
                istate.lastFileNrMaster    = varargin{3};
                istate.nrOfFiles = maxFileNr - minFileNr + 1;
                fractionReady = (istate.lastFileNrMaster - minFileNr + 1) / istate.nrOfFiles;
                timeLeft = etime(clock, clockStart1) * (1 - fractionReady) / fractionReady;
                waitbarMessage = sprintf('Removing parameter files.\n%.1f%% done.\nest. time left: %s\n', ...
                    100 * fractionReady, formatTime(round(timeLeft), 'short'));
            end

        case 'delete'
            % delete waitbar
            if waitbarExisting
                delete(waitbarHandle);
            end
            initialized = 0;
            return

        otherwise
            error('Command "%s" unknown.', command);
    end

    if ~initialized && ~isempty(waitbarMessage) && (...
            (fractionReady  > 0.05 && etime(clock, clockInit0) >  5.0) || ...
            etime(clock, clockInit0) > 10.0 )

        % remove all waitbars generated before
        delete(findobj('Tag', tag));

        % generate new waitbar
        waitbarHandle = waitbar(fractionReady, waitbarMessage, 'Name', 'Multicore progress', 'Tag', tag, ...
            'CreateCancelBtn', {@cancelCallback, multicoreCancelHandle});
        set(waitbarHandle, 'HandleVisibility', 'on', 'CloseRequestFcn', 'closereq');
        initialized = 1;
        clockUpdate = clock;
    end

    if initialized && waitbarExisting && updateNow
        % update waitbar
        waitbar(fractionReady, waitbarHandle, waitbarMessage);
        clockUpdate = clock;
    end

end % function