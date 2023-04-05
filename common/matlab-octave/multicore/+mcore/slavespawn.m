function spawnstate = slavespawn (spawnopts, spawnstate)
% launches matlab or octave multicore slave processes and monitors them
% them
%
% Syntax
% 
% spawnstate = mcore.slavespawn ()
% spawnstate = mcore.slavespawn (spawnopts)
% spawnstate = mcore.slavespawn (spawnopts, spawnstate)
%
% Description
%
% slavespawn launches a predetermined number of slaves and keeps them
% active, starting new processes to get the desired number. 
% 
% 
%
% Inputs
% 
%  spawnopts - structure containing options determining the spwning
%   process. The following fields must be present in the structure:
%
%   PauseTime : A pause is included between launching slaves, as lauhching
%     them all simultaneously occasionally causes problems with file system
%     access as they all almost simultaneously attempt to access the same
%     files when parsing the Matlab path.
%
%   InitialWait : if SpawnCount is one, an initial delay in seconds to wait
%     before checking if additional slaves are needed and launching them. 
%
%   StartTime : 
%
%   MulticoreSharedDir : shared directory to use for the multicore
%     communication.
%
%   MaxSlaves : target number of slaves to reach
%
%   MatOrOct : character vector containing 'm' or 'o' indicating whether
%     Matlab of Octave workers are to be created.
%
%   SlaveStartDir : the startup directory 
%
%   IgnoreParamFiles : true/false flag indicating whether to ignore the
%     precence of parameter files when determining whether to launch new
%     slave proceses. If this is true, slaves will be launched regardless
%     of whether there are any files to processed waiting in the shared
%     directory. If false, slaves will only be launched if there are files
%     waiting to be processed.
%
%   OutputFilePrefix : 
%
%  spawnstate - structure containing information about the spawning process
%   and state. The following fields must be present in the structure:
%
%   SlavesSubmitted : The number of slaves submitted in the previous call
%     to mcore.slavespawn. Default is zero.
%
%   SlavesSubmissionTime : Date number with the time of the previous slave
%     submissions being started.
%
%   SpawnCount : Count of the number of times the spawning process has been
%     attempted.
%
%
% Output
%
%  spawnstate - same as input, but updated with latest information
%

    logfcn = @(fid, msg) fprintf ( fid, [ 'MCORE SLAVE SPAWN: [', datestr(now()), '] ', msg, sprintf('\n')] );

    if nargin < 2
        spawnstate = [];
    end

    if isempty (spawnstate)
        
        spawnstate = struct ('SlavesSubmitted', 0, ...
                             'SlavesSubmissionTime', now (), ...
                             ... 'CurrentSlaves', 0, ...
                             ... 'SlavesDeletionTime', now (), ...
                             ... 'DeletedSlaves', 0, ...
                             'SpawnCount', 0 );
                        
        return;
         
    end

    % get a date vector representing the current time
    currenttime = datevec (now ());
    currenttime = currenttime(4:end);

    submit_time = dateadd ( spawnstate.SlavesSubmissionTime, ...
                           [0,0,0,0,0,spawnstate.SlavesSubmitted*4*spawnopts.PauseTime] );

    if spawnstate.SpawnCount == 1
        % after the first spawning event wait for InitialWait
        % before spawning again
        submit_time = dateadd ( submit_time, [0,0,0,0,0,spawnopts.InitialWait] );
    end

    if (isempty (spawnopts.StartTime) || time2seconds (currenttime) > time2seconds (spawnopts.StartTime)) ...
            && (spawnstate.SlavesSubmitted == 0 || submit_time < now)

        spawnstate.SlavesSubmitted = 0;

        % count the number of parameter files in the multicore
        % directory
        nparamfiles = mcore.countparameterfiles (spawnopts.MulticoreSharedDir);

        nactiveslaves = mcore.countslaveIDfiles (spawnopts.MulticoreSharedDir);

        if (nparamfiles > 0 || spawnopts.IgnoreParamFiles) && (nactiveslaves < spawnopts.MaxSlaves) 
            % start some slaves
            
            % choose the number of slaves to launch, but no more
            % than a hard limit specified in spawnopts.MaxSlaves
            nslaves = spawnopts.MaxSlaves - nactiveslaves;

            % start the slaves
            logfcn (1, sprintf ('Attempting to launch %d slaves.', nslaves));

            mcore.startnslaves ( nslaves, ...
                          'MulticoreSharedDir', spawnopts.MulticoreSharedDir, ...
                          'SlaveType', spawnopts.MatOrOct, ...
                          'PauseTime', spawnopts.PauseTime, ...
                          'StartDir', spawnopts.SlaveStartDir, ...
                          'CountExisting', false, ...
                          'OutputFilePrefix', spawnopts.OutputFilePrefix);

            % increment the count of spawning events (mainly so we can
            % optionally wait for a while after the first spawning event
            % before doing another)
            spawnstate.SpawnCount = spawnstate.SpawnCount + 1;

            spawnstate.SlavesSubmitted = nslaves;
            spawnstate.SlavesSubmissionTime = now;

            logfcn (1, sprintf ('Requested %d new matlab slaves', nslaves));

        end

    end

end

function s = time2seconds (time)
% timetoseconds: takes a 3 elemnt vector representing a time of day and
% returns the number second since midnight this represents
    
    s = (time(1) * 3600) + (time(2) * 60) + time(3);

end
