function spawnstate = slavespawn (spawnopts, spawnstate)
% launches matlab or octave multicore slave processes and monitors them
% them
%
% Syntax
% 
% spawnstate = mcore.slavespawn (spawnopts, spawnstate)
%
% Description
%
% slavespawn launches a predetermined number of slaves and keeps them
% active, starting new processes to get the desired number. 
% 
% Inputs
% 
%  spawnopts - structure containing options determining the spwning
%   process. The following fields must be present in the structure:
%
%   pausetime : time between launching slaves
%
%   initialwait : 
%
%   starttime : 
%
%   sharedir : 
%
%   maxslaves : 
%
%   matoroct : character vector containing 'm' or 'o' indicating whether
%     Matlab of Octave workers are to be created.
%
%   slavestartdir : 
%
%   ignoreparamfiles : 
% 
%  spawnstate - 
%
% Output
%
%  spawnstate - 
%

    logfcn = @(fid, msg) fprintf ( fid, [ 'MCORE SLAVE SPAWN: [', datestr(now()), '] ', msg, sprintf('\n')] );

    if nargin < 2 || isempty (spawnstate)
        
        spawnstate = struct ('SlavesSubmitted', 0, ...
                             'SlavesSubmissionTime', now (), ...
                             'CurrentSlaves', 0, ...
                             'SlavesDeletionTime', now (), ...
                             'DeletedSlaves', 0, ...
                             'SpawnCount', 0 );
                        
        return;
         
    end

    % get a date vector representing the current time
    currenttime = datevec (now ());
    currenttime = currenttime(4:end);

    submit_time = dateadd ( spawnstate.SlavesSubmissionTime, ...
                           [0,0,0,0,0,spawnstate.SlavesSubmitted*4*spawnopts.pausetime] );

    if spawnstate.SpawnCount == 1
        % after the first spawning event wait for initialwait
        % before spawning again
        submit_time = dateadd ( submit_time, [0,0,0,0,0,spawnopts.initialwait] );
    end

    if (isempty (spawnopts.starttime) || time2seconds (currenttime) > time2seconds (spawnopts.starttime)) ...
            && (spawnstate.SlavesSubmitted == 0 || submit_time < now)

        spawnstate.SlavesSubmitted = 0;

        % count the number of parameter files in the multicore
        % directory
        nparamfiles = mcore.countparameterfiles (spawnopts.sharedir);

        nactiveslaves = mcore.countslaveIDfiles (spawnopts.sharedir);

        if (nparamfiles > 0 || spawnopts.ignoreparamfiles) && (nactiveslaves < spawnopts.maxslaves) 
            % start some slaves
            
            % choose the number of slaves to launch, but no more
            % than a hard limit specified in spawnopts.maxslaves
            nslaves = spawnopts.maxslaves - nactiveslaves;

            % start the slaves
            logfcn (1, sprintf ('Attempting to launch %d slaves.', nslaves));

            mcore.startnslaves ( nslaves, ...
                          'MulticoreSharedDir', spawnopts.sharedir, ...
                          'SlaveType', spawnopts.matoroct, ...
                          'PauseTime', spawnopts.pausetime, ...
                          'StartDir', spawnopts.slavestartdir, ...
                          'CountExisting', false);

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
