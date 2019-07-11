function spawnstate = mcoreslavespawn (spawnopts, spawnstate)
% launches matlab or octave multicore slave processes and monitors them
% them
%
% Syntax
% 
% spawnstate = mcoreslavespawn (spawnopts, spawnstate)
%
% Description
%
% mcoreslavespawn launches a predetermined number of slaves and keeps them
% active, starting new processes to get the desired number. 
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
        nparamfiles = countparameterfiles (spawnopts.sharedir);

        nactiveslaves = countactiveslaves (spawnopts.sharedir);

        if (nparamfiles > 0) && (nactiveslaves < spawnopts.maxslaves) 
            % start some slaves
            
            % choose the number of slaves to launch, but no more
            % than a hard limit specified in spawnopts.maxslaves
            nslaves = spawnopts.maxslaves - nactiveslaves;

            % start the slaves
            logfcn (1, sprintf ('Attempting to launch %d slaves.', nslaves));

            startslaves ( nslaves, ...
                          'MulticoreSharedDir', spawnopts.sharedir, ...
                          'SlaveType', spawnopts.matoroct, ...
                          'PauseTime', spawnopts.pausetime, ...
                          'StartDir', spawnopts.slavestartdir );

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

function nslavefiles = countactiveslaves (dirpath)

    nslavefiles = countmatfilesstartingwith ('slaveID_', dirpath);

end

function nparamfiles = countparameterfiles (dirpath)

    nparamfiles = countmatfilesstartingwith ('parameters_', dirpath);
    
end

function nfiles = countmatfilesstartingwith (startstr, dirpath)

    nfiles = numel (dir (fullfile (dirpath, [startstr, '*.mat'])));

end

function startslaves (nslaves, varargin)
% launches slave processes for multicore, on same machine as master process

    Inputs.MulticoreSharedDir = fullfile (tempdir (), 'multicore_slaves');
    Inputs.SlaveType = 'm';
    Inputs.PauseTime = 5;
    Inputs.StartDir = '~/Documents/MATLAB/';
    % default end date is a very long time in the future (at the time of
    % writing) slaves will quit after this time
    Inputs.EndDate = [2101,2,3,4,5,6]; 
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if ~exist (Inputs.MulticoreSharedDir, 'dir')
        mkdir (Inputs.MulticoreSharedDir);
    end
    
    switch Inputs.SlaveType
        case 'm'
            launchscript = 'matlabmulticoreslaves';
        case 'o'
            launchscript = 'octavemulticoreslaves';
        otherwise
            error ('slavetype should be ''m'' for matlab of ''o'' for octave')
    end
    
    sleeptime = 0;
    for n = 1:nslaves
        
        if isunix ()
            launchargs = sprintf ( ' ''%s'' [%d,%d,%d,%d,%d,%d] ''%s'' %d > /dev/null &', ...
                Inputs.MulticoreSharedDir, ...
                Inputs.EndDate(1), Inputs.EndDate(2), Inputs.EndDate(3), Inputs.EndDate(4), Inputs.EndDate(5), Inputs.EndDate(6), ...
                Inputs.StartDir, ...
                sleeptime );
            
            fulllaunchscript = fullfile (getmfilepath (mfilename ()), [launchscript, '.sh']);
        else
            error ('Windows slave launch not yet implemented');
            fulllaunchscript = fullfile (getmfilepath (mfilename ()), [launchscript, '.cmd']);
        end
    
        system (['"', fulllaunchscript , '"' launchargs]);
        
        sleeptime = sleeptime + Inputs.PauseTime;

    end
    
end
