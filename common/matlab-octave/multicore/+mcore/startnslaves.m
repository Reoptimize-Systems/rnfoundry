function spawnstate = startnslaves (nslaves, varargin)
% launches slave processes for multicore, on same machine as master process
%
% Syntax
%
% startnslaves (nslaves)
% startnslaves (nslaves, 'Param', value)
%
% Description
%
% Starts one or more slave processes on the local machine. Slaves can be
% stopped by deleting individual slave ID file (files created by slaves
% named slave_ID*.mat in the shared multicore directory), or by calling
% mcore.quitallslaves.
%
% Input
%
%  nslaves - The number of slaves to launch on the local machine.
%
% Additional optional parameters can be supplied as Parameter-value pairs
%
%  'MulticoreSharedDir' - shared directory to use for the multicore
%    communication. By default the value returned by the function
%
%    mcore.defaultmulticoredir ()
%
%    is used.
%
%  'CountExisting' - true/false flag indicating whether existing slaves
%    should be counted in the total number to be launched. Default is true,
%    so only the number of slaves will be launched which are necessary to
%    take the total number of active slaves to the value of nslaves.
%
%  'SlaveType' - either 'm' or 'o' to indicated if Matlab (m) or Octave (o)
%    slaves should be launched. Default is 'm'.
%
%  'PauseTime' - A pause is included between launching slaves, as lauhching
%    them all simultaneously occasionally causes problems with file system
%    access as they all almost simultaneously attempt to access the same
%    files when parsing the Matlab path. Default is 5s between slave
%    launches.
%
%  'StartDir' - The startup directory of the slaves when starting.
%
%  'EndDate' - datetime object or date vector (e.g. [2101,2,3,4,5,6])
%    indicating when the slaves should quit. Allows the user to set an
%    expiry date for slave so it will self-terminate after a given time. By
%    default this is set to the date [2101,2,3,4,5,6], far in the future so
%    there is effectively no end date.
%
%  'OutputFilePrefix' - character vector containing file prefix to use for 
%    piping the output of the script which launches the slave process.
%    Default is /dev/null so no output is saved.
%
%
% See Also: mcore.startslave, mcore.quitallslaves
%

    Inputs.MulticoreSharedDir = mcore.defaultmulticoredir ();
    if mcore.isoctave ()
        Inputs.SlaveType = 'o';
    else
        Inputs.SlaveType = 'm';
    end
    Inputs.PauseTime = 5;
    Inputs.StartDir = userpath ();
    % default end date is a very long time in the future (at the time of
    % writing) slaves will quit after this time
    Inputs.EndDate = [2101,2,3,4,5,6];
%     Inputs.OutputFilePrefix = '/dev/null';
    Inputs.OutputFilePrefix = 'slave';
    Inputs.CountExisting = true;
    Inputs.SlaveNiceProcessPriority = 5;
    Inputs.CopyWorkspace = false;

    Inputs = parse_pv_pairs (Inputs, varargin);

    check.isLogicalScalar (Inputs.CountExisting, true, 'CountExisting');
    check.isScalarInteger (Inputs.SlaveNiceProcessPriority, true, 'SlaveNiceProcessPriority');
    
    if (Inputs.SlaveNiceProcessPriority < -20) || (Inputs.SlaveNiceProcessPriority > 19)
        error ('SlaveNiceProcessPriority must be between -20 and 19');
    end

    if isdatetime (Inputs.EndDate)

        Inputs.EndDate = datevec (Inputs.EndDate);

    elseif ~(isreal (Inputs.EndDate) && numel(Inputs.EndDate) == 6)

        error ('EndDate should be a datetime object or a 6 element date vector');

    end

    if ispc
        error ('Starting multiple slaves is currently only supported on linux');
    end

    if ~exist (Inputs.MulticoreSharedDir, 'dir')
        mkdir (Inputs.MulticoreSharedDir);
    end

    if Inputs.CopyWorkspace
        workspace_file = [ tempname(Inputs.MulticoreSharedDir), '.mat' ];
        evalin('base', sprintf('save(''%s'');', workspace_file));
    else
        workspace_file = '';
    end

    switch Inputs.SlaveType
        case 'm'
            launchscript = 'matlabmulticoreslaves';
        case 'o'
            launchscript = 'octavemulticoreslaves';
        otherwise
            error ('slavetype should be ''m'' for matlab of ''o'' for octave')
    end

    if Inputs.CountExisting

        spawnopts.MaxSlaves = nslaves;
        spawnopts.MulticoreSharedDir = Inputs.MulticoreSharedDir;
        spawnopts.MatOrOct = Inputs.SlaveType;
        spawnopts.PauseTime = Inputs.PauseTime;
        spawnopts.SlaveStartDir = Inputs.StartDir;
        spawnopts.InitialWait = 0;
        spawnopts.PauseTime = Inputs.PauseTime;
        spawnopts.StartTime = [];
        spawnopts.IgnoreParamFiles = true;
        spawnopts.OutputFilePrefix = Inputs.OutputFilePrefix;
        spawnopts.SlaveNiceProcessPriority = Inputs.SlaveNiceProcessPriority;

        spawnstate = mcore.slavespawn (spawnopts);

        spawnstate = mcore.slavespawn (spawnopts, spawnstate);

    else

        % get the full path to the launch script
        launchscript = fullfile (getmfilepath ('mcore.startnslaves'), '..', launchscript);

        spawnstate = [];

        sleeptime = 0;

        for n = 1:nslaves

            slaveID = randi(99999);

            if strcmp (Inputs.OutputFilePrefix, '/dev/null')
                outfile = '/dev/null';
            else
                outfile = sprintf('%s_%d.txt', fullfile(Inputs.MulticoreSharedDir, Inputs.OutputFilePrefix), slaveID);
            end
            
            if isunix ()
                launchargs = sprintf ( ' ''%s'' [%d,%d,%d,%d,%d,%d] ''%s'' %d %d %d %s > "%s" &', ...
                    Inputs.MulticoreSharedDir, ...
                    Inputs.EndDate(1), Inputs.EndDate(2), Inputs.EndDate(3), Inputs.EndDate(4), Inputs.EndDate(5), Inputs.EndDate(6), ...
                    Inputs.StartDir, ...
                    sleeptime, ...
                    slaveID, ...
                    Inputs.SlaveNiceProcessPriority, ...
                    workspace_file, ...
                    outfile ...
                );
                
                fulllaunchscript = [launchscript, '.sh'];
            else
                fulllaunchscript  = [launchscript, '.cmd'];
            end
        
            cleansystem (['"', fulllaunchscript , '"' launchargs]);
            
            sleeptime = sleeptime + Inputs.PauseTime;
 
        end

    end
    
end
