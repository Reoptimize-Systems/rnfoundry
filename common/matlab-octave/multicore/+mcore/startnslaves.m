function startnslaves (nslaves, varargin)
% launches slave processes for multicore, on same machine as master process

    Inputs.MulticoreSharedDir = mcore.defaultmulticoredir ();
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
            fulllaunchscript  = fullfile (getmfilepath (mfilename ()), [launchscript, '.cmd']);
        end
    
        system (['"', fulllaunchscript , '"' launchargs]);
        
        sleeptime = sleeptime + Inputs.PauseTime;

    end
    
end