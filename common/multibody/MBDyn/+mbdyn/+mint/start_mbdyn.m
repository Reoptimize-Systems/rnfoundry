function [status, cmdout] = start_mbdyn (inputfile, varargin)
% run mbdyn with the appropriate commands

    options.Verbosity = 0;
    options.StartWaitTime = 2;
    options.MBDynExecutable = mbdyn.mint.find_mbdyn (false);
    options.MBDynOutputFile = '';
    options.OutputPrefix = '';
    
    options = parse_pv_pairs (options, varargin);
    
    if ~exist (options.MBDynExecutable, 'file')
        error ('MBDyn executable was not found in the specified location.');
    end
    
    if isempty (options.OutputPrefix)
         [pathstr,name] = fileparts (inputfile);
         options.OutputPrefix = fullfile (pathstr,name);
    end
    
    if isempty (options.MBDynOutputFile)
        options.MBDynOutputFile = [options.OutputPrefix, '.txt'];
    end
    
    if options.Verbosity > 0
        Pcmds = ['-', repmat('P', [1, options.Verbosity])];
    else
        Pcmds = '';
    end

    % start mbdyn
    cmdline = sprintf ('%s %s -f "%s" -o "%s" > "%s" 2>&1 &', ...
                        options.MBDynExecutable, ...
                        Pcmds, ...
                        inputfile, ...
                        options.OutputPrefix, ...
                        options.MBDynOutputFile  ...
                                         );

    if options.Verbosity > 0
        fprintf (1, 'Starting MBDyn with command:\n%s\n', cmdline);
    end

    [status, cmdout] = cleansystem ( cmdline );

    pause (options.StartWaitTime);

end