function [mbstatus, cmdout, pid] = start_mbdyn (inputfile, varargin)
% runs mbdyn with the appropriate commands
%
% Syntax
%
% [status, cmdout] = start_mbdyn (inputfile)
% [status, cmdout] = start_mbdyn (..., 'Parameter', value)
%
% Description
%
% Finds and rund the MBDyn executable with the provided file as the MBDyn
% input file.
%
% Input
%
%  inputfile - character vector containing the path to the MBDyn input file
%   containing the problem to be solved by MBDyn.
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'StartWaitTime' - optional scalar value with a time in seconds for which
%    to pause after starting MBDyn. This allows MBDyn to initialise and/or
%    run the problem before the user's other code loads the solution or
%    begins comunicating with MBDyn. Default is 2 if not supplied.
%
%  'MBDynExecutable' - optional character vector containing the full path
%    to the MBDyn executable to use. If not supplied, start_mbdyn will look
%    in the various normal places for this on typical systems. To learn
%    more about where start_mbdyn will look for MBDyn if this option is not
%    supplied, see the function mbdyn.mint.find_mbdyn.
%
%  'OutputPrefix' - optional output file prefix for MBDyn to store the
%    solution in. This is a file path without any extension, e.g.
%    /home/jbloggs/my_mbdyn_sim_output. MBDyn will create multiple files
%    with various extensions like /home/jbloggs/my_mbdyn_sim_output.log and
%    /home/jbloggs/my_mbdyn_sim_output.mov etc. If no OutputPrefix is
%    supplied, the input file is used with any file extension stripped.
%
%  'MBDynOutputFile' - optional file path to which MBDyn command line
%    output will be redirected for later inspection. If not supplied, the
%    input file path is used, with any file name extension stripped, and
%    the extension '.txt' added.
%
%  'Verbosity' - optional scalar integer indicating the amount of output
%    MBDyn should produce at the command line (and therefore in the
%    MBDynOutputFile). The value of Verbosity indicates how many -P flags
%    are passed to MBDyn on the command line when starting. The more -P's
%    the output MBDyn produces. Default is 0, for no additional output.
%
%  'Block' - optional true/false flag indicating whether start_mbdyn will
%    wait for MBDyn to complete the requested operation, or run it
%    asynchronously in the background. Deafault is true, so start_mbdyn
%    will wait while MBDyn runs the command.
%
%
% Output
%
%  status - the status returned by the operating system after attempting to
%   run the MBDyn command.
%
%  cmdout - string containing the std output produced by the MBDyn command
%
%
%
% See Also: 
%

    options.Verbosity = 0;
    options.StartWaitTime = 0.1;
    options.MBDynExecutable = mbdyn.mint.find_mbdyn (false);
    options.MBDynOutputFile = '';
    options.OutputPrefix = '';
    options.Block = true;
    
    options = parse_pv_pairs (options, varargin);
    
    mbdyn.pre.base.checkScalarInteger (options.Verbosity, true, 'Verbosity');
    mbdyn.pre.base.checkNumericScalar (options.StartWaitTime, true, 'StartWaitTime');
    assert (ischar (options.MBDynExecutable), 'MBDynExecutable must be a character vector');
    assert (ischar (options.MBDynOutputFile), 'MBDynOutputFile must be a character vector');
    assert (ischar (options.OutputPrefix), 'OutputPrefix must be a character vector');
    mbdyn.pre.base.checkLogicalScalar (options.Block, true, 'Block');
    
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
    if ispc
        
        if options.Block
            blockstr = '/wait';
        else
            blockstr = '';
        end
        
        % see https://superuser.com/questions/338277/windows-cmd-batch-start-and-output-redirection
        % and https://ss64.com/nt/start.html
        cmdline = sprintf ( 'start "mbdyn" /B %s "%s" %s -f "%s" -o "%s" ^1^> "%s" ^2^>^&^1', ...
                            blockstr, ...
                            options.MBDynExecutable, ...
                            Pcmds, ...
                            inputfile, ...
                            options.OutputPrefix, ...
                            options.MBDynOutputFile  ...
                          );
    else
        
        if options.Block
            blockstr = '';
        else
            blockstr = '&';
        end
        
        % see https://serverfault.com/questions/205498/how-to-get-pid-of-just-started-process
        cmdline = sprintf ( 'sh -c ''echo $$; exec "%s" %s -f "%s" -o "%s" > "%s" 2>&1 %s''', ...
                            options.MBDynExecutable, ...
                            Pcmds, ...
                            inputfile, ...
                            options.OutputPrefix, ...
                            options.MBDynOutputFile,  ...
                            blockstr ...
                          ); 
    end

    if options.Verbosity > 0
        fprintf (1, 'Starting MBDyn with command:\n%s\n', cmdline);
    end

    [mbstatus, cmdout] = mbdyn.mint.cleansystem ( cmdline );
    
    pid = [];
    
    if mbstatus == 0
        
        if ispc
            pid = [];
        else
            pid = str2double (cmdout);
            
            if ~options.Block
                
                % check the PID is actually mbdyn
                [status, cmdout] = mbdyn.mint.cleansystem ( sprintf ('ps -p %s -o comm=', int2str (pid)) );

                if strcmpi (cmdout, 'mbdyn')
                    % do nothing, we've got the mbdyn pid
                else
                    % try pid + 1
                    pid = pid + 1;
                    
                    [status, cmdout] = mbdyn.mint.cleansystem ( sprintf ('ps -p %s -o comm=', int2str (pid)) );

                    if ~strncmpi (cmdout, 'mbdyn', 5)
                        pid = [];
                    end
                    
                end
                
            end
            
        end

    end

    pause (options.StartWaitTime);

end