function [mbstatus, cmdout, pid, inputfile, cmdline] = start_mbdyn (input, varargin)
% runs mbdyn with the appropriate commands
%
% Syntax
%
% [status, cmdout] = mbdyn.mint.start_mbdyn (inputfile)
% [status, cmdout] = mbdyn.mint.start_mbdyn (input_system)
% [status, cmdout] = mbdyn.mint.start_mbdyn (..., 'Parameter', value)
%
% Description
%
% Finds and runs the MBDyn executable with the provided file as the MBDyn
% input file, or generates an input file if an mbdyn.pre.system object is
% provided.
%
% Input
%
%  inputfile - character vector containing the path to the MBDyn input file
%   containing the problem to be solved by MBDyn.
%
%  input_system - an mbdyn.pre.system object, the generateMBDynInputFile
%   will be called to create the MBDyn input file in the system temporary
%   directory. The location of the created file will be returned in
%   the 'inputfile' output argument.
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
%    asynchronously in the background. Default is true, so start_mbdyn
%    will wait while MBDyn runs the command.
%
%  'ThrowExceptions' - optional true/false flag indicating whether MBDyn
%    will throw exceptions (can be useful for debugging). Default is false. 
%
%
% Output
%
%  status - the status returned by the operating system after attempting to
%   run the MBDyn command.
%
%  cmdout - character vector containing the std output produced by the
%    MBDyn command
%
%  pid - 
%
%  inputfile - character vector containing the path of the input file used
%    when starting MBDyn
%
%  cmdline - character vector containing the command used to launch mbdyn
%
% See Also: mbdyn.pre.system.generateMBDynInputFile
%

    options.Verbosity = 0;
    options.StartWaitTime = 0.1;
    options.MBDynExecutable = mbdyn.mint.find_mbdyn (false);
    options.MBDynOutputFile = '';
    options.OutputPrefix = '';
    options.Block = true;
    options.ThrowExceptions = false;
    
    options = parse_pv_pairs (options, varargin);
    
    mbdyn.pre.base.checkScalarInteger (options.Verbosity, true, 'Verbosity');
    mbdyn.pre.base.checkNumericScalar (options.StartWaitTime, true, 'StartWaitTime');
    assert (ischar (options.MBDynExecutable), 'MBDynExecutable must be a character vector');
    assert (ischar (options.MBDynOutputFile), 'MBDynOutputFile must be a character vector');
    assert (ischar (options.OutputPrefix), 'OutputPrefix must be a character vector');
    mbdyn.pre.base.checkLogicalScalar (options.Block, true, 'Block');
    
    if isa (input, 'mbdyn.pre.system')
        inputfile = [ tempname, '.mbd'];
        input.generateMBDynInputFile (inputfile);
    elseif ischar (input)
        if exist (input, 'file') == 2
            inputfile = input;
        else
            error ('The input file %s does exist.', input);
        end
    else
        error ('Input must be a file name or a mbdyn.pre.system objet');
    end
    
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
    
    if options.ThrowExceptions
        excepstr = '-e';
    else
        excepstr = '';
    end
    
    if isoctave ()
        options.OutputPrefix = make_absolute_filename (options.OutputPrefix);
        options.MBDynOutputFile = make_absolute_filename (options.MBDynOutputFile);
        inputfile = make_absolute_filename (inputfile);
    end

    % start mbdyn
    if ispc
        
        if options.Block
            blockstr = '/wait';
        else
            blockstr = '';
        end
        
        if isoctave ()
          
            % see https://superuser.com/questions/338277/windows-cmd-batch-start-and-output-redirection
            % and https://ss64.com/nt/start.html
            cmdline = sprintf ( '"%s" %s %s -f "%s" -o "%s" ^1^> "%s" ^2^>^&^1', ...
                                options.MBDynExecutable, ...
                                Pcmds, ...
                                excepstr, ...
                                inputfile, ...
                                options.OutputPrefix, ...
                                options.MBDynOutputFile  ...
                              );

        else
        
            % see https://superuser.com/questions/338277/windows-cmd-batch-start-and-output-redirection
            % and https://ss64.com/nt/start.html
            %cmdline = sprintf ( 'start "mbdyn" /B %s "%s" %s %s -f "%s" -o "%s" ^1^> "%s" ^2^>^&^1', ...
            cmdline = sprintf ( 'start "mbdyn" %s "%s" %s %s -f "%s" -o "%s" ^1^> "%s" ^2^>^&^1', ...
                                blockstr, ...
                                options.MBDynExecutable, ...
                                Pcmds, ...
                                excepstr, ...
                                inputfile, ...
                                options.OutputPrefix, ...
                                options.MBDynOutputFile  ...
                              );
        end
    else
        
        if options.Block
            blockstr = '';
        else
            blockstr = '&';
        end
        
        if isoctave ()
            % see https://serverfault.com/questions/205498/how-to-get-pid-of-just-started-process
            cmdline = sprintf ( 'sh -c ''echo $$; exec "%s" %s %s -f "%s" -o "%s" > "%s" 2>&1''', ...
                                options.MBDynExecutable, ...
                                Pcmds, ...
                                excepstr, ...
                                inputfile, ...
                                options.OutputPrefix, ...
                                options.MBDynOutputFile ...
                              ); 
        else
            % see https://serverfault.com/questions/205498/how-to-get-pid-of-just-started-process
            cmdline = sprintf ( 'sh -c ''echo $$; exec "%s" %s %s -f "%s" -o "%s" > "%s" 2>&1 %s''', ...
                                options.MBDynExecutable, ...
                                Pcmds, ...
                                excepstr, ...
                                inputfile, ...
                                options.OutputPrefix, ...
                                options.MBDynOutputFile,  ...
                                blockstr ...
                              ); 
        end
    end

    if options.Verbosity > 0
        fprintf (1, 'Starting MBDyn with command:\n%s\n', cmdline);
    end

    if isoctave ()
        if options.Block
            % Ah! the delights of Octave!
            [mbstatus, cmdout] = system ( cmdline, true, 'sync' );
            if mbstatus > -1
                pid = mbstatus;
                mbstatus = 0;
            else
                pid = [];
            end
        else
            % Ah! the delights of Octave!
            mbstatus = system ( cmdline, false, 'async' );
            cmdout = '';
            if mbstatus > -1
                pid = mbstatus;
                mbstatus = 0;
            else
                pid = [];
            end
        end
    else
        [mbstatus, cmdout] = mbdyn.mint.cleansystem ( cmdline );
    end
    
    pid = [];
    
    if mbstatus == 0
        
        if ispc
            if ~isoctave ()
                pid = [];
            end
        else
            pid = str2double (cmdout);
            
%             if ~options.Block
%                 
%                 % check the PID is actually mbdyn
%                 [status, cmdout] = mbdyn.mint.cleansystem ( sprintf ('ps -p %s -o comm=', int2str (pid)) );
% 
%                 if strcmpi (cmdout, 'mbdyn')
%                     % do nothing, we've got the mbdyn pid
%                 else
%                     % try pid + 1
%                     pid = pid + 1;
%                     
%                     [status, cmdout] = mbdyn.mint.cleansystem ( sprintf ('ps -p %s -o comm=', int2str (pid)) );
% 
%                     if ~strncmpi (cmdout, 'mbdyn', 5)
%                         pid = [];
%                     end
%                     
%                 end
%                 
%             end
            
        end

    end

    pause (options.StartWaitTime);

end