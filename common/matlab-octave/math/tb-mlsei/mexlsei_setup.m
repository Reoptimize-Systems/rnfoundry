function mexlsei_setup(varargin)
% This file contains the procedure for building the mexlsei function from
% scratch. Not every step in this file may be necessary as some will likely
% have been performed already, such as steps 2 and 3. 
% 

    rootdir = fullfile(fileparts(which('mexlsei_setup')));
    
    
    options.F2CLibFilePath = '';
    options.IncludePath = fullfile(rootdir, 'c');
    options.ForceF2CLibRecompile = false;
    options.ForceLSEIWrite = false;
    options.Verbose = false;
    options.ExtraMexArgs = {};
    options.MexExtension = mexext ();
    options.W64CrossBuild = false;
    options.ThrowBuildErrors = false;
    
    options = parse_pv_pairs (options, varargin);
    
    fprintf (1, 'Setting up mexlsei\n');
    
    origdir = pwd;
    
    % return to original directory on exit or failure
    CC = onCleanup (@() cd(origdir));
    
    %% 1. Change to the installation directory (where this file is)
    % you may need to do this manually if mexlsei_setup is not on the path

    cd (rootdir);

    %% 2. Build the f2c lib which must be included if necessary

    % get the local machine architecture
    machine = computer ('arch');
    

    if isempty (options.F2CLibFilePath)
    
        % Put the library file to the directory above
        libfiledir = fullfile(rootdir, 'c');
        
        if options.W64CrossBuild

                % will be using mxe to cross-compile from linux
                libfilename = 'libf2c';
                libfiledir = '';
                
        elseif strcmp(machine,'win32')

            % To compile use Visual C++ Express Edition 2008
            % Open the visual studio 2008 Command prompt and run the
            % following
            %
            % > nmake -f makefile.vc32

            libfilename = 'vcf2c.lib';
            if ~exist(libfilename, 'file') || options.ForceF2CLibRecompile
                warning ('Manual steps are required for win32 f2c lib compilation, see comments in mexlsei_setup.m. Exiting mexlsei setup now.');
                return;
            end

        elseif strcmp(machine , 'mingw32-i686')
            % is probably Octave on 32 bit windows
            libfilename = 'libf2cmingw32x86.a';
            if ~exist(libfilename, 'file') || options.ForceF2CLibRecompile
                cd(fullfile(libfiledir, 'libf77'));
                % compile using gcc
                system('mingw32-make -f makefile.u');
                movefile(libfilename, libfiledir);
            end

        elseif strcmp(machine , 'mingw32-x86_64')
            % is probably Octave on 64 bit windows
            libfilename = 'libf2cmingw32x86_64.a';
            if ~exist(libfilename, 'file') || options.ForceF2CLibRecompile
                cd(fullfile(libfiledir, 'libf77'));
                % compile using gcc
                system('mingw32-make -f makefile.u');
                movefile(libfilename, libfiledir);
            end

        elseif strcmp(machine, 'win64')

            % To compile use Visual C++ Express Edition 2008
            % Open the visual studio Cross-Tools x64 Command prompt and run the
            % following
            %
            % > nmake -f makefile.vc64
            libfilename = 'vc64f2c.lib';
            if ~exist(libfilename, 'file') || options.ForceF2CLibRecompile
                warning ('Manual steps are required for win64 f2c lib compilation, see comments in mexlsei_setup.m. Exiting mexlsei setup now.');
                return;
            end

        elseif strcmp(machine,'glnxa64') || strcmp(machine, 'gnu-linux-x86_64')

            % see if libf2c is installed as a system library 
            [~,ldconfout] = system('ldconfig -p | grep libf2c');

            if isempty (ldconfout)
                res = f2cdocontinue ();
                if ~res, return; end
            end

            libfilename = 'libf2cx64.a';
            if (isempty(ldconfout) && ~exist(libfilename, 'file')) || options.ForceF2CLibRecompile
                cd(fullfile(libfiledir, 'libf77'));
                % compile using gcc on 64 bit platform
                system('make -f makefile.lx64');
                movefile(libfilename, libfiledir);
            elseif ~isempty(ldconfout)
                libfilename = 'libf2c';
                libfiledir = '';
            end

        elseif strcmp(machine,'glnx86')

            % see if libf2c is installed as a system library 
            [~,ldconfout] = system('ldconfig -p | grep libf2c');

            if isempty (ldconfout)
                res = f2cdocontinue ();
                if ~res, return; end
            end

            if ~res, return; end

            libfilename = 'libf2cx86.a';
            if (isempty(ldconfout) && ~exist(libfilename, 'file'))  || options.ForceF2CLibRecompile
                cd(fullfile(libfiledir, 'libf77'));
                % compile using gcc
                system('make -f makefile.lx86');
                movefile(libfilename, libfiledir);
            elseif ~isempty(ldconfout)
                libfilename = 'libf2c';
                libfiledir = '';
            end

        end

        % change back to mexlsei_setup directory
        cd(rootdir);
        
        includepath = libfiledir;
    else
        % use provided library file details
        [libfiledir, libfilename] = fileparts (options.F2CLibFilePath);
        
        if isempty (options.IncludePath)
            includepath = libfiledir;
        else
            includepath = options.IncludePath;
        end
    end

    %% 3. Convert the fortran code files into a single C file dlsei.c using
    %% f2c

    if ~exist(fullfile(rootdir, 'c', 'dlsei.c'), 'file') || options.ForceLSEIWrite
        
        if ispc
            
            % use compiled f2c on windows
            system('type .\fortran\*.f | ".\c\f2c.exe" -A -R > .\c\dlsei.c');

        else

            [~,out] = system('which f2c');
            if isempty(out)
                error('The program f2c must be available.')
            end
            system('cat ./fortran/*.f | f2c -A -R > ./c/dlsei.c');

        end
        
        % for some reason, the F2C library on recent (1996+) Linux systems
        % includes a reference to the MAIN__ entry point so any program
        % compiled with the f2c library needs such an entry point, even
        % though it's not used.  (ESF 1997.04.07), so we make sure the
        % symbol is there
        fid = fopen(fullfile(rootdir, 'c', 'dlsei.c'), 'a');
            
        fprintf(fid, ...
                [ '\n\n', ...
                'void MAIN__(void)\n', ...
                '{\n', ...
                '/*fprintf(stderr,"If this appears, F2C is in conflict with C\\n");\n', ...
                'exit(10);*/\n'...
                'return;', ...
                '}', ...
                ] );
            
        fclose(fid);
        
    end

    %% 4. Now mex the files, including the appropriate lib

    % the c functions d1mach.c and i1mach.c replace the original fortran
    % versions from the slatec lib. Their source was contained in the comments
    % of fortran versions of the same functions available from the BLAS
    % colection on netlib here:
    %
    % http://www.netlib.org/blas/
    %
    % or more directly:
    % 
    % http://www.netlib.org/blas/i1mach.f
    % http://www.netlib.org/blas/r1mach.f
    % http://www.netlib.org/blas/d1mach.f
    %

    if isunix || isoctave
        ismscompiler = false;
    else
        cc = mex.getCompilerConfigurations ('C');
        if strncmpi (cc.Manufacturer, 'Microsoft', 9)
            ismscompiler = true;
        else
            ismscompiler = false;
        end
    end
    
    if strcmpi (libfilename(1:3), 'lib')
        
        % we force linking to the static library. This fixes the
        % "libf2c.so: undefined reference to `MAIN__' " error
        libfilemexcallname = [':', libfilename, '.a'];

    elseif strcmpi (libfilename(end-4:end), '.lib')
        libfilemexcallname = libfilename(1:3);
    else
        libfilemexcallname = libfilename;
    end
    
    mexinputs = { './c/mexlsei.c', ...
                  './c/dlsei.c', ...
                  './c/i1mach.c', ...
                  './c/d1mach.c', ...
                  sprintf('-l%s', libfilemexcallname) };
              
	if ~isoctave ()
        mexinputs = [mexinputs, {['EXE="existfile.', options.MexExtension, '"']}];
    end
              
    if ~isempty (includepath)
        mexinputs = [mexinputs, { sprintf('-I"%s"', includepath) }];
    end
    
    if ~isempty (libfiledir)
        mexinputs = [mexinputs, { sprintf('-L"%s"', libfiledir) }];
    end
    
    if options.Verbose
        mexinputs = [mexinputs, { '-v' }];
    end
    
    if ismscompiler
        
        % To compile using microsoft compiler
        % I had to move the library file (vcf2c.lib) to a directory with no
        % spaces, 'C:\libraries' on my system, you might need to change
        % this directory for your windows machine
        mexinputs = [mexinputs, { 'LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:LIBCMT.lib"'}];
    end
    
    % call mex
    try
        mex(mexinputs{:}, options.ExtraMexArgs{:});
        
        fprintf (1, 'Finished setting up mexlsei\n');
    catch err
        if options.ThrowBuildErrors == true
            rethrow (err);
        else
            warning ( 'mexlsei compilation falied, error message was:\n%s', ...
                      err.message );
        end
    end
    
end


function res = f2cdocontinue ()

    S = input ( sprintf( ...
['f2c is not installed in your system, it is advisable to quit and install it\n', ...
 'using your packaging system (if on linux, or find and download it for windows).\n', ...
 'Otherwise mexlsei_setup will now attempt to build f2c and libf2c from source\n', ...
 'before continuing.\n', ...
 '\n', ...
 'Do you wish to continue? (Y/n)'] ), 's');
 
    if strcmpi (S, 'n')
        res = false;
    else
        res = true;
    end

end
