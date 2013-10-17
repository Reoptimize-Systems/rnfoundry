function mexlsei_setup(forcef2clibrecompile, forcedlseiwrite)
% This file contains the procedure for building the mexlsei function from
% scratch. Not every step in this file may be necessary as some will likely
% have been performed already, such as steps 2 and 3. 
% 

    % optional settings
    if nargin < 1
        forcef2clibrecompile = false;
    end
    
    if nargin < 2
        forcedlseiwrite = false;
    end
    
    origdir = pwd;
    
    %% 1. Change to the installation directory (where this file is)
    % you may need to do this manually if mexlsei_setup is not on the path

    rootdir = fullfile(fileparts(which('mexlsei_setup')));
    cd(rootdir);

    %% 2. Build the f2c lib which must be included if necessary

    cd(fullfile(pwd, 'c', 'libf2c'));

    % get the local machine architecture
    machine = computer('arch');

    % Put the library file to the directory above
    libfiledir = fullfile(rootdir, 'c');

    if strcmp(machine,'win32')

        % To compile use Visual C++ Express Edition 2008
        % Open the visual studio 2008 Command prompt and run the
        % following
        %
        % > nmake -f makefile.vc32

        libfilename = 'vcf2c.lib';
        if ~exist(libfilename, 'file') || forcef2clibrecompile
            fprintf(1, 'Manual steps are required for win32 f2c lib compilation');
            return;
        end

    elseif strcmp(machine , 'mingw32-i686')

        libfilename = 'libf2cmingw32x86.a';
        if ~exist(libfilename, 'file') || forcef2clibrecompile
            % compile using gcc
            system('make -f makefile.u');
            movefile(libfilename, libfiledir);
        end

    elseif strcmp(machine, 'win64')

        % To compile use Visual C++ Express Edition 2008
        % Open the visual studio Cross-Tools x64 Command prompt and run the
        % following
        %
        % > nmake -f makefile.vc64
        libfilename = 'vc64f2c.lib';
        if ~exist(libfilename, 'file') || forcef2clibrecompile
            fprintf(1, 'Manual steps are required for win64 f2c lib compilation');
            return;
        end

    elseif strcmp(machine,'glnxa64') || strcmp(machine, 'gnu-linux-x86_64')

        % see if libf2c is installed as a system library 
        [~,ldconfout] = system('ldconfig -p | grep libf2c');

        libfilename = 'libf2cx64.a';
        if (isempty(ldconfout) && ~exist(libfilename, 'file')) || forcef2clibrecompile
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

        libfilename = 'libf2cx86.a';
        if (isempty(ldconfout) && ~exist(libfilename, 'file'))  || forcef2clibrecompile
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

    %% 3. Convert the fortran code files into a single C file dlsei.c using
    %% f2c

    if ~exist(fullfile(rootdir, 'c', 'dlsei.c'), 'file') || forcedlseiwrite
        
        if ispc
            
            % use compiled f2c on windows
            sysytem('type .\fortran\*.f | ".\c\f2c.exe" -A -R > .\c\dlsei.c');

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
                [ '\\n\\n', ...
                'void MAIN__()\n', ...
                '{\n', ...
                'fprintf(stderr,"If this appears, F2C is in conflict with C\\n");\n', ...
                'exit(10);\n'...
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

    if strcmp(machine,'win32')

        % To compile use Visual C++ Express Edition 2008
        % I had to move the library file (vcf2c.lib) to a directory with no
        % spaces, 'C:\libraries' on my system, you will probably need to change
        % this directory for your windows machine
        mex .\c\mexlsei.c .\c\dlsei.c .\c\i1mach.c .\c\d1mach.c -outdir . -L"C:\libraries" -lvcf2c LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:LIBCMT.lib" -v

    elseif strcmp(machine,'win64')
        % To compile use Visual C++ Express Edition 2008
        % I had to move the library file (vc64f2c.lib) to a directory with no
        % spaces, 'C:\libraries' on my system, you will probably need to change
        % this directory for your windows machine
        mex .\c\mexlsei.c .\c\dlsei.c .\c\i1mach.c .\c\d1mach.c -outdir . -L"C:\libraries" -lvc64f2c LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:LIBCMT.lib" -v

    else

        % is probably linux, or possibly octave on windows
        % compile using gcc
        if isempty(libfiledir)

            mex('./c/mexlsei.c', ...
                './c/dlsei.c', ...
                './c/i1mach.c', ...
                './c/d1mach.c', ...
                ['-l', fullfile(libfiledir, libfilename)], ...
                '-v');
        else
             mex('./c/mexlsei.c', ...
                './c/dlsei.c', ...
                './c/i1mach.c', ...
                './c/d1mach.c', ...
                ['-L"', libfiledir, '"'], ...
                ['-l', fullfile(libfiledir, libfilename)], ...
                '-v');
            
        end
    end

    % restore the current directory
    cd(origdir);
    
%%

%load Test_lsei.mat
%
%copyfile('~/Postgrad_Research/fortran/mexlsei/mexlsei.mexw64', '/home/s0237326/Postgrad_Research/MATLAB_Scripts/subversion/matlab/Useful_Functions/mlsei/')
%
%[x, rnorme, rnorml, mode] = mlsei(A, b, Mineq, rhsineq, Meq, rhseq)

end