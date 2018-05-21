function rnfoundry_setup (varargin)
% setup function for the RenewNet Foundry. Sets paths and attempts to
% compile mexfunctions for best performance of the tools
%
% Syntax
%
% rnfoundry_setup ()
% rnfoundry_setup ('Parameter', value)
%
% Description
%
% rnfoundry_setup sets up the RenewNet Foundry matlab code on a new
% installation. At it's simplest, it will just add some directories to your
% Matlab/Octave path. At it's most complicated it will download missing
% packages required for optimum performance and feature sets, and compile
% from scratch any matlab mex files which need compiled. What it does, or
% attempts to do, will depend on what has already been done before, or what
% has been supplied, or what you tell it to do using various options.
%
% If using Matlab (rather than Octave), to do any of the compilation you
% must have set up a C++ compiler for your Matlab installation in order to
% use the high performance functions. See Matlab's documentation for
% further information on setting up C++ compiler. Good starting points for
% this are
%
% http://uk.mathworks.com/help/matlab/matlab_external/choose-c-or-c-compilers.html
%
% and the 'doc mex' command
%
% Input
%
% rnfoundry_setup can be called with no arguments, or for finer control,
% over the install process, it can be called with additional optional
% arguemtns supplied as Parameter-Value pairs. The avaialble options are:
%
%  'RunTests' : Flag determining whether tor runs some scripts to test the 
%    setup after installation is complete. Default is false.
%
%  'ForceExistfileSetup' : Forces the recompilation of the existfile mex
%    function. A non-mex version will be used if not present. Default is
%    false.
%
%  'SkipExistfileSetup' : Skips compilation of the existfile mex
%    function, even if it is not on the path. A slower non-mex version will
%    be used if not present. Default is false.
%
%  'ForceMexLseiSetup' : Forces the recompilation of the mexlsei mex 
%    function even if it already on the path. mexlsei is not required if
%    your system has the 'quadprog' function. Default is false.
%
%  'SkipMexLseiSetup' : Skips the compilation of the mexlsei mex function.
%    Features requiring this will not be available. mexlsei is not required
%    if your system has the 'quadprog' function.
%
%  'ForceMexLseiF2cLibRecompile' : Forces the recompilation of the f2c
%    library which must be linked to by the mexlsei mex function. Default
%    is false.
%
%  'ForceMexLseiCFileCreation' : Forces the creation of the C language file
%    dlsei from the original fortran sources of dlsei using f2c. Default is
%    false in which case a presupplied version is used.
%
%  'ForceMexSLMSetup' : Forces the recompilation of the mexslmeval mex 
%    function even if it already on the path. Compilation of this function
%    requires the GNU scientific library (libgsl and libgslblas). Default
%    is false.
%
%  'SkipMexSLMSetup' : Skips the recompilation of the mexslmeval mex 
%    function even if it is not on the path. A slower non-compiled version
%    will be used if the compiled version is not present. Default is false.
%
%  'ForceMexPPValSetup' : Forces the recompilation of the mexppval mex 
%    function even if it already on the path. Default is false.
%
%  'SkipMexPPValSetup' : Skips compilation of the mexppval mex function
%    even if it is not on the path. Default is false.
%
%  'ForceMexmPhaseWLSetup' : Forces the recompilation of the mexmPhaseWL
%    mex function even if it already on the path. Default is false.
%
%  'SkipMexmPhaseWLSetup' : Skips the recompilation of the mexmPhaseWL
%    mex function even if it is not on the path. There is no non-mex
%    alternative so winding design funcitons will not work. Default is
%    false.
%
%  'ForceMBDynSetup' : Forces the recompilation of the mbdyn related
%    mex functions even if they are already on the path. Default is false.
%
%  'SkipMBDynSetup' : Skips the recompilation of the mbdyn related
%    mex functions even if it they are not on the path. There is no non-mex
%    alternative multibody dynamics modelling functions will not work.
%    Default is false.
%
%  'ForceAllMex' : Equivalent to setting all the Force* options above to
%    true. Default is false.
%
%  'PreventXFemmCheck' :  Many functions in the renewnet foundry require 
%    the 'xfemm' finite element analysis package. This option determines
%    whether rnfoundry_setup checks to see if xfemm is already installed
%    (by looking for xfemm functions in the path). Default is false, so
%    rnfoundry_setup WILL check to see if xfemm is installed and display a
%    notice if it is not.
%


% the following options no longer apply as https breaks the sourceforge
% download process, might restore some day
%
%  'XFemmInstallPrefix' : Many functions in the renewnet foundry require 
%    the 'xfemm' finite element analysis package. rnfoundry_setup can 
%    download and install this package if desired. By defualt the package
%    will be installed in the same directory as the one containing
%    rnfoundry_setup.m, you can use this option to set this to a different
%    directory.
%
%  'XFemmDownloadSource' : Many functions in the renewnet foundry require 
%    the 'xfemm' finite element analysis package. rnfoundry_setup can
%    download and install this package if desired. To change the default
%    download location (i.e. the remote url pointing to the package on the
%    internet) you can set this option. The default is a location on
%    Sourceforge.net, and depends on your machine architecture. It's fairly
%    unlikely you'll ever want to change this option.
%
%

    % set up the matlab path first to get access to a load of utility
    % functions we can then use
    thisfilepath = fileparts (which ('rnfoundry_setup'));
    addpath(genpath (thisfilepath));
    
    
    workdir = pwd ();
    % restore working directory on error or exit
    CC = onCleanup (@() cd (workdir));
    
    Inputs.RunTests = false;
    Inputs.Verbose = false;
    Inputs.ForceAllMex = false;
    % mexlsei related
    Inputs.ForceMexLseiSetup = false;
    Inputs.SkipMexLseiSetup = false;
    Inputs.ForceMexLseiF2cLibRecompile = false;
    Inputs.ForceMexLseiCFileCreation = false;
    % slm fitting tool related
    Inputs.ForceMexSLMSetup = false;
    Inputs.SkipMexSLMSetup = false;
    % mbdyn
    Inputs.ForceMBDynSetup = false;
    Inputs.SkipMBDynSetup = false;
    [Inputs.MBCLibDir, Inputs.MBCIncludeDir] = mbdyn.mint.find_libmbc ();
    if isunix
        Inputs.GSLLibDir = ''; % for mexslmeval
        Inputs.GSLIncludeDir = ''; % for mexslmeval
        Inputs.F2CLibPath = ''; % for mlse
    else
        Inputs.GSLLibDir = fullfile (thisfilepath, 'x86_64-w64-mingw32', 'lib'); % for mexslmeval
        Inputs.GSLIncludeDir = fullfile (thisfilepath, 'x86_64-w64-mingw32', 'include'); % for mexslmeval
        if isoctave
            Inputs.F2CLibPath = fullfile (thisfilepath, 'x86_64-w64-mingw32', 'lib', 'libf2c.a'); % for mlse
        else
            Inputs.F2CLibPath = fullfile (thisfilepath, 'x86_64-w64-mingw32', 'lib', 'f2c.lib'); % for mlse
        end
    end
    % mex ppval related
    Inputs.ForceMexPPValSetup = false;
    Inputs.SkipMexPPValSetup = false;
    % force setting up mexmPhaseWL
    Inputs.ForceMexmPhaseWLSetup = false;
    Inputs.SkipMexmPhaseWLSetup = false;
    % xfemm related
    Inputs.PreventXFemmCheck = false;
%     if ispc 
%         Inputs.XFemmDownloadSource = 'https://downloads.sourceforge.net/project/xfemm/Release/Release%201.8/xfemm_v1_8_mingw_win64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fxfemm%2Ffiles%2FRelease%2FRelease%25201.8%2F&ts=1488547864&use_mirror=netix';
%     elseif isunix
%         Inputs.XFemmDownloadSource = 'https://downloads.sourceforge.net/project/xfemm/Release/Release%201.8/xfemm_v1_8_linux64.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fxfemm%2Ffiles%2FRelease%2FRelease%25201.8%2F&ts=1488547898&use_mirror=netcologne';
%     else
%         Inputs.XFemmDownloadSource = '';
%     end
    Inputs.XFemmInstallPrefix = fullfile (thisfilepath, 'common');
    Inputs.ForceExistfileSetup = false;
    Inputs.SkipExistfileSetup = false;
    
    % now parse the pv pairs
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if Inputs.ForceAllMex
        Inputs.ForceMexLseiSetup = true;
        Inputs.ForceMexSLMSetup = true;
        Inputs.ForceMBDynSetup = true;
        Inputs.ForceMexPPValSetup = true;
    end
    
    didcompwarn = false;
    
    %% existfile
    if Inputs.ForceExistfileSetup && Inputs.SkipExistfileSetup
        error ('The options ForceExistfileSetup and SkipExistfileSetup are both set to true');
    end
    
    if ~Inputs.SkipExistfileSetup
        if Inputs.ForceExistfileSetup || (exist (['existfile.', mexext], 'file') ~= 3)
            didcompwarn = compilerwarning (didcompwarn);
            % set up existfile mex function
            mexcompile_existfile ('Verbose', Inputs.Verbose);
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as it already exists\n', 'existfile')
            end
        end
    end
    
    %% mexlsei
    if Inputs.ForceMexLseiSetup && Inputs.SkipMexLseiSetup
        error ('The options ForceMexLseiSetup and SkipMexLseiSetup are both set to true');
    end
    
    if ~Inputs.SkipMexLseiSetup
        if Inputs.ForceMexLseiSetup || (exist (['mexlsei.', mexext], 'file') ~= 3)
            didcompwarn = compilerwarning (didcompwarn);
            mexlsei_setup ( 'ForceF2CLibRecompile', Inputs.ForceMexLseiF2cLibRecompile, ...
                            'ForceLSEIWrite', Inputs.ForceMexLseiCFileCreation, ...
                            'F2CLibFilePath', Inputs.F2CLibPath, ...
                            'Verbose', Inputs.Verbose );
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as it already exists\n', 'lsei')
            end
        end
    end
    
    %% mexslmeval
    if Inputs.ForceMexSLMSetup && Inputs.SkipMexSLMSetup
        error ('The options ForceMexSLMSetup and SkipMexSLMSetup are both set to true');
    end
    
    if ~Inputs.SkipMexSLMSetup
        if Inputs.ForceMexSLMSetup || (exist (['mexslmeval.', mexext], 'file') ~= 3)
            didcompwarn = compilerwarning (didcompwarn);
            mexslmeval_setup ( 'Verbose', Inputs.Verbose, ...
                               'GSLLibDir', Inputs.GSLLibDir, ...
                               'GSLIncludeDir', Inputs.GSLIncludeDir);
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as it already exists\n', 'slmeval')
            end
        end
    end
    
    %% mbdyn
    if Inputs.ForceMBDynSetup && Inputs.SkipMBDynSetup
        error ('The options ForceMBDynSetup and SkipMBDynSetup are both set to true');
    end
    
    % octave's which function doesn't currently work for functions in
    % package directories, so wrap the call to getmfilepath in a try catch
    % with a fallback. If Octave is fixed this will then start working
    % again
    try
        mexMBCNodalexists = exist (fullfile (getmfilepath ('mbdyn.mint.MBCNodal'), ['mexMBCNodal.', mexext()]), 'file') == 3;
    catch
        % this version is more fragile to files/directories being moved
        % around
        mexMBCNodalexists = exist ( fullfile ( mbdyn_rootdir (), ...
                                              '+mbdyn', ...
                                              '+mint', ...
                                              ['mexMBCNodal.', mexext()]), ...
                                    'file' ) == 3;
    end
    
    if ~Inputs.SkipMBDynSetup
        if Inputs.ForceMBDynSetup || ~mexMBCNodalexists
            didcompwarn = compilerwarning (didcompwarn);
            mexmbdyn_setup ( 'Verbose', Inputs.Verbose, ...
                             'MBCLibDir', Inputs.MBCLibDir, ...
                             'MBCIncludeDir', Inputs.MBCIncludeDir);
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as it already exists\n', 'mbdyn')
            end
        end
    end
    
    %% ppmval ppuval
    if Inputs.ForceMexPPValSetup && Inputs.SkipMexPPValSetup
        error ('The options ForceMexPPValSetup and SkipMexPPValSetup are both set to true');
    end
    
    if ~Inputs.SkipMexPPValSetup
        if Inputs.ForceMexPPValSetup ...
                || (exist (['ppmval.', mexext], 'file') ~= 3)...
                || (exist (['ppuval.', mexext], 'file') ~= 3)
            didcompwarn = compilerwarning (didcompwarn);
            mexppval_setup ('Verbose', Inputs.Verbose);
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as they already exists', 'ppmval and ppuval')
            end
        end
    end
    
    %% mexmPhaseWL
    if Inputs.ForceMexmPhaseWLSetup && Inputs.SkipMexmPhaseWLSetup
        error ('The options ForceMexmPhaseWLSetup and SkipMexmPhaseWLSetup are both set to true');
    end
    
    if ~Inputs.SkipMexmPhaseWLSetup
        if Inputs.ForceMexmPhaseWLSetup || (exist (['mexmPhaseWL.', mexext], 'file') ~= 3)
            didcompwarn = compilerwarning (didcompwarn);
            mmake.make ('', fullfile (pm_machines_tools_rootdir (), 'common', 'winding-layout', 'MMakefile.m'));
            mmake.make ('tidy', fullfile (pm_machines_tools_rootdir (), 'common', 'winding-layout', 'MMakefile.m'));
        else
            if Inputs.Verbose
                fprintf (1, 'Not compiling %s mex as it already exists\n', 'mPhaseWL')
            end
        end
    end
    
    xfemm_main_page_url = 'https://sourceforge.net/projects/xfemm/';
    
    % check for the existence of xfemm package
    if ~Inputs.PreventXFemmCheck
        
        if exist ('mexfmesher', 'file') ~= 3
            
            fprintf (1, [ ...
                    '\n\n', ...
                    '**************************    XFEMM  NOTICE    **************************\n\n', ...
                    'You do not appear to have the xfemm package for performing electromagnetic\n', ...
                    'simulations which is required for many functions in the foundry\n', ...
                    'to work. You can obtain this package from Sourceforce here:\n\n', ...
                    xfemm_main_page_url, ...
                    '\n\n', ...
                    '**************************    XFEMM NOTICE    **************************\n\n' ]);
                
%             response = '';
%             while ~(strcmpi (response, 'Y') || strcmpi (response, 'N') )
%                 response = input( sprintf ([ ...
%                     '\n\n', ...
%                     '**************************     NOTICE    **************************\n', ...
%                     'You do not appear to have the xfemm package for performing electromagnetic\n', ...
%                     'simulations which is required for many functions in the foundry\n', ...
%                     'to work. Do you want to try to download and install it? (y/n): ']), 's');
%             end
%             
%             if upper(response) == 'Y'
%                 
%                 doxfemm = true;
%                 if ispc
%                     xfemm_prefix = 'xfemm_mingw_win64';
%                     xfemmfile = fullfile (tempdir, [xfemm_prefix, '.zip']);
%                     urlwrite ( Inputs.XFemmDownloadSource, xfemmfile );
%                     unpackfcn = @unzip;
%                 elseif isunix
%                     xfemm_prefix = 'xfemm_linux64';
%                     xfemmfile = fullfile (tempdir, [xfemm_prefix, '.tar.gz']);
%                     urlwrite ( Inputs.XFemmDownloadSource, xfemmfile );
%                     unpackfcn = @untar;
%                 else
%                     fprintf ('No xfemm compiled package is currently available for mac, skipping.\n');
%                     doxfemm = false;
%                 end
% 
%                 if doxfemm
%                     
%                     fprintf (1, 'Downloading xfemm package from:\n%s\n', Inputs.XFemmDownloadSource);
%                     
%                     % unpack the download to the appropriate location
%                     unpackfcn (xfemmfile, Inputs.XFemmInstallPrefix);
% 
%                     % remove the downloaded package
%                     delete (xfemmfile);
% 
%                     % add mfemm setup function location to path and run it
%                     addpath (fullfile (Inputs.XFemmInstallPrefix, xfemm_prefix, 'mfemm'));
%                     mfemm_setup ();
%                     
%                     fprintf (1, 'xfemm package installation complete\n');
% 
%                 end
%             else
%                 fprintf ('Skipping xfemm install. You can obtain this package from Sourceforge later if you wish.\n');
%             end
            
        end
    end
    
    cd (workdir);
    
    
    if ~isoctave
        % print message about adding files to the path (octave currently
        % doesn't support displaying help in subfunctions)
        fprintf(1, '\n\n');
        help rnfoundry_setup>addedtopathmsg

        fprintf (1, '\n');
    
    end
    
    if Inputs.RunTests
        runtests ();
    end
    
end


function runtests ()
    % run some examples to test the installation
    
    run_one_test ('example_basic_heaving_buoy_simulation');
    close all;
    
    run_one_test ('example_buoy_sim_with_ACTM');
    close all;
    
    run_one_test ('example_radial_flux_permanent_magnet_machine_sim');
    close all;
    
    run_one_test ('example_radial_flux_pm_with_ratio_specification');
    close all;
    
    if isoctave
        run_one_test ('example_rm3.run_wecsim_in_octave');
    else
        run_one_test ('example_rm3.run_wecsim');
    end
    close all;
    
end


function run_one_test (script_name)

    fprintf (1, 'Running test script: "%s"\n', script_name);
    
    try
        eval (script_name);
    catch err
        fprintf (1, 'The following test script: \n"%s"\nresulted in the following error:\n%s\n', err.message);
    end

end


function addedtopathmsg ()
% rnfoundry_setup has added the directories containing mfiles to the search
% path. This change will not be saved when you close Matlab/Octave. You
% will need to have these directories on the path in future for the
% RenewNet Foundry code to work when you restart Matlab/Octave. You can add
% them again by using the 'Set path' dialog in Matlab, or using the
% 'addpath' command, or by rerunning rnfoundry_setup again, which calls
% addpath internally to do this. rnfoundry_setup will not recompile the mex
% functions if they already exist and will just modify your path.
%
% If you open the 'Set path' dialog now, you could also save the changes
% that have just been made.
%
% Run the command 'help path' for more information on the path. In general
% you can save changes to the Matlab path that persist across sessions
% using the 'Set path' dialog, or you can make a startup.m file which runs
% every time Matlab starts (.octaverc in Octave) containing the addpath
% commands.
%

end

function didcompwarn = compilerwarning (didcompwarn)

    if (didcompwarn == false) && (isoctave == false)
        cc = mex.getCompilerConfigurations('C++');
        
        if numel (cc) == 0
           warning ( ['The renewnet foundry code will have best performance if ' ...
                      'you have set up a C++ compiler for matlab using "mex -setup". ' ...
                      'No C++ compiler seems to have been set up on your system yet. ' ... 
                      'You may wish to set this up and re-run rnfoundry_setup.']) 
        end
    end
    
    didcompwarn = true;
    
end
