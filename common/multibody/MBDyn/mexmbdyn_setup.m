function mexmbdyn_setup (varargin)
% compiles mexMBCNodal functions for interaction with MBDyn multibody
% solver
%
% Syntax
%
% mexmbdyn_setup ()
% mexmbdyn_setup ('Parameter', Value)
% 
% Inputs
%
% Inputs may be supplied as parameter-value pairs, the available options
% are:
%
% 'Verbose' - logical (true/false) flag determining whether verbose output
%   from the compilation is provided. Default is false if not supplied.
%
% 'Debug' - logical (true/false) flag determining whether to perform a
%   debug build (with debugging symbols and -DDEBUG). Default is false if
%   not supplied.
%
% 'MBCLibDir' - directory in which to look for the MBDyn mbc library header
%   file. The default value of this arguement is platform-dependant.
%
% 'MBCIncludeDir' -  directory in which to look for the MBDyn mbc library 
%   file. The default value of this arguement is platform-dependant.
%
% 'PreventMBDynCheck' - Prevent checking for the existance of MBDyn on the
%   system and attempt to build the interface regardless. Default is false
%   if not supplied.
%
% 'ForceMexMBCNodalSharedMem' - Force building of the shared memory
%   communication method version of the MBDyn interface. By default this is
%   not built on Windows (but *is* built on other systems). Requires the
%   BOOST C++ library to be available. Default is false if not supplied.
%
% Example
%
%  mexmbdyn_setup ('Verbose', true)
%
%
%

    mbdyn_download_page_url = 'https://www.mbdyn.org/?Software_Download';

    % make sure all files are on matlab path
    addpath (genpath (fileparts (which ('mexmbdyn_setup'))));
    
    options.Verbose = false;
    options.Debug = false;
    options.ThrowBuildErrors = false;
    options.PreventMBDynCheck = false;
    options.ForceMexMBCNodalSharedMem = false;
    options.MBCNodalExtraMexArgs = {};
    options.MBCNodalSharedMemExtraMexArgs = {};
    options.MexExtension = mexext ();
    options.W64CrossBuild = false;

    [options.MBCLibDir, options.MBCIncludeDir, libwasfound, headerwasfound] = mbdyn.mint.find_libmbc ();
    
    options = parse_pv_pairs (options, varargin);
    
    if isoctave
        % matlab cannot display links in command output if you add html
        % tags
        mbdyn_download_page_url_with_link = mbdyn_download_page_url;
    else
        % matlab can display links in command output if you add html tags
        mbdyn_download_page_url_with_link = sprintf ('<a href="%s">%s</a>', mbdyn_download_page_url, mbdyn_download_page_url);
    end
    
    % check for the existence of MBDyn package
    if ~options.PreventMBDynCheck
        
        if ~(libwasfound || headerwasfound)
            
            fprintf ( 1, ...
                      [ ...
                        '\n\n', ...
                        '**************************    MBDyn  NOTICE    **************************\n\n', ...
                        TextWrapper.wraptext( [ ...
                        'You do not appear to have the MBDyn package installed, or at least the mbc library and ', ...
                        'header could not be found in the normal (or specified) locations. This package is needed for performing multibody dynamics ', ...
                        'simulations which is required for some functions in the foundry to work ', ...
                        '(particularly the wecSim code). You can obtain this package from here:' ] ), ...
                        '\n\n', ...
                        mbdyn_download_page_url_with_link, ...
                        '\n\n', ...
                        '**************************    MBDyn NOTICE    **************************\n\n' ...
                      ], ...
                      mbdyn_download_page_url_with_link, ...
                      mbdyn_download_page_url );
            
        elseif libwasfound && ~headerwasfound
            
            fprintf ( 1, ...
                      [ ...
                        '\n\n', ...
                        '**************************    MBDyn  NOTICE    **************************\n\n', ...
                        TextWrapper.wraptext( [ ...
                        'The mbc library from the MBDyn package was found, but not the required header files ', ...
                        'in the normal (or specified) locations. Check your MBDyn installation, it may be broken. ', ...
                        'The MBDyn package is needed for performing multibody dynamics ', ...
                        'simulations which is required for some functions in the foundry to work ', ...
                        '(particularly the wecSim code). You can obtain this package from here:' ] ), ...
                        '\n\n', ...
                        mbdyn_download_page_url_with_link, ...
                        '\n\n', ...
                        '**************************    MBDyn NOTICE    **************************\n\n' ...
                      ], ...
                      mbdyn_download_page_url_with_link, ...
                      mbdyn_download_page_url );
                
        elseif headerwasfound && ~libwasfound
            
            fprintf ( 1, ...
                      [ ...
                        '\n\n', ...
                        '**************************    MBDyn  NOTICE    **************************\n\n', ...
                        TextWrapper.wraptext( [ ...
                        'The header files for the MBDyn package were found, but not the mbc library ', ...
                        'in the normal (or specified) locations. Check your MBDyn installation, it may be broken. ', ...
                        'The MBDyn package is needed for performing multibody dynamics ', ...
                        'simulations which is required for some functions in the foundry to work ', ...
                        '(particularly the wecSim code). You can obtain this package from here:' ] ), ...
                        '\n\n', ...
                        mbdyn_download_page_url_with_link, ...
                        '\n\n', ...
                        '**************************    MBDyn NOTICE    **************************\n\n' ...
                      ], ...
                      mbdyn_download_page_url_with_link, ...
                      mbdyn_download_page_url );
            
        end
        
        if ~(libwasfound && headerwasfound)
            
            fprintf (1, 'Exiting MBDyn setup.\n');
            return;
            
        end
        
    end
    

    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mexMBCNodal and mexMBCNodalSharedMem.\n');
    
    cd(fullfile(getmfilepath (mfilename), '+mbdyn', '+mint'));

    mexMBCNodal_mexargs = {'mexMBCNodal.cpp'};
    mexMBCNodalSharedMem_mexargs = {'mexMBCNodalSharedMem.cpp'};
    
    if ~isoctave ()
        mexMBCNodal_mexargs = [ mexMBCNodal_mexargs, ...
                                { ['EXE="mexMBCNodal.', options.MexExtension, '"'], ...
                                  'CXXFLAGS="$CXXFLAGS -std=c++11"' } ...
                              ];
                            
        mexMBCNodalSharedMem_mexargs = [ mexMBCNodalSharedMem_mexargs, ...
                                         { ['EXE="mexMBCNodalSharedMem.', options.MexExtension, '"'], ...
                                           'CXXFLAGS="$CXXFLAGS -std=c++11"' } ...
                                       ];
    end
    
    if ~isempty (options.MBCIncludeDir)
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {['-I"', options.MBCIncludeDir, '"']}];
        mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {['-I"', options.MBCIncludeDir, '"']}];
    end
    
    if ~isempty (options.MBCLibDir)
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {['-L"', options.MBCLibDir, '"']}];
        mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {['-L"', options.MBCLibDir, '"']}];
    end
    
    if options.Verbose
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-v'}];
        mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {'-v'}];
    end
    
    if options.Debug
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-DDEBUG'}];
        mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {'-DDEBUG'}];
    end
    
    mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-lmbc'}];
    mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {'-lmbc'}];
    %mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {'-lmbc', 'LDFLAGS="$LDFLAGS -Wl,-rpath,"/opt/lib""'}];
    
    if ispc () || options.W64CrossBuild
        % note that this library *must* appear after -lmbc or there will be
        % linking errors on windows
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-lws2_32'}];
    end
    
    % compiling mexMBCNodal
    success = false;
    try
        mex (mexMBCNodal_mexargs{:}, options.MBCNodalExtraMexArgs{:});
        success = true; 
    catch err
        
        if ispc && ~isempty (strfind (err.message, '__imp___acrt_iob_func'))
            % FUCK! this is a pain work around missing symbols for certain
            % mingw versions
            workaround_c_prefix = 'mingw_missing___imp___acrt_iob_func_workaround';
            
            mingwroot = getenv ('MW_MINGW64_LOC');
            
            workaround_c = help ('mexmbdyn_setup>create_mingw_missing___imp___acrt_iob_func_workaround_c');
            
            output_workaround_c = fullfile (options.MBCLibDir, [workaround_c_prefix, '.c']);
            
            str2txtfile ( output_workaround_c, workaround_c );
            
            cd (fullfile (mingwroot, 'bin'));
            
            system (sprintf ('"%s" -fPIC -O2 -c "%s" -o "%s"', 'gcc.exe', output_workaround_c, [output_workaround_c(1:end-2), '.o']));

            cd (options.MBCLibDir);
            
            system (sprintf ('"%s" cr "%s" "%s"', fullfile (mingwroot, 'bin', 'ar.exe'), ['lib', workaround_c_prefix, '.a'], [workaround_c_prefix, '.o']));
            
            copyfile (fullfile (options.MBCLibDir, ['lib', workaround_c_prefix, '.a']), fullfile (options.MBCLibDir, ['lib', workaround_c_prefix, '.lib']));
            
            % add the library to the link commands
            mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {['-l', workaround_c_prefix]}];
            mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {['-l', workaround_c_prefix]}];
            
            cd(fullfile(getmfilepath (mfilename), '+mbdyn', '+mint'));
            
            try
                mex (mexMBCNodal_mexargs{:}, options.MBCNodalExtraMexArgs{:});
                success = true; 
            catch err
                if options.ThrowBuildErrors
                    rethrow (err);
                else
                    warning ('MEXMBDYN:compilefailed', ...
                             'Unable to compile mex function mexMBCNodal. Error reported was:\n%s', ...
                             err.message);
                end
            end
            
        else
            if options.ThrowBuildErrors
                rethrow (err);
            else
                warning ('MEXMBDYN:compilefailed', ...
                         'Unable to compile mex function mexMBCNodal. Error reported was:\n%s', ...
                         err.message);
            end
            
        end

    end
    
    if success
        fprintf (1, 'Successfully compiled mexMBCNodal.\n');
    else
        fprintf (1, 'mexMBCNodal was not able to be compiled.\n');
    end
    
    if options.ForceMexMBCNodalSharedMem
%         
%     if (ispc () || ~isempty (strfind (options.MBCIncludeDir, 'mingw'))) &&  ~options.ForceMexMBCNodalSharedMem
%          fprintf (1, 'Not compiling mexMBCNodalSharedMem as we are on Windows.\n');
%     else
        % compiling mexMBCNodalSharedMem
        success = false;
        try
            mex (mexMBCNodalSharedMem_mexargs{:}, options.MBCNodalSharedMemExtraMexArgs{:});
            success = true;
        catch err
            
            if options.ThrowBuildErrors
                rethrow (err);
            else
                warning ('MEXMBDYN:compilefailed', ...
                    'Unable to compile mex function mexMBCNodalSharedMem. Error reported was:\n%s', ...
                    err.message);
            end
            
        end

        if success
            fprintf (1, 'Successfully compiled mexMBCNodalSharedMem.\n');
        else
            fprintf (1, 'mexMBCNodalSharedMem was not able to be compiled.\n');
        end
    
    end
    
    fprintf (1, 'Exiting MBDyn setup.\n');
    
end


function create_mingw_missing___imp___acrt_iob_func_workaround_c ()
% /* mingw_missing___imp___acrt_iob_func_workaround.c */
% #define _CRTBLD
% #include <stdio.h>
% 
% FILE *__cdecl __acrt_iob_func(unsigned index)
% {
%     return &(__iob_func()[index]);
% }
% 
% typedef FILE *__cdecl (*_f__acrt_iob_func)(unsigned index);
% _f__acrt_iob_func __MINGW_IMP_SYMBOL(__acrt_iob_func) = __acrt_iob_func;

end