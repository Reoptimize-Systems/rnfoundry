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
% Example
%
%  mexmbdyn_setup ('Verbose', true)
%
% 
%
%

    options.Verbose = false;
    options.Debug = false;
    options.ThrowErrors = false;
    
    if ispc
        switch computer ('arch')
            case 'win64'
                options.MBCLibDir = fullfile (getmfilepath ('mexmbdyn_setup'), 'mbdyn_win64', 'lib');
                options.MBCIncludeDir = fullfile (getmfilepath ('mexmbdyn_setup'), 'mbdyn_win64', 'include');
            case 'win32'
                options.MBCLibDir = fullfile (getmfilepath ('mexmbdyn_setup'), 'mbdyn_win32', 'lib');
                options.MBCIncludeDir = fullfile (getmfilepath ('mexmbdyn_setup'), 'mbdyn_win32', 'include');
            otherwise
                error ('Arch not supported.')
        end
    else
        options.MBCLibDir = '';
        options.MBCIncludeDir = '';
    end
    
    options = parse_pv_pairs (options, varargin);

    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mexMBCNodal and mexMBCNodalSharedMem.\n');
    
    cd(fullfile(getmfilepath (mfilename), '+mbdyn', '+mint'));

    mexMBCNodal_mexargs = {'mexMBCNodal.cpp'};
    mexMBCNodalSharedMem_mexargs = {'mexMBCNodalSharedMem.cpp'};
    
    if ~isoctave
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'CXXFLAGS="$CXXFLAGS -std=c++11"'}];
        mexMBCNodalSharedMem_mexargs = [mexMBCNodalSharedMem_mexargs, {'CXXFLAGS="$CXXFLAGS -std=c++11"'}];
    end
    
    if ispc
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-lws2_32'}];
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
    
    % compiling mexMBCNodal
    try
        mex (mexMBCNodal_mexargs{:});
        fprintf (1, 'Finished setting up mmexMBCNodal.\n');
    catch err
        if options.ThrowErrors
            rethrow (err);
        else
            warning ('MEXMBDYN:compilefailed', ...
                'Unable to compile mex function mexMBCNodal. Error reported was:\n%s', ...
                err.message);
        end
    end
    
    % compiling mexMBCNodalSharedMem
    try
        mex (mexMBCNodalSharedMem_mexargs{:});
        fprintf (1, 'Finished setting up mexMBCNodalSharedMem.\n');
    catch err
        if options.ThrowErrors
            rethrow (err);
        else
            warning ('MEXMBDYN:compilefailed', ...
                'Unable to compile mex function mexMBCNodalSharedMem. Error reported was:\n%s', ...
                err.message);
        end
    end
    
end