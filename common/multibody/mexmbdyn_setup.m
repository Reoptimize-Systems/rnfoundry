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
% Inputs may be supplied as parameter-value paris, currently only one
% option is supposted: 'Verbose' which is a flag determining whether
% verbose output from the compilation is provided. e.g.
%
%  mexmbdyn_setup ('Verbose', true)
%
% prints the additional information.
%
%

    options.Verbose = false;
    options.Debug = false;
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
    
    fprintf (1, 'Setting up mexMBCNodal.\n');
    
    cd(fullfile(getmfilepath (mfilename), '+mbdyn', '+mint'));

    mexMBCNodal_mexargs = {'mexMBCNodal.cpp', '-lmbc', 'CXXFLAGS="$CXXFLAGS -std=c++11"'};
    
    if ispc
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-lws2_32.lib'}];
    end
    
    if ~isempty (options.MBCIncludeDir)
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {['-I"', options.MBCIncludeDir, '"']}];
    end
    
    if ~isempty (options.MBCLibDir)
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {['-L"', options.MBCLibDir, '"']}];
    end
    
    if options.Verbose
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-v'}];
    end
    
    if options.Debug
        mexMBCNodal_mexargs = [mexMBCNodal_mexargs, {'-DDEBUG'}];
    end
    
    % compiling mexMBCNodal
    try
        mex (mexMBCNodal_mexargs{:});
    catch err
        warning ('MEXMBDYN:compilefailed', ...
            'Unable to compile mex functions mexMBCNodal. Error reported was:\n%s', ...
            err.message);
    end
    
    fprintf (1, 'Finished setting up mmexMBCNodal.\n');
    
end