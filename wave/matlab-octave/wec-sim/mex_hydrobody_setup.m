function mex_hydrobody_setup (varargin)
% compiles mex_hydrobody function from EWST
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
% 'MexExtension' - 
%
% 'W64CrossBuild' -  
%
% 'ThrowBuildErrors' - 
%
% Example
%
%  mexmbdyn_setup ('Verbose', true)
%
%
%

    % make sure all files are on matlab path
    addpath (genpath (fileparts (which ('mex_hydrobody_setup'))));
    
    options.Verbose = false;
    options.Debug = false;
    options.DefineDEBUG = [];
    options.ThrowBuildErrors = false;
    options.MexExtension = mexext ();
    options.W64CrossBuild = false;
    options.MBCNodalExtraMexArgs = {};

    options = parse_pv_pairs (options, varargin);
    
    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mex_hydrobody.\n');
    
    cd(fullfile(getmfilepath (mfilename), '+wsim'));

    mex_hydrobody_mexargs = {'mex_hydrobody.cpp', 'hydrobody.cpp', '-DMATLAB_MEX'};
    
    if ~isoctave ()
        mex_hydrobody_mexargs = [ mex_hydrobody_mexargs, ...
                                { ['EXE="mex_hydrobody.', options.MexExtension, '"'], ...
                                  'CXXFLAGS="$CXXFLAGS -std=c++11 -Wfatal-errors"' } ...
                              ];
    end
    
    
    if options.Verbose
        mex_hydrobody_mexargs = [mex_hydrobody_mexargs, {'-v'}];
    end
    
    if options.Debug
        
        mex_hydrobody_mexargs = [mex_hydrobody_mexargs, {'-g'}];
        
        if isempty (options.DefineDEBUG)
            options.DefineDEBUG = true;
        else
            check.isLogicalScalar (options.DefineDEBUG, true, 'DefineDEBUG');
        end
        
        if options.DefineDEBUG
            mex_hydrobody_mexargs = [mex_hydrobody_mexargs, {'-DDEBUG'}];
        end
        
    end
    
    % compiling mex_hydrobody
    success = false;
    try
        mex (mex_hydrobody_mexargs{:}, options.MBCNodalExtraMexArgs{:});
        success = true; 
    catch err
        if options.ThrowBuildErrors
            rethrow (err);
        else
            warning ('MEXHYDRO:compilefailed', ...
                'Unable to compile mex function mex_hydrobody. Error reported was:\n%s', ...
                err.message);
        end
    end
    
    if success
        fprintf (1, 'Successfully compiled mex_hydrobody.\n');
    else
        fprintf (1, 'mex_hydrobody was not able to be compiled.\n');
    end
    
    fprintf (1, 'Exiting mex_hydrobody setup.\n');
    
end