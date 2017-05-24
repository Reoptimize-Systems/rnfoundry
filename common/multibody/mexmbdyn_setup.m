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
    
    options = parse_pv_pairs (options, varargin);

    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mexMBCNodal.\n');
    
    cd(fullfile(getmfilepath (mfilename), '+mbdyn', '+mint'));

    mexMBCNodal_mexargs = {'mexMBCNodal.cpp', '-lmbc'};
    
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