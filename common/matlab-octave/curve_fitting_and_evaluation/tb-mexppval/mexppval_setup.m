function mexppval_setup (varargin)
% compiles ppuva and ppmval functions for fast piecewise polynomial
% evaluation
%
% Syntax
%
% mexppval_setup ()
% mexppval_setup ('Parameter', Value)
% 
% Inputs
%
% Inputs may be supplied as parameter-value paris, currently only one
% option is supposted: 'Verbose' which is a flag determining whether
% verbose output from the compilation is provided. e.g.
%
%  mexppval_setup ('Verbose', true)
%
% prints the additional information.
%
%

    options.Verbose = false;
    
    options = parse_pv_pairs (options, varargin);

    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mex ppval and ppuval.\n');
    
    cd(getmfilepath (mfilename));

    ppmval_mexargs = {'ppmval.cpp', 'interpUtil.cpp'};
    
    ppuval_mexargs = {'ppuval.cpp', 'interpUtil.cpp'};
    
    if options.Verbose
        ppmval_mexargs = [ppmval_mexargs, {'-v'}];
        ppuval_mexargs = [ppuval_mexargs, {'-v'}];
    end
    
    % compiling ppmval
    try
        mex (ppmval_mexargs{:});
        % compiling ppuval
        mex (ppuval_mexargs{:});
    catch
        warning ('Unable to compile mex functions ppmval and ppuval. Do you have a compiler setup?');
    end
    
    fprintf (1, 'Exiting mex ppval and ppuval setup.\n');
    
end