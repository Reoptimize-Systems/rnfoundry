function mexslmeval_setup (varargin)
% builds the mexslmeval mex function from source using mex
%
% Syntax
% 
% mexslmeval_setup ()
%
% Description
%
% Builds the mexslm mex funcion from the C++ source. Requires a working mex
% compiler to be set up in Matlab. mexslmeval uses the GNU scientific
% library and you must therefore have this on your system and on your
% compiler's path. To be specific mexslmeval links to libgsl and
% libgslcblas for its histogram functions.
%
%
% See also: slmeval, slmengine

    options.MexOpts = '';
    options.DoCrossBuildWin64 = false;
    options.Verbose = false;
    options.GSLLibDir = '';
    options.GSLIncludeDir = '';
    
    options = parse_pv_pairs (options, varargin);
    
    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mexslmeval.\n');
    
    cd(getmfilepath (mfilename));
    
    mexargs = {'mexslmeval.cpp', '-lgsl', '-lgslcblas'};
    
    if options.Verbose
        mexargs = [ mexargs, {'-v'}];
    end
    
    if ~isempty (options.GSLLibDir)
        % add the location of the gsl libs files
        mexargs = [ mexargs, ...
            { ['-L"' options.GSLLibDir, '"'] } ];
    end
    
    if ~isempty (options.GSLIncludeDir)
        % add the location of the gsl header files
        mexargs = [ mexargs, ...
            { ['-I"' options.GSLIncludeDir, '"']} ];
    end
    
    if ~isempty(options.MexOpts)
        if exist (options.MexOpts, 'file')
            mexargs = [mexargs, {['-f "', options.MexOpts, '"']}];
        else
            error ('Specified alternative mex options file does not exist.')
        end
    end

    try
        % note the order of the linking commands seams to matter here
        mex(mexargs{:})
    catch
        warning ('mexslmeval compilation failed, you may be missing required libraries, gsl and gslcblas');
    end
    
    fprintf (1, 'Finished building mexslmeval.\n');

end