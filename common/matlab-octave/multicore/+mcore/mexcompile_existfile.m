function mexcompile_existfile (varargin)
% compile the existfile mex function

    options.Verbose = false;
    options.ExtraMexArgs = {};
    options.MexExtension = mexext ();
    options.ThrowBuildErrors = false;
    
    options = parse_pv_pairs (options, varargin);
    
    filepath = getmfilepath (['mcore.', mfilename]);
    
    cdir = pwd;
    CC = onCleanup (@() cd (cdir));
    
    cd (filepath)
    
    mexargs = {'existfile.c'};
    
    if ~mcore.isoctave ()
        mexargs = [mexargs, { ['EXE="existfile.', options.MexExtension, '"']}];
    end
    
    if options.Verbose
        mexargs = [mexargs, {'-v'}];
    end
    
    try
        fprintf (1, 'Building mcore.existfile mex file\n');
        mex (mexargs{:}, options.ExtraMexArgs{:});
    catch err
        if options.ThrowBuildErrors == true
            rethrow (err);
        else
            warning ('existfile mex compilation falied with message:\n%s', err.message);
        end
    end

end