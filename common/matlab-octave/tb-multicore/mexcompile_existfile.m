function mexcompile_existfile (varargin)
% compile the existfile mex function

    options.Verbose = false;
    
    options = parse_pv_pairs (options, varargin);
    

    filepath = getmfilepath (mfilename);
    
    cdir = pwd;
    CC = onCleanup (@() cd (cdir));
    
    cd (filepath)
    
    mexargs = {'existfile.c'};
    
    if options.Verbose
        mexargs = [mexargs, {'-v'}];
    end
    
    try
        mex (mexargs{:});
    catch
        warning ('existfile mex compilation falied');
    end

end