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

    options.Verbose = false;
    options.GSLLibDir = '';
    options.GSLIncludeDir = '';
    options.ExtraMexArgs = {};
    options.MexExtension = mexext ();
    options.ThrowBuildErrors = false;
    
    options = parse_pv_pairs (options, varargin);
    
    CC = onCleanup (@() cd(pwd));
    
    fprintf (1, 'Setting up mexslmeval.\n');
    
    cd(getmfilepath (mfilename));
    
    mexargs = {'mexslmeval.cpp', '-lgsl', '-lgslcblas'};
    
    if ~isoctave ()
        mexargs = [mexargs, { ['EXE="mexslmeval.', options.MexExtension, '"']}];
    end
    
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

    try
        % note the order of the linking commands seams to matter here
        mex(mexargs{:}, options.ExtraMexArgs{:})
    catch err
        
        if ispc && ~isempty (strfind (err.message, '__imp___acrt_iob_func'))
            % FUCK! this is a pain work around missing symbols for certain
            % mingw versions
            workaround_c_prefix = 'mingw_missing___imp___acrt_iob_func_workaround';
            
            mingwroot = getenv ('MW_MINGW64_LOC');
            
            workaround_c = help ('mexmbdyn_setup>create_mingw_missing___imp___acrt_iob_func_workaround_c');
            
            output_workaround_c = fullfile (options.GSLLibDir, [workaround_c_prefix, '.c']);
            
            str2txtfile ( output_workaround_c, workaround_c );
            
            cd (fullfile (mingwroot, 'bin'));
            
            system (sprintf ('"%s" -fPIC -O2 -c "%s" -o "%s"', 'gcc.exe', output_workaround_c, [output_workaround_c(1:end-2), '.o']));

            cd (options.GSLLibDir);
            
            system (sprintf ('"%s" cr "%s" "%s"', fullfile (mingwroot, 'bin', 'ar.exe'), ['lib', workaround_c_prefix, '.a'], [workaround_c_prefix, '.o']));
            
            copyfile (fullfile (options.GSLLibDir, ['lib', workaround_c_prefix, '.a']), fullfile (options.GSLLibDir, ['lib', workaround_c_prefix, '.lib']));
            
            % add the library to the link commands
            mexargs = [mexargs, {['-l', workaround_c_prefix]}];
            
            cd(fullfile(getmfilepath (mfilename)));
            
            try
                mex(mexargs{:}, options.ExtraMexArgs{:})
                success = true; 
            catch err
                if options.ThrowBuildErrors
                    rethrow (err);
                else
                    warning ('mexslmeval compilation failed with the following error message:\n%s\nYou may be missing required libraries, gsl and gslcblas', err.message);
                end
            end
            
        else
            
            if options.ThrowBuildErrors
                rethrow (err);
            else
                warning ('mexslmeval compilation failed with the following error message:\n%s\nYou may be missing required libraries, gsl and gslcblas', err.message);
            end
        
        end
    end
    
    fprintf (1, 'Finished building mexslmeval.\n');

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