function mexslmeval_setup ()
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

    CC = onCleanup (@() cd(pwd));
    
    cd(getmfilepath (mfilename));

    try
        % note the order of the linking commands seams to matter here
        mex mexslmeval.cpp -lgsl -lgslcblas
    catch
        warning ('mexslmeval compilation failed, you may be missing required libraries, gsl and gslcblas');
    end

end