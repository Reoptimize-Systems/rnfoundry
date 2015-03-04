function mexslmeval_setup ()
% creates the mexslmeval mex file

    CC = onCleanup (@() cd(pwd));
    
    cd(getmfilepath (mfilename));

    try
        % note the order of the linking commands seams to matter here
        mex mexslmeval.cpp -lgsl -lgslcblas
    catch
        warning ('mexslmeval compilation failed, you may be missing required libraries, gsl and gslcblas');
    end

end