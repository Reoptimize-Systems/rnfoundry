function mexslmeval_setup ()

    CC = onCleanup (@() cd(pwd));
    
    cd(getmfilepath (mfilename));

    mex mexslmeval.cpp -lgsl -lgslcblas

end