function mexppval_setup()

    CC = onCleanup (@() cd(pwd));
    
    cd(getmfilepath (mfilename));

    % compiling ppmval
    mex ppmval.cpp interpUtil.cpp
    % compiling ppuval
    mex ppuval.cpp interpUtil.cpp
    
end