function mexppval_setup()

    CC = onCleanup (@() cd(pwd));
    
    cd(getmfilepath (mfilename));

    % compiling ppmval
    try
        mex ppmval.cpp interpUtil.cpp;
        % compiling ppuval
        mex ppuval.cpp interpUtil.cpp;
    catch
        warning ('Unable to compile mex functions ppmval and ppuval. Do you have a compiler setup?');
    end
end