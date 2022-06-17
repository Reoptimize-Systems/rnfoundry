function default_dir = defaultmulticoredir ()
% returns the default 


    default_dir = fullfile (fileparts(which('mcore.startslave')), '..', 'mcore_comms_directory');

    if exist (default_dir, "dir") ~=  7

        mkdir (default_dir);

    end


end