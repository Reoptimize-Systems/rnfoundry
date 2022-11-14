function nfiles = countmatfilesstartingwith (startstr, dirpath)

    nfiles = numel (dir (fullfile (dirpath, [startstr, '*.mat'])));

end
