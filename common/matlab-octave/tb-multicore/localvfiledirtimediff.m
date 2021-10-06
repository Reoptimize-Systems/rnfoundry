function tdiff = localvfiledirtimediff(remotedir, slaveID)
% localvfiledirtimediff: determines difference between local computer time
% and a remote directory time by examining file information
%
%
    tempfile = fullfile(remotedir, ['temp_', int2str(slaveID), '.mat']);

    comptime = now;
    
    generateemptyfile(tempfile);
    
    finfo = dir(tempfile);
    
    mbdelete2(tempfile);
    
    try
        fdate = datenum(finfo.date, 'dd-mmm-yyyy hh:mm:ss');
    catch
        fdate = comptime;
    end
    
    if isempty(fdate)
        fdate = comptime;
    end
    
    tdiff = comptime - fdate;

end