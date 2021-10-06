function buoydat = buoynum2buoydata(libdir, buoynum, buoydat)
% buoynum2buoydata: converts a number to buoy data where the number
% corresponds to one of the buoy library directories when numbered
% alphabetically

    if isempty(libdir)
        % use default directory if not supplied
        libdir = fullfile(fileparts(which('buoynum2buoydata')), 'buoylib');
    end
    
    temp = dir(fullfile(libdir, 'cyl_*'));
    
    if nargin == 1
        % return the max number of available buoys in the directory
        buoydat = numel(temp);
    else
        
        check.isScalarInteger(buoynum, true, 'buoynum');
        
        buoynum = round(buoynum);
        
        if nargin < 3 || isempty(buoydat)
            buoydat.dummy = [];
        end

        if buoynum <= length(temp)
            buoydir = temp(buoynum).name;
        else
            error('BUOYSIM:invalidbuoynum', ...
                'buoy number not valid, buoy %d does not exist in %s\nOnly %d buoy directories were found\nbuoynum2buoydata directory is %s', buoynum, libdir, length(temp), which('buoynum2buoydata'));
        end

        buoydat = buoydatafromdir(libdir, buoydir, buoydat);

        if isfield(buoydat, 'dummy')
            buoydat = rmfield(buoydat, 'dummy');
        end
        
    end
end

