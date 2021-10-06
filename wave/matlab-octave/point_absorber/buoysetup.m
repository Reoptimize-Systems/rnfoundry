function buoydat = buoysetup(buoy, buoylibdirectory, buoydat)
% buoysetup: load the appropriate buoy details from a buoy library
% directory based on either the buoy directory name, or on the number of
% the buoy in the directory
%
% Syntax
%
% buoydat = buoysetup(buoy, buoylibdir, buoydat)
%
% Input
%
% buoy - either a scalar number which is the desired buoy directory number,
%        if the buoy subdirectories in the buoy library directory were
%        ordered alphabetically, or a string containing the desired buoy
%        subdirectory name, e.g. 'cyl_3dia_1dr' for the 3m, diameter, 1m
%        draft buoy.
%
% buoylibdir - optional location of the buoy library directory, if not
%        supplied this will be assumed to be a subdirectory 'buoylib' in 
%        the same directory as this function. The function must be on the 
%        matlab path for this to work.
%
% buoydat - optional structure to which the buoy information is to be 
%        added. If not supplied a new structure is created.
%
% Output
%
% buoydat - A structure with the following fields added or overwritten:
%
%       HeaveFile - string containing the location of the file containing
%       the heave radiation coefficients relative to the buoy directory
%
%       SurgeFile - string containing the location of the file containing
%       the surge radiation coefficients relative to the buoy directory
%
%       ExcitationFile - string containing the location of the file
%       containing the excitation coefficients relative to the buoy
%       directory
%
%       HydroCoeffsFile - string containing the location of the file
%       containing the buoy hydrodynamic coefficients relative to the buoy
%       directory
%
%       BuoyParameters - a structure containing various other important
%       physical characteristics of the buoy and it's capabilities in
%       simulation, such as it's radius, draft, max and min possible
%       simulation frequencies for the coefficients in the other files etc.
%

    if nargin < 3 || isempty(buoydat)
        buoydat = struct();
    end

    if nargin < 2 || isempty(buoylibdirectory)
        buoylibdirectory = buoylibdir ();
    end

    if isnumeric(buoy) && isscalar(buoy) 
        
        % make sure the buoy number is an integer
        buoy = ceil(buoy);
        
        % get the total number of buoys in the buoy library directory
        numbuoys = buoynum2buoydata(buoylibdirectory);
        
        % if the desired buoy number is sensible, look up the buoy info
        if buoy <= numbuoys

            buoydat = buoynum2buoydata(buoylibdirectory, buoy, buoydat);

        else
            error('There are less buoys in the library than the requested buoy number.');
        end

    elseif ischar(buoy)

        % the buoy directory name has been supplied directly, so look up
        % the information from the buoy library directory
        buoydat = buoydatafromdir(buoylibdirectory, buoy, buoydat);

    else

        error('BUOYSIM:incorrectinput', 'Buoy must be a string or scalar number.')

    end
        

end