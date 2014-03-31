function buoydat = buoydatafromdir(libdir, buoydir, buoydat)
% buoydatafromdir: loads the appropriate data from the buoy directory
% assuming certain naming conventions for the files
%
% Syntax
%
% buoydat = buoydatafromdir(libdir, buoydir)
% 
% buoydat = buoydatafromdir(libdir, buoydir, buoydat)
%
% Input
%
%   libdir - The location of the buoy library directory, this directory
%            should contain the directories which have the buoy files 
%            inside, with each buoy having it's own directory, see 
%            'buoydir'.
% 
%   buoydir - The directory of the buoy files, this is the path relative to
%             the buoy library directory
% 
%   buoydat - (optional) A structure which is to be filled out with
%                fields containing the appropriate buoy file locations and
%                buoy data. The following fields will be added to the
%                structure, or overwritten if already present.
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
%
% Output
%
%   buoydat - the completed buoydat structure as described above
%
%

    if isempty(libdir)
        libdir = getbuoylibdir();
    end
    
    fullbuoydir = fullfile(libdir, buoydir);
    
    % Check if the directory exists then proceed to extract the required
    % information
    if exist(fullbuoydir, 'dir')

        % First get the HeaveFile
        tempfile = dir(fullfile(fullbuoydir, '*heave_radiation*'));

        % All directories are specified relative to the buoy library
        % directory so it is easier to manually construct the fiel
        % locations using appropriate separators for the platform
        if ispc
            dirslash = '\';
        else
            dirslash = '/';
        end

        % Assign the HeaveFile variable in buoydat
        buoydat.HeaveFile = [buoydir, dirslash, tempfile.name];

        % Next get the SurgeFile
        tempfile = dir(fullfile(fullbuoydir, '*surge_radiation*'));

        buoydat.SurgeFile = [buoydir, dirslash, tempfile.name];

        % Next get the ExcitationFile File, the .1 file, previously the
        % HydroCoeffsFile
        tempfile = dir(fullfile(fullbuoydir, '*excitation_force*'));

        buoydat.ExcitationFile = [buoydir, dirslash, tempfile.name];

        % Next get the HydroCoeffsFile file, the .2 file, previously the
        % BuoyVarsFile
        tempfile = dir(fullfile(fullbuoydir, '*hydro_coefficients*'));

        buoydat.HydroCoeffsFile = [buoydir, dirslash, tempfile.name];

        % Finally get the buoyparams file
        tempfile = dir(fullfile(fullbuoydir, '*buoyparams*'));

        buoydat.BuoyParameters = load(fullfile(fullbuoydir, tempfile.name));
        
        buoydat.buoy = buoydir;
    
    else
        error('BUOYSIM:dirnotfound', 'Buoy directory:\n%s\n\ndoes not exist.', fullbuoydir)
    end

end