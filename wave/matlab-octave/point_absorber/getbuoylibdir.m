function libdir = getbuoylibdir()
% getbuoylibdir: returns the directory that this function resides in. This
% should be placed in the same directory as the directory containing the
% buoy library files directories

    libdir = fullfile(fileparts(which('getbuoylibdir')), 'buoylib');

end