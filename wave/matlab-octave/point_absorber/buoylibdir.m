function libdir = buoylibdir()
% returns the directory containing the buoy library files directories, the
% top level buoy library directory
%
% Syntax
%
% libdir = buoylibdir()
%

    libdir = fullfile(getmfilepath (mfilename), 'buoylib');

end