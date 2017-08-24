% publish examples

example_dir = getmfilepath (mfilename);

% make sure directories exist
pubdir = fullfile (example_dir, 'published');
mkdir (pubdir);

htmlpubdir = fullfile (pubdir, 'html');
mkdir (fullfile (htmlpubdir, 'html'));


% example_radial_flux_permanent_magnet_machine_sim
outdir = fullfile (htmlpubdir, 'example_radial_flux_permanent_magnet_machine_sim');

mkdir (outdir);

publish ( fullfile (example_dir, 'example_radial_flux_permanent_magnet_machine_sim'), ...
          'outputDir', outdir);