function all_hydro_data = write_hydrobody_mat_files (hydro, outdir, varargin)
% Converts a hydro structure to .mat files for wsim.hydroBody object(s)
%
% Syntax
%
% wsim.bemio.write_hydrobody_mat_files (hydro, outdir)
% wsim.bemio.write_hydrobody_mat_files (..., 'Parameter', value);
%
% Description
%
% wsim.bemio.write_hydrobody_mat_files takes in a hydro structure as would
% be produced by the BEMIO functions such as Read_NEMOH, Read_WAMIT,
% Combine_BEM etc. and generates a set of .mat files, one for each body in
% the system, which can be loaded directly by wsim.hydroBody objects rather
% than using the intermediate H5 format files. This is particularly useful
% when using Octave as it does not support the h5read function.
%
% Each mat file will be given a name taken from the contents of hydro.body,
% i.e. the mat file will be named <name>.mat where <name> is taken from
% hydro.body{bodyind}.
%
% Input
%
%  hydro - scalar structure containing the BEMIO hydrodynamic data for all
%   bodies in a hydrodynamic system
%
%  outdir - directory in which to create the .mat files containing the
%   hydrodynamic data for each body
%
% Additional arguments may be supplied using parameter-value pairs. The
% available options are:
%
%  'AppendToFileNames' - optional character vector to append to every body
%    file name. Can be useful to denote different versions of the data for
%    the same bodies.
%
% Output
%
%  all_hydro_data - array of structures, one for each body, containing the same
%   information as was written to the .mat files.
%
% See also:  Write_H5.m
%
%

    options.AppendToFileNames = '';
    
    options = parse_pv_pairs (options, varargin);
    
    % some input checking
    assert (isstruct (hydro) && isscalar (hydro), ...
        'hydro must be a scalar structure representing a hydrodynamic system' );
    assert (ischar (outdir), 'outdir must be a character vector' );
    assert (exist (outdir, 'dir') == 7, ...
        'The output directory outdir does not appear to exist' );
    
    assert (ischar (options.AppendToFileNames), 'AppendToFileNames must be a character vector');
    
    n = 0;
    all_hydro_data = struct ();

    for bodyind = 1:hydro.Nb

        m = hydro.dof(bodyind);

        name = hydro.body{bodyind};

%         all_hydro_data(bodyind) = struct ();

        all_hydro_data(bodyind).simulation_parameters.scaled = 0;
        all_hydro_data(bodyind).simulation_parameters.wave_dir = hydro.beta;
        all_hydro_data(bodyind).simulation_parameters.water_depth = hydro.h;
        all_hydro_data(bodyind).simulation_parameters.w = hydro.w;
        all_hydro_data(bodyind).simulation_parameters.T = hydro.T;
        all_hydro_data(bodyind).properties.name = name;
        all_hydro_data(bodyind).properties.body_number = bodyind - 1;
        all_hydro_data(bodyind).properties.cg = hydro.cg(:,bodyind)';
        all_hydro_data(bodyind).properties.cb = hydro.cb(:,bodyind)';
        all_hydro_data(bodyind).properties.disp_vol = hydro.Vo(bodyind);
        all_hydro_data(bodyind).hydro_coeffs.linear_restoring_stiffness = hydro.C(:,:,bodyind);

        all_hydro_data(bodyind).hydro_coeffs.excitation.re = permute(hydro.ex_re((n+1):(n+m),:,:),[3 2 1]);
        all_hydro_data(bodyind).hydro_coeffs.excitation.im = permute(hydro.ex_im((n+1):(n+m),:,:),[3 2 1]);

        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.f = permute(hydro.ex_K((n+1):(n+m),:,:),[3 2 1]);
        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.t = hydro.ex_t;
        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.w = hydro.ex_w;
        all_hydro_data(bodyind).hydro_coeffs.added_mass.all = permute(hydro.A((n+1):(n+m),:,:),[3 2 1]);
        all_hydro_data(bodyind).hydro_coeffs.added_mass.inf_freq = permute(hydro.Ainf((n+1):(n+m),:),[2 1]);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.all = permute(hydro.B((n+1):(n+m),:,:),[3 2 1]);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.impulse_response_fun.K = permute(hydro.ra_K((n+1):(n+m),:,:),[3 2 1]);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.impulse_response_fun.t = hydro.ra_t;

        if isfield(hydro,'ss_A')

            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.it = permute(hydro.ss_O((n+1):(n+m),:),[2 1]);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.A.all = permute(hydro.ss_A((n+1):(n+m),:,:,:),[4 3 2 1]);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.B.all = permute(hydro.ss_B((n+1):(n+m),:,:,:),[4 3 2 1]);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.C.all = permute(hydro.ss_C((n+1):(n+m),:,:,:),[4 3 2 1]);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.D.all = permute(hydro.ss_D((n+1):(n+m),:),[2 1]);

        end

        % write the .mat file for this body
        outfile = fullfile (outdir, [name, options.AppendToFileNames, '.mat']);
        hydroData = all_hydro_data(bodyind);
        save (outfile, '-struct', 'hydroData');

        n = n + m;

    end

end
