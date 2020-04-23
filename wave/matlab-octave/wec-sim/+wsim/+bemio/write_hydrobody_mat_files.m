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
    options.PrependToFileNames = '';
    
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
    m_add = 0;
    
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
        all_hydro_data(bodyind).properties.dof = hydro.dof(bodyind);
        all_hydro_data(bodyind).properties.dof_start = m_add + 1;
        all_hydro_data(bodyind).properties.dof_end = m_add + m;
        
        if isfield (hydro,'gbm') 
            % Only if generalized body modes have been used
            tmp = hydro.gbm((n+1):(n+m),:,1);
            all_hydro_data(bodyind).gbm.mass = tmp(dof_start+6:dof_end, dof_start+6:dof_end);
            
            tmp = hydro.gbm((n+1):(n+m),:,2);
            all_hydro_data(bodyind).gbm.damping = tmp(obj.dof_start+6:obj.dof_end,obj.dof_start+6:obj.dof_end);
            
            tmp = hydro.gbm((n+1):(n+m),:,3);
            all_hydro_data(bodyind).gbm.stiffness = tmp(obj.dof_start+6:obj.dof_end,obj.dof_start+6:obj.dof_end);
            
            % the below permutes were too confusing for me to unpick
            % quickly, so just replicating the first premute from Write_H5,
            % then the second permute which would be done in h5load
            tmp = permute(hydro.gbm((n+1):(n+m),:,4),[3 2 1]);
            all_hydro_data(bodyind).hydro_coeffs.linear_restoring_stiffness = permute (tmp(1,m_add + 1:m_add + m,:), [3,2,1]);
        else
            all_hydro_data(bodyind).hydro_coeffs.linear_restoring_stiffness = hydro.C(:,:,bodyind);
        end
        
        m_add = m_add + m;
        
        if isfield (hydro,'md_mc')
            % Only if mean drift variables (momentum conservation) have been calculated in BEM
             all_hydro_data(bodyind).mean_drift_momentum_conservation = hydro.md_mc((n+1):(n+m),:,:);
%             H5_Create_Write_Att(filename,['/body' num2str(i) '/hydro_coeffs/mean_drift/momentum_conservation/val/'],permute(hydro.md_mc((n+1):(n+m),:,:),[3 2 1]),'Value of mean drift force (momentum conservation)','');
        end
    
        if isfield (hydro,'md_cs')
            % Only if mean drift variables (control surface approach) have been calculated in BEM
            all_hydro_data(bodyind).mean_drift_control_surface = hydro.md_cs((n+1):(n+m),:,:);
%             H5_Create_Write_Att(filename,['/body' num2str(i) '/hydro_coeffs/mean_drift/control_surface/val/'],permute(hydro.md_cs((n+1):(n+m),:,:),[3 2 1]),'Value of mean drift force (control surface)','');
        end
    
        all_hydro_data(bodyind).hydro_coeffs.excitation.re = hydro.ex_re((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.excitation.im = hydro.ex_im((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.excitation.mag = hydro.ex_ma((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.excitation.phase = hydro.ex_ph((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.f = hydro.ex_K((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.t = hydro.ex_t;
        all_hydro_data(bodyind).hydro_coeffs.excitation.impulse_response_fun.w = hydro.ex_w;
        all_hydro_data(bodyind).hydro_coeffs.added_mass.all = hydro.A((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.added_mass.inf_freq = hydro.Ainf((n+1):(n+m),:);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.all = hydro.B((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.impulse_response_fun.K = hydro.ra_K((n+1):(n+m),:,:);
        all_hydro_data(bodyind).hydro_coeffs.radiation_damping.impulse_response_fun.t = hydro.ra_t;

        if isfield(hydro,'ss_A')

            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.it = hydro.ss_O((n+1):(n+m),:);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.A.all = hydro.ss_A((n+1):(n+m),:,:,:);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.B.all = hydro.ss_B((n+1):(n+m),:,:,:);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.C.all = hydro.ss_C((n+1):(n+m),:,:,:);
            all_hydro_data(bodyind).hydro_coeffs.radiation_damping.state_space.D.all = hydro.ss_D((n+1):(n+m),:);

        end
        
%         try tmp = wsim.bemio.h5load(filename, [body_name '/properties/mass']);
%             obj.hydroData.gbm.mass      = tmp(obj.dof_start+6:obj.dof_end,obj.dof_start+6:obj.dof_end); clear tmp; end;
%         try tmp = wsim.bemio.h5load(filename, [body_name '/properties/stiffness']);
%             obj.hydroData.gbm.stiffness = tmp(obj.dof_start+6:obj.dof_end,obj.dof_start+6:obj.dof_end); clear tmp; end;
%         try tmp = wsim.bemio.h5load(filename, [body_name '/properties/damping']);
%             obj.hydroData.gbm.damping   = tmp(obj.dof_start+6:obj.dof_end,obj.dof_start+6:obj.dof_end); clear tmp;end;
%
%         if (obj.dof_gbm>0)
%             obj.linearDamping = [obj.linearDamping(1:6) zeros(1,obj.dof_gbm)];
%         end
%         
%         if obj.meanDriftForce == 0
%             obj.hydroData.hydro_coeffs.mean_drift = 0.*obj.hydroData.hydro_coeffs.excitation.re;
%         elseif obj.meanDriftForce == 1
%             obj.hydroData.hydro_coeffs.mean_drift =  wsim.bemio.h5load(filename, [body_name '/hydro_coeffs/mean_drift/control_surface/val']);
%         elseif obj.meanDriftForce == 2
%             obj.hydroData.hydro_coeffs.mean_drift =  wsim.bemio.h5load(filename, [body_name '/hydro_coeffs/mean_drift/momentum_conservation/val']);
%         else
%             error('Wrong flag for mean drift force.')
%         end

        % write the .mat file for this body
        outfile = fullfile (outdir, [options.PrependToFileNames, name, options.AppendToFileNames, '.mat']);
        hydroData = all_hydro_data(bodyind);
        save (outfile, '-struct', 'hydroData');

        n = n + m;

    end

end
