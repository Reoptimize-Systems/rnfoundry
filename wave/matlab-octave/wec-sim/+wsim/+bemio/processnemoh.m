function hydro = processnemoh (filedir, varargin)
% Reads data from a NEMOH working folder and generates wsim hydro structure 
%
% Syntax
%
% hydro = processnemoh (filedir)
% hydro = processnemoh (..., 'Parameter' value)
%
% Input
% 
%  filedir - NEMOH working folder, must include:
%   - Nemoh.cal
%   - Mesh/Hydrostatics.dat (or Hydrostatiscs_0.dat, Hydrostatics_1.dat,
%     etc. for multiple bodies)
%   - Mesh/KH.dat (or KH_0.dat, KH_1.dat, etc. for multiple bodies)
%   - Results/RadiationCoefficients.tec
%   - Results/ExcitationForce.tec
%   - Results/DiffractionForce.tec - If simu.nlHydro = 3 will be used
%   - Results/FKForce.tec - If simu.nlHydro = 3 will be used
%
% Additional arguments may be supplied using parameter-value pairs. The
% available options are:
%
% 'HydroStructure' - An empty structure, or an array of one or more
%   structures containing existing hydro data loaded from BEM simulation
%   results. The structure must have the following fields: rho, g, h, Nf,
%   w, T, body, Nh, beta, dof, cg, cb, Vo, C, A, B, ex_ma, ex_ph, ex_re,
%   ex_im, sc_ma, sc_ph, sc_re, sc_im, fk_ma, fk_ph, fk_re, fk_im. It may
%   also optionally have the field Nb, which is the number of bodies
%
%   It may also be an empty structure, in which case it has the new data
%   fields added. If the option is not supplied, this is the default value
%   of the option.
%
% Output
%
%  hydro - array of one or more structures (one for each body) containing
%    the following fields:
%
%      Nb : number of bodies in the problem
%
%      rho : density of the fluid used for the problem (kg/m^3)
% 
%      g : acceleration due to gravity for problem (m/s^2)
% 
%      h : water depth (0 for infinite)
%
%      Nf : number of frequencies the data has been generated for
% 
%      w : array of wave frequencies the data has been generated for
% 
%      T : array of wave periods the data has been generated for
% 
%      Nh : number of wave headings the data has been generated for
% 
%      body : name of each body
% 
%      dof : degrees of freedom for each body (default is 6)
%
%      cg : centre of gravity of the body
% 
%      cb : centre of buoancy of the body
% 
%      Vo : displacement volume of the body
%
%      beta : array of wave headings the data has been generated for
% 
%      A : radiation added mass coefficients
% 
%      B : radiation added damping coefficients
%
%      C : linear restoring stiffness
% 
%      ex_ma : magnitude of excitation forces
% 
%      ex_ph : phase of excitation forces
% 
%      ex_re : real part of excitation forces
% 
%      ex_im : imaginary part of excitation forces
% 
%      sc_ma : magnitude of diffraction forces
% 
%      sc_ph : phase of diffraction forces
% 
%      sc_re : real part of diffraction forces
% 
%      sc_im : imaginary part of diffraction forces
% 
%      fk_ma : magnitude of Froude-Krylov forces
% 
%      fk_ph : phase of Froude-Krylov forces
% 
%      fk_re : real part of Froude-Krylov forces
% 
%      fk_im : imaginary part of Froude-Krylov forces
%

    options.HydroStructure = struct ();
    options.ForceCoG = [];
    options.ForceCoB = [];
%     options.DoRadiationIRF = true;
%     options.IRFDuration = [];
%     options.IRFNSteps = [];
%     options.IRFNOmega = [];
%     options.IRFOmegaMin = [];
%     options.IRFOmegaMax = [];
    
    options = parse_pv_pairs (options, varargin);
    
    hydro = options.HydroStructure;
    
    % Check if we need to expand or replace the hydro structure array
    nhydro = numel (hydro);  
    
    if nhydro == 1

        if isfield (hydro(nhydro), 'Nb')
            if hydro(nhydro).Nb == 1
                % there's only one body, we replace the data in the
                % structure
                hydroind = 1;
            else
                % there are multiple bodies, so we're not replacing
                % existing data in a single structure, but are appending a
                % new structure to an array of hydro structures
                hydroind = 2;
            end
        else
            % there's only one body, we replace the data in the structure
            hydroind = 1;
        end

    elseif nhydro > 1  
        % there are multiple bodies already loaded, so we're not replacing
        % existing data, but are appending a new structure to the array
        hydroind = nhydro + 1;

    end

%     p = waitbar(0,'Reading NEMOH output file...');  % Progress bar

    hydro(hydroind).code = 'NEMOH';

    [hydro(hydroind).filedir, hydro(hydroind).file, ~] = fileparts (filedir); % Base name
    
    if exist ( fullfile (filedir, 'mesh'), 'dir') == 7
        
        meshdir = fullfile (filedir, 'mesh');
        
    elseif exist ( fullfile (filedir, 'Mesh'), 'dir') == 7
        
        meshdir = fullfile (filedir, 'Mesh');
        
    else
        error ('mesh (or Mesh) directory not found in working folder');
    end
    
    if exist ( fullfile (filedir, 'results'), 'dir') == 7
        
        resultsdir = fullfile (filedir, 'results');
        
    elseif exist ( fullfile (filedir, 'Results'), 'dir') == 7
        
        resultsdir = fullfile (filedir, 'Results');
        
    else
        error ('results (or Results) directory not found in working folder');
    end

    %% nemoh.cal file
    
    if exist (fullfile (filedir, 'nemoh.cal'), 'file')
        nemohfile = fullfile (filedir, 'nemoh.cal');
    elseif exist (fullfile (filedir, 'Nemoh.cal'), 'file')
        nemohfile = fullfile (filedir, 'Nemoh.cal');
    else
        error ('Nemoh.cal (or nemoh.cal) file was not found in the working folder')
    end

    fileID = fopen (nemohfile);
    
    raw = textscan (fileID, '%[^\n\r]');  %Read nemoh.cal
    
    raw = raw{:};
    
    fclose(fileID);

    bodyind = 0;
    for n = 1:numel (raw)

        if isempty ( strfind (raw{n}, 'Fluid specific volume')) == false
            tmp = textscan(raw{n},'%f');
            hydro(hydroind).rho = tmp{1};  % Density
        end

        if isempty ( strfind (raw{n}, 'Gravity')) == false
            tmp = textscan(raw{n},'%f');
            hydro(hydroind).g = tmp{1};  % Gravity
        end

        if isempty ( strfind (raw{n}, 'Water depth')) == false
            tmp = textscan(raw{n},'%f');
            if tmp{1} == 0 
                hydro(hydroind).h = Inf;
            else
                hydro(hydroind).h = tmp{1};  % Water depth
            end
        end

        if isempty ( strfind (raw{n}, 'Number of bodies')) == false
            tmp = textscan(raw{n},'%f');
            hydro(hydroind).Nb = tmp{1};  % Number of bodies
        end

        if isempty ( strfind (raw{n}, 'Name of mesh file')) == false
            bodyind = bodyind + 1; % increment the body index
            tmp = strsplit (raw{n},{'.', filesep()});
            hydro(hydroind).body{bodyind} = tmp{length(tmp) - 1};  % Body names
        end

        if isempty(strfind(raw{n}, 'Number of wave frequencies')) == false
            tmp = textscan(raw{n},'%f %f %f');
            hydro(hydroind).Nf = tmp{1};  % Number of wave frequencies
            hydro(hydroind).w = linspace (tmp{2}, tmp{3}, tmp{1});  % Wave frequencies
            hydro(hydroind).T = 2*pi./hydro(hydroind).w;  % Wave periods
        end

        if isempty(strfind(raw{n}, 'Number of wave directions')) == false
            tmp = textscan (raw{n},'%f %f %f');
            hydro(hydroind).Nh = tmp{1};  % Number of wave headings
            hydro(hydroind).beta = linspace (tmp{2}, tmp{3}, tmp{1});  % Wave headings
        end

    end
    % waitbar(1/7);

    %% Hydrostatics file(s)
    for bodyind = 1:hydro(hydroind).Nb

        hydro(hydroind).dof(bodyind) = 6;  % Default degrees of freedom for each body is 6

        if hydro(hydroind).Nb == 1
            fileID = fopen(fullfile (meshdir, 'Hydrostatics.dat'));
        else
            fileID = fopen(fullfile (meshdir, sprintf ('Hydrostatics_%d.dat', bodyind-1)));
        end

        raw = textscan(fileID,'%[^\n\r]');  % Read Hydrostatics.dat

        raw = raw{:};

        fclose(fileID);

        for i = 1:3
            tmp = textscan (raw{i},'%s %s %f %s %s %s %f');
            if ~isempty(options.ForceCoG{bodyind})
                hydro(hydroind).cg(i,bodyind) = options.ForceCoG{bodyind}(i);
            else
                hydro(hydroind).cg(i,bodyind) = tmp{7};  % Center of gravity
            end
            
            if ~isempty(options.ForceCoB{bodyind})
                hydro(hydroind).cb(i,bodyind) = options.ForceCoB{bodyind}(i);
            else
                hydro(hydroind).cb(i,bodyind) = tmp{3};  % Center of buoyancy
            end
        end

        tmp = textscan(raw{4},'%s %s %f');

        hydro(hydroind).Vo(bodyind) = tmp{3};  % Displacement volume

    end
    % waitbar(2/7);

    %% KH file(s)
    for bodyind = 1:hydro(hydroind).Nb

        if hydro(hydroind).Nb == 1
            fileID = fopen (fullfile (meshdir, 'KH.dat'));
        else
            fileID = fopen ( fullfile (meshdir, sprintf ('KH_%d.dat', bodyind-1)) );
        end

        raw = textscan(fileID,'%[^\n\r]');

        raw = raw{:};

        fclose(fileID);

        for i = 1:6
            tmp = textscan(raw{i},'%f');
            hydro(hydroind).C(i,:,bodyind) = tmp{1,1}(1:6);  % Linear restoring stiffness
        end

    end
    % waitbar(3/7);

    %% Radiation Coefficient file
    fileID = fopen (fullfile (resultsdir, 'RadiationCoefficients.tec'));

    raw = textscan(fileID,'%[^\n\r]');
    raw = raw{:};
    fclose(fileID);

    i = 0;

    for n = 1:numel (raw)

        if isempty(strfind(raw{n},'Motion of body')) == false

            i = i+1;

            for k = 1:hydro(hydroind).Nf
                tmp = textscan(raw{n+k},'%f');
                hydro(hydroind).A(i,:,k) = tmp{1,1}(2:2:end);  % Added mass
                hydro(hydroind).B(i,:,k) = tmp{1,1}(3:2:end);  % Radiation damping
            end

        end

    end
    % waitbar(4/7);

    %% Excitation Force file

    fileID = fopen (fullfile (resultsdir, 'ExcitationForce.tec'));

    raw = textscan (fileID, '%[^\n\r]');
    raw = raw{:};
    fclose (fileID);

    i = 0;
    for n = 1:numel (raw)

        if isempty (strfind (raw{n}, 'Diffraction force')) == false

            i = i+1;

            for k = 1:hydro(hydroind).Nf
                tmp = textscan (raw{n+k}, '%f');
                hydro(hydroind).ex_ma(:,i,k) = tmp{1,1}(2:2:end);  % Magnitude of exciting force
                hydro(hydroind).ex_ph(:,i,k) = -tmp{1,1}(3:2:end);  % Phase of exciting force (-ph, since NEMOH's x-dir is flipped)
            end

        end

    end

    hydro(hydroind).ex_re = hydro(hydroind).ex_ma .* cos(hydro(hydroind).ex_ph);  % Real part of exciting force
    hydro(hydroind).ex_im = hydro(hydroind).ex_ma .* sin(hydro(hydroind).ex_ph);  % Imaginary part of exciting force

    % waitbar(5/7);

    %% Diffraction Force file (scattering)

    DiffractionForceFile = fullfile (resultsdir, 'DiffractionForce.tec');

    if exist (DiffractionForceFile, 'file') == 2

        fileID = fopen (DiffractionForceFile);

        raw = textscan (fileID, '%[^\n\r]');
        raw = raw{:};
        fclose (fileID);

        i = 0;

        for n = 1:numel (raw)

            if isempty (strfind (raw{n}, 'Diffraction force')) == false

                i = i+1;

                for k = 1:hydro(hydroind).Nf
                    tmp = textscan (raw{n+k}, '%f');
                    hydro(hydroind).sc_ma(:,i,k) = tmp{1,1}(2:2:end);  % Magnitude of diffraction force
                    hydro(hydroind).sc_ph(:,i,k) = -tmp{1,1}(3:2:end);  % Phase of diffraction force 
                end

            end

        end

        hydro(hydroind).sc_re = hydro(hydroind).sc_ma .* cos(hydro(hydroind).sc_ph);  % Real part of diffraction force
        hydro(hydroind).sc_im = hydro(hydroind).sc_ma .* sin(hydro(hydroind).sc_ph);  % Imaginary part of diffraction force

    end
    % waitbar(6/7);

    %% Froude-Krylov force file

    FKForceFile = fullfile (resultsdir, 'FKForce.tec');

    if exist (FKForceFile, 'file') == 2

        fileID = fopen (FKForceFile);

        raw = textscan (fileID, '%[^\n\r]');
        raw = raw{:};
        fclose(fileID);

        i = 0;
        for n = 1:numel (raw)

            if isempty (strfind (raw{n}, 'FKforce')) == false

                i = i+1;

                for k = 1:hydro(hydroind).Nf
                    tmp = textscan(raw{n+k},'%f');
                    hydro(hydroind).fk_ma(:,i,k) = tmp{1,1}(2:2:end);  % Magnitude of Froude-Krylov force
                    hydro(hydroind).fk_ph(:,i,k) = -tmp{1,1}(3:2:end);  % Phase of Froude-Krylov force 
                end

            end

        end

        hydro(hydroind).fk_re = hydro(hydroind).fk_ma .* cos(hydro(hydroind).fk_ph);  % Real part of Froude-Krylov force
        hydro(hydroind).fk_im = hydro(hydroind).fk_ma .* sin(hydro(hydroind).fk_ph);  % Imaginary part of Froude-Krylov force

    end
    % waitbar(7/7);

    hydro = Normalize (hydro);  % Normalize the data according the WAMIT convention

%     close(p); % close the waitbar

%     if options.DoRadiationIRF
%         hydro = Radiation_IRF ( hydro, ...
%                                 options.RadIRFDuration, ...;
%                                 options.RadIRFNSteps, ...
%                                 options.RadIRFNOmega, ...
%                                 options.RadIRFOmegaMin, ...
%                                 options.RadIRFOmegaMax );
%     end
%     
%     if options.DoRadiationIRF
%     hydro = Radiation_IRF_SS (hydro, [], []);
% 
%     hydro = Excitation_IRF (hydro,157, [], [], [], 1.9);

end
