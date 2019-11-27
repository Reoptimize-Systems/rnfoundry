function hydro = processnemoh (filedir, varargin)
% Reads data from a NEMOH working folder and generates wsim hydro structure 
%
% Syntax
%
% hydro = wsim.bemio.processnemoh (filedir)
% hydro = wsim.bemio.processnemoh (..., 'Parameter', value)
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
% 'ForceCoG' - option to replace all or some of the bodies' Centres of
%   Gravity (CoG). Should contain a cell array with the same number of
%   elements as bodies in the Nemoh solution being loaded. For each body,
%   the corresponding cell of the array will be checked for a replacement 3
%   element centre of gravity vector. Each cell should either contain an
%   empty matrix, or a replacement 3 element centre of gravity vector for
%   the corresponding body. If is empty the CoG is not replaced, if it is a
%   vector, the CoG is replaced with this vector. e.g.
%
%   hydro = wsim.bemio.processnemoh ( 'my_nemoh_ouput_dir', ...
%                                     'ForceCoG', { [0, 0, 0.5], [], [0, 0,1.4] });
%
%   loads Nemoh output in a folder 'my_nemoh_ouput_dir' with three bodies,
%   and replaces the CoG of the first and third body, but leaves the second
%   body as it is in the Nemoh files.
%
%
% 'ForceCoB' - option to replace all or some of the bodies' Centres of
%   Buoyancy (CoB). Should contain a cell array with the same number of
%   elements as bodies in the Nemoh solution being loaded. For each body,
%   the corresponding cell of the array will be checked for a replacement 3
%   element centre of buoyancy vector. Each cell should either contain an
%   empty matrix, or a replacement 3 element centre of buoyancy vector for
%   the corresponding body. If is empty the CoG is not replaced, if it is a
%   vector, the CoB is replaced with this vector. e.g.
%
%   hydro = wsim.bemio.processnemoh ( 'my_nemoh_ouput_dir', ...
%                                     'ForceCoB', { [0, 0, 0.5], [], [0, 0,1.4] });
%
%   loads Nemoh output in a folder 'my_nemoh_ouput_dir' with three bodies,
%   and replaces the CoB of the first and third body, but leaves the second
%   body as it is in the Nemoh files.
%
%
% 'ForceKH' - option to replace all or some of the bodies' linear restoring 
%   matrices. Should contain a cell array with the same number of elements
%   as bodies in the Nemoh solution being loaded. For each body, the
%   corresponding cell of the array will be checked for a replacement 6
%   element vector representing the diagonal of the linear restoring
%   matrix. Each cell should either contain an empty matrix, or a
%   replacement 6 element vector for the corresponding body. If is empty
%   the linear restoring matrix is not replaced, if it is a vector, the
%   linear restoring matrix is replaced with this vector. 
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
    options.ForceCoG = {};
    options.ForceCoB = {};
    options.ForceVolume = {};
    options.ForceKH = {};
    options.ProcessKochin = true;
%     options.DoRadiationIRF = true;
%     options.IRFDuration = [];
%     options.IRFNSteps = [];
%     options.IRFNOmega = [];
%     options.IRFOmegaMin = [];
%     options.IRFOmegaMax = [];
    
    options = parse_pv_pairs (options, varargin);
    
    assert (iscell (options.ForceCoG), 'ForceCoG must be a cell array');
    assert (iscell (options.ForceCoB), 'ForceCoB must be a cell array');
    assert (iscell (options.ForceVolume), 'ForceVolume must be a cell array');
    assert (iscell (options.ForceKH), 'ForceKH must be a cell array');
    
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
    
    % if no replacement CoG or CoB was supplied, ensure the default is the
    % right size.
    if isempty (options.ForceCoG)
        options.ForceCoG = repmat ({[]}, 1, hydro(hydroind).Nb);
    end
    if isempty (options.ForceCoB)
        options.ForceCoB = repmat ({[]}, 1, hydro(hydroind).Nb);
    end
    if isempty (options.ForceVolume)
        options.ForceVolume = repmat ({[]}, 1, hydro(hydroind).Nb);
    end
    if isempty (options.ForceKH)
        options.ForceKH = repmat ({[]}, 1, hydro(hydroind).Nb);
    end

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

        if ~isempty(options.ForceVolume{bodyind})
            hydro(hydroind).Vo(bodyind) = options.ForceVolume{bodyind};
        else
            hydro(hydroind).Vo(bodyind) = tmp{3};  % Displacement volume
        end

    end
    % waitbar(2/7);

    %% KH file(s)
    for bodyind = 1:hydro(hydroind).Nb

        if isempty(options.ForceKH{bodyind})
            
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
            
        else
            hydro(hydroind).C(1:6,:,bodyind) = options.ForceKH{bodyind};
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
    
    %================= READING KOCHIN FILES ===================%
    
    %clear Kochin_BVP x theta i H
    if exist(fullfile(resultsdir, 'Kochin.    1.dat'),'file') == 2 ...
            && options.ProcessKochin == true
        
        nb_DOF=size(hydro(hydroind).ex_ma,1);
        nBodies=hydro(hydroind).Nb;
        nw=hydro(hydroind).Nf;
        for j=1:nw
            for i=1:(nb_DOF*nBodies+1)
                clear Kochin
                x=(nb_DOF*nBodies+1)*(j-1)+i;
                switch  numel(num2str(x))
                    case 1
                        filename=['Kochin.    ',num2str(x),'.dat'];
                    case 2
                        filename=['Kochin.   ',num2str(x),'.dat'];
                    case 3
                        filename=['Kochin.  ',num2str(x),'.dat'];
                    case 4
                        filename=['Kochin. ',num2str(x),'.dat'];
                    case 5
                        filename=['Kochin.',num2str(x),'.dat'];
                end
                
                kochin_file_path = fullfile (resultsdir, filename);

                if exist(kochin_file_path, 'file') ~= 2
                    error ('Looked for Kochin file %s, but it does not appear to exist.', kochin_file_path);
                end
                
                fileID = fopen (kochin_file_path, 'r');
                
                Kochin = fscanf(fileID,'%f');

                for ntheta=1:size(Kochin,1)/3
                    theta(ntheta) = Kochin(3*(ntheta-1)+1);
                    Kochin_BVP(ntheta,1,x) = Kochin(3*(ntheta-1)+2);
                    Kochin_BVP(ntheta,2,x) = Kochin(3*(ntheta-1)+3);
                end
                
                fclose(fileID);
                
            end
        end

    %------Calculate RAO-------
        w=hydro(hydroind).w;
        i=sqrt(-1);
        Fe = hydro(hydroind).ex_re - hydro(hydroind).ex_im*i; % Added by Toan
        A= hydro(hydroind).A;
        B = hydro(hydroind).B ;
        f=w/(2*pi);
        M=zeros(6,6);
        for nn =1:3
            M(nn,nn) = hydro(hydroind).Vo(bodyind)*hydro(hydroind).rho;
        end

        %% Read Inertia
        if hydro(hydroind).Nb == 1
            fileID = fopen (fullfile (meshdir, 'Inertia_hull.dat'));
        else
            fileID = fopen ( [ fullfile(meshdir, 'Inertia_'),num2str(bodyind-1), '.dat' ] );
        end
        
        if fileID ~= -1
        
            for i=1:3
                ligne=fscanf(fileID,'%g %g',3);
                M(i+3,4:6)=ligne;
            end

            KHyd = squeeze(hydro(hydroind).C);
            for k=1:length(w)
                for j=1:nb_DOF
        %           RAO1(k,j)=(Fe(j,k))/(-(M+A(j,j,k))*(w(k))^2-1i*w(k)*(B(j,j,k))+KHyd(j,j)); % No coupling between the DoF
                    RAO1(k,j)=(Fe(j,k))/(-(M(j,j)+A(j,j,k))*(w(k))^2-1i*w(k)*(B(j,j,k))+KHyd(j,j)); % No coupling between the DoF
                end
            end

            RAO = RAO1;
        
        %------Initialisation-----
            first_constant= zeros(1,nw);
            second_constant= zeros(1,nw);
            Fdrift_x=zeros(1,nw);
            Fdrift_y=zeros(1,nw);
            H=zeros(ntheta,nw);
            ampl_wave = 1;

        %--------------- CALCULATION-----------------------%
            Kochin_BVP_complex(:,:)=Kochin_BVP(:,1,:).*exp(1i*Kochin_BVP(:,2,:)); % H complex
            w=hydro(hydroind).w;
            depth=hydro(hydroind).h;
            for j=1:nw
                m0(j)=wave_number(w(j),depth);
            %       m0(j)=wave_number(w(j)/(2*pi),depth);
                k0(j)=(w(j)^2)/9.81; % wave number at an infinite water depth

                local1=zeros(ntheta,nb_DOF*nBodies+1);
                local1(:,1)=ampl_wave*Kochin_BVP_complex(:,(nb_DOF*nBodies+1)*(j-1)+1)*exp(1i*pi/2);%

                for i=2:(nb_DOF*nBodies+1)% sum of the radiation terms times velocities RAOs
                    x=(nb_DOF*nBodies+1)*(j-1)+i;
                    local1(:,i)=ampl_wave*(RAO(j,i-1))*(Kochin_BVP_complex(:,x))*exp(1i*pi/2)*(-1i*w(j));
                end
                H(:,j)=sum(local1,2); % H= Kochin function per frequency
                H_real(:,j)=real(H(:,j));
                H_imag(:,j)=imag(H(:,j));
            end
            rad = pi/180*hydro(hydroind).beta; % conversion degrees to radians
            ind_beta=find(abs(theta-rad)==min(abs(theta-rad))); % ind_beta used for determining the appropriate angle of H(dir)
            ind_beta=min(ind_beta); % in case of 2 min found
            for j=1:nw
                % FORMULA (2.170) in Delhommeau Thesis
                first_constant(j)=-2*pi*ampl_wave*hydro(hydroind).rho*w(j);
                second_constant(j)=-(8*pi*hydro(hydroind).rho*m0(j)*(k0(j)*depth)^2)/(depth*(m0(j)^2*depth^2-k0(j)^2*depth^2+k0(j)*depth));
                Fdrift_x(j)=first_constant(j)*cos(rad)*imag(H(ind_beta,j)) + second_constant(j)*imag(trapz(theta,H_real(:,j).*imag(conj(H(:,j))).*cos(theta')));
                Fdrift_y(j)=first_constant(j)*sin(rad)*imag(H(ind_beta,j)) + second_constant(j)*imag(trapz(theta,H_real(:,j).*imag(conj(H(:,j))).*sin(theta')));
            end
            hydro(hydroind).md_mc=hydro(hydroind).ex_ma.*0;
            hydro(hydroind).md_mc(1,1,:) = Fdrift_x./hydro(hydroind).rho./9.81;
            hydro(hydroind).md_mc(3,1,:) = Fdrift_y./hydro(hydroind).rho./9.81;
        
        end
    
    end

    hydro = Normalize(hydro);  % Normalize the data according the WAMIT convention

end
