%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef hydrobody < handle

    properties (SetAccess = 'private', GetAccess = 'public') % hdf5 file
        hydroData         = struct()                                            % Hydrodynamic data from BEM or user defined.
    end

    properties (SetAccess = 'public', GetAccess = 'public') % input file
        name              = []                                                  % Body name. For WEC bodies this is given in the h5 file.
        mass              = []                                                  % Mass in kg or specify 'equilibrium' to have mass= dis vol * density
        momOfInertia      = []                                                  % Moment of inertia [Ixx Iyy Izz] in kg*m^2
        cg                = []                                                  % Center of gravity [x y z] in meters. For WEC bodies this is given in the h5 file.
        dispVol           = []                                                  % Displaced volume at equilibrium position in meters cubed. For WEC bodies this is given in the h5 file.
        geometryFile      = 'NONE'                                              % Location of geomtry stl files
        viscDrag          = struct(...                                          % Structure defining the viscous (quadratic) drag
                                   'cd',                   [0 0 0 0 0 0], ...       % Viscous (quadratic) drag cd, vector length 6
                                   'characteristicArea',   [0 0 0 0 0 0])           % Characteristic area for viscous drag, vector length 6
        initDisp          = struct(...                                          % Structure defining the initial displacement
                                   'initLinDisp',          [0 0 0], ...             % Initial displacement of center fo gravity - used for decay tests (format: [displacment in m], default = [0 0 0])
                                   'initAngularDispAxis',  [0 1 0], ...             % Initial displacement of cog - axis of rotation - used for decay tests (format: [x y z], default = [1 0 0])
                                   'initAngularDispAngle', 0)                       % Initial displacement of cog - Angle of rotation - used for decay tests (format: [radians], default = 0)
        linearDamping     = [0 0 0 0 0 0]                                       % Linear drag coefficient, vector length 6
        userDefinedExcIRF = []                                                  % Excitation IRF from BEMIO used for User-Defined Time-Series
        viz               = struct(...                                          % Structur defining visualization properties
                                   'color', [1 1 0], ...                            % Visualization color for either SimMechanics Explorer or Paraview.
                                   'opacity', 1)                                    % Visualization opacity for either SimMechanics Explorer or Paraview.
        morrisonElement   = struct(...                                          % Structure defining the Morrison Elements
                                   'cd',                 [0 0 0], ...               % Viscous (quadratic) drag cd, vector length 3
                                   'ca',                 [0 0 0], ...               % Added mass coefficent for Morrison Element (format [Ca_x Ca_y Ca_z], default = [0 0 0])
                                   'characteristicArea', [0 0 0], ...               % Characteristic area for Morrison Elements calculations (format [Area_x Area_y Area_z], default = [0 0 0])
                                   'VME',                 0     , ...               % Characteristic volume for Morrison Element (default = 0)
                                   'rgME',               [0 0 0])                   % Vector from center of gravity to point of application for Morrison Element (format [X Y Z], default = [0 0 0]).
%         nhBody            = 0                                                   % Flag for non-hydro body
    end

    properties (SetAccess = 'public', GetAccess = 'public') % body geometry stl file
        bodyGeometry      = struct(...                                          % Structure defining body's mesh
                                   'numFace', [], ...                               % Number of faces
                                   'numVertex', [], ...                             % Number of vertices
                                   'vertex', [], ...                                % List of vertices
                                   'face', [], ...                                  % List of faces
                                   'norm', [], ...                                  % List of normal vectors
                                   'area', [], ...                                  % List of cell areas
                                   'center', [])                                    % List of cell centers
    end

    properties (SetAccess = 'public', GetAccess = 'public') %internal
        hydroForce        = struct()                                            % Hydrodynamic forces and coefficients used during simulation.
        h5File            = ''                                                  % hdf5 file containing the hydrodynamic data
        hydroDataBodyNum  = []                                                  % Body number within the hdf5 file.
        massCalcMethod    = []                                                  % Method used to obtain mass: 'user', 'fixed', 'equilibrium'
        bodyNumber        = []                                                  % bodyNumber in WEC-Sim as defined in the input file. Can be different from the BEM body number.
        bodyTotal         = []                                                  % Total number of WEC-Sim bodies (body block iterations)
        lenJ              = []                                                  % Matrices length. 6 for no body-to-body interactions. 6*numBodies if body-to-body interactions.
    end


    properties (SetAccess = 'private', GetAccess = 'public') % internal

        excitationMethod;
        doNonLinearFKExcitation;
        radiationMethod;
        hydroRestoringForceMethod;
        freeSurfaceMethod;
        bodyToBodyInteraction;
        doMorrisonElementViscousDrag;
        caseDir;

    end

    properties (SetAccess = 'private', GetAccess = 'private') % internal

        waves;
        simu;

        excitationMethodNum;
        radiationMethodNum;
        hydroRestoringForceMethodNum;
        freeSurfaceMethodNum;

        % properties used to calculate nonFKForce at reduced sample time
        oldForce    = [];
        oldWp       = [];
        oldWpMeanFs = [];

        % nonlinear buoyancy
        oldNonLinBuoyancyF = [];
        oldNonLinBuoyancyP = [];
        
        % history of last few time steps
        timeStepHist;
            
        % accel history store (used by added mass calc transport delay)
        accelHist;
        
        stepCount;

        % wave radiation force convolution integral states
        CIdt;
        radForceVelocity;
        radForceOldTime;
        radForceOldF_FM;
        radForce_IRKB_interp;

        % wave radiation forces state-space system object
        radForceSS;
        
        % wave elevation
        oldElev;

    end

    % public pre-processing related methods (and constructor)
    methods (Access = 'public') %modify object = T; output = F

        function obj = hydrobody (filename, varargin)
            % constructor for the hydrobody class
            %
            % Syntax
            %
            % hb = hydrobody (filename)
            % hb = hydrobody (filename, caseDir)
            %
            % Input
            %
            %  filename - string containing the h5 file containing the
            %    hydrodynamic data for the body
            %
            %  CaseDirectory - optional string containng the path to the 
            %    case directory. If not supplied, the case directory is
            %    assumed to be the current working directory.
            %
            % Output
            %
            %  hb - a hydrobody object
            %
            %
            
%             options.Waves = [];
%             options.Simu = [];
            options.CaseDirectory = pwd();
            
            options = parse_pv_pairs (options, varargin);
            
            if exist (fullfile (options.CaseDirectory, filename), 'file') ~= 2
                if exist (options.CaseDirectory, 'dir') ~= 7
                    error ('The case directory %s does not appear to exist', options.CaseDirectory);
                else
                    error ('The h5 file was not found in the case directory:\n%s', options.CaseDirectory);
                end
            end

            % hydro data file
            obj.h5File = filename;
            obj.caseDir = options.CaseDirectory;

        end

        function readH5File (obj)
            % Reads the HDF5 file containing the hydrodynamic data for the
            % body
            %
            % Syntax
            %
            % readH5File(hb)
            %
            % Description
            %
            % readH5File reads the HDF5 file containing the hydrodynamic
            % data for the body and stores it in the body properties. The
            % file location is expected to be in the case directory
            % provided on construction of the object.
            %
            % Generating the hydrodynamic data:
            %
            % The hydrobody requires frequency-domain hydrodynamic
            % coefficients (added mass, radiation damping, and wave
            % excitation). Typically, these hydrodynamic coefficients for
            % each body of the WEC device are generated using a boundary
            % element method (BEM) code (e.g., WAMIT, NEMOH or AQWA). The
            % HDF5 file must then be generated from the output of these
            % codes.
            %
            % Create HDF5 file:
            %
            % readH5File reads the hydrodynamic data in HDF5 format from
            % the (<hydrodata_file_name>.h5) file provided when the object
            % was constructed. A helper tool, BEMIO, is available to parse
            % BEM solutions (from WAMIT, NEMOH and AQWA) into the required
            % HDF5 data structure. 
            %
            % Input
            %
            %  hb - hydrobody object
            %
            %
            
            filename = fullfile (obj.caseDir, obj.h5File);
            name = ['/body' num2str(obj.bodyNumber)];
            obj.cg = h5read(filename,[name '/properties/cg']);
            obj.cg = obj.cg';
            obj.dispVol = h5read(filename,[name '/properties/disp_vol']);
            obj.name = h5read(filename,[name '/properties/name']);
            try obj.name = obj.name{1}; end
            obj.hydroData.simulation_parameters.scaled = h5read(filename,'/simulation_parameters/scaled');
            obj.hydroData.simulation_parameters.wave_dir = h5read(filename,'/simulation_parameters/wave_dir');
            obj.hydroData.simulation_parameters.water_depth = h5read(filename,'/simulation_parameters/water_depth');
            obj.hydroData.simulation_parameters.w = h5read(filename,'/simulation_parameters/w');
            obj.hydroData.simulation_parameters.T = h5read(filename,'/simulation_parameters/T');
            obj.hydroData.properties.name = h5read(filename,[name '/properties/name']);
            try obj.hydroData.properties.name = obj.hydroData.properties.name{1}; end
            obj.hydroData.properties.body_number = h5read(filename,[name '/properties/body_number']);
            obj.hydroData.properties.cg = h5read(filename,[name '/properties/cg']);
            obj.hydroData.properties.disp_vol = h5read(filename,[name '/properties/disp_vol']);
            obj.hydroData.hydro_coeffs.linear_restoring_stiffness = h5load(filename, [name '/hydro_coeffs/linear_restoring_stiffness']);
            obj.hydroData.hydro_coeffs.excitation.re = h5load(filename, [name '/hydro_coeffs/excitation/re']);
            obj.hydroData.hydro_coeffs.excitation.im = h5load(filename, [name '/hydro_coeffs/excitation/im']);
            try obj.hydroData.hydro_coeffs.excitation.impulse_response_fun.f = h5load(filename, [name '/hydro_coeffs/excitation/impulse_response_fun/f']); end
            try obj.hydroData.hydro_coeffs.excitation.impulse_response_fun.t = h5load(filename, [name '/hydro_coeffs/excitation/impulse_response_fun/t']); end
            obj.hydroData.hydro_coeffs.added_mass.all = h5load(filename, [name '/hydro_coeffs/added_mass/all']);
            obj.hydroData.hydro_coeffs.added_mass.inf_freq = h5load(filename, [name '/hydro_coeffs/added_mass/inf_freq']);
            obj.hydroData.hydro_coeffs.radiation_damping.all = h5load(filename, [name '/hydro_coeffs/radiation_damping/all']);
            try obj.hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.K = h5load(filename, [name '/hydro_coeffs/radiation_damping/impulse_response_fun/K']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.t = h5load(filename, [name '/hydro_coeffs/radiation_damping/impulse_response_fun/t']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.state_space.it = h5load(filename, [name '/hydro_coeffs/radiation_damping/state_space/it']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.state_space.A.all = h5load(filename, [name '/hydro_coeffs/radiation_damping/state_space/A/all']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.state_space.B.all = h5load(filename, [name '/hydro_coeffs/radiation_damping/state_space/B/all']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.state_space.C.all = h5load(filename, [name '/hydro_coeffs/radiation_damping/state_space/C/all']); end
            try obj.hydroData.hydro_coeffs.radiation_damping.state_space.D.all = h5load(filename, [name '/hydro_coeffs/radiation_damping/state_space/D/all']); end
        end

        function loadHydroData (obj, hydroData)
            % Loads hydroData structure from matlab variable as alternative
            % to reading the h5 file. Used when doing multiple condition
            % runs on the same system
            %
            
            obj.hydroData = hydroData;
            obj.cg        = hydroData.properties.cg';
            obj.dispVol   = hydroData.properties.disp_vol;
            obj.name      = hydroData.properties.name;
        end

        function hydroForcePre(obj)
            % HydroForce Pre-processing calculations
            % 1. Set the linear hydrodynamic restoring coefficient, viscous
            %    drag, and linear damping matrices
            % 2. Set the wave excitation force
                  
            rho = obj.simu.rho;
            g = obj.simu.g;
            
            obj.setMassMatrix(rho, obj.simu.nlHydro)
            
            k = obj.hydroData.hydro_coeffs.linear_restoring_stiffness;
            
            obj.hydroForce.linearHydroRestCoef = k .* rho .* g;
            
            obj.hydroForce.visDrag = diag (0.5 * rho .* obj.viscDrag.cd .* obj.viscDrag.characteristicArea);
            
            obj.hydroForce.linearDamping = diag (obj.linearDamping);
            
            obj.hydroForce.userDefinedFe = zeros (length (obj.waves.waveAmpTime(:,2)),6);   %initializing userDefinedFe for non user-defined cases
            
            switch obj.waves.type
                case {'noWave'}
                    obj.noExcitation ()
                    obj.constAddedMassAndDamping ();
                    obj.excitationMethodNum = 0;
                    
                case {'noWaveCIC'}
                    obj.noExcitation ()
                    obj.irfInfAddedMassAndDamping ();
                    obj.excitationMethodNum = 0;
                    
                case {'regular'}
                    obj.regExcitation ();
                    obj.constAddedMassAndDamping ();
                    obj.excitationMethodNum = 1;
                    
                case {'regularCIC'}
                    obj.regExcitation ();
                    obj.irfInfAddedMassAndDamping ();
                    obj.excitationMethodNum = 1;
                    
                case {'irregular','irregularImport'}
                    obj.irrExcitation ();
                    obj.irfInfAddedMassAndDamping ();
                    obj.excitationMethodNum = 2;
                    
                case {'userDefined'}
                    obj.userDefinedExcitation (obj.waves.waveAmpTime);
                    obj.irfInfAddedMassAndDamping ();
                    obj.excitationMethodNum = 3;
                    
            end

            % store the description of the wave type for information later
            obj.excitationMethod = obj.waves.type;

        end

        function adjustMassMatrix(obj)
            % Merge diagonal term of added mass matrix to the mass matrix
            % 1. Store the original mass and added-mass properties
            % 2. Add diagonal added-mass inertia to moment of inertia
            % 3. Add the maximum diagonal translational added-mass to body mass
            iBod = obj.bodyNumber;
            obj.hydroForce.storage.mass = obj.mass;
            obj.hydroForce.storage.momOfInertia = obj.momOfInertia;
            obj.hydroForce.storage.fAddedMass = obj.hydroForce.fAddedMass;
            
            if obj.simu.b2b == true
                
                tmp.fadm = diag (obj.hydroForce.fAddedMass(:,1+(iBod-1)*6:6+(iBod-1)*6));
                tmp.adjmass = sum(tmp.fadm(1:3)) * obj.simu.adjMassWeightFun;
                obj.mass = obj.mass + tmp.adjmass;
                obj.momOfInertia = obj.momOfInertia+tmp.fadm(4:6)';
                obj.hydroForce.fAddedMass(1,1+(iBod-1)*6) = obj.hydroForce.fAddedMass(1,1+(iBod-1)*6) - tmp.adjmass;
                obj.hydroForce.fAddedMass(2,2+(iBod-1)*6) = obj.hydroForce.fAddedMass(2,2+(iBod-1)*6) - tmp.adjmass;
                obj.hydroForce.fAddedMass(3,3+(iBod-1)*6) = obj.hydroForce.fAddedMass(3,3+(iBod-1)*6) - tmp.adjmass;
                obj.hydroForce.fAddedMass(4,4+(iBod-1)*6) = 0;
                obj.hydroForce.fAddedMass(5,5+(iBod-1)*6) = 0;
                obj.hydroForce.fAddedMass(6,6+(iBod-1)*6) = 0;
                
            else
                
                tmp.fadm = diag(obj.hydroForce.fAddedMass);
                tmp.adjmass = sum (tmp.fadm(1:3)) * obj.simu.adjMassWeightFun;
                obj.mass = obj.mass + tmp.adjmass;
                obj.momOfInertia = obj.momOfInertia + tmp.fadm(4:6)';
                obj.hydroForce.fAddedMass(1,1) = obj.hydroForce.fAddedMass(1,1) - tmp.adjmass;
                obj.hydroForce.fAddedMass(2,2) = obj.hydroForce.fAddedMass(2,2) - tmp.adjmass;
                obj.hydroForce.fAddedMass(3,3) = obj.hydroForce.fAddedMass(3,3) - tmp.adjmass;
                obj.hydroForce.fAddedMass(4,4) = 0;
                obj.hydroForce.fAddedMass(5,5) = 0;
                obj.hydroForce.fAddedMass(6,6) = 0;
                
            end
        end

        function restoreMassMatrix(obj)
            % Restore the mass and added-mass matrix back to the original value
            tmp = struct;
            tmp.mass = obj.mass;
            tmp.momOfInertia = obj.momOfInertia;
            tmp.hydroForce_fAddedMass = obj.hydroForce.fAddedMass;
            obj.mass = obj.hydroForce.storage.mass;
            obj.momOfInertia = obj.hydroForce.storage.momOfInertia;
            obj.hydroForce.fAddedMass = obj.hydroForce.storage.fAddedMass;
            obj.hydroForce.storage = tmp; clear tmp
        end

        function storeForceAddedMass(obj,am_mod,ft_mod)
            % Store the modified added mass and total forces history (inputs)
            obj.hydroForce.storage.output_forceAddedMass = am_mod;
            obj.hydroForce.storage.output_forceTotal = ft_mod;
        end

        function setInitDisp(obj, x_rot, ax_rot, ang_rot, addLinDisp)
            % Function to set the initial displacement when having initial rotation
            % x_rot: rotation point
            % ax_rot: axis about which to rotate (must be a normal vector)
            % ang_rot: rotation angle in radians
            % addLinDisp: initial linear displacement (in addition to the displacement caused by rotation)

            relCoord = obj.cg - x_rot;

            rotatedRelCoord = hydrobody.rotateXYZ(relCoord, ax_rot, ang_rot);

            newCoord = rotatedRelCoord + x_rot;

            linDisp = newCoord - obj.cg;

            obj.initDisp.initLinDisp = linDisp + addLinDisp;

            obj.initDisp.initAngularDispAxis = ax_rot;

            obj.initDisp.initAngularDispAngle = ang_rot;

        end

        function listInfo(obj)
            % List body info
            fprintf('\n\t***** Body Number %G, Name: %s *****\n',obj.hydroData.properties.body_number,obj.hydroData.properties.name)
            fprintf('\tBody CG                          (m) = [%G,%G,%G]\n',obj.hydroData.properties.cg)
            fprintf('\tBody Mass                       (kg) = %G \n',obj.mass);
            fprintf('\tBody Diagonal MOI              (kgm2)= [%G,%G,%G]\n',obj.momOfInertia)
        end

        function bodyGeo(obj,fname)
            % Reads mesh file and calculates areas and centroids
            try
                [obj.bodyGeometry.vertex, obj.bodyGeometry.face, obj.bodyGeometry.norm] = import_stl_fast (fname,1,1);
            catch
                [obj.bodyGeometry.vertex, obj.bodyGeometry.face, obj.bodyGeometry.norm] = import_stl_fast (fname,1,2);
            end
            obj.bodyGeometry.numFace = length (obj.bodyGeometry.face);
            obj.bodyGeometry.numVertex = length (obj.bodyGeometry.vertex);
            obj.checkStl ();
            obj.triArea ();
            obj.triCenter ();
        end

        function triArea(obj)
            % Function to calculate the area of a triangle
            points = obj.bodyGeometry.vertex;
            faces = obj.bodyGeometry.face;
            v1 = points(faces(:,3),:) - points(faces(:,1),:);
            v2 = points(faces(:,2),:) - points(faces(:,1),:);
            av_tmp =  1/2 .* (cross(v1,v2));
            area_mag = sqrt (av_tmp(:,1).^2 + av_tmp(:,2).^2 + av_tmp(:,3).^2);
            obj.bodyGeometry.area = area_mag;
        end

        function checkStl(obj)
            % Function to check STL file
            tnorm = obj.bodyGeometry.norm;
            %av = zeros(length(area_mag),3);
            %av(:,1) = area_mag.*tnorm(:,1);
            %av(:,2) = area_mag.*tnorm(:,2);
            %av(:,3) = area_mag.*tnorm(:,3);
            %if sum(sum(sign(av_tmp))) ~= sum(sum(sign(av)))
            %    warning(['The order of triangle vertices in ' obj.geometryFile ' do not follow the right hand rule. ' ...
            %        'This will causes visualization errors in the SimMechanics Explorer'])
            %end
            norm_mag = sqrt (tnorm(:,1).^2 + tnorm(:,2).^2 + tnorm(:,3).^2);
            check = sum(norm_mag)/length(norm_mag);
            if check > 1.01 || check < 0.99
                error(['length of normal vectors in ' obj.geometryFile ' is not equal to one.'])
            end
        end

        function triCenter(obj)
            % Calculate the center coordinate of the body geometries'
            % triangles
            %
            points = obj.bodyGeometry.vertex;
            faces = obj.bodyGeometry.face;
            c = zeros (length (faces), 3);
            c(:,1) = ( points(faces(:,1),1) + points(faces(:,2),1) + points(faces(:,3),1) ) ./ 3;
            c(:,2) = ( points(faces(:,1),2) + points(faces(:,2),2) + points(faces(:,3),2) ) ./ 3;
            c(:,3) = ( points(faces(:,1),3) + points(faces(:,2),3) + points(faces(:,3),3) ) ./ 3;
            obj.bodyGeometry.center = c;
        end

        function plotStl(obj)
            % Plots the body's mesh and normal vectors
            c = obj.bodyGeometry.center;
            tri = obj.bodyGeometry.face;
            p = obj.bodyGeometry.vertex;
            n = obj.bodyGeometry.norm;
            figure()
            hold on
            trimesh (tri,p(:,1),p(:,2),p(:,3))
            quiver3 (c(:,1),c(:,2),c(:,3),n(:,1),n(:,2),n(:,3))
        end

        function checkinputs(obj)
            % Checks the user inputs
            % hydro data file
            if exist(fullfile (obj.caseDir, obj.h5File),'file')==0 % && obj.nhBody==0
                error('The hdf5 file %s does not exist', ...
                    fullfile (obj.caseDir, obj.h5File))
            end
            % geometry file
            if exist(fullfile (obj.caseDir, obj.geometryFile),'file') == 0
                error('Could not locate and open geometry file %s', ...
                    fullfile (obj.caseDir, obj.geometryFile))
            end
        end
        
        function [node, body] = makeMBDynComponents (obj)
            
            gref = mbdyn.pre.globalref;
            
            ref_hydrobody = mbdyn.pre.reference ( obj.cg, [], [], [], 'Parent', gref);

            node = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', ref_hydrobody.pos);

            body = mbdyn.pre.body ( obj.mass,  ...
                                  [0;0;0], ...
                                  diag (obj.momOfInertia), ...
                                  node, ...
                                  'STLFile', fullfile (obj.caseDir, obj.geometryFile) );

            
        end
        
        
    end

    % non-public pre-processing methods
    methods (Access = 'protected') %modify object = T; output = F


        function noExcitation(obj)
            % Set excitation force for no excitation case
            obj.hydroForce.fExt.re=zeros(1,6);
            obj.hydroForce.fExt.im=zeros(1,6);
        end

        function regExcitation(obj)
            % Regular wave excitation force
            % Used by hydroForcePre
            re = obj.hydroData.hydro_coeffs.excitation.re(:,:,:) .* obj.simu.rho .* obj.simu.g;
            im = obj.hydroData.hydro_coeffs.excitation.im(:,:,:) .* obj.simu.rho .* obj.simu.g;
            obj.hydroForce.fExt.re = zeros(1,6);
            obj.hydroForce.fExt.im = zeros(1,6);
            
            for ii=1:6
                
                if length(obj.hydroData.simulation_parameters.wave_dir) > 1
                    
                    [X,Y] = meshgrid ( obj.hydroData.simulation_parameters.w, ...
                                       obj.hydroData.simulation_parameters.wave_dir);
                    
                    obj.hydroForce.fExt.re(ii) = interp2 (X, Y, squeeze(re(ii,:,:)), obj.waves.w, obj.waves.waveDir);
                    obj.hydroForce.fExt.im(ii) = interp2 (X, Y, squeeze(im(ii,:,:)), obj.waves.w, obj.waves.waveDir);
                    
                elseif obj.hydroData.simulation_parameters.wave_dir == obj.waves.waveDir
                    
                    obj.hydroForce.fExt.re(ii) = interp1 (obj.hydroData.simulation_parameters.w, squeeze(re(ii,1,:)), obj.waves.w, 'spline');
                    obj.hydroForce.fExt.im(ii) = interp1 (obj.hydroData.simulation_parameters.w, squeeze(im(ii,1,:)), obj.waves.w, 'spline');
                    
                end
                
            end
            
        end

        function irrExcitation(obj)
            % Irregular wave excitation force
            % Used by hydroForcePre
            re = obj.hydroData.hydro_coeffs.excitation.re(:,:,:) .* obj.simu.rho .* obj.simu.g;
            im = obj.hydroData.hydro_coeffs.excitation.im(:,:,:) .* obj.simu.rho .* obj.simu.g;
            
            obj.hydroForce.fExt.re = zeros (obj.waves.numFreq, 6);
            obj.hydroForce.fExt.im = zeros (obj.waves.numFreq, 6);
            
            for ii=1:6
                if length (obj.hydroData.simulation_parameters.wave_dir) > 1
                    
                    [X,Y] = meshgrid (obj.hydroData.simulation_parameters.w, ...
                                      obj.hydroData.simulation_parameters.wave_dir);
                                  
                    obj.hydroForce.fExt.re(:,ii) = interp2 (X, Y, squeeze(re(ii,:,:)), obj.waves.w, obj.waves.waveDir);
                    
                    obj.hydroForce.fExt.im(:,ii) = interp2 (X, Y, squeeze(im(ii,:,:)), obj.waves.w, obj.waves.waveDir);
                    
                elseif obj.hydroData.simulation_parameters.wave_dir == obj.waves.waveDir
                    
                    obj.hydroForce.fExt.re(:,ii) = interp1 (obj.hydroData.simulation_parameters.w, squeeze(re(ii,1,:)), obj.waves.w, 'spline');
                    
                    obj.hydroForce.fExt.im(:,ii) = interp1 (obj.hydroData.simulation_parameters.w, squeeze(im(ii,1,:)), obj.waves.w, 'spline');
                    
                end
            end
        end

        function userDefinedExcitation(obj, waveAmpTime)
            % Calculated User-Defined wave excitation force with non-causal
            % convolution Used by hydroForcePre
            kf = obj.hydroData.hydro_coeffs.excitation.impulse_response_fun.f .* obj.simu.rho .* obj.simu.g;
            kt = obj.hydroData.hydro_coeffs.excitation.impulse_response_fun.t;
            t =  min(kt):obj.simu.dt:max(kt);
            
            for ii = 1:6
                if length (obj.hydroData.simulation_parameters.wave_dir) > 1
                    
                    [X,Y] = meshgrid (kt, obj.hydroData.simulation_parameters.wave_dir);
                    
                    kernel = squeeze (kf(ii,:,:));
                    
                    obj.userDefinedExcIRF = interp2 (X, Y, kernel, t, obj.waves.waveDir);
                    
                elseif obj.hydroData.simulation_parameters.wave_dir == obj.waves.waveDir
                    
                    kernel = squeeze (kf(ii,1,:));
                    
                    obj.userDefinedExcIRF = interp1 (kt,kernel,min(kt):obj.simu.dt:max(kt));
                    
                end
                obj.hydroForce.userDefinedFe(:,ii) = conv (waveAmpTime(:,2), obj.userDefinedExcIRF, 'same') * obj.simu.dt;
                
            end
            
            obj.hydroForce.fExt.re = zeros(1,6);
            obj.hydroForce.fExt.im = zeros(1,6);
            
        end

        function constAddedMassAndDamping(obj)
            % Set added mass and damping for a specific frequency
            % Used by hydroForcePre
            
            am = obj.hydroData.hydro_coeffs.added_mass.all .* obj.simu.rho;
            
            rd = obj.hydroData.hydro_coeffs.radiation_damping.all .* obj.simu.rho;
            
            for i = 1:length (obj.hydroData.simulation_parameters.w)
                rd(:,:,i) = rd(:,:,i) .* obj.hydroData.simulation_parameters.w(i);
            end
            % Change matrix size: B2B [6x6n], noB2B [6x6]
            if obj.simu.b2b == true
                
                lenJ = 6*obj.bodyTotal;
                obj.hydroForce.fAddedMass = zeros(6,lenJ);
                obj.hydroForce.fDamping = zeros(6,lenJ);
                obj.hydroForce.totDOF  =zeros(6,lenJ);
                for ii=1:6
                    for jj=1:lenJ
                        obj.hydroForce.fAddedMass(ii,jj) = interp1 (obj.hydroData.simulation_parameters.w,squeeze(am(ii,jj,:)), obj.waves.w, 'spline');
                        obj.hydroForce.fDamping  (ii,jj) = interp1 (obj.hydroData.simulation_parameters.w,squeeze(rd(ii,jj,:)), obj.waves.w, 'spline');
                    end
                end
                    
            else
                lenJ = 6;
                obj.hydroForce.fAddedMass = zeros(6,lenJ);
                obj.hydroForce.fDamping = zeros(6,lenJ);
                obj.hydroForce.totDOF  =zeros(6,lenJ);

                for ii=1:6
                    for jj=1:lenJ
                        jjj = 6*(obj.bodyNumber-1)+jj;
                        obj.hydroForce.fAddedMass(ii,jj) = interp1 (obj.hydroData.simulation_parameters.w,squeeze(am(ii,jjj,:)), obj.waves.w, 'spline');
                        obj.hydroForce.fDamping  (ii,jj) = interp1 (obj.hydroData.simulation_parameters.w,squeeze(rd(ii,jjj,:)), obj.waves.w, 'spline');
                    end
                end
                    
            end
        end

        function irfInfAddedMassAndDamping(obj)
            % Set radiation force properties using impulse response function
            % Used by hydroForcePre
            % Added mass at infinite frequency
            % Convolution integral raditation damping
            % State space formulation
            
            iBod = obj.simu.numWecBodies;
            
            if obj.simu.b2b == true
                lenJ = obj.bodyTotal*6;
            else
                lenJ = 6;
            end
            
            % Convolution integral formulation
            if obj.simu.b2b == true
                obj.hydroForce.fAddedMass = obj.hydroData.hydro_coeffs.added_mass.inf_freq .* obj.simu.rho;
            else
                obj.hydroForce.fAddedMass = obj.hydroData.hydro_coeffs.added_mass.inf_freq(1:6,(iBod-1)*6+1:(iBod-1)*6+6) .* obj.simu.rho;
            end
            
            % Radition IRF
            obj.hydroForce.fDamping=zeros(6,lenJ);
            irfk = obj.hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.K .* obj.simu.rho;
            irft = obj.hydroData.hydro_coeffs.radiation_damping.impulse_response_fun.t;
            %obj.hydroForce.irkb=zeros(CIkt,6,lenJ);
            if obj.simu.b2b == true
                for ii=1:6
                    for jj=1:lenJ
                        obj.hydroForce.irkb(:,ii,jj) = interp1 (irft,squeeze(irfk(ii,jj,:)), obj.simu.CTTime, 'spline');
                    end
                end
            else
                for ii=1:6
                    for jj=1:lenJ
                        jjj = (iBod-1)*6+jj;
                        obj.hydroForce.irkb(:,ii,jj) = interp1 (irft,squeeze(irfk(ii,jjj,:)), obj.simu.CTTime, 'spline');
                    end
                end
            end
            % State Space Formulation
            if obj.simu.ssCalc > 0
                
                if obj.simu.b2b == 1
                    
                    for ii = 1:6
                        
                        for jj = 1:lenJ
                            
                            arraySize = obj.hydroData.hydro_coeffs.radiation_damping.state_space.it(ii,jj);
                            
                            if ii == 1 && jj == 1 % Begin construction of combined state, input, and output matrices
                                
                                Af(1:arraySize,1:arraySize) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                                
                                Bf(1:arraySize,jj)          = obj.hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                                
                                Cf(ii,1:arraySize)          = obj.hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                                
                            else
                                
                                Af(size(Af,1)+1:size(Af,1)+arraySize,size(Af,2)+1:size(Af,2)+arraySize) ...
                                    = obj.hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                                
                                Bf(size(Bf,1)+1:size(Bf,1)+arraySize,jj) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                                
                                Cf(ii,size(Cf,2)+1:size(Cf,2)+arraySize) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                                
                            end
                            
                        end
                        
                    end
                    
                    obj.hydroForce.ssRadf.D = zeros (6,lenJ);
                    
                else
                    
                    for ii = 1:6
                        
                        for jj = (iBod-1)*6+1:(iBod-1)*6+6
                            
                            jInd = jj-(iBod-1)*6;
                            
                            arraySize = obj.hydroData.hydro_coeffs.radiation_damping.state_space.it(ii,jj);
                            
                            if ii == 1 && jInd == 1 % Begin construction of combined state, input, and output matrices
                                
                                Af(1:arraySize,1:arraySize) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                                Bf(1:arraySize,jInd)        = obj.hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                                Cf(ii,1:arraySize)          = obj.hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                                
                            else
                                
                                Af(size(Af,1)+1:size(Af,1)+arraySize,size(Af,2)+1:size(Af,2)+arraySize) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.A.all(ii,jj,1:arraySize,1:arraySize);
                                Bf(size(Bf,1)+1:size(Bf,1)+arraySize,jInd) = obj.hydroData.hydro_coeffs.radiation_damping.state_space.B.all(ii,jj,1:arraySize,1);
                                Cf(ii,size(Cf,2)+1:size(Cf,2)+arraySize)   = obj.hydroData.hydro_coeffs.radiation_damping.state_space.C.all(ii,jj,1,1:arraySize);
                                
                            end
                            
                        end
                        
                    end
                    
                    obj.hydroForce.ssRadf.D = zeros (6,6);
                end
                
                obj.hydroForce.ssRadf.A = Af;
                obj.hydroForce.ssRadf.B = Bf;
                obj.hydroForce.ssRadf.C = Cf .* obj.simu.rho;
                
            end
        end

        function setMassMatrix(obj, rho, nlHydro)
            % Sets mass for the special cases of body at equilibrium or fixed
            % Used by hydroForcePre
            if strcmp(obj.mass, 'equilibrium')
                obj.massCalcMethod = obj.mass;
                if nlHydro == 0
                    obj.mass = obj.hydroData.properties.disp_vol * rho;
                else
                    cg_tmp = obj.hydroData.properties.cg;
                    z = obj.bodyGeometry.center(:,3) + cg_tmp(3);
                    z(z>0) = 0;
                    area = obj.bodyGeometry.area;
                    av = [area area area] .* -obj.bodyGeometry.norm;
                    tmp = rho*[z z z].*-av;
                    obj.mass = sum(tmp(:,3));
                end
            elseif strcmp(obj.mass, 'fixed')
                obj.massCalcMethod = obj.mass;
                obj.mass = 999;
                obj.momOfInertia = [999 999 999];
            else
                obj.massCalcMethod = 'user';
            end
        end

    end

    % public transient simulation methods
    methods (Access = 'public')
        
        
        function odeSimSetup (obj, waves, simu, bodynum)
            % sets up the body in preparation for a transient simulation
            %
            % Syntax
            %
            % odeSimSetup (hb, waves, simu, bodynum)
            %
            % Desciription
            %
            % odeSimSetup initialises various parmaters and settings in
            % preparation for performing a transient simulation based on
            % the ODE solver routines. 
            %
            % Input
            %
            %   hb - hydrobody object
            %
            %   waves - waveClass object with the desired wave parameters
            %     to be used in the simulation
            %
            %   simu - simulationClass Object with the desired simulation
            %     parameters to be used in the simulation
            %
            %   bodynum - the number associated with this body. Typically
            %     this is generated by a parent hydrosys object
            %
            % 
            
            assert (isa (simu, 'wsim.simsettings'), 'waves must be a wsim.simsettings object')
            assert (isa (waves, 'wsim.wavesettings'), 'waves must be a wsim.wavesettings object');
            assert (isscalar (bodynum) &&  isint2eps (bodynum) , 'bodynum must be a scalar integer')

            % store waves and simu for later access
            obj.waves = waves;
            obj.simu = simu;
            obj.bodyNumber = bodynum;

            % Morrison Element
            if obj.simu.morrisonElement == 0
                obj.doMorrisonElementViscousDrag = false;
            elseif obj.simu.morrisonElement == 1
                obj.doMorrisonElementViscousDrag = true;
            end

            % Wave type

            % linear excitation type
            if obj.waves.typeNum < 10
                obj.excitationMethod = 'no waves';
                obj.excitationMethodNum = 0;
            elseif obj.waves.typeNum >= 10 && obj.waves.typeNum < 20
                obj.excitationMethod = 'regular waves';
                obj.excitationMethodNum = 1;
            elseif obj.waves.typeNum >= 20 && obj.waves.typeNum < 30
                obj.excitationMethod = 'irregular waves';
                obj.excitationMethodNum = 2;
            elseif obj.waves.typeNum >= 30
                obj.excitationMethod = 'user defined waves';
                obj.excitationMethodNum = 3;
            end

            % nonlinear excitation type
            if obj.simu.nlHydro == 0
                obj.doNonLinearFKExcitation = false;
            elseif obj.simu.nlHydro > 0
                obj.doNonLinearFKExcitation = true;
            end

            if obj.simu.nlHydro < 2
                obj.freeSurfaceMethod = 'mean';
                obj.freeSurfaceMethodNum = 0;
            elseif obj.simu.nlHydro == 2
                obj.freeSurfaceMethod = 'instantaneous';
                obj.freeSurfaceMethodNum = 1;
            end

            % Radiation Damping
            if obj.waves.typeNum == 0 || obj.waves.typeNum == 10 %'noWave' & 'regular'
                % constant radiation coefficients
                obj.radiationMethod = 'constant radiation coefficients';
                obj.radiationMethodNum = 0;
            elseif obj.simu.ssCalc == 1
                % state space radiation forces
                obj.radiationMethod = 'state space representation';
                obj.radiationMethodNum = 2;
            elseif obj.simu.ssCalc == 2
                % state space radiation forces
                obj.radiationMethod = 'state space representation using external solver';
                obj.radiationMethodNum = 3;
            else
                % convolution integral radiation forces
                obj.radiationMethod = 'convolution integral';
                obj.radiationMethodNum = 1;
            end

            % Body2Body
            if obj.simu.b2b == 0
                obj.bodyToBodyInteraction = false;
            else
                obj.bodyToBodyInteraction = true;
            end

            % first reset everything
            odeSimReset (obj);

            % then do the setup again
            obj.hydroForcePre ();

            % Radiation Damping
            if obj.radiationMethodNum == 1 && ~isempty (fieldnames(obj.hydroForce))

                % reset the radiation force convolution integral related states
                obj.CIdt = obj.simu.CTTime(2) - obj.simu.CTTime(1);

                obj.radForceVelocity = zeros(length(obj.lenJ),length(obj.simu.CTTime));

                obj.radForceOldTime = 0;

                obj.radForceOldF_FM = zeros(6,1);

                IRKB_reordered = permute(obj.hydroForce.irkb, [3 1 2]);

                interp_factor = 1;

                obj.radForce_IRKB_interp = IRKB_reordered(:,1:interp_factor:end, :);

            elseif obj.radiationMethodNum == 2

                % initialise the radiation force state space solver object
                obj.radForceSS = stateSpace ( obj.hydroForce.ssRadf.A, ...
                                              obj.hydroForce.ssRadf.B, ...
                                              obj.hydroForce.ssRadf.C, ...
                                              obj.hydroForce.ssRadf.D, ...
                                              zeros (size (obj.hydroForce.ssRadf.A,2), 1) );
                                          
                % initialise the fixed step integration
                ufcn = @(arg1, arg2) -K * arg2;

                obj.radForceSS.initIntegration (0, ufcn);
                
            end

        end
        
        
        function odeSimReset (obj)
            % resets the body in readiness for a transient simulation
            %
            % Syntax
            %
            % odeSimReset (hb)
            %
            % Desciription
            %
            % odeSimReset resets various internal storeage parameters and
            % settings in preparation for performing a transient simulation
            % based on the ODE solver routines, returning the hydrobody to
            % the state it is in just after calling odeSimSetup. This
            % should be called before re-running a transient simulation
            % with the same parameters.
            %
            % Input
            %
            %   hb - hydrobody object
            %
            %
            
            obj.stepCount = 0;
            
            nsteps = 3;
            
            % reset time history store
            obj.timeStepHist = zeros (1,nsteps);
            
            % reset accel history store (used by added mass calc transport
            % delay)
            if obj.simu.b2b
                obj.accelHist = zeros (nsteps, 6 * obj.simu.numWecBodies);
            else
                obj.accelHist = zeros (nsteps, 6);
            end
            
            % reset wave elevation
            obj.oldElev = [];
            
            obj.oldForce    = [];
            obj.oldWp       = [];
            obj.oldWpMeanFs = [];

            % reset non-linear buoyancy calc stuff
            obj.oldNonLinBuoyancyF = [];
            obj.oldNonLinBuoyancyP = [];

            % reset radiation force related stuff
            if isa (obj.radForceSS, 'stateSpace')
                obj.radForceSS.reset ();
            end

        end
        
        function advanceStep (obj, t, accel, vel)
            % advance to the next time step, accepting the current time
            % step and data into stored solution histories
            
            obj.timeStepHist = circshift (obj.timeStepHist, [0,-1]);
            
            obj.timeStepHist(end) = t;
            
            % update the acceleration history
            obj.accelHist = circshift (obj.accelHist, [-1,0]);
            
            obj.accelHist(end,:) = accel';
            
            obj.stepCount = obj.stepCount + 1;
            
        end
        
        function [forces, breakdown] = hydroForces (obj, t, pos, vel, accel)
            % hydroForces calculates the hydrodynamic forces acting on a
            % body
            %
            % Syntax
            %
            %  [forces, breakdown] = hydroForces (hb, t, x, vel, accel, elv)
            %
            % Input
            %
            %  hb - hydrobody object
            %
            %  t - global simulation time
            %
            %  pos - (6 x 1) displcement of this body in x,y and z and
            %    rotatation around the x, y and z axes
            %
            %  vel - (6 x n) velocities of all bodies in system
            %
            %  accel - (6 x n) acceleration of all bodies in system
            %
            %  elv - wave elevation
            %
            % Output
            %
            %  forces - (6 x 1) forces and moments acting on the body
            %
            %  breakdown - structure containing more detailed
            %    breakdown of the forces acting on the body. The fields
            %    present depend on the simulation settings and can include
            %    the following:
            %
            %    F_ExcitLin : 
            %
            %    F_Excit : 
            %
            %    F_ExcitRamp : 
            %
            %    F_ViscousDamping : 
            %
            %    F_addedmass : 
            %
            %    F_RadiationDamping : 
            %
            %    F_Restoring : 
            %
            %    BodyHSPressure : 
            %
            %    F_ExcitLinNonLin : 
            %
            %    WaveNonLinearPressure :
            %
            %    WaveLinearPressure : 
            %
            %    F_MorrisonElement : 
            %
            

            
            % always do linear excitation forces
            breakdown.F_ExcitLin = linearExcitationForces (obj, t);

            % always do viscous damping
            breakdown.F_ViscousDamping = viscousDamping (obj, vel(:,obj.bodyNumber));

            % always do radiation forces
            [breakdown.F_addedmass, breakdown.F_RadiationDamping] = radiationForces (obj, t, vel, accel);

            % calculate the wave elevation
            elv = waveElevation(obj, pos, t);
            
            % hydrostatic restoring forces
            [breakdown.F_Restoring, breakdown.BodyHSPressure ] =  hydrostaticForces (obj, t, pos, elv);

            if obj.doNonLinearFKExcitation

                [ breakdown.F_ExcitLinNonLin, ...
                  breakdown.WaveNonLinearPressure, ...
                  breakdown.WaveLinearPressure ] = nonlinearExcitationForces (obj, t, pos, elv);

            else
                breakdown.F_ExcitLinNonLin = zeros (6, 1);
            end

            if obj.doMorrisonElementViscousDrag

                breakdown.F_MorrisonElement = morrisonElementForce ( obj, t, pos, ...
                                                                     vel(:,obj.bodyNumber), ...
                                                                     accel(:,obj.bodyNumber) );

            else

                breakdown.F_MorrisonElement = zeros (6, 1);

            end


            breakdown.F_Excit = breakdown.F_ExcitLin + breakdown.F_ExcitLinNonLin;

            breakdown.F_ExcitRamp = applyRamp (obj, t, breakdown.F_Excit);

            forces = breakdown.F_ExcitRamp ...
                     - breakdown.F_ViscousDamping ...
                     - breakdown.F_addedmass ...
                     - breakdown.F_Restoring ...
                     - breakdown.F_RadiationDamping ...
                     - breakdown.F_MorrisonElement;

        end

        function forces = viscousDamping (obj, vel)

            forces = obj.hydroForce.visDrag * ( vel .* abs (vel));

        end

        function forces = morrisonElementForce (obj, t, pos, vel, accel)

            % TODO: convert morrison element simulink models
            switch obj.excitationMethodNum

                case 0
                    % no wave
                    forces = [0 0 0 0 0 0];

                case 1
                    % regular wave
                    forces = [0 0 0 0 0 0];

                case 2
                    % irregular wave
                    forces = [0 0 0 0 0 0];

            end

        end

        function forces = linearExcitationForces (obj, t)
            % calculates linear wave excitation forces during transient
            % simulation
            %
            % Syntax
            %
            % forces = linearExcitationForces (obj, t)
            %
            % Input

            switch obj.excitationMethodNum

                case 0
                    % no wave
                    forces = [0; 0; 0; 0; 0; 0];

                case 1
                    % regular wave

                    % Calculates the wave force, F_wave, for the case of Regular Waves.

                    % F_wave =   A * cos(w * t) * Re{Fext}
                    %            -  A * sin(w * t) * Im{Fext}

                    wt = obj.waves.w(1,:) .* t;

                    forces = obj.waves.A(1,:) .* ( ...
                                cos (wt) .* obj.hydroForce.fExt.re(1,:) ...
                                - sin (wt) .* obj.hydroForce.fExt.im(1,:) ...
                                             ).';

                case 2
                    % irregular wave

                    % Calculates the wave force, F_wave, for the case of Irregular Waves.
                    %
                    % F_wave = sum( F_wave(i))
                    %
                    % where i = each frequency bin.

                    % TOD: check correct dimension/orientation of force output
                    A1 = bsxfun (@plus, obj.waves.w * t, pi/2);

                    B1 = sin (bsxfun (@plus, A1, obj.waves.phaseRand));

                    B11 = sin (bsxfun (@plus, obj.waves.w * t, obj.waves.phaseRand));

                    C1 = sqrt (bsxfun (@times, obj.waves.A, obj.waves.dw));

                    D1 = bsxfun (@times, obj.hydroForce.fExt.re, C1);

                    D11 = bsxfun (@times, obj.hydroForce.fExt.im, C1);

                    E1 = bsxfun (@times, B1, D1);

                    E11 = bsxfun (@times, B11, D11);

                    forces = sum (bsxfun (@minus, E1, E11))';


                case 3
                    % user defined

                    % Calculates the wave force, F_wave, for the case of User Defined Waves.
                    %
                    % F_wave = convolution calculation [1x6]

                    error ('not yet implemented')
                    % TODO: make interpolation function for user defined waves, using ppval (C++ version)

            end

        end

        function [forces, wavenonlinearpressure, wavelinearpressure] = nonlinearExcitationForces (obj, t, pos, elv)

                pos = pos - [ obj.hydroData.properties.cg, 0, 0, 0];

                [forces, wavenonlinearpressure, wavelinearpressure]  = nonFKForce (obj, t, pos, elv);

        end

        function [F_addedmass, F_RadiationDamping] = radiationForces (obj, t, vel, accel)
            % calculates the wave radiation forces
            %
            % Syntax
            %
            % [F_addedmass, F_RadiationDamping] = radiationForces (obj, t, vel, accel)
            %
            % Input
            %
            %  t - current simulation time
            %
            %  vel - (6 x 1) translational and angular velocity of the
            %    body, or, if body to body interactions are being
            %    considered, a (6 x n) vector of all the body
            %    velocities.
            %
            %  accel - (6 x 1) translational and angular acceleration of
            %    the body, or, if body to body interactions are being
            %    considered, a (6 x n) vector of all the body
            %    accelerations.
            %
            %
            % Output
            %
            % F_addedmass - force due to added mass
            %
            % F_RadiationDamping - force due to wave radiation damping
            %
            %

            % matrix multiplication with acceleration
            delay = 10e-8;
            if t > (obj.simu.startTime + delay)
                if obj.stepCount > 2
                    % extrapolate acceleration from previous time steps
%                     thisaccel = interp1 (obj.timeStepHist', obj.accelHist, t-delay, 'linear', 'extrap');

%                     thisaccel = zeros (size (obj.accelHist,2), 1);
%                     for ind = 1:numel (thisaccel)
% %                         thisaccel(ind) = obj.lagrangeinterp (obj.timeStepHist',obj.accelHist(:,ind),t-delay);
% 
%                         p = polyfit (obj.timeStepHist(end-1:end)', obj.accelHist(end-1:end,ind), 1);
%                         thisaccel(ind) = polyval (p, t-delay);
%                     end

                thisaccel = obj.linearInterp ( obj.timeStepHist(end-1), ...
                                               obj.timeStepHist(end), ...
                                               obj.accelHist(end-1,:), ...
                                               obj.accelHist(end,:), ...
                                               t-delay );
                    
                elseif obj.stepCount == 2
                    thisaccel = interp1 (obj.timeStepHist(end-obj.stepCount+1:end)', obj.accelHist(end-obj.stepCount+1:end,:), t-delay, 'linear', 'extrap');
                elseif obj.stepCount == 1
                    thisaccel = obj.accelHist (end,:)';
                end
                F_addedmass = obj.hydroForce.fAddedMass * thisaccel(:);
%                 F_addedmass = obj.hydroForce.fAddedMass * accel(:);
            else
                F_addedmass = [0;0;0;0;0;0];
            end

            switch obj.radiationMethodNum

                case 0
                    % simple static coefficients
                    F_RadiationDamping = obj.hydroForce.fDamping * vel(:);

                case 1
                    % convolution
                    F_RadiationDamping = convolutionIntegral (obj, vel(:), t);

                case 2
                    % state space
                    F_RadiationDamping = obj.radForceSS.outputs (vel(:));
                    
                case 3
                    % state space, but calculated by some external solver,
                    % e.g. by supplying MBDyn with the state-space matrices
                    F_RadiationDamping = [0;0;0;0;0;0];

            end

        end
        
        function statederivs = radForceSSDerivatives (obj, u)
            
            statederivs = obj.radForceSS.derivatives (u);
            
        end
        
        function status = radForceODEOutputfcn (obj, t, x, flag, varargin)
            % OutputFcn to be called after every completed ode time step
            % when using the state-space representation of the radiation
            % forces
            %
            % Syntax
            %
            % radForceODEOutputfcn (hb, t, x, flag)
            %
            %
            % Input
            %
            %  hb - hydrobody oject
            %
            %  t - current time step
            %
            %  x - state variables at the current time step
            
            status = obj.radForceSS.outputfcn (t, x, flag);
            
        end

        function [forces, body_hspressure_out] =  hydrostaticForces (obj, t, pos, waveElv)
            % calculates the hydrostatic forces acting on the body
            
            pos = pos - [ obj.cg; 0; 0; 0 ];

            switch obj.freeSurfaceMethodNum

                case 0

                    body_hspressure_out = [];

                    forces = obj.hydroForce.linearHydroRestCoef * pos;

                    % Add Net Bouyancy Force to Z-Direction
                    forces(3) = forces(3) + ...
                        ((obj.simu.g .* obj.hydroForce.storage.mass) - (obj.simu.rho .* obj.simu.g .*  obj.dispVol));

                case 1

                    [forces, body_hspressure_out]  = nonLinearBuoyancy( obj, pos, waveElv, t );

                    % Add Net Bouyancy Force to Z-Direction
                    forces = -forces + [ 0, 0, (obj.simu.g .* obj.hydroForce.storage.mass), 0, 0, 0 ];

            end

%             forces = forces;

        end
        
        function f = waveElevation(obj, pos, t)
            % Function to calculate the wave elevation at the ceontroids of triangulated surface
            % NOTE: This function assumes that the STL file is imported with its CG at 0,0,0

            if obj.freeSurfaceMethodNum == 1
                
                if isempty(obj.oldElev)
                    f = calc_elev(pos,t);
                    obj.oldElev = f;
                else
                    if mod(t, obj.simu.dtFeNonlin) < obj.simu.dt/2
                        f = calc_elev (pos,t);
                        obj.oldElev = f;
                    else
                        f = obj.oldElev;
                    end
                end
                
            else
                f = 0;
            end
        end
        
    end

    % non-public transient simulation methods
    methods (Access = 'protected')

        function ramped = applyRamp (obj, t, nominal)
            % apply a time based ramp to the input
            %
            %

            if t < obj.simu.rampT
                % (3 * pi/2) == 4.712388980384690
                ramped  = nominal .* 0.5 .* (1 + sin( pi .* (t ./ obj.simu.rampT) + 4.712388980384690));
            else
                ramped = nominal;
            end

        end

        function [f, wp, wpMeanFS]  = nonFKForce (obj, t, pos, elv)
            % Function to calculate wave exitation force and moment on a
            % triangulated surface
            %
            % NOTE: This function assumes that the STL file is imported
            % with its CG at 0,0,0

            % Logic to calculate nonFKForce at reduced sample time
            if isempty(obj.oldForce) || (mod (t, obj.simu.dtFeNonlin) < obj.simu.dt/2)

                % TODO: reduce calc_nonFKForce method inputs
                % no need to have all these inputs as they could be
                % accessed directly from the class properties in the
                % subfunction
                [f, wp, wpMeanFS] = calc_nonFKForce ( obj, pos, elv, t );
                obj.oldForce = f;
                obj.oldWp = wp;
                obj.oldWpMeanFs = wpMeanFS;

            else

                f = obj.oldForce;
                wp = obj.oldWp;
                wpMeanFS = obj.oldWpMeanFs;

            end

        end

        function [f, wp, wpMeanFS]  = calc_nonFKForce (obj, pos, elv, t)
            % Function to apply translation and rotation, and calculate
            % forces.

            % Compute new tri coords after cog rotation and translation
            centerMeanFS = hydrobody.offsetXYZ (obj.bodyGeometry.center, obj.hydroData.properties.cg);
            avMeanFS     = obj.bodyGeometry.tnorm .* [obj.bodyGeometry.area, obj.bodyGeometry.area, obj.bodyGeometry.area];

            % Compute new tri coords after cog rotation and translation
            center = hydrobody.rotateXYZ (obj.bodyGeometry.center, [1 0 0], pos(4));
            center = hydrobody.rotateXYZ (center, [0 1 0], pos(5));
            center = hydrobody.rotateXYZ (center, [0 0 1], pos(6));
            center = hydrobody.offsetXYZ (center, pos);
            center = hydrobody.offsetXYZ (center, obj.hydroData.properties.cg);

            % Compute new normal vectors coords after cog rotation and translation
            tnorm = hydrobody.rotateXYZ (obj.bodyGeometry.tnorm, [1 0 0], pos(4));
            tnorm = hydrobody.rotateXYZ (tnorm, [0 1 0], pos(5));
            tnorm = hydrobody.rotateXYZ (tnorm, [0 0 1], pos(6));

            % Compute area vectors
            av = tnorm .* [obj.bodyGeometry.area, obj.bodyGeometry.area, obj.bodyGeometry.area];

            % Calculate the free surface
            wpMeanFS = pDis (obj, centerMeanFS, 0, t);

            wp = pDis (obj, center, elv, t);

            % Calculate forces
            f_linear    = FK ( centerMeanFS,          obj.hydroData.properties.cg, avMeanFS, wpMeanFS );
            f_nonLinear = FK ( center      , pos(1:3) + obj.hydroData.properties.cg,       av,       wp );
            f = f_nonLinear - f_linear;

        end

        function f = pDis (obj, center, elv, t)
            % Function to calculate pressure distribution

            f = zeros (length (center(:,3)), 1);
            z = zeros (length (center(:,1)), 1);

            if obj.waves.typeNum < 10

            elseif obj.waves.typeNum < 20

                f = obj.simu.rho .* obj.simu.g .* obj.waves.AH(1) .* cos(obj.waves.k(1) .* center(:,1) - obj.waves.w(1) * t);

                if obj.waves.deepWaterWave == 0

                    z = (center(:,3) - elv) .* obj.waves.wDepth ./ (obj.waves.wDepth + elv);

                    f = f .* (cosh(obj.waves.k(1) .* (z + obj.waves.wDepth)) ./ cosh(obj.waves.k(1) * obj.waves.wDepth));

                else

                    z = (center(:,3) - elv);

                    f = f .* exp(obj.waves.k(1) .* z);

                end

            elseif obj.waves.typeNum < 30

                for i=1:length(obj.waves.AH)

                    if obj.waves.deepWaterWave == 0 && obj.waves.wDepth <= 0.5*pi/obj.waves.k(i)

                        z = (center(:,3) - elv) .* obj.waves.wDepth ./ (obj.waves.wDepth + elv);

                        f_tmp = obj.simu.rho .* obj.simu.g .* sqrt(obj.waves.AH(i) * obj.waves.dw) .* cos(obj.waves.k(i) .* center(:,1) - obj.waves.w(i) * t - obj.waves.phaseRand(i));

                        f = f + f_tmp .* (cosh(obj.waves.k(i) .* (z + obj.waves.wDepth)) ./ cosh(obj.waves.k(i) .* obj.waves.wDepth));

                    else

                        z = (center(:,3) - elv);

                        f_tmp = obj.simu.rho .* obj.simu.g .* sqrt (obj.waves.AH(i) * obj.waves.dw) .* cos (obj.waves.k(i) .* center(:,1) - obj.waves.w(i) * t - obj.waves.phaseRand(i));

                        f = f + f_tmp .* exp(obj.waves.k(i) .* z);

                    end
                end

            end

            f(z > 0) = 0;

        end

        function f = FK(center, instcg, av, wp)
            % Function to calculate the force and moment about the cog due
            % to Froude-Krylov pressure

            f = zeros(6,1);

            % Calculate the hydrostatic pressure at each triangle center
            pressureVect = [wp, wp, wp] .* -av;

            % Compute force about cog
            f(1:3) = sum (pressureVect);

            % Compute moment about cog
            tmp1 = ones (length (center(:,1)), 1);
            tmp2 = tmp1 * instcg';
            center2cgVec = center - tmp2;

            f(4:6)= sum (cross (center2cgVec, pressureVect));
        end

        function F_FM = convolutionIntegral(obj, vel, t)
            % Function to calculate convolution integral

            if abs(t - obj.radForceOldTime - obj.CIdt) < 1e-8

                obj.radForceVelocity      = circshift(obj.radForceVelocity, 1, 2);

                obj.radForceVelocity(:,1) = vel(:);

                % integrate
                time_series = bsxfun(@times, obj.radForce_IRKB_interp, obj.radForceVelocity);

                F_FM = squeeze(trapz(obj.simu.CTTime, sum(time_series, 1)));

                obj.radForceOldF_FM = F_FM;

                obj.radForceOldTime = t;
            else
                % use the old value which at the start of a simulation
                % is always zeros
                F_FM = obj.radForceOldF_FM;
            end

        end

        function [f,p]  = nonLinearBuoyancy (obj, pos, elv, t)
            % Function to calculate buoyancy force and moment on a
            % triangulated surface NOTE: This function assumes that the STL
            % file is imported with its CG at 0,0,0

            if isempty(obj.oldNonLinBuoyancyF)

                [f,p] = calc_nonLinearBuoyancy (obj, pos, elv);

                obj.oldNonLinBuoyancyF  = f;

                obj.oldNonLinBuoyancyP  = p;

            else

                if mod(t, obj.simu.dtFeNonlin) < obj.simu.dt/2

                    [f,p] = calc_nonLinearBuoyancy (obj, pos,elv);

                    obj.oldNonLinBuoyancyF  = f;

                    obj.oldNonLinBuoyancyP  = p;

                else

                    f  = obj.oldNonLinBuoyancyF;

                    p  = obj.oldNonLinBuoyancyP;

                end
            end
        end

        function [f,p] = calc_nonLinearBuoyancy (obj, pos, elv)
            % Function to apply translation and rotation and calculate forces

            % Compute new tri coords after cog rotation and translation
            center = hydrobody.rotateXYZ(obj.bodyGeometry.center, [1 0 0], pos(4));
            center = hydrobody.rotateXYZ(center, [0 1 0], pos(5));
            center = hydrobody.rotateXYZ(center, [0 0 1], pos(6));
            center = hydrobody.offsetXYZ(center, pos);
            center = hydrobody.offsetXYZ(center, obj.cg);

            % Compute new normal vectors coords after cog rotation
            tnorm = hydrobody.rotateXYZ(obj.bodyGeometry.norm, [1 0 0], pos(4));
            tnorm = hydrobody.rotateXYZ(tnorm, [0 1 0], pos(5));
            tnorm = hydrobody.rotateXYZ(tnorm, [0 0 1], pos(6));

            % Calculate the hydrostatic forces
            av = tnorm .* [obj.bodyGeometry.area, obj.bodyGeometry.area, obj.bodyGeometry.area];

            [f,p] = fHydrostatic (obj, center, elv, pos(1:3) + obj.cg, av);

        end

        function [f,p] = fHydrostatic(obj, center, elv, instcg, av)
            % Function to calculate the force and moment about the cog due
            % to hydrostatic pressure
            f = zeros(6,1);

            % Zero out regions above the mean free surface
            z = center(:,3); z((z-elv)>0)=0;

            % Calculate the hydrostatic pressure at each triangle center
            pressureVect = obj.simu.rho * obj.simu.g .* [-z -z -z] .* -av;
            p = obj.simu.rho * obj.simu.g .* -z;

            % Compute force about cog
            f(1:3) = sum(pressureVect);

            tmp1 = ones(length(center(:,1)),1);
            tmp2 = tmp1*instcg';
            center2cgVec = center - tmp2;

            f(4:6)= sum(cross(center2cgVec,pressureVect));
        end
        
        function f = calc_elev(obj, pos, t)
            % Function to rotate and translate body and call wave elevation function at new locations 

            % Compute new tri center coords after cog rotation and translation
            center = obj.rotateXYZ (obj.bodyGeometry.center, [1, 0, 0], pos(4));
            center = obj.rotateXYZ (center, [0, 1, 0], pos(5));
            center = obj.rotateXYZ (center, [0, 0, 1], pos(6));
            center = obj.offsetXYZ (center, pos);
            center = obj.offsetXYZ (center, obj.cg);
            
            % Calculate the free surface
            f = waveElev (center,t);
        end
        
        function f = waveElev (obj, center, t)
            % Function to calculate the wave elevation at an array of points
            
            f = zeros(length(center(:,3)),1);
            
            cx = center(:,1);
            cy = center(:,2);
            
            X = cx * cos (obj.waves.waveDir * pi/180) ...
                + cy * sin (obj.waves.waveDir * pi/180);
            
            if obj.waves.typeNum <10
                
            elseif obj.waves.typeNum <20
                
                f = obj.waves.AH(1) .* cos(obj.waves.k(1) .* X - obj.waves.w(1) * t);
                
            elseif obj.waves.typeNum <30
                
                tmp = sqrt (obj.waves.AH .* obj.waves.dw);
                
                tmp1 = ones (1, length (center(:,1)));
                
                tmp2 = (obj.waves.w .* t + obj.waves.phaseRand) * tmp1;
                
                tmp3 = cos (obj.waves.k * X'- tmp2);
                
                f(:,1) = tmp3' * tmp;
                
            end
            
            % apply ramp if we are not past the initial ramp time
            if t <= obj.simu.rampT
                rampF = (1 + cos (pi + pi * t / rampT)) / 2;
                f = f .* rampF;
            end
        end

    end

	% public post-processing related methods
    methods (Access = 'public') %modify object = F; output = T

        function fam = forceAddedMass(obj,acc,B2B)
            % Recomputes the real added mass force time history for the
            % body
            %
            % Syntax
            %
            % fam = forceAddedMass(hb,acc,B2B)
            %
            % Input
            %
            %  hb - hydrobody object
            %
            %  acc - (n x 6) time history of body accelerations for which
            %   the added mass is to be recalculated
            %
            %  B2B - flag indicating whether body-to-body interactions are
            %    present
            %
            % Output
            %
            %  fam - added mass recalculated from time history of
            %    accelerations and body added mass.
            %
            %
            
            iBod = obj.bodyNumber;
            fam = zeros(size(acc));
            for i =1:6
                tmp = zeros(length(acc(:,i)),1);
                for j = 1:6
                    if B2B == 1
                        jj = (iBod-1)*6+j;
                    else
                        jj = j;
                    end
                    iam = obj.hydroForce.fAddedMass(i,jj);
                    tmp = tmp + acc(:,j) .* iam;
                end
                fam(:,i) = tmp;
            end
        end

        function write_paraview_vtp (obj, t, pos_all, bodyname, model, simdate, hspressure, wavenonlinearpressure, wavelinearpressure)
            % Writes vtp files for visualization with ParaView
            numVertex = obj.bodyGeometry.numVertex;
            numFace = obj.bodyGeometry.numFace;
            vertex = obj.bodyGeometry.vertex;
            face = obj.bodyGeometry.face;
            cellareas = obj.bodyGeometry.area;
            for it = 1:length(t)
                % calculate new position
                pos = pos_all(it,:);
                vertex_mod = hydrobody.rotateXYZ(vertex,[1 0 0],pos(4));
                vertex_mod = hydrobody.rotateXYZ(vertex_mod,[0 1 0],pos(5));
                vertex_mod = hydrobody.rotateXYZ(vertex_mod,[0 0 1],pos(6));
                vertex_mod = hydrobody.offsetXYZ(vertex_mod,pos(1:3));
                % open file
                filename = ['vtk' filesep 'body' num2str(obj.bodyNumber) '_' bodyname filesep bodyname '_' num2str(it) '.vtp'];
                fid = fopen(filename, 'w');
                % write header
                fprintf(fid, '<?xml version="1.0"?>\n');
                fprintf(fid, ['<!-- WEC-Sim Visualization using ParaView -->\n']);
                fprintf(fid, ['<!--   model: ' model ' - ran on ' simdate ' -->\n']);
                fprintf(fid, ['<!--   body:  ' bodyname ' -->\n']);
                fprintf(fid, ['<!--   time:  ' num2str(t(it)) ' -->\n']);
                fprintf(fid, '<VTKFile type="PolyData" version="0.1">\n');
                fprintf(fid, '  <PolyData>\n');
                % write body info
                fprintf(fid,['    <Piece NumberOfPoints="' num2str(numVertex) '" NumberOfPolys="' num2str(numFace) '">\n']);
                % write points
                fprintf(fid,'      <Points>\n');
                fprintf(fid,'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
                for ii = 1:numVertex
                    fprintf(fid, '          %5.5f %5.5f %5.5f\n', vertex_mod(ii,:));
                end;
                clear vertex_mod
                fprintf(fid,'        </DataArray>\n');
                fprintf(fid,'      </Points>\n');
                % write tirangles connectivity
                fprintf(fid,'      <Polys>\n');
                fprintf(fid,'        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
                for ii = 1:numFace
                    fprintf(fid, '          %i %i %i\n', face(ii,:)-1);
                end;
                fprintf(fid,'        </DataArray>\n');
                fprintf(fid,'        <DataArray type="Int32" Name="offsets" format="ascii">\n');
                fprintf(fid, '         ');
                for ii = 1:numFace
                    n = ii * 3;
                    fprintf(fid, ' %i', n);
                end;
                fprintf(fid, '\n');
                fprintf(fid,'        </DataArray>\n');
                fprintf(fid, '      </Polys>\n');
                % write cell data
                fprintf(fid,'      <CellData>\n');
                % Cell Areas
                fprintf(fid,'        <DataArray type="Float32" Name="Cell Area" NumberOfComponents="1" format="ascii">\n');
                for ii = 1:numFace
                    fprintf(fid, '          %i', cellareas(ii));
                end;
                fprintf(fid, '\n');
                fprintf(fid,'        </DataArray>\n');
                % Hydrostatic Pressure
                if ~isempty(hspressure)
                    fprintf(fid,'        <DataArray type="Float32" Name="Hydrostatic Pressure" NumberOfComponents="1" format="ascii">\n');
                    for ii = 1:numFace
                        fprintf(fid, '          %i', hspressure.signals.values(it,ii));
                    end;
                    fprintf(fid, '\n');
                    fprintf(fid,'        </DataArray>\n');
                end
                % Non-Linear Froude-Krylov Wave Pressure
                if ~isempty(wavenonlinearpressure)
                    fprintf(fid,'        <DataArray type="Float32" Name="Wave Pressure NonLinear" NumberOfComponents="1" format="ascii">\n');
                    for ii = 1:numFace
                        fprintf(fid, '          %i', wavenonlinearpressure.signals.values(it,ii));
                    end;
                    fprintf(fid, '\n');
                    fprintf(fid,'        </DataArray>\n');
                end
                % Linear Froude-Krylov Wave Pressure
                if ~isempty(wavelinearpressure)
                    fprintf(fid,'        <DataArray type="Float32" Name="Wave Pressure Linear" NumberOfComponents="1" format="ascii">\n');
                    for ii = 1:numFace
                        fprintf(fid, '          %i', wavelinearpressure.signals.values(it,ii));
                    end;
                    fprintf(fid, '\n');
                    fprintf(fid,'        </DataArray>\n');
                end
                fprintf(fid,'      </CellData>\n');
                % end file
                fprintf(fid, '    </Piece>\n');
                fprintf(fid, '  </PolyData>\n');
                fprintf(fid, '</VTKFile>');
                % close file
                fclose(fid);
            end
        end

    end

    % Static methods
    methods (Static)

        function xn = rotateXYZ (x, ax, theta)
            % Function to rotate a point about an arbitrary axis
            % x: 3-componenet coordinates
            % ax: axis about which to rotate (must be a normal vector)
            % theta: rotation angle
            % xn: new coordinates after rotation
            rotMat = zeros(3);
            rotMat(1,1) = ax(1)*ax(1)*(1-cos(theta))    + cos(theta);
            rotMat(1,2) = ax(2)*ax(1)*(1-cos(theta))    + ax(3)*sin(theta);
            rotMat(1,3) = ax(3)*ax(1)*(1-cos(theta))    - ax(2)*sin(theta);
            rotMat(2,1) = ax(1)*ax(2)*(1-cos(theta))    - ax(3)*sin(theta);
            rotMat(2,2) = ax(2)*ax(2)*(1-cos(theta))    + cos(theta);
            rotMat(2,3) = ax(3)*ax(2)*(1-cos(theta))    + ax(1)*sin(theta);
            rotMat(3,1) = ax(1)*ax(3)*(1-cos(theta))    + ax(2)*sin(theta);
            rotMat(3,2) = ax(2)*ax(3)*(1-cos(theta))    - ax(1)*sin(theta);
            rotMat(3,3) = ax(3)*ax(3)*(1-cos(theta))    + cos(theta);
            xn = x * rotMat;
        end

        function verts_out = offsetXYZ (verts, x)
            % Function to move the position vertices
            verts_out(:,1) = verts(:,1) + x(1);
            verts_out(:,2) = verts(:,2) + x(2);
            verts_out(:,3) = verts(:,3) + x(3);
        end
        
        function v = lagrangeInterp (x,y,u)
            
            n = length(x);
            v = zeros(size(u));
            
            for k = 1:n
                w = ones(size(u));
                for j = [1:k-1 k+1:n]
                    w = (u-x(j))./(x(k)-x(j)).*w;
                end
                v = v + w*y(k);
            end
            
        end
        
        function v = linearInterp (x1, x2, y1, y2, u)
            
            m = (y2 - y1) ./ (x2 - x1);
            c = y1 - m.*x1;
            
            v = m.*u + c;
            
        end

    end

end
