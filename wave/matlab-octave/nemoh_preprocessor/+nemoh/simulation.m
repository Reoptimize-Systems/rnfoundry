classdef simulation < nemoh.base
    % class for running NEMOH hydrodynamic BEM solver simulations
    %
    % Syntax
    %
    % ns = nemoh.simulation (inputdir, installdir)
    % ns = nemoh.simulation (..., 'Parameter', value)
    %
    % Description
    %
    % nemoh.simulation is a class used to generate input files for the
    % NEMOH hydrodynamic BEM solver programs and run the solver on them.
    %
    % simulation Methods:
    %
    %   simulation - constructor
    %   addBody - add a nemoh body to the simulation
    %   writeMesherInputFiles - writes the input files for the NEMOH mesh
    %     program (.dat files)
    %   processMeshes - performs pre-simulation mesh processing for all
    %     bodies (generates KH.dat, Hydrostatics.dat etc.)
    %   writeNemoh - writes the NEMOH input file (Nemoh.cal)
    %   run - runs the NEMOH programs on the generated input files
    %
    %
    % See also: example_nemoh_cylinder.m
    %
    
    properties (GetAccess = public, SetAccess = private)
        
        meshProgPath; % nemoh mesh program file path
        preProcProgPath; % nemoh preProc program file path
        solverProgPath; % nemoh solver program file path
        postProcProgPath; % nemoh postProc program file path
        
        inputDataDirectory; % directory where Nemoh input files will be created
        meshDirectory; % directory where some mesh related files for NEMO are found

        rho; % density of fluid (kg/m^3)
        g; % accleration due to gravity (m/s^2)
        
        bodies;
        simReady;
        
    end
    
    
    methods
        
        function self = simulation (inputdir, varargin)
            % constructor for the nemoh.simulation class
            %
            % Syntax
            %
            % ns = nemoh.simulation (inputdir, installdir)
            % ns = nemoh.simulation (..., 'Parameter', value)
            %
            % Input
            %
            %  inputdir - directory where Nemoh input files will be created
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            % 'InstallDir' - directory where the Nemoh executables can be
            %   found. If not supplied, the executables are assumed to be
            %   on the computer's path so they can be invoked without the
            %   full file path by just using their name.
            %
            % 'rho' - Density of water in the simulation (kg/m^3), Default
            %   is 1025 if not supplied.
            %
            % 'g' - Graviational acceleration used in the simulation
            %   (m/s^2). 9.81 is used if not supplied.
            %
            % 'Bodies' - an array of one or more nemoh.body objects to be
            %   added to the simulation. Bodies can also be added later
            %   using the addBody method.
            %
            % Output
            %
            %  ns - a nemoh.system object
            %
            
            options.rho = 1025;
            options.g = 9.81;
            options.Bodies = [];
            options.InstallDir = '';
            
            options = parse_pv_pairs (options, varargin);
            
            if ispc
                prognames = {'Mesh', 'preProcessor', 'Solver', 'postProcessor'};
                meshdir = 'Mesh';
            else
                prognames = {'mesh', 'preProc', 'solver', 'postProc'};
                meshdir = 'mesh';
            end
            
            if ~isempty (options.InstallDir)
                % check files actually exist
                for ind = 1:numel (prognames)
                    if ~exist (fullfile (options.InstallDir, prognames{ind}), 'file')
                        error ('%s binary was not found in provided install directory:\n%s', ...
                            prognames{ind}, options.InstallDir);
                    end
                end
            end
            
            self.checkNumericScalar (options.rho, true, 'rho');
            self.checkNumericScalar (options.g, true, 'g');
            
            if ~exist (inputdir, 'dir')
                
                [status, msg] = mkdir (inputdir);
                
                if status == false
                    error ('inputdir:\n%s\ndoes not exist, and an attempt to create the directory failed with the following error:\n%s', ...
                        msg);
                end
                
            end
            
            self.meshProgPath = fullfile (options.InstallDir, prognames{1});
            self.preProcProgPath = fullfile (options.InstallDir, prognames{2});
            self.solverProgPath = fullfile (options.InstallDir, prognames{3});
            self.postProcProgPath = fullfile (options.InstallDir, prognames{4});
            self.inputDataDirectory = inputdir;
            self.meshDirectory = fullfile (self.inputDataDirectory, meshdir);
            self.rho = options.rho;
            self.g = options.g;
            
            if ~isempty (options.Bodies)
                for ind = 1:numel (options.Bodies)
                    self.addBody (options.Bodies(ind));
                end
            end
            
            % mark sim as not ready to run
            self.simReady = false;
            
        end
        
        function run (self, varargin)
            % run Nemoh on the simulation input files
            %
            % Syntax
            %
            % run (ns)
            % run (ns, 'Parameter', value)
            %
            % Description
            %
            % run invokes the Nemo programs on the generated input files.
            % Before running it also first generates the ID.dat file and
            % input.txt file for the problem. See the Nemoh documentation
            % for more information about these files. After generating the
            % files, the preprocessor, solver and postprocessor are all run
            % on the input files in sequence.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            % 'Solver' - integer or string. Used to select the solver for
            %   the problem, the following values are possible:
            %
            %   0 or 'direct gauss' or 'directgauss'
            %   1 or 'GMRES'
            %   2 or 'GMRES-FHM' or 'GMRES with FHM'
            %
            %   The string forms are case-insensitive (e.g. 'gmres' works
            %   just as well as 'GMRES'.
            %
            % 'SavePotential' - logical (true/false) flag indicating whther
            %   the potential should be saved for visualisation.
            %
            % 'GMRESRestart' - scalar integer. Restart criteria for the
            %   gmres solver types. Ignored for other solvers. Default is
            %   20.
            %
            % 'GMRESTol' - scalar value. Tolerance for the gmres solver
            %   types. Ignored for other solvers. Default is 5e-7.
            %
            % 'GMRESMaxIterations' - scalar integer. Maximum number of
            %   iterations allowed for the gmres solver types. Ignored for
            %   other solvers. Default is 100.
            %
            % 'Verbose' - logical (true/false) flag determning whether to
            %   display some text on the command line informing on the
            %   progress of the simulation.
            %
            %
            % See also: nemoh.simulation.writeNemoh
            %
            
            options.Verbose = false;
            options.Solver = 0;
            options.GMRESRestart = 20;
            options.GMRESTol = 5e-7;
            options.GMRESMaxIterations = 100;
            options.SavePotential = true;
            
            options = parse_pv_pairs (options, varargin);

            self.checkLogicalScalar (options.Verbose, true, 'Verbose');
            self.checkLogicalScalar (options.SavePotential, true, 'SavePotential');
            
            if self.simReady
                
                % make the results directory if it does not exist
                mkdir ( fullfile (self.inputDataDirectory, 'results'));
                
                CC = onCleanup (@() cd (pwd));
                cd (self.inputDataDirectory);
                
                % write ID file
                fid = fopen ('ID.dat', 'w');
                CC = onCleanup (@() fclose (fid));
                
                fprintf (fid, '1\n.');
                
                clear CC;
                
                % write input.txt
                fid = fopen ('input.txt', 'w');
                CC = onCleanup (@() fclose (fid));
                
                fprintf (fid, '--- Calculation parameters ------------------------------------------------------------------------------------\n');
                
                if ischar (options.Solver)
                    % change to lower case
                    options.Solver = lower (options.Solver);
                elseif self.checkScalarInteger (options.Solver, false) ~= true
                    error ('Solver must be a string or scalar integer');
                end
                
                switch options.Solver
                    
                    case {0, 'direct gauss', 'directgauss'}
                        
                        fprintf (fid, '0\t! direct gauss solver\n');
                        
                    case {1, 'GMRES', 'gmres'}
                        
                        fprintf (fid, '1\t! GMRES solver\n');
                        
                        self.checkNumericScalar (options.GMRESRestart, true, 'GMRESRestart');
                        self.checkNumericScalar (options.GMRESTol, true, 'GMRESTol');
                        self.checkNumericScalar (options.GMRESMaxIterations, true, 'GMRESMaxIterations');
                        
                        fprintf (fid, '%d\t! IRES, Restart parameter for GMRES\n', options.GMRESRestart);
                        fprintf (fid, '%s\t! TOL_GMRES, Stopping criterion for GMRES\n', self.formatNumber (options.GMRESTol));
                        fprintf (fid, '%d\t! MAXIT, Maximum iterations for GMRES\n', options.GMRESMaxIterations);
                    
                    case {2, 'gmres-fhm', 'gmres with fhm'}
                        
                        fprintf (fid, '2    ! GMRES with FMM acceleration\n');
                        
                        self.checkNumericScalar (options.GMRESRestart, true, 'GMRESRestart');
                        self.checkNumericScalar (options.GMRESTol, true, 'GMRESTol');
                        self.checkNumericScalar (options.GMRESMaxIterations, true, 'GMRESMaxIterations');
                        
                        fprintf (fid, '%d\t! IRES, Restart parameter for GMRES\n', options.GMRESRestart);
                        fprintf (fid, '%s\t! TOL_GMRES, Stopping criterion for GMRES\n', self.formatNumber (options.GMRESTol));
                        fprintf (fid, '%d\t! MAXIT, Maximum iterations for GMRES\n', options.GMRESMaxIterations);
                        
                    otherwise
                        error ('Unrecognised Solver option');
                        
                end
                
                fprintf (fid, '%s\t! Sav_potential, Save potential for visualization \n', ...
                    self.formatInteger ( double (options.SavePotential)));
                
                clear CC;
                
                if options.Verbose, fprintf('------ Starting NEMOH pre-processor----------- \n'); end
                
                status = system (sprintf ('"%s"', self.preProcProgPath));
                
                if status ~= 0
                    error ('pre-processing failed with error code %d', status);
                end
                
                if options.Verbose, fprintf('------ Solving BVPs ------------- \n'); end
                
                status = system (sprintf ('"%s"', self.solverProgPath));
                
                if status ~= 0
                    error ('solving failed with error code %d', status);
                end
                
                if options.Verbose, fprintf('------ Postprocessing results --- \n'); end
                
                status = system (sprintf ('"%s"', self.postProcProgPath));
                
                if status ~= 0
                    error ('post-processing failed with error code %d', status);
                end
                
                if options.Verbose, fprintf('------ Analysis complete --- \n'); end
                
            else
                
                error ('Nemoh simulation is not ready to be run.');
                
            end
            
        end
        
        function addBody (self, body)
            % add one or more bodies to the simulation
            %
            % Syntax
            %
            % addBody (ns, body)
            %
            % Description
            %
            % Adds one or more nemoh bodies to the simulation. The bodies'
            % data are used to generate the bodies section of the Nemoh.cal
            % input file. When added, the id of the body is set to a new
            % value, the values for rho and g are set to the simulation
            % values and the MeshProgPath value is set to the same value as
            % the simulation property. Adding a body means that the
            % simulation is reset, so wrteMeshes, processMeshes and
            % writeNemoh must all be re-run before the simulation can be
            % run again (i.e. by using the 'run' method).
            %
            % Input
            %
            %  ns - nemoh.simulation object.
            %
            %  body - nemoh.body object (or array of nemoh.body objects)
            %
            
            if ~isa (body, 'nemoh.body')
                error ('body must be a nemoh.body object');
            end
            
            for ind = 1:numel (body)
                
                % work around for the fact that 
                %   self.bodies = [self.bodies, body];
                % won't work in octave due to bug 
                n = numel (self.bodies);
                if n == 0
                    self.bodies = body;
                else
                    self.bodies(end+1) = body;
                end

                % set ids starting from zero (we got 'n' before adding a
                % body)
                self.bodies(end).setID (n);

                self.bodies(end).setMeshProgPath (self.meshProgPath);

                self.bodies(end).setRho (self.rho);

                self.bodies(end).setg (self.g);
            
            end
            
            self.simReady = false;
            
        end
        
        function [hax, hfig] = drawMesh (self, varargin)
            
            options.PlotForces = true;
            options.Axes = [];
            options.AddTitle = true;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.PlotForces, true, 'PlotForces');
            self.checkLogicalScalar (options.AddTitle, true, 'AddTitle');
            
            if isempty (options.Axes)
                hfig = figure ();
                hax = axes (hfig);
                view (hax, 3);
                axis equal;
            else
                hax = options.Axes;
                hfig = get (hax, 'parent');
            end
            
            hold on;
            CC = onCleanup (@() hold ('off'));
            
            colormap (hax, summer (numel (self.bodies)));
            
            ColOrd = get (hax,'ColorOrder');
            % Determine the number of colors in
            % the matrix
            [m,n] = size(ColOrd);

            for ind = 1:numel (self.bodies)
                
                % Determine which row to use in the
                % Color Order Matrix
                ColRow = rem(ind,m);
                if ColRow == 0
                  ColRow = m;
                end
                % Get the color
                Col = ColOrd(ColRow,:);
    
                self.bodies(ind).drawMesh ( 'Axes', hax, ...
                                            'PlotForces', options.PlotForces, ...
                                            'AddTitle', false, ...
                                            'EdgeColor', Col );
            end
            
            clear CC
            
            if options.AddTitle
                title ('NEMOH Simulation, All Bodies Mesh');
            end
            
        end
        
        function writeMesherInputFiles (self, varargin)
            % writes nemoh mesh input files
            %
            % Syntax
            %
            % writeMeshes (ns)
            %
            % Description
            %
            % writeMeshes calls the writeMesh method for each body in the
            % simulation. This method writes mesh related files to disk.
            % What files are written depends on the type of the mesh for
            % that body. e.g. for a axisymmetric mesh specified using a 2D
            % profile it will write the course mesh files to disk ready for
            % processing by the Nemoh meshing program.
            %
            % After running writeMeshes, processMeshes should be called to
            % do any post-processing. e.g. for the example just given,
            % processMeshes runs the Nemoh meshing program on the course
            % input mesh files for the body to generate the refined mesh
            % file in the format Nemoh can read and loads some initial
            % calulations performed by the meshing function.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            % Output
            %
            %  none
            %
            % See Also: nemoh.simulation.processMeshes
            %
            
%             options.TargetNumberOfPanels = [];
%             
%             options = parse_pv_pairs (options, varargin);
%             
%             if ~isempty (options.TargetNumberOfPanels)
%                 assert ( isnumeric (options.TargetNumberOfPanels) ...
%                            && numel(options.TargetNumberOfPanels) == numel (self.bodies), ...
%                          'TargetNumberOfPanels must be a vector of integers the same length as the number of bodies in the simulation' ) ;              
%             end
            
            % make the Mesh directory if it does not exist
            mkdir ( self.meshDirectory );
            
            for ind = 1:numel (self.bodies)
                
%                 if ~isempty (options.TargetNumberOfPanels)
%                     self.checkScalarInteger ( options.TargetNumberOfPanels(ind), ...
%                         true, ...
%                         sprintf('TargetNumberOfPanels (for body %d)', ind) );
%                     
%                     targetpanels = options.TargetNumberOfPanels(ind);
%                 else
%                     targetpanels = [];
%                 end
                    
                self.bodies(ind).writeMesh ();
            end
            
        end
        
        function processMeshes (self, varargin)
            % performs post-processing on mesh files in prep for simulation
            %
            % Syntax
            %
            % processMeshes (ns)
            %
            % Description
            %
            % After running writeMeshes, processMeshes should be called to
            % do any post-processing. e.g. for a axisymmetric mesh
            % specified using a 2D profile, processMeshes runs the Nemoh
            % meshing program on the course input mesh files for the body
            % (created by writeMeshes) to generate the refined mesh file in
            % the format Nemoh can read and loads some initial calulations
            % performed by the meshing function.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            
            for ind = 1:numel (self.bodies)
                if self.bodies(ind).meshProcessed == false
                    
                    % process the mesh file for the body
                    self.bodies(ind).processMesh ();
                    
                    hydrostaticsfile = fullfile (self.bodies(ind).meshDirectory, 'Hydrostatics.dat');
                    khfile = fullfile (self.bodies(ind).meshDirectory, 'KH.dat');
                    
                    if numel (self.bodies) == 1
                        % copy the hydrostatics and KH file to top level mesh
                        % directory with same name
                        if exist (hydrostaticsfile, 'file')
                            copyfile ( hydrostaticsfile, ...
                                       fullfile (self.meshDirectory, 'Hydrostatics.dat') );
                        end
                        
                        if exist (hydrostaticsfile, 'file')
                            copyfile ( khfile, ...
                                       fullfile (self.meshDirectory, 'KH.dat') );
                        end
                        
                    else
                        % copy the hydrostatics and KH file to top level mesh
                        % directory with new name
                        if exist (hydrostaticsfile, 'file')
                            copyfile ( hydrostaticsfile, ...
                                       fullfile (self.meshDirectory, sprintf('Hydrostatics_%d.dat', ind-1)) );
                        end
                        
                        if exist (hydrostaticsfile, 'file')
                            copyfile ( khfile, ...
                                       fullfile (self.meshDirectory, sprintf('KH_%d.dat', ind-1)) );
                        end
                        
                    end
            
                end
            end
            
        end
        
        function loadProcessedMeshes (self)
            % loads processed nemoh mesh files
            %
            % Syntax
            %
            % loadProcessedMeshes (ns)
            %
            % Description
            %
            % After running processMeshes, loadProcessedMeshes can be
            % called to load the processed mesh files into each body. This
            % allows the processed mesh files to be viewed using drawMesh.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            % Output
            %
            %  none
            %
            % See Also: 
            %
            
            for ind = 1:numel (self.bodies)
                self.bodies(ind).loadProcessedMesh ();
            end
            
        end
        
        function loadProcessedMeshData (self)
            % loads processed nemoh mesh data such as displacement,
            % buoyancy center, hydrostatic stiffness etc.
            %
            % Syntax
            %
            % loadMeshData (ns)
            %
            % Description
            %
            % After running processMeshes, loadProcessedMeshData can be
            % called to load the processed mesh data into each body. The
            % data loaded is the displacement, buoyancy center, hydrostatic
            % stiffness, an estimate of the masses and the inertia matrix.
            % The results are then accessible via the body properties kH,
            % hydrostaticForceX, hydrostaticForceY, hydrostaticForceZ, xB,
            % yB, zB, mass, WPA and inertia.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            % Output
            %
            %  none
            %
            % See Also: 
            %
            
            for ind = 1:numel (self.bodies)
                self.bodies(ind).loadProcessedMeshData ();
            end
            
        end
        
        function writeNemoh (self, varargin)
            % generates the Nemoh.cal input file for NEMOH
            %
            % Syntax
            %
            % writeNemoh (ns)
            % writeNemoh (ns, 'Parameter', value)
            %
            % Description
            %
            % writeNemoh generates the Nemoh.cal input file for the NEMOH
            % BEM solver programs.
            %
            % Input
            %
            %  ns - nemoh.simulation object
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            %  'DoIRFCalculation' - logical (true/false) flag determining
            %    whether the impulse response function for the radiation
            %    force is calculated. Default is true if not supplied.
            %
            %  'IRFTimeStep' - Scalar time step in seconds to be used for
            %    the IRF calcualtion. Default is 0.1 if not supplied.
            %
            %  'IRFDuration' - Scalar duration in seconds for the IRF
            %    calculation. Default is 10 if not supplied.
            %
            %  'NWaveFreqs' - Scalar integer number of wave frequencies for
            %    which to calculate the hydrodynamic date. Default is 1 if
            %    not supplied.
            %
            %  'MinMaxWaveFreqs' - If NWaveFreqs is set to 1, this is a
            %    scalar value, the wave frequency at which to do the
            %    simulation, in this case a two element vector with both
            %    values the same is also acceptable. If NWaveFreqs > 1,
            %    this must be a two element vector containing the min and
            %    max values of the frequency range to be computed. Default
            %    is 0.8 if not supplied.
            %
            %  'NDirectionAngles' - Number of directions to be used to
            %    calculate the Kochin function. Default is 0 if not
            %    supplied, meaning no Kochin function calculation is
            %    performed.
            %
            %  'MinMaxDirectionAngles' - Two element vector containing the
            %    min and max angles for the Kochin function calculation.
            %    Default is 0 if not supplied.
            %
            %  'FreeSurfacePoints' - Two element vector containing the
            %    number of points at which to calcualte the free surface
            %    elevation. The first element is the number of points in
            %    the x direction (0 for no calculations) and y direction
            %    By default no calculations are performed (the default
            %    value used to set the input in Nemoh.cal is [0,50]).
            %
            %  'DomainSize' - Two element vector containing the domain size
            %    in the x and y direction. Default is [400, 400] if not
            %    supplied.
            %
            %  'WaterDepth' - Scalar value of the water depth, a value of 0
            %    means infinite depth. Default is 0 if not supplied.
            %
            %  'WaveMeasurementPoint' - Two element vector containing the
            %    wave measurement point [XEFF YEFF]. Default is [0, 0] if
            %    not supplied.
            %
            %
            % Output
            %
            %  The Nemoh.cal file is generated in the specified input data
            %  directory.
            %
            % See Also: nemoh.simulation/run
            %
            
            options.DoIRFCalculation = true;
            options.IRFTimeStep = 0.1;
            options.IRFDuration = 10;
            options.NWaveFreqs = 1;
            options.MinMaxWaveFreqs = 0.8;
            options.NDirectionAngles = 1;
            options.MinMaxDirectionAngles = [0, 0];
            options.FreeSurfacePoints = [0, 50];
            options.DomainSize = [400, 400];
            options.WaterDepth = 0;
            options.WaveMeasurementPoint = [0, 0];
            
            options = parse_pv_pairs (options, varargin);
            
            assert (islogical (options.DoIRFCalculation), ...
                'DoIRFCalculation must be a logical value (true or false)');
            
            assert (numel (options.NWaveFreqs) == 1 && self.checkScalarInteger (options.NWaveFreqs, false), ...
                'NWaveFreqs must be a scalar integer');
            
            assert (isscalar (options.IRFTimeStep), ...
                'IRFTimeStep must be a scalar value.' );
            
            assert (isscalar (options.IRFDuration), ...
                'IRFDuration must be a scalar value.' );
            
            if options.NWaveFreqs == 1
                if numel (options.MinMaxWaveFreqs) == 2
                    if options.MinMaxWaveFreqs(1) ~= options.MinMaxWaveFreqs(2)
                        error ('NWaveFreqs is 1 but numel(MinMaxWaveFreqs) > 1 and MinMaxWaveFreqs(1) is not equal to MinMaxWaveFreqs(2)');
                    end
                elseif numel (options.MinMaxWaveFreqs) == 1
                    options.MinMaxWaveFreqs(2) = options.MinMaxWaveFreqs(1);
                else
                    error ('MinMaxWaveFreqs must be a 1 or 2 element vector when NWaveFreqs is 1');
                end
            elseif options.NWaveFreqs < 1
                error ('NWaveFreqs must be greater than or equal to 1');
            else
                if numel (options.MinMaxWaveFreqs) ~= 2
                    error ('If NWaveFreqs is greater than 1, MinMaxWaveFreqs must be a 2 element vector');
                end
            end
            
            assert ( numel (options.NDirectionAngles) == 1 ...
                     && self.checkScalarInteger (options.NDirectionAngles, false), ...
                'NDirectionAngles must be a scalar integer');
            
            if options.NDirectionAngles == 1
                if numel (options.MinMaxDirectionAngles) == 2
                    if options.MinMaxDirectionAngles(1) ~= options.MinMaxDirectionAngles(2)
                        error ('NDirectionAngles is 1 but numel(MinMaxDirectionAngles) > 1 and MinMaxDirectionAngles(1) is not equal to MinMaxDirectionAngles(2)');
                    end
                elseif numel (options.MinMaxDirectionAngles) == 1
                    options.MinMaxDirectionAngles(2) = options.MinMaxDirectionAngles(1);
                else
                    error ('MinMaxDirectionAngles must be a 1 or 2 element vector when NDirectionAngles is 1');
                end
            elseif options.NDirectionAngles < 0
                error ('NDirectionAngles must be greater than or equal to 0');
            else
                if numel (options.MinMaxWaveFreqs) ~= 2
                    error ('If NDirectionAngles is greater than 1, MinMaxDirectionAngles must be a 2 element vector');
                end
            end
            
            assert (numel (options.FreeSurfacePoints) == 2, ...
                'FreeSurfacePoints must be a two element vector.');
            
            assert (self.checkScalarInteger (options.FreeSurfacePoints(1), false) ...
                    && self.checkScalarInteger (options.FreeSurfacePoints(2), false), ...
                'FreeSurfacePoints elements must be integers');
                
            assert (numel (options.DomainSize) == 2, ...
                'DomainSize must be a two element vector.');
            
            assert (isscalar (options.WaterDepth), ...
                'WaterDepth must be a scalar value.' );
            
            assert (numel (options.WaveMeasurementPoint) == 2, ...
                'WaveMeasurementPoint must be a two element vector.');
            
            if options.DoIRFCalculation
                doirf = '1';
            else
                doirf = '0';
            end
            
            % Write Nemoh input file
            fid = fopen (fullfile (self.inputDataDirectory, 'Nemoh.cal'), 'w');
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
            fprintf (fid,'%s\t\t\t\t! RHO \t\t\t! KG/M**3 \t! Fluid specific volume \n', self.formatNumber (self.rho));
            fprintf (fid,'%s\t\t\t\t! G\t\t\t! M/S**2\t! Gravity \n', self.formatNumber (self.g));
            fprintf (fid,'%s                 ! DEPTH\t\t\t! M\t\t! Water depth (m), 0 for deep water.\n', self.formatNumber (options.WaterDepth));
            fprintf (fid,'%s\t%s              ! XEFF YEFF\t\t! M\t\t! Wave measurement point\n', self.formatNumber (options.WaveMeasurementPoint(1)), self.formatNumber (options.WaveMeasurementPoint(2)));
            fprintf (fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%d\t\t\t\t! Number of bodies\n', numel (self.bodies));
            
            writeBodies (self, fid);
            
            fprintf (fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%d\t%s\t%s\t\t! Number of wave frequencies, Min, and Max (rad/s)\n', options.NWaveFreqs, self.formatNumber (options.MinMaxWaveFreqs(1)), self.formatNumber (options.MinMaxWaveFreqs(2)));
            fprintf (fid,'%d\t%s\t%s\t\t! Number of wave directions, Min and Max (degrees)\n', options.NDirectionAngles, self.formatNumber (options.MinMaxDirectionAngles(1)), self.formatNumber (options.MinMaxDirectionAngles(2)));
            fprintf (fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%s\t%s\t%s\t\t! IRF \t\t\t\t! IRF calculation (0 for no calculation), time step and duration\n', doirf, self.formatNumber (options.IRFTimeStep), self.formatNumber (options.IRFDuration));
            fprintf (fid,'0\t\t\t\t! Show pressure\n');
            fprintf (fid,'%d\t%s\t%s\t\t! Kochin function \t\t! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n', ...
                options.NDirectionAngles, self.formatNumber (options.MinMaxDirectionAngles(1)), self.formatNumber (options.MinMaxDirectionAngles(2)));
            fprintf (fid,'%d\t%d\t%s\t%s\t! Free surface elevation \t! Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction\n', ...
                options.FreeSurfacePoints(1), options.FreeSurfacePoints(2), self.formatNumber (options.DomainSize(1)), self.formatNumber (options.DomainSize(2)) );
            fprintf (fid,'---');
            
            self.simReady = true;
            
        end

    end
    
    methods (Access = private)
        
        function writeBodies (self, fid)
            % writes all body sections of Nemoh.cal to the specified file
            %
            
            str = generateBodiesStr(self);
            
            fprintf (fid, '%s', str);
            
        end

        function str = generateBodiesStr (self)
            
            str = '';
            
            for bodyind = 1:numel (self.bodies)
                
                str = [str, sprintf('--- Body %d -----------------------------------------------------------------------------------------------------------------------\n', bodyind)];
            
                str = [str, self.bodies(bodyind).generateBodyStr()];
                
            end
            
        end
        
    end
    
end