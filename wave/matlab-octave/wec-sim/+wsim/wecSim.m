classdef wecSim < handle

    properties
        wavePlotNumPointsX = 50;
        wavePlotNumPointsY = 50;
        wavePlotDomainSize = [];
    end

    properties (GetAccess = public, SetAccess = private)

        % sim control variables
        readyToRun; % true/false flag indicating whether the simulation is ready to be run
        simStarted; % true/false flag indicating whether the simulation has started (i.e. simStart had been called)
        simComplete; % true/false flag indicating whether the simulation is finished (i.e. simFinish has been called)
        showProgress; % true/false flag indicating whether a progress bar will be shown
        
        powerTakeOffs;
        wecController;
        mBDynOutputFile;
        mBDynInputFile;
        caseDirectory;
        
        % simInfo - structure containing information about the simulation
        %  The structure is populated when simStart is called (or the run
        %  method, which ultimately calls simStart). It will contain the
        %  following fields:
        %
        %  TStart : the simulation initial time
        %
        %  TEnd : the simulation end time
        %
        %  TStep : the magnitude of the time step of the simulation
        %
        %  MBDynSystem : the mbdyn.pre.system object representing the
        %    multibody system 
        %
        %  HydroSystem : wsim.hydroSystem object representing the
        %    hydrodynamic interaction system
        %
        %  HydroMotionSyncSteps : 
        %
        %  OutputDirectory : the location of the directory where results
        %    will be output
        %
        %  CaseDirectory : the directory containing the case information
        %    for the simulation
        %
        simInfo;

        loggingSettings;

        % logging variables
        lastTime;                   % last value of simulation time
        lastPositions;              % last calculated cartesian positions
        lastAngularPositions;       % last calculated angular positions (euler angles)
        lastVelocities;             % last calculated velocites
        lastAngularVelocities;      % last calculated angular velocities
        lastAccelerations;          % last calculated accelerations
        lastAngularAccelerations;   % last calculated angular accelerations
        lastForceHydro;             % last calculated total of all hydrodynamic forces
        lastForceExcitation;        % last calculated value of all hydrodynamic excitation forces
        lastForceExcitationRamp;    % last calculated value of all hydrodynamic excitation forces with applied ramp
        lastForceExcitationLin;     % last calculated value of hydrodynamic linear excitation forces
        lastForceExcitationNonLin;  % last calculated value of hydrodynamic nonlinear excitation forces
        lastForceRadiationDamping;  % last calculated value of hydrodynamic radiation damping forces
        lastForceRestoring;         % last calculated value of hydrodynamic restoring (buoyancy) forces
        lastForceMorrison;          % last calculated value of hydrodynamic morrison forces
        lastForceViscousDamping;    % last calculated value of hydrodynamic viscous damping forces
        lastForceAddedMassUncorrected; % last calculated value of hydrodynamic forces due to added mass (uncorrected)
        lastNodeForcesUncorrected; % last calculated value ofcombined forces but with uncorrected added mass forces
        lastMomentHydro;             % last calculated total of all hydrodynamic moments
        lastMomentExcitation;        % last calculated value of all hydrodynamic excitation moments
        lastMomentExcitationRamp;    % last calculated value of all hydrodynamic excitation moments with applied ramp
        lastMomentExcitationLin;     % last calculated value of hydrodynamic linear excitation moments
        lastMomentExcitationNonLin;  % last calculated value of hydrodynamic nonlinear excitation moments
        lastMomentRadiationDamping;  % last calculated value of hydrodynamic radiation damping moments
        lastMomentRestoring;         % last calculated value of hydrodynamic restoring (buoyancy) moments
        lastMomentMorrison;          % last calculated value of hydrodynamic morrison moments
        lastMomentViscousDamping;    % last calculated value of hydrodynamic viscous damping moments
        lastMomentAddedMassUncorrected; % last calculated value of hydrodynamic moments due to added mass (uncorrected)
        lastNodeMomentsUncorrected; % last calculated value ofcombined moments but with uncorrected added mass forces
        lastForceAddedMass; % last calculated value of hydrodynamic forces due to added mass (corrected)
        lastMomentAddedMass; % last calculated value of hydrodynamic moments due to added mass (corrected)
        lastForceIterations % Number of iterations in the last time step

        logger; % wsim.logger object containing the logged simulation data from the current or most recent sim

        minForceIterations; % minimum number of iterations which will be performed at each time step
        maxForceIterations; % maximum number of iteration allowed before convergence
        absForceTolerance;  % absolute tolerance on the forces for convergence
        absMomentTolerance; % absolute tolerance on the moments for convergence
        relForceTolerance;  % relative  tolerance on the forces for convergence
        
        mBDynSystem;
        hydroSystem;
    end

    properties (GetAccess = private, SetAccess = private)

        ptoIndexMap;
        hydroNodeIndexMap;
        
        hydroMotionSyncStep = 1;
        hydroMotionSyncStepCount = 1;
        simStepCount;
        mBDynMBCNodal;
        outputFilePrefix;


        nMBDynNodes;
        mBDynPostProc;

        wavePlotX;
        wavePlotY;

        drawAxesH;
        waveSurfH;

        defaultAbsForceTolerance;
        defaultAbsMomentTolerance;
        
        progressBar;

    end

    methods

        function self = wecSim (hsys, mbsys, varargin)
            % wsim.wecSim constructor
            %
            % Syntax
            %
            % wsobj = wsim.wecSim (hsys, mbsys)
            % wsobj = wsim.wecSim (..., 'Parameter', Value)
            %
            % Description
            %
            % wsim.wecSim is a class which manages wave energy converter
            % simulations.
            %
            % Input
            %
            %  hsys - a wsim.hydroSystem object
            %
            %  mbsys - a mbdyn.pre.sysem object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'PTOs' - optional scalar wsim.powerTakeOff object or a cell
            %    array of multiple wsim.powerTakeOff objects (or derived
            %    objects). The wsim.powerTakeOff class is used to simplify
            %    the calculation of forces and moments for the multibody
            %    simulation and make it easier to apply them in the
            %    multibody simulation with the correct directions etc.
            %    Rather than using the basic wsim.powerTakeOff class,
            %    usually a derived class such as wsim.rotaryPowerTakeOff or
            %    wsim.linearPowerTakeOff is used. See the help for these
            %    classes for more information.
            %
            %  'LoggingSettings' - optional wsim.loggingSettings object.
            %    wsim.loggingSetting is a small class which holds the
            %    settings determining what is logged during the simulation.
            %    The simulation data is stored in a wsim.logger object. The
            %    logger object for a simulation is returned by the 'run'
            %    method. See the help for both wsim.loggingSettings and
            %    wsim.logger for more information. If the LoggingSettings
            %    option is not used a new loggingSettings object is created
            %    internally with all logging turned off.
            %
            %  'Controller' - optional wsim.wecController or derived
            %    object. This is a controller, usually also used in
            %    conjuction with a PTO which also accepts a controller
            %    object as input (e.g. wsim.linearPTOWithController). The
            %    same controller must be given to both the PTO and the
            %    wsim.wecSim object. The controller will be given access to
            %    the wsim.wecSim object in which is resides, so it will
            %    have access to the data log via the logger property of the
            %    wsim.wecSim object.
            %
            %
            % Output
            %
            %  wsobj - wsim.wecSim
            %
            %
            %
            % See Also: wsim.hydroSystem, wsim.logger,
            %           wsim.loggingSettings, mbdyn.pre.system
            %

            options.PTOs = {};
            options.NEMOHSim = [];
            options.LoggingSettings = wsim.loggingSettings ();
            options.Controller = [];

            options = parse_pv_pairs (options, varargin);

            assert (isa (hsys, 'wsim.hydroSystem'), ...
                'hsys must be an wsim.hydroSystem object');

            assert (isa (mbsys, 'mbdyn.pre.system'), ...
                'mbsys must be an mbdyn.pre.system object');

            if ~isempty (options.NEMOHSim)
                assert (isa (options.NEMOHSim, 'nemoh.simulation'), ...
                    'NEMOHSim must be a nemoh.simulation object' );
            end

            if ~isempty (options.Controller)
                assert (isa (options.Controller, 'wsim.wecController'), ...
                    'Controller must be a wsim.wecController object' );
            end

            self.caseDirectory = hsys.simu.caseDir;
            self.loggingSettings = options.LoggingSettings;
            self.mBDynSystem = mbsys;
            self.hydroSystem = hsys;
            self.readyToRun = false;
            self.simComplete = false;
            self.nMBDynNodes = nan;
            self.wecController = options.Controller;

            % add the spplied PTOs
            if ~isempty (options.PTOs)
                if ~iscell (options.PTOs)
                    options.PTOs = {options.PTOs};
                end
                for ptoind = 1:numel (options.PTOs)
                    self.addPTO (options.PTOs{ptoind});
                end
            end

        end

        function addPTO (self, pto)
            % add one or more PTO objects to a simulation
            %
            % Syntax
            %
            % addPTO (wsobj, pto)
            %
            % Description
            %
            % addPTO adds one or more wsim.powerTakeOff (or derived)
            % objects to a wsim.wecSim simulation. The wsim.powerTakeOff
            % class is used to simplify the calculation of forces and
            % moments for the multibody simulation and make it easier to
            % apply them in the multibody simulation with the correct
            % directions etc. Rather than using the basic wsim.powerTakeOff
            % class, usually a derived class such as
            % wsim.rotaryPowerTakeOff or wsim.linearPowerTakeOff is used.
            % See the help for these classes for more information.
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            %  pto - wsim.powerTakeOff object (or derived object), or a
            %   cell array of wsim.powerTakeOff objects.
            %
            %
            % See Also: wsim.powerTakeOff, wsim.rotaryPowerTakeOff,
            %           wsim.linearPowerTakeOff
            %
            %

            assert (isa (pto, 'wsim.powerTakeOff'), ...
                'pto must be a wsim.powerTakeOff object (or derived class)');

            % check that the PTO nodes are in the MBDyn system external
            % structual force element
            strinfo = self.mBDynSystem.externalStructuralInfo ();

            for ptoind = 1:numel(self.powerTakeOffs)

                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)

                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.referenceNode

                        foundptonode = true;

                        break;

                    end

                end

                if foundptonode == false
                    error ('Could not find reference node for new PTO in MBDyn system structural external force element');
                end

                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)

                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.otherNode

                        foundptonode = true;

                        break;

                    end

                end

                if foundptonode == false
                    error ('Could not find non-reference node for new PTO in MBDyn system structural external force element');
                end

            end

            % append it to the existing PTO objects
            self.powerTakeOffs = [ self.powerTakeOffs, {pto}];

            % set the id of the pto
            self.powerTakeOffs{end}.id = numel (self.powerTakeOffs);

            % mark ready to run false, as the PTO index map will need to be
            % updated before proceeding to run a simulation
            self.readyToRun = false;

        end

        function prepare (self)
            % prepare the wsim.wecSim simulation so it is ready to run
            %
            % Syntax
            %
            % prepare (wsobj)
            %
            % Description
            %
            % prepare performs various preprocessing taks to get the
            % simulation into a state where it is ready to be run. prepare
            % must be called before the run method is called.
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            %
            % See Also: wsim.wecSim.run
            %

            self.mapPTOForceInds ();
            self.mapHydroForceInds ();
            self.initDataStructures ();

            if ~isempty (self.wecController)
                % give the controller access to this object so it can get
                % the data
                self.wecController.wecSimObj = self;
            end

            % TODO: check if NEMOH data needs updated

            % ensure hydro system transient simulation is ready
%             self.hydroSystem.timeDomainSimSetup ();

            % clear any mbdyn.postproc object
            self.mBDynPostProc = [];

            % calculate the default absolute force tolerance based on an
            % acceleration for each body
            absforcetols = zeros (3, self.hydroSystem.nHydroBodies);
            absmomenttols = zeros (3, self.hydroSystem.nHydroBodies);
            minaccel = repmat (1e-3, 3, 1);
            minomegap = repmat (1e-3, 3, 1);
            for hbind = 1:self.hydroSystem.nHydroBodies

                absforcetols(1:3,hbind) = min ([1; 1; 1], self.hydroSystem.getBodyProperty (hbind, 'mass') * minaccel);
                absmomenttols(1:3,hbind) = min ([1; 1; 1], diag (self.hydroSystem.getBodyProperty (hbind, 'momOfInertia')) * minomegap);

            end
            self.defaultAbsForceTolerance = absforcetols;
            self.defaultAbsMomentTolerance = absmomenttols;

            self.readyToRun = true;
            self.simComplete = false;

        end

        function [ datalog, mbdyn_postproc ] = run (self, varargin)
            % Run entire the WEC simulation
            %
            % Syntax
            %
            % [datalog, mbdyn_postproc] = run (wsobj)
            % [datalog, mbdyn_postproc] = run (..., 'Parameter', Value)
            %
            % Description
            %
            % run runs the wsim.wecSim simulation to completion without
            % stopping.
            %
            % This is a shortcut for calling simStart, then calling simStep
            % or simSteps until the simulation is finished and then calling
            % simFinish. The run method is used when no interaction is
            % required during the simulation, and is the most common way of
            % running a simulation.
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'MBDynInputFile' - optional character vector containing the
            %    full path which should be used for the MBDyn input file.
            %    If not supplied, a default file name will generated with
            %    the extension '.mbd' and placed in the outut directory.
            %
            %  'Verbosity' - flag controlling how much text output is
            %    produced at the command line during simulation, and in the
            %    MBDyn output file during the simulation. Verbosity should
            %    be a scalar integer. The higher the value, the more output
            %    will be produced. Default is 0, indicating the minimum
            %    amount of output. The value of verbosity adds that umber
            %    of -p input arguments to the command line used to launch
            %    MBDyn (more -p inputs means more verbose output, see the
            %    documentation for MBDYn for more information).
            %
            %  'AbsForceTolerance' - positive scalar value representing the
            %    absolute tolerance on the forces applied. The system is
            %    solved iteratively with forces and motion being
            %    recalculated repeatedly until convergence is achieved.
            %    This is the smallest absolute change possible before the
            %    force calculation is considered to be converged. Default
            %    is 100N if not supplied.
            %
            %  'RelForceTolerance' - positive scalar value representing the
            %    relative tolerance on the forces applied. The system is
            %    solved iteratively with forces and motion being
            %    recalculated repeatedly until convergence is achieved.
            %    This is the smallest change possible in the force relative
            %    to it's absolute value before the force calculation is
            %    considered to be converged. Default is 1e-5 if not
            %    supplied.
            %
            %  'MinIterations' - scalar integer indicating the minimum
            %    number of iterations of the force/motion calculation which
            %    should be performed. Default is 0.
            %
            %  'TimeExecution' - true/false flag indicating whether to time
            %    the execution of the simulation and report it at the end
            %    of the sim.
            %
            %  'HydroMotionSyncSteps' -
            %
            %  'ForceMBDynNetCDF' - true/false flag indicating whether to
            %    ensure MBDyn will generate a netcdf format output file by
            %    modifying the supplied mbdyn.pre.system object
            %    OutputResults setting in the control data (if necessary).
            %    Default is true.
            %
            %  'OutputFilePrefix' - character vector containing the output
            %    path prefix. This
            %    is the name of output files, but without their file
            %    extension. e.g. /home/jbloggs/my_mbdyn_sim will create the
            %    files:
            %
            %    /home/jbloggs/my_mbdyn_sim.frc
            %    /home/jbloggs/my_mbdyn_sim.ine
            %    /home/jbloggs/my_mbdyn_sim.out
            %    /home/jbloggs/my_mbdyn_sim.mov
            %    /home/jbloggs/my_mbdyn_sim.jnt
            %    /home/jbloggs/my_mbdyn_sim.log
            %
            %    and/or possibly a netcdf format file:
            %
            %    /home/jbloggs/my_mbdyn_sim.nc
            %
            %    A windows example might look like
            %    C:\Users\IEUser\Documents\my_mbdyn_sim
            %    producing the files:
            %
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.frc
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.ine
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.out
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.mov
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.jnt
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.log
            %
            %    and/or:
            %
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.nc
            %
            %    The netcdf format is preferred for the mbdyn.postproc
            %    class output in the mbdy_postproc variable. By default
            %    wsim.wecSim will force the output of a netcdf format file
            %    (see the ForceMBDynNetCDF output above).
            %
            %  'MaxIterations' - optionaL scalar integer indicating the
            %    maximum number of interations which should be performed in
            %    any simulation step. If this number of iterations is
            %    exceeded, the simulation will simply advance anyway using
            %    the last set of forces calculated, regardless of meeting
            %    tolerances etc. so this option should be used with care.
            %    Default is Inf, so there is no limit to the number of
            %    iterations that will be performed.
            %
            %  'TimeExecution' - optional true/false flag indicating
            %    whether to time the execution of the simulation. Default
            %    is false.
            %
            %  'ShowProgress' - optional true/false flag indicating
            %    whether to show the percentage simulation complete using a
            %    text based progress bar in the command window.
            %
            % Output
            %
            %  datalog - wsim.logger object containing the logged data (if
            %   any) from the simulation.
            %
            %  mbdyn_postproc - mbdyn.postproc object which can be used to
            %   access the full output for all joints and elements from
            %   MBDyn produced during the simulation. See the help for
            %   mbdyn.postproc for more information.
            %
            %
            % See Also: wsim.logger, mbdyn.postproc
            %

            % options passed to simStart
            options.MBDynInputFile = '';
            thedate = datestr(now (), 'yyyy-mm-dd_HH-MM-SS-FFF');
            options.OutputFilePrefix = fullfile (self.caseDirectory, ['output_', thedate], ['mbdyn_sim_results_', thedate]);
            options.Verbosity = 0;
            options.AbsForceTolerance = self.defaultAbsForceTolerance;
            options.AbsMomentTolerance = self.defaultAbsMomentTolerance;
            options.RelForceTolerance = 1e-6;
            options.MinIterations = 0;
            options.MaxIterations = self.mBDynSystem.problems{1}.maxIterations;
            options.HydroMotionSyncSteps = 1;
            options.ForceMBDynNetCDF = true;
            options.SyncMBDynNetCDF = false;
            options.ShowProgress = false;

            % specific to run
            options.TimeExecution = false;

            options = parse_pv_pairs (options, varargin);

            check.isLogicalScalar (options.TimeExecution, true, 'TimeExecution');

            [status, datalog] = self.simStart ( 'MBDynInputFile', options.MBDynInputFile, ...
                            'OutputFilePrefix', options.OutputFilePrefix, ...
                            'Verbosity', options.Verbosity, ...
                            'AbsForceTolerance', options.AbsForceTolerance, ...
                            'AbsMomentTolerance', options.AbsMomentTolerance, ...
                            'RelForceTolerance', options.RelForceTolerance, ...
                            'MinIterations', options.MinIterations, ...
                            'MaxIterations', options.MaxIterations, ...
                            'HydroMotionSyncSteps', options.HydroMotionSyncSteps, ...
                            'ForceMBDynNetCDF', options.ForceMBDynNetCDF, ...
                            'SyncMBDynNetCDF', options.SyncMBDynNetCDF, ...
                            'ShowProgress', options.ShowProgress );

            if options.TimeExecution, tic; end

            while status == 0

                % step through the simulation
                status = self.simStep ();

            end

            if nargout > 0
                mbdyn_postproc = self.simFinish ();
            else
                self.simFinish ();
            end

            if options.TimeExecution, toc; end

        end

        function [status, datalog] = simStart (self, varargin)
            % Start a WEC simulation
            %
            % Syntax
            %
            % [datalog, mbdyn_postproc] = simStart (wsobj)
            % [datalog, mbdyn_postproc] = simStart (..., 'Parameter', Value)
            %
            % Description
            %
            % simStart starts wsim.wecSim simulation, performing the
            % initial time step.
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'MBDynInputFile' - optional character vector containing the
            %    full path which should be used for the MBDyn input file.
            %    If not supplied, a default file name will generated with
            %    the extension '.mbd' and placed in the outut directory.
            %
            %  'Verbosity' - flag controlling how much text output is
            %    produced at the command line during simulation, and in the
            %    MBDyn output file during the simulation. Verbosity should
            %    be a scalar integer. The higher the value, the more output
            %    will be produced. Default is 0, indicating the minimum
            %    amount of output. The value of verbosity adds that umber
            %    of -p input arguments to the command line used to launch
            %    MBDyn (more -p inputs means more verbose output, see the
            %    documentation for MBDYn for more information).
            %
            %  'AbsForceTolerance' - positive scalar value representing the
            %    absolute tolerance on the forces applied. The system is
            %    solved iteratively with forces and motion being
            %    recalculated repeatedly until convergence is achieved.
            %    This is the smallest absolute change possible before the
            %    force calculation is considered to be converged. Default
            %    is 100N if not supplied.
            %
            %  'RelForceTolerance' - positive scalar value representing the
            %    relative tolerance on the forces applied. The system is
            %    solved iteratively with forces and motion being
            %    recalculated repeatedly until convergence is achieved.
            %    This is the smallest change possible in the force relative
            %    to it's absolute value before the force calculation is
            %    considered to be converged. Default is 1e-5 if not
            %    supplied.
            %
            %  'MinIterations' - scalar integer indicating the minimum
            %    number of iterations of the force/motion calculation which
            %    should be performed. Default is 0.
            %
            %  'TimeExecution' - true/false flag indicating whether to time
            %    the execution of the simulation and report it at the end
            %    of the sim.
            %
            %  'HydroMotionSyncSteps' -
            %
            %  'ForceMBDynNetCDF' - true/false flag indicating whether to
            %    ensure MBDyn will generate a netcdf format output file by
            %    modifying the supplied mbdyn.pre.system object
            %    OutputResults setting in the control data (if necessary).
            %    Default is true.
            %
            %  'OutputFilePrefix' - character vector containing the output
            %    path prefix. This
            %    is the name of output files, but without their file
            %    extension. e.g. /home/jbloggs/my_mbdyn_sim will create the
            %    files:
            %
            %    /home/jbloggs/my_mbdyn_sim.frc
            %    /home/jbloggs/my_mbdyn_sim.ine
            %    /home/jbloggs/my_mbdyn_sim.out
            %    /home/jbloggs/my_mbdyn_sim.mov
            %    /home/jbloggs/my_mbdyn_sim.jnt
            %    /home/jbloggs/my_mbdyn_sim.log
            %
            %    and/or possibly a netcdf format file:
            %
            %    /home/jbloggs/my_mbdyn_sim.nc
            %
            %    A windows example might look like
            %    C:\Users\IEUser\Documents\my_mbdyn_sim
            %    producing the files:
            %
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.frc
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.ine
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.out
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.mov
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.jnt
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.log
            %
            %    and/or:
            %
            %    C:\Users\JBloggs\Documents\my_mbdyn_sim.nc
            %
            %    The netcdf format is preferred for the mbdyn.postproc
            %    class output in the mbdy_postproc variable. By default
            %    wsim.wecSim will force the output of a netcdf format file
            %    (see the ForceMBDynNetCDF output above).
            %
            %  'MaxIterations' - optionaL scalar integer indicating the
            %    maximum number of interations which should be performed in
            %    any simulation step. If this number of iterations is
            %    exceeded, the simulation will simply advance anyway using
            %    the last set of forces calculated, regardless of meeting
            %    tolerances etc. so this option should be used with care.
            %    Default is Inf, so there is no limit to the number of
            %    iterations that will be performed.
            %


            options.MBDynInputFile = '';
            thedate = datestr(now (), 'yyyy-mm-dd_HH-MM-SS-FFF');
            options.OutputFilePrefix = fullfile (self.caseDirectory, ['output_', thedate], ['mbdyn_sim_results_', thedate]);
            options.Verbosity = 0;
            options.AbsForceTolerance = self.defaultAbsForceTolerance;
            options.AbsMomentTolerance = self.defaultAbsMomentTolerance;
            options.RelForceTolerance = 1e-5;
            options.MinIterations = 0;
            options.MaxIterations = self.mBDynSystem.problems{1}.maxIterations;
            options.HydroMotionSyncSteps = 1;
            options.ForceMBDynNetCDF = true;
            options.SyncMBDynNetCDF = false;
            options.ShowProgress = false;

            options = parse_pv_pairs (options, varargin);

            if isempty (options.MBDynInputFile)
                [ pathstr, ~ ] = fileparts (options.OutputFilePrefix);
                self.mBDynInputFile = fullfile (pathstr, ['mbdyn_input_file_', thedate, '.mbd']);
            else
                assert (ischar (options.MBDynInputFile), ...
                    'MBDynInputFile must be string containing the file path where the MBDyn input file should be generated');
            end

            check.isPositiveScalarInteger (options.Verbosity, true, 'Verbosity', true);

            assert (ischar (options.OutputFilePrefix), ...
                'OutputFilePrefix must be a string');

            if exist (options.OutputFilePrefix, 'dir') == 0
                pathstr = fileparts (options.OutputFilePrefix);
                mkdir (pathstr);
            end

            % check for positive numeric scalar
            if check.isNumericScalar (options.AbsForceTolerance, false, 'AbsForceTolerance', 1)

                options.AbsForceTolerance = repmat (options.AbsForceTolerance, 3, self.hydroSystem.nHydroBodies);

            elseif size (options.AbsForceTolerance, 1) == 3 ...
                   	&& size (options.AbsForceTolerance, 2) == 1

                options.AbsForceTolerance = repmat (options.AbsForceTolerance, 1, self.hydroSystem.nHydroBodies);

            else
                assert ( ( size (options.AbsForceTolerance, 1) == 3 ...
                                && size (options.AbsForceTolerance, 2) == self.hydroSystem.nHydroBodies ...
                         ), ...
                         'AbsForceTolerance must be a scalar value, or a three element vector, or a (3 x nHydroBodies) matrix' );

            end

            assert (all (options.AbsForceTolerance(:) > 0), 'All AbsForceTolerance values must be > 0');

            % check for positive numeric scalar
            if check.isNumericScalar (options.AbsMomentTolerance, false, 'AbsMomentTolerance', 1)

                options.AbsMomentTolerance = repmat (options.AbsMomentTolerance, 3, self.hydroSystem.nHydroBodies);

            elseif size (options.AbsMomentTolerance, 1) == 3 ...
                   	&& size (options.AbsMomentTolerance, 2) == 1

                options.AbsMomentTolerance = repmat (options.AbsMomentTolerance, 1, self.hydroSystem.nHydroBodies);

            else
                assert ( ( size (options.AbsMomentTolerance, 1) == 3 ...
                                && size (options.AbsMomentTolerance, 2) == self.hydroSystem.nHydroBodies ...
                         ), ...
                         'AbsMomentTolerance must be a scalar value, or a three element vector, or a (3 x nHydroBodies) matrix' );

            end

            assert (all (options.AbsMomentTolerance(:) > 0), 'All AbsMomentTolerance values must be > 0');


%             check.isNumericScalar (options.AbsForceTolerance, true, 'AbsForceTolerance', 1);
            check.isNumericScalar (options.RelForceTolerance, true, 'RelForceTolerance', 1);

            check.isPositiveScalarInteger (options.MaxIterations, true, 'Verbosity');

            assert (options.MinIterations <= options.MaxIterations, ...
                'MinIterations (%d) is not less than or equal to MaxIterations (%d)', ...
                options.MinIterations,  options.MaxIterations)

            check.isPositiveScalarInteger (options.HydroMotionSyncSteps, true, 'HydroMotionSyncSteps');

            check.isLogicalScalar (options.ForceMBDynNetCDF, true, 'ForceMBDynNetCDF');
            check.isLogicalScalar (options.SyncMBDynNetCDF, true, 'SyncMBDynNetCDF');
            check.isLogicalScalar (options.ShowProgress, true, 'ShowProgress');

            % -----------------   input checking finished

            self.minForceIterations = options.MinIterations;
            self.absForceTolerance = options.AbsForceTolerance;
            self.absMomentTolerance = options.AbsMomentTolerance;
            self.relForceTolerance = options.RelForceTolerance;
            self.maxForceIterations = options.MaxIterations;
            self.outputFilePrefix = options.OutputFilePrefix;
            self.showProgress = options.ShowProgress;

            if self.readyToRun == false
                error ('Simulation is not ready to run, have you run ''prepare'' yet?');
            end

            if options.ForceMBDynNetCDF
                % ensure mbdyn will output a netcdf file so we can load the
                % results from this after the sim
                if isempty (self.mBDynSystem.controlData.OutputResults)
                    self.mBDynSystem.controlData.OutputResults = {'netcdf'};
                    if options.SyncMBDynNetCDF
                        self.mBDynSystem.controlData.OutputResults = [self.mBDynSystem.controlData.OutputResults, {'sync'}];
                    end
                    self.mBDynSystem.controlData.OutputResults = [self.mBDynSystem.controlData.OutputResults, {'no text'}];
                end
            end

            self.hydroMotionSyncStep = options.HydroMotionSyncSteps;
            % set the hydro--mbdyn motion sync step count equal to the
            % number of steps to sync on so that forces are calculated the
            % first time applyHydroForces is called (where this step count
            % is used)
            self.hydroMotionSyncStepCount = self.hydroMotionSyncStep;

            self.simInfo.TStart = self.mBDynSystem.problems{1}.initialTime;
            self.simInfo.TEnd = self.mBDynSystem.problems{1}.finalTime;
            self.simInfo.TStep = self.mBDynSystem.problems{1}.timeStep;
            self.simInfo.MBDynSystem = self.mBDynSystem;
            self.simInfo.HydroSystem = self.hydroSystem;
            self.simInfo.HydroMotionSyncSteps = options.HydroMotionSyncSteps;
            self.simInfo.OutputDirectory = fileparts (self.outputFilePrefix);
            self.simInfo.CaseDirectory = self.caseDirectory;

            % start the controller if there is one
            if ~isempty (self.wecController)
                self.wecController.start ();
            end

            % start the PTO components
            self.startPTOs ();

            % generate input file and start mbdyn

            % create the communicator object. As an mbdyn system object is
            % supplied, the mbdyn input file will be generated
            % automatically
           try
                self.mBDynMBCNodal = mbdyn.mint.MBCNodal ('MBDynPreProc', self.mBDynSystem, ...
                    'UseMoments', true, ...
                    'MBDynInputFile', self.mBDynInputFile, ...
                    'OverwriteInputFile', true, ...
                    'OutputPrefix', self.outputFilePrefix, ...
                    'NodeOrientationType', 'euler 123' ...
                    );
           catch err
               self.readyToRun = false;
               if ~isempty (self.mBDynOutputFile)
                   if exist (self.mBDynOutputFile, 'file')
                       self.displayLastNLinesOfFile (self.mBDynOutputFile, 50);
                   end
               end
               self.cleanUpAfterError ();
               disp(err)
               if ~isempty (self.mBDynOutputFile)
                   error ('Starting MBDyn communication falied, aborting sim, some output might have been sent to the following file:\n%s\nIf so, this may help diagnose the error.', ...
                           self.mBDynOutputFile)
               else
                   error ('Starting MBDyn communication falied, aborting sim, unfortunately no MBDyn Output file is available to view, which means MBDyn probably never even started running.');
               end
           end

            % copy over the input file location to make it easier to
            % examine later if required
            self.mBDynOutputFile = self.mBDynMBCNodal.MBDynOutputFile;

            % ensure MBCNodal is destroyed in the event of a problem
            % (closes communication to MBDyn and tells it to quit so
            % sockets and so on are also cleaned up)
%             if isoctave ()
%                 % FIXME: octave calls the delete method multiple times, so
%                 % calling it has been disabled in the development sources,
%                 % so we make sure it is called for the MBCNodal when wecSim
%                 % finishes. This makes sure the sockets are closed and
%                 % memory is freed. See Octave bug #46497
%                 CC = onCleanup (@() self.mBDynMBCNodal.delete ());
%             end

            self.mBDynMBCNodal.start ('Verbosity', options.Verbosity);

            % get the number of nodes in the problem
            self.nMBDynNodes = self.mBDynMBCNodal.GetNodes ();

            % preallocate variables to hold total forces and moments at
            % each time step
            forces_and_moments = zeros (6, self.nMBDynNodes);

            % take initial time from the MBDyn system problem
            self.lastTime = self.simInfo.TStart;

            % fetch the initial configuration of the nodes from MBDyn and
            % check for errors
            status = self.mBDynMBCNodal.GetMotion ();

            if status ~= 0
                self.readyToRun = false;
                if exist (self.mBDynOutputFile, 'file')
                    self.displayLastNLinesOfFile (self.mBDynOutputFile, 50);
                end
                self.cleanUpAfterError ();
                error ('mbdyn returned %d, aborting sim, check output file:\n%s\nfor clues at to why this happened.', ...
                        status, self.mBDynOutputFile)
            end
            
            if self.showProgress
                self.progressBar = ui.progressbar ();
                self.progressBar.init ('TextLeader', sprintf('Sim Progress (of %ds): ', self.simInfo.TEnd - self.simInfo.TStart));
                self.progressBar.dispProgress (0);
            end

            % get the current motion of the multibody system which was sent
            % by MBDyn
            [pos, vel, accel] = self.getMotion (self.mBDynMBCNodal);

            % store most recently calculated motion so it can be used by
            % the PTO or controller if necessary
            self.lastPositions = pos(1:3,:);
            self.lastAngularPositions = pos(4:6,:);
            self.lastVelocities = vel(1:3,:);
            self.lastAngularVelocities = vel(4:6,:);
            self.lastAccelerations = accel(1:3,:);
            self.lastAngularAccelerations = accel(4:6,:);
            
            % calculate and apply hydrodynamic forces based on the motion
            forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

            % get and apply the PTO forces and moments, in this case the
            % position and velocities etc are obtained from the mbsys
            % object stored in the mb object, hence no need to supply them.
            % The node positions etc. are updated by the call to
            % mb.GetMotion (). This is possible because all PTO objects
            % also have access to the same mbsys object
            forces_and_moments = self.applyPTOForces (forces_and_moments);

            % set the new forces an moments ready to be sent to MBDyn
            self.mBDynMBCNodal.F (forces_and_moments(1:3,:));
            self.mBDynMBCNodal.M (forces_and_moments(4:6,:));

            % send the forces and moments to MBDyn, but noting that we have
            % not yet converged (so MBDyn will send new motion based on
            % these forces, and not advance the self.lastTime step)
            mbconv = self.mBDynMBCNodal.applyForcesAndMoments (false);

            self.lastForceIterations =  1;

            self.lastNodeForcesUncorrected = forces_and_moments(1:3,:);
            self.lastNodeMomentsUncorrected = forces_and_moments(4:6,:);

            self.advanceStep ();

            % now begin the simulation loop, beginning from the next time
            % index
            self.simStepCount = 2;

            if nargout > 1
                datalog = self.logger;
            end

            self.simStarted = true;

        end


        function status = simStep (self)
            % manually advance one simulation step

            self.lastForceIterations = 0;

            assert (self.simStarted, 'simStart must be called before calling simStep');

            status = self.mBDynMBCNodal.GetMotion ();

            if status ~= 0
                self.readyToRun = false;
                return;
            end

            self.lastTime = self.lastTime + self.mBDynSystem.problems{1}.timeStep;

            % get the current motion of the multibody system
            [pos, vel, accel] = self.getMotion (self.mBDynMBCNodal);

            forces_and_moments = zeros (6, self.nMBDynNodes);

            % calculate new hydrodynamic interaction forces (also
            % updates self.lastForceHydro with new values)
            forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

            % Add PTO forces
            forces_and_moments = self.applyPTOForces (forces_and_moments);

            self.mBDynMBCNodal.F (forces_and_moments(1:3,:));
            self.mBDynMBCNodal.M (forces_and_moments(4:6,:));

            mbconv = self.mBDynMBCNodal.applyForcesAndMoments (false);

            self.lastForceIterations = self.lastForceIterations + 1;

            status = self.mBDynMBCNodal.GetMotion ();

            if status ~= 0
                % quit, MBDyn has finished simulation (or some error
                % occured)
                self.readyToRun = false;
                self.simStepCount = self.simStepCount + 1;

                return;
            end

            % get the current motion of the multibody system
            [newpos, newvel, newaccel] = self.getMotion (self.mBDynMBCNodal);

            repeat_force = true;
            if mbconv == 0

                % check if new motion is within tolerances or if we need to
                % recalculate forces
                posdiff = abs (newpos - pos);
                veldiff = abs (newvel - vel);

                if all ([posdiff(:); veldiff(:)] < 1e-8)
                    repeat_force = false;
                end

            end

            if repeat_force

                pos = newpos;
                vel = newvel;
                accel = newaccel;

                % repeat the force calulation to test convergence

                % clear out the previous forces and moments
                forces_and_moments = zeros (6, self.nMBDynNodes);

                % get the previous calculated values of hydrodynamic forces
                prev_hydro_forces_and_moments = [ self.lastForceHydro; self.lastMomentHydro ];

                % calculate new hydrodynamic interaction forces (also
                % updates self.lastForceHydro with new values)
                forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

                % PTO forces
                forces_and_moments = self.applyPTOForces (forces_and_moments);

                self.mBDynMBCNodal.F (forces_and_moments(1:3,:));
                self.mBDynMBCNodal.M (forces_and_moments(4:6,:));

                mbconv = self.mBDynMBCNodal.applyForcesAndMoments (false);

                self.lastForceIterations = self.lastForceIterations + 1;

                % check for force convergence

                hydroforcediff = abs (prev_hydro_forces_and_moments(1:3,:) - self.lastForceHydro);
                hydromomentdiff = abs (prev_hydro_forces_and_moments(4:6,:) - self.lastMomentHydro);

                do_iter = false;

                if mbconv ~= 0

                    % MBDyn hasn't converged, so iterate
                    do_iter = true;

                else
                    % find any components for which the differences are
                    % bigger than the absolute tolerances
                    fabs_exceeded = hydroforcediff > self.absForceTolerance;
                    mabs_exceeded = hydromomentdiff > self.absMomentTolerance;

                    % check these components (if any) against the relative
                    % tolerances
                    if ( anyalldims ( hydroforcediff(fabs_exceeded) > (self.relForceTolerance * self.lastForceHydro(fabs_exceeded)) ) ...
                             || anyalldims ( hydromomentdiff(mabs_exceeded) > (self.relForceTolerance * self.lastMomentHydro(mabs_exceeded)) ) ...
                       )

                        do_iter = true;

                    end

                end

                % iterate until force/motion convergence (or max
                % iterations)
                while do_iter == true

                    status = self.mBDynMBCNodal.GetMotion ();

                    if status ~= 0
                        self.readyToRun = false;
                        break;
                    end

                    % clear out the previous forces and moments
                    forces_and_moments = zeros (6, self.nMBDynNodes);

                    % get the current motion of the multibody system
                    [pos, vel, accel] = self.getMotion (self.mBDynMBCNodal);

                    % get the previous calculated values of hydrodynamic forces
                    prev_hydro_forces_and_moments = [ self.lastForceHydro; self.lastMomentHydro ];

                    % calculate new hydrodynamic interaction forces (also
                    % updates self.lastForceHydro with new values)
                    forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

                    % PTO forces
                    forces_and_moments = self.applyPTOForces (forces_and_moments);

                    self.mBDynMBCNodal.F (forces_and_moments(1:3,:));
                    self.mBDynMBCNodal.M (forces_and_moments(4:6,:));

                    mbconv = self.mBDynMBCNodal.applyForcesAndMoments (false);

                    self.lastForceIterations = self.lastForceIterations + 1;

                    if self.lastForceIterations > self.maxForceIterations
                        warning ('mbdyn force iterations exceeded max allowed');
                        self.readyToRun = false;
                        status = -2;
                        return;
                    end

                    hydroforcediff = abs (prev_hydro_forces_and_moments(1:3,:) - self.lastForceHydro);
                    hydromomentdiff = abs (prev_hydro_forces_and_moments(4:6,:) - self.lastMomentHydro);

                    do_iter = false;

                    if mbconv ~= 0

                        % MBDyn hasn't converged, so iterate
                        do_iter = true;

                    else
                        % find any components for which the differences are
                        % bigger than the absolute tolerances
                        fabs_exceeded = hydroforcediff > self.absForceTolerance;
                        mabs_exceeded = hydromomentdiff > self.absMomentTolerance;

                        % check these components (if any) against the relative
                        % tolerances
                        if ( anyalldims ( hydroforcediff(fabs_exceeded) > (self.relForceTolerance * self.lastForceHydro(fabs_exceeded)) ) ...
                               || anyalldims ( hydromomentdiff(mabs_exceeded) > (self.relForceTolerance * self.lastMomentHydro(mabs_exceeded)) ) ...
                           )

                            do_iter = true;

                        end

                    end

                end

                % get latest motion from MBDyn
                status = self.mBDynMBCNodal.GetMotion ();

            end

            if status ~= 0

                % quit, MBDyn sim is finished or there's an error
                self.readyToRun = false;

                self.simStepCount = self.simStepCount + 1;

                return;
            end

            % apply the last calculated forces from the iteration to
            % the system
            self.mBDynMBCNodal.F (forces_and_moments(1:3,:));
            self.mBDynMBCNodal.M (forces_and_moments(4:6,:));

            mbconv = self.mBDynMBCNodal.applyForcesAndMoments (true);

            % store most recently calculated motion so it can be logged
            self.lastPositions = pos(1:3,:);
            self.lastAngularPositions = pos(4:6,:);
            self.lastVelocities = vel(1:3,:);
            self.lastAngularVelocities = vel(4:6,:);
            self.lastAccelerations = accel(1:3,:);
            self.lastAngularAccelerations = accel(4:6,:);
            self.lastNodeForcesUncorrected = forces_and_moments(1:3,:);
            self.lastNodeMomentsUncorrected = forces_and_moments(4:6,:);

            self.advanceStep ();

            self.simStepCount = self.simStepCount + 1;
            
            if self.showProgress
                
                percent_progress = 100 * self.lastTime / self.simInfo.TEnd;
                
                if percent_progress >= self.progressBar.lastProgress + 5
                    self.progressBar.dispProgress (percent_progress);
                end
                
            end

        end

        function [status, stepcount] = simSteps (self, nsteps)
            % manually advance multiple simulation steps

            status = 0;
            stepcount = 0;

            while status == 0 && stepcount < nsteps

                status = simStep (self);

                stepcount = stepcount + 1;

            end

        end


        function mbdyn_postproc = simFinish (self)
            % finalaise a simulation and clean up

            assert ( self.simStarted, 'You must call simStart before calling simFinish' );

            % clear the MBCNodal object, hich should trigger its delete
            % method
            self.mBDynMBCNodal = [];

            % finish off
            fprintf ( 1, 'Reached time %fs (%fs of an intended %fs), in %d steps, postprocessing ...\n', ...
                      self.lastTime, ...
                      self.lastTime - self.simInfo.TStart, ...
                      self.simInfo.TEnd - self.simInfo.TStart, ...
                      self.simStepCount-1 );

            if self.loggingSettings.nodeForces ...
                    || self.loggingSettings.forceAddedMass ...
                    || self.loggingSettings.nodeMoments ...
                    || self.loggingSettings.momentAddedMass

                if self.loggingSettings.nodeForces
                    self.logger.setSeries ('NodeForces', self.logger.data.NodeForcesUncorrected);
                end

                if self.loggingSettings.nodeMoments
                    self.logger.setSeries ('NodeMoments', self.logger.data.NodeMomentsUncorrected);
                end

                % the total forces on hydrodynamic bodies must be corrected
                [ corrected_node_forces_and_moments, ...
                  F_and_M_added_mass ] = ...
                    correctAddedMassForce ( self.hydroSystem ...
                                            , [ self.logger.data.NodeForcesUncorrected(:,self.hydroNodeIndexMap(:,1),:); ...
                                                self.logger.data.NodeMomentsUncorrected(:,self.hydroNodeIndexMap(:,1),:) ] ...
                                            , [ self.logger.data.ForceAddedMassUncorrected; ...
                                                self.logger.data.MomentAddedMassUncorrected ] ...
                                            , [ self.logger.data.Accelerations; ...
                                                self.logger.data.AngularAccelerations ] ...
                                          );

                if self.loggingSettings.nodeForces
                    self.logger.data.NodeForces(:,self.hydroNodeIndexMap(:,1),:) = corrected_node_forces_and_moments(1:3,:,:);
                end

                if self.loggingSettings.nodeMoments
                    self.logger.data.NodeMoments(:,self.hydroNodeIndexMap(:,1),:) = corrected_node_forces_and_moments(4:6,:,:);
                end

                if self.loggingSettings.forceAddedMass
                    self.logger.setSeries ('ForceAddedMass', F_and_M_added_mass(1:3,:,:));
                end

                if self.loggingSettings.momentAddedMass
                    self.logger.setSeries ('MomentAddedMass', F_and_M_added_mass(4:6,:,:));
                end

            end


            % tell the PTOs we are done
            self.finishPTOs ();

            % tell the controller we are done (if there is one)
            if ~isempty (self.wecController)
                self.wecController.start ();
            end
            
            if self.showProgress
                self.progressBar.done ();
            end

            fprintf (1, 'Simulation complete\n');

            self.logger.truncateAllVariables ();

            self.simComplete = true;
            self.simStarted = false;
            self.readyToRun = false;

            if nargout > 0
                mbdyn_postproc = mbdyn.postproc (self.outputFilePrefix, self.mBDynSystem);
                self.mBDynPostProc = mbdyn_postproc;
            end
            
            if self.simStepCount < 5
                if ~isempty (self.mBDynOutputFile)
                   if exist (self.mBDynOutputFile, 'file')
                       fprintf (1, 'Simulation achieved less than 5 steps before stopping, printing the MBDyn output:\n');
                       self.displayLastNLinesOfFile (self.mBDynOutputFile, 50);
                   end
                end
            end 

        end


        function plotdata = drawStep (self, tind, varargin)
            % draw the system at the given time step index
            %
            % Syntax
            %
            % plotdata = drawStep (wsobj, tind, 'Parameter', value)
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            %  Addtional optional arguments are supplied as parameter-value
            %  pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is true.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %

            if ~self.simComplete
                error ('Simulation data is not available for plotting')
            end

            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Light = false;
            options.Title = true;
            options.OnlyNodes = [];
            options.ForceRedraw = false;
            options.DrawWaves = true;
            options.View = [];
            options.FigPositionAndSize = [];

            options = parse_pv_pairs (options, varargin);

            self.checkAxes (options.PlotAxes);

            if isempty (self.mBDynPostProc)
                self.mBDynPostProc = mbdyn.postproc (self.mBDynInputFile, self.mBDynSystem);
            end

            if isempty (options.OnlyNodes)
                options.OnlyNodes = self.mBDynPostProc.nNodes;
            end

            if options.DrawWaves

                wavedrawfcn = @(hax, tind) wavePlot3D ( self, ...
                                                        [], ...
                                                        'PlotAxes', hax, ...
                                                        'DomainSize', 'fromaxes', ...
                                                        'TimeIndex', tind );

                plotdata = self.mBDynPostProc.drawStep ( tind, ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              'ForceRedraw', options.ForceRedraw, ...
                                              'ExternalDrawFcn', wavedrawfcn, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize );

            else

                plotdata = self.mBDynPostProc.drawStep ( tind, ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              'ForceRedraw', options.ForceRedraw, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize  );

            end

        end


        function animate (self, varargin)
            % animate the wave energy system
            %
            % Syntax
            %
            % animate (wsobj, 'Parameter', value)
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            % Addtional optional arguments are supplied as parameter-value
            % pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is false.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            %  'ExternalDrawFcn' - function handle or string with function
            %    which will be called after drawing the scene is complete.
            %    Intended to be used by external programs to add their own
            %    features to the scene. The funciton must take two
            %    arguments, the first is the handle to the axes contining
            %    the plot, the second is the time step index of the
            %    simulation.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
            %
            %  'VideoFile' - optional file name into which the aimation
            %    will be saved as a video. Default is empty, in which case
            %    no video file will be produced unless a VideoWriter object
            %    is supplied instead through the 'VideoWriter' option (see
            %    below).
            %
            %  'VideoWriter' - optional VideoWriter object which will be
            %    used to create a video instead of creating one internally.
            %    This option cannot be used at the same time as the
            %    'VideoFile' option. If this option is used, the
            %    VideoProfile option (see below) is ignored. However,
            %    animate will still modify the FrameRate property orf the
            %    supplied VideoWriter object internally, and the VideoSpeed
            %    option (see below) will also be applied. Default is empty,
            %    in which case no video file will be produced unless a
            %    video file path is supplied instead through the
            %    'VideoFile' option (see above).
            %
            %  'VideoSpeed' - optional speed multiplier for the video when
            %    using the 'VideoFile'or 'VideoWriter' option. It must be a
            %    scalar value greater than 0. The animation will play a
            %    speed multiplied by this factor (by changing the video
            %    frame rate) rather than real time. Default is 1 (no
            %    speedup).
            %
            %  'VideoProfile' - optional character vector indicating what
            %    type of video compression to use when using the
            %    'VideoFile' option. This option coresponds to the profile
            %    option for the VideoWrite class, soo see the help for
            %    VideoWriter to see what the possible options are. Default
            %    is 'Motion JPEG AVI'.
            %

            if ~self.simComplete
                error ('No simulation results are available for plotting')
            end

            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.PlotTrajectories = false;
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Skip = 1;
            options.Light = false;
            options.VideoFile = [];
            options.VideoSpeed = 1;
            options.VideoProfile = 'Motion JPEG AVI';
            options.VideoWriter = [];
            options.VideoQuality = 75;
            options.OnlyNodes = [];
            options.ExternalDrawFcn = [];
            options.View = [];
            options.DrawWaves = true;
            options.FigPositionAndSize  = [];
            options.Title = false;

            options = parse_pv_pairs (options, varargin);

            self.checkAxes (options.PlotAxes);

            if isempty (self.mBDynPostProc)
                self.mBDynPostProc = mbdyn.postproc (self.mBDynInputFile, self.mBDynSystem);
            end

            if isempty (options.OnlyNodes)
                options.OnlyNodes = 1:self.mBDynPostProc.nNodes;
            end

            if options.DrawWaves

                wavedrawfcn = @(hax, tind) wavePlot3D ( self, ...
                                                        [], ...
                                                        'PlotAxes', hax, ...
                                                        'DomainSize', 'fromaxes', ...
                                                        'TimeIndex', tind );

                self.mBDynPostProc.animate ( ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Skip', options.Skip, ...
                                              'Light', options.Light, ...
                                              'VideoFile', options.VideoFile, ...
                                              'VideoSpeed', options.VideoSpeed, ...
                                              'VideoProfile', options.VideoProfile, ...
                                              'VideoWriter', options.VideoWriter, ...
                                              'VideoQuality', options.VideoQuality, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              ...'ForceRedraw', options.ForceRedraw, ...
                                              'ExternalDrawFcn', wavedrawfcn, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize );

            else

                self.mBDynPostProc.animate (  ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'VideoFile', options.VideoFile, ...
                                              'VideoSpeed', options.VideoSpeed, ...
                                              'VideoProfile', VideoProfile, ...
                                              'VideoWriter', options.VideoWriter, ...
                                              'VideoQuality', options.VideoQuality, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              ...'ForceRedraw', options.ForceRedraw, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize  );

            end

        end


        function [hsurf, hax] = wavePlot3D (self, t, varargin)
            % plot the wave elevation as a 3D surface

            options.PlotAxes = [];
            options.NumPointsX = self.wavePlotNumPointsX;
            options.NumPointsY = self.wavePlotNumPointsY;
            options.DomainSize = self.wavePlotDomainSize;
            options.DomainCorner = [];
            options.TimeIndex = [];
            options.ForceRedraw = false;

            options = parse_pv_pairs (options, varargin);

            check.isLogicalScalar (options.ForceRedraw, true, 'ForceRedraw');

            self.checkAxes (options.PlotAxes);

            hax = self.drawAxesH;

            if ~isempty (options.TimeIndex)
                % override the t value with logged time from index
                check.isPositiveScalarInteger (options.TimeIndex, true, 'TimeIndex', false);
                t = self.logger.data.Time (options.TimeIndex);
            end

%             xlim = get (hax, 'Xlim');
%             ylim = get (hax, 'Ylim');

            if ischar (options.DomainSize)
                switch options.DomainSize

                    case 'fromaxes'

                        xlim = get (hax, 'Xlim');
                        ylim = get (hax, 'Ylim');
                        options.DomainSize = [ xlim(2) - xlim(1), ylim(2) - ylim(1) ];
                        options.DomainCorner = [ xlim(1), ylim(1) ];

                end
            end

            % the following test is supposed to capture when we are doing a
            % brand new plot, or changing the domain of the existing plot
            if isempty (self.wavePlotDomainSize) && isempty (options.DomainSize)

                options.DomainSize = [ 10, 10 ];

            elseif isscalar (options.DomainSize)

                options.DomainSize = [options.DomainSize, options.DomainSize];

            end

            if isempty (options.DomainCorner)
                options.DomainCorner = [ -options.DomainSize(1)/2, -options.DomainSize(2)/2 ];
            end

            % check if the desired plot parameters have changed at all,
            % requiring the surface to be re-plotted
            if isempty (self.waveSurfH) ...
                || isempty (self.wavePlotDomainSize) ...
                || options.NumPointsX ~= self.wavePlotNumPointsX ...
                || options.NumPointsY ~= self.wavePlotNumPointsY ...
                || ~all (options.DomainSize == self.wavePlotDomainSize) ...
                || options.ForceRedraw ...
                || ~ishghandle (self.waveSurfH)

                self.wavePlotNumPointsX = options.NumPointsX;
                self.wavePlotNumPointsY = options.NumPointsY;
                self.wavePlotDomainSize = options.DomainSize;

                % make the grid
                x = linspace ( options.DomainCorner(1), ...
                               options.DomainCorner(1) + options.DomainSize(1), ...
                               self.wavePlotNumPointsX );

                y = linspace ( options.DomainCorner(2), ...
                               options.DomainCorner(2) + options.DomainSize(2), ...
                               self.wavePlotNumPointsY );

                [self.wavePlotX, self.wavePlotY] = meshgrid (x, y);

                wavecolour = [0.3010, 0.7450, 0.9330];

                if ~isempty (self.waveSurfH) && ishghandle (self.waveSurfH)
                    delete (self.waveSurfH);
                end

                hold on
                self.waveSurfH = surf ( hax, ...
                                        self.wavePlotX, ...
                                        self.wavePlotY, ...
                                        self.hydroSystem.waves.waveElevationGrid (t, self.wavePlotX, self.wavePlotY), ...
                                        'EdgeColor', wavecolour, ...
                                        'FaceColor', wavecolour, ...
                                        'EdgeAlpha', 0.25, ...
                                        'FaceAlpha', 0.25 );
                hold off

            else
                % there is an existing plot

                % recalculate Z
                Z = self.hydroSystem.waves.waveElevationGrid (t, self.wavePlotX, self.wavePlotY);

                % update the z data for the surface plot
                set ( self.waveSurfH, 'Zdata', Z );

                % drawnow () % don't know if we need a drawnow call here
            end

            if nargout > 0
                hsurf = self.waveSurfH;
            end

        end

    end

    methods
        % getter/setter methods go here

        function set.loggingSettings (self, newsettings)

            assert (isa (newsettings, 'wsim.loggingSettings'), ...
                'loggingSettings property must be an wsim.loggingSettings object');

            self.loggingSettings = newsettings;

        end

        function set.wavePlotNumPointsX (self, n)

            check.isPositiveScalarInteger (n, true, 'wavePlotNumPointsX', true);

            self.wavePlotNumPointsX = n;

        end

        function set.wavePlotNumPointsY (self, n)

            check.isPositiveScalarInteger (n, true, 'wavePlotNumPointsY', true);

            self.wavePlotNumPointsY = n;

        end

    end

    methods (Access = private)

        function [pos, vel, accel] = getMotion (self, mb)

            eul = zeros (3, self.nMBDynNodes);

            R = mb.GetRot();
            for Rind = 1:size (R,3)
%                 om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
%                 eul(1:3,Rind) = om.euler123 ();

                eul(1:3,Rind) = self.euler123(R(:,:,Rind));

            end

            pos = [ mb.NodePositions();
                    eul ];

%             pos = [ mb.NodePositions();
%                     mb.GetRot() ];

            vel = [ mb.NodeVelocities();
                    mb.NodeOmegas() ];

            accel = [ mb.NodeAccelerations();
                      mb.NodeAngularAccels() ];

        end

        function eul = euler123 (self, R)
            % returns the extrinsic euler123 angles corresponding to the
            % orientation matrix

            alpha = -atan2 ( R(2,3) , R(3,3) );

            beta = atan2 ( R(1,3) ...
                            , ( cos (alpha)*R(3,3) - sin (alpha)*R(2,3) ) );

            gamma = atan2 ( ( cos (alpha)*R(2,1) - sin (alpha)*R(3,1) ) ...
                           , ( cos (alpha)*R(2,2) - sin (alpha)*R(3,2) ) );

            eul = [ alpha;
                    beta;
                    gamma ];

        end

        function forces_and_moments = applyPTOForces (self, forces_and_moments)

            % get the PTO forces and moments and put them in the corect
            % places in the forces_and_moments matrix, based on an index
            % mapping done by the method mapPTOForceInds(), which is called
            % by the prepare method before simulation can be started.
            for ptoind = 1:size (self.ptoIndexMap, 1)

                % Each PTO has a reference node and another node, the PTO
                % forces/moments are applied between the two nodes. The
                % returned PTO force is applied to the reference node, and
                % the opposite force is applied to the other node

                ptoForceAndTorque = self.powerTakeOffs{ptoind}.forceAndMoment (self.lastTime);

                forces_and_moments (:,self.ptoIndexMap(ptoind,1)) = ...
                    forces_and_moments (:,self.ptoIndexMap(ptoind,1)) + ptoForceAndTorque(:,1);

                forces_and_moments (:,self.ptoIndexMap(ptoind,2)) = ...
                    forces_and_moments (:,self.ptoIndexMap(ptoind,2)) + ptoForceAndTorque(:,2);

            end


        end

        function startPTOs (self)

            for ptoind = 1:size (self.ptoIndexMap, 1)

                self.powerTakeOffs{ptoind}.start (self.simInfo);

            end

        end


        function forces_and_moments = applyHydroForces (self, forces_and_moments, time, pos, vel, accel)

            if self.hydroMotionSyncStepCount == self.hydroMotionSyncStep

                % calculate hydrodynamic forces based on the motion
                [hydro_forces_and_moments, out] = self.hydroSystem.hydroForces ( ...
                                                    time, ...
                                                    pos(:,self.hydroNodeIndexMap), ...
                                                    vel(:,self.hydroNodeIndexMap), ...
                                                    accel(:,self.hydroNodeIndexMap) ...
                                                                               );

                % store the last generated forces for later use (e.g. logging)
                self.lastForceHydro = hydro_forces_and_moments(1:3,:);
                self.lastForceExcitation = out.F_Excit(1:3,:);
                self.lastForceExcitationRamp = out.F_ExcitRamp(1:3,:);
                self.lastForceExcitationLin = out.F_ExcitLin(1:3,:);
                self.lastForceExcitationNonLin = out.F_ExcitNonLin(1:3,:);
                self.lastForceRadiationDamping = out.F_RadiationDamping(1:3,:);
                self.lastForceRestoring = out.F_Restoring(1:3,:);
                self.lastForceMorrison = out.F_MorrisonElement(1:3,:);
                self.lastForceViscousDamping = out.F_ViscousDamping(1:3,:);
                self.lastForceAddedMassUncorrected = out.F_AddedMass(1:3,:);

                self.lastMomentHydro = hydro_forces_and_moments(4:6,:);
                self.lastMomentExcitation = out.F_Excit(4:6,:);
                self.lastMomentExcitationRamp = out.F_ExcitRamp(4:6,:);
                self.lastMomentExcitationLin = out.F_ExcitLin(4:6,:);
                self.lastMomentExcitationNonLin = out.F_ExcitNonLin(4:6,:);
                self.lastMomentRadiationDamping = out.F_RadiationDamping(4:6,:);
                self.lastMomentRestoring = out.F_Restoring(4:6,:);
                self.lastMomentMorrison = out.F_MorrisonElement(4:6,:);
                self.lastMomentViscousDamping = out.F_ViscousDamping(4:6,:);
                self.lastMomentAddedMassUncorrected = out.F_AddedMass(4:6,:);
                
                % calculate the 'real' added mass forces, as they will be
                % reported at the end of the simulation
                F_added_mass = dynamicRealAddedMassForce (self.hydroSystem, accel(:,self.hydroNodeIndexMap));

                self.lastForceAddedMass = F_added_mass(1:3,:);
                self.lastMomentAddedMass = F_added_mass(4:6,:);
                
                % reset step count
                self.hydroMotionSyncStepCount = 1;

            else

                % just return the last set of hydro forces calculated
                hydro_forces_and_moments = [ self.lastForceHydro; self.lastMomentHydro ];

                self.hydroMotionSyncStepCount = self.hydroMotionSyncStepCount + 1;

            end

            % add hydrodynamic forces to the correct nodes (which are the
            % nodes attached to bodies with hydrodynamic interaction). This
            % is done using a mapping of the indexes created by the
            % mapHydroForceInds method, which is called by the 'prepare'
            % method
            for ind = 1:size (self.hydroNodeIndexMap, 1)
                forces_and_moments (:,self.hydroNodeIndexMap(ind,1)) = ...
                    forces_and_moments (:,self.hydroNodeIndexMap(ind,1)) + hydro_forces_and_moments(:,ind);
            end

        end

        function mapPTOForceInds (self)
            % gets the node indices for application of power take-off forces
            %
            % Syntax
            %
            % wsim.mapPTOForceInds (ws)
            %
            % Description
            %
            % mapPTOForceInds creates a matrix whose values give the
            % appropriate indices of the force and moment matrix for the
            % whole system to which to add the forces from each PTO. Not
            % every node in the system is necessarily accessible via the
            % interface, only those associated with the MBDyn external
            % structural forces element, so this is actually a mapping of
            % the index of the structural external force nodes to the PTO
            % nodes. The mapping is done by comparing the mbdyn structural
            % nodes associated with the PTO objects with all the nodes in
            % the external sturctural force element.
            %
            % mapPTOForceInds updates the ptoIndexMap of the wsim.wecSim
            % object to contain a two column matrix. Each row corresponds
            % to each of the PTO objects. The first column is the index of
            % the column of the force_and_moment matrix corresponding to
            % the reference node of each PTO, the second column is the
            % index of the column of the force_and_moment matrix
            % corresponding to the other node of each PTO.
            %
            % Input
            %
            %  ws - wsim.wecSim object
            %
            %

            strinfo = self.mBDynSystem.externalStructuralInfo ();

            % preallocate the index map for all the nodes
            self.ptoIndexMap = zeros (numel (self.powerTakeOffs), 2);

            for ptoind = 1:numel(self.powerTakeOffs)
                % loop through the power take-off objects, for each one
                % match up the PTO nodes with the nodes in the whole
                % system. Each PTO has a reference node and another node,
                % the PTO force/moment are applied between the two nodes.
                % The returned PTO force is applied to the reference node,
                % and the opposite force is applied to the other node
                %

                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)

                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.referenceNode

                        foundptonode = true;

                        self.ptoIndexMap(ptoind, 1) = nodeind;

                        break;

                    end

                end

                if foundptonode == false
                    error ('Could not find reference node for PTO %d in MBDyn system structural external force element', ptoind);
                end

                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)

                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.otherNode

                        foundptonode = true;

                        self.ptoIndexMap(ptoind, 2) = nodeind;

                        break;

                    end

                end

                if foundptonode == false
                    error ('Could not find non-reference node for PTO %d in MBDyn system structural external force element', ptoind);
                end

            end

        end

        function mapHydroForceInds (self)
            % gets the node indices for application of hydrodynamic forces
            %
            % Syntax
            %
            % wsim.mapHydroForceInds (ws)
            %
            % Description
            %
            % mapPTOForceInds creates a matrix whose values give the
            % appropriate indices of the force and moment matrix for the
            % whole system to which to add the forces due to hydrodynamic
            % interaction. Not every node in the system is necessarily
            % accessible via the interface, only those associated with the
            % MBDyn external structural forces element, so this is actually
            % a mapping of the index of the structural external force nodes
            % to the hydro body nodes. The mapping is done by comparing the
            % mbdyn structural nodes associated with the hydrobody objects
            % with all the nodes in the external structural force element.
            %
            % mapHydroForceInds updates the hydroNodeIndexMap of the
            % wsim.wecSim object to contain a column vector containing the
            % index of the column of the force_and_moment matrix
            % corresponding to each hydrobody in the system.
            %
            % Input
            %
            %  ws - wsim.wecSim object
            %


            strinfo = self.mBDynSystem.externalStructuralInfo ();

            self.hydroNodeIndexMap = zeros (numel (self.hydroSystem.bodyMBDynNodes), 1);

            for hydronodeind = 1:numel(self.hydroSystem.bodyMBDynNodes)

                foundnode = false;
                for nodeind = 1:numel (strinfo.Nodes)

                    if strinfo.Nodes{nodeind} == self.hydroSystem.bodyMBDynNodes{hydronodeind}

                        foundnode = true;

                        self.hydroNodeIndexMap(hydronodeind, 1) = nodeind;

                        break;

                    end

                end

                if foundnode == false
                    error ('Could not find corresponding external structural element node for hydro node %d', hydronodeind);
                end

            end

        end


        function initDataStructures (self)
            % intialise logging data structures

            extforceinfo = self.mBDynSystem.externalStructuralInfo ();

            % make a clean logger object to store simulation data
            self.logger = wsim.logger ();

            % preallocate data vectors for the predicted number of time
             % steps
            nsteps = ceil ((self.hydroSystem.simu.endTime - self.hydroSystem.simu.startTime) ./ self.hydroSystem.simu.dt + 1);

            if self.loggingSettings.windowed == true
                nsteps = min ([self.loggingSettings.windowSize, nsteps]);
            end

            % always store time, as is is independent variable for other
            % vars in logger object
            self.logger.addVariable ( 'Time', [1, 1], ...
                                      'Desc', 'main hydrodynamic/multibody time step', ...
                                      'AxisLabel', 'Time [s]', ...
                                      'PreallocateStorage', nsteps, ...
                                      'Windowed', self.loggingSettings.windowed );

            if self.loggingSettings.positions

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['X node ', int2str(noden)]; ...
                                          ['Y node ', int2str(noden)]; ...
                                          ['Z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'Positions', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian positions of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Positions [m]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs);
            end

            if self.loggingSettings.angularPositions

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['\theta_1 node ', int2str(noden)]; ...
                                          ['\theta_2 node ', int2str(noden)]; ...
                                          ['\theta_3 node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'AngularPositions', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular positions of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Euler Angles [rad]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.velocities

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['XP node ', int2str(noden)]; ...
                                          ['YP node ', int2str(noden)]; ...
                                          ['ZP node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'Velocities', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian velocities of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Velocities [ms^{-1}]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.angularVelocities

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['\omega_1 node ', int2str(noden)]; ...
                                          ['\omega_2 node ', int2str(noden)]; ...
                                          ['\omega_3 node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'AngularVelocities', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular velocities of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Angular Velocities [rads^{-1}]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.accelerations ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['XPP node ', int2str(noden)]; ...
                                          ['YPP node ', int2str(noden)]; ...
                                          ['ZPP node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'Accelerations', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian accelerations of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Accelerations [ms^{-2}]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.angularAccelerations ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['\omegaP_1 node ', int2str(noden)]; ...
                                          ['\omegaP_2 node ', int2str(noden)]; ...
                                          ['\omegaP_3 node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'AngularAccelerations', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular accelerations of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Angular Accelerations [rad s^{-2}]', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end


            % all the hydro forces


            if self.loggingSettings.nodeForces || self.loggingSettings.forceAddedMass

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'NodeForces', [3, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces for all nodes with external forces', ...
                                          'AxisLabel', 'Total Forces [N] on External Struct Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.nodeForcesUncorrected ...
                    || self.loggingSettings.forceAddedMass ...
                    || self.loggingSettings.momentAddedMass ...
                    || self.loggingSettings.nodeMoments ...
                    || self.loggingSettings.momentAddedMass ...
                    || self.loggingSettings.nodeForces

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'NodeForcesUncorrected', [3, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces (with uncorrected added mass forces) for all external structural nodes with external forces', ...
                                          'AxisLabel', 'Uncorrected Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end


            if self.loggingSettings.forceHydro ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceHydro', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of all hydrodynamic forces for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceExcitation

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceExcitation', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Excitation Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceExcitationRamp

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceExcitationRamp', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes, but with a ramp function applied', ...
                                          'AxisLabel', 'Ramped Hydro Excitation Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceExcitationLin

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceExcitationLin', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'linear hydrodynamic excitation forces for all hydro nodes', ...
                                          'AxisLabel', 'Linear Hydro Excitation Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceExcitationNonLin

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceExcitationNonLin', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'nonlinear hydrodynamic excitation forces for all hydro nodes', ...
                                          'AxisLabel', 'Nonlinear Hydro Excitation Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceRadiationDamping

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceRadiationDamping', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic radiation and damping forces for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Radiation Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceRestoring

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceRestoring', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic restoring forces for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Restoring Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceMorrison

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceMorrison', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'morrison forces for all hydro nodes', ...
                                          'AxisLabel', 'Morrison Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceViscousDamping

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceViscousDamping', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic viscous damping forces for all hydro nodes', ...
                                          'AxisLabel', 'Viscous Damping Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceAddedMass

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceAddedMass', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic added mass forces for all hydro nodes', ...
                                          'AxisLabel', 'Added Mass Forces [N] on Nodes ', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end

            if self.loggingSettings.forceAddedMassUncorrected

                legstrs = {};
                for noden = 1:extforceinfo.NNodes
                    legstrs = [legstrs, { ['F_x node ', int2str(noden)]; ...
                                          ['F_y node ', int2str(noden)]; ...
                                          ['F_z node ', int2str(noden)] }];
                end

                self.logger.addVariable ( 'ForceAddedMassUncorrected', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'uncorrected hydrodynamic added mass forces for all hydro nodes', ...
                                          'AxisLabel', 'Uncorrected Added Mass Forces [N] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3, ...
                                          'Legends', legstrs );
            end



            % all the hydro moments


            if self.loggingSettings.nodeMoments || self.loggingSettings.momentAddedMass
                self.logger.addVariable ( 'NodeMoments', [3, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all moments for all nodes with external moments', ...
                                          'AxisLabel', 'Total Moments [Nm] on External Struct Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.nodeMomentsUncorrected ...
                    || self.loggingSettings.momentAddedMass ...
                    || self.loggingSettings.nodeMoments ...
                    || self.loggingSettings.momentAddedMass ...
                    || self.loggingSettings.nodeForces ...
                    || self.loggingSettings.forceAddedMass

                self.logger.addVariable ( 'NodeMomentsUncorrected', [3, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all moments (with uncorrected added mass moments) for all external structural nodes with external moments', ...
                                          'AxisLabel', 'Uncorrected Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end


            if self.loggingSettings.momentHydro ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass
                self.logger.addVariable ( 'MomentHydro', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of all hydrodynamic moments for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentExcitation
                self.logger.addVariable ( 'MomentExcitation', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation moments for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Excitation Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentExcitationRamp
                self.logger.addVariable ( 'MomentExcitationRamp', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation moments for all hydro nodes, but with a ramp function applied', ...
                                          'AxisLabel', 'Ramped Hydro Excitation Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentExcitationLin
                self.logger.addVariable ( 'MomentExcitationLin', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'linear hydrodynamic excitation moments for all hydro nodes', ...
                                          'AxisLabel', 'Linear Hydro Excitation Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentExcitationNonLin
                self.logger.addVariable ( 'MomentExcitationNonLin', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'nonlinear hydrodynamic excitation moments for all hydro nodes', ...
                                          'AxisLabel', 'Nonlinear Hydro Excitation Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentRadiationDamping
                self.logger.addVariable ( 'MomentRadiationDamping', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic radiation and damping moments for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Radiation Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentRestoring
                self.logger.addVariable ( 'MomentRestoring', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic restoring moments for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Restoring Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentMorrison
                self.logger.addVariable ( 'MomentMorrison', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'morrison moments for all hydro nodes', ...
                                          'AxisLabel', 'Morrison Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentViscousDamping
                self.logger.addVariable ( 'MomentViscousDamping', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic viscous damping moments for all hydro nodes', ...
                                          'AxisLabel', 'Viscous Damping Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentAddedMass
                self.logger.addVariable ( 'MomentAddedMass', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic added mass moments for all hydro nodes', ...
                                          'AxisLabel', 'Added Mass Moments [Nm] on Nodes ', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.momentAddedMassUncorrected
                self.logger.addVariable ( 'MomentAddedMassUncorrected', [3, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'uncorrected hydrodynamic added mass moments for all hydro nodes', ...
                                          'AxisLabel', 'Uncorrected Added Mass Moments [Nm] on Nodes', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time', ...
                                          'ForceLogDimension', 3 );
            end

            if self.loggingSettings.forceIterations
                self.logger.addVariable ( 'ForceIterations', [1, 1], ...
                                          'Desc', 'Number of iterations performed for force convergence on each time step', ...
                                          'AxisLabel', 'Force Iterations', ...
                                          'PreallocateStorage', nsteps, ...
                                          'Windowed', self.loggingSettings.windowed, ...
                                          'Indep', 'Time' );
            end

            % initialise PTO internal logging

            if self.loggingSettings.powerTakeOffInternal

                for ptoind = 1:numel (self.powerTakeOffs)

                    self.powerTakeOffs{ptoind}.loggingSetup (self.logger);

                end

            else

                for ptoind = 1:numel (self.powerTakeOffs)

                    self.powerTakeOffs{ptoind}.loggingOn = false;

                end

            end

        end


        function logData (self)
            % log data after a time step has completed

            % always store time, as is is independent variable for other
            % vars in logger object
            self.logger.logVal ( 'Time', self.lastTime, false, false );

            if self.loggingSettings.positions
                self.logger.logVal ( 'Positions', self.lastPositions, false, false );
            end

            if self.loggingSettings.angularPositions
                self.logger.logVal ( 'AngularPositions', self.lastAngularPositions, false, false );
            end

            if self.loggingSettings.velocities
                self.logger.logVal ( 'Velocities', self.lastVelocities, false, false );
            end

            if self.loggingSettings.angularVelocities
                self.logger.logVal ( 'AngularVelocities', self.lastAngularVelocities, false, false );
            end

            if self.loggingSettings.angularAccelerations ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                self.logger.logVal ( 'AngularAccelerations', self.lastAngularAccelerations, false, false );

            end

            if self.loggingSettings.accelerations ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                self.logger.logVal ( 'Accelerations', self.lastAccelerations, false, false );

            end

            % log the forces

            if self.loggingSettings.nodeForcesUncorrected ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass

                self.logger.logVal ( 'NodeForcesUncorrected', self.lastNodeForcesUncorrected, false, false );
            end

            if self.loggingSettings.forceHydro ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass
                self.logger.logVal ( 'ForceHydro', self.lastForceHydro, false, false );
            end

            if self.loggingSettings.forceExcitation
                self.logger.logVal ( 'ForceExcitation', self.lastForceExcitation, false, false );
            end

            if self.loggingSettings.forceExcitationRamp
                self.logger.logVal ( 'ForceExcitationRamp', self.lastForceExcitationRamp, false, false );
            end

            if self.loggingSettings.forceExcitationLin
                self.logger.logVal ( 'ForceExcitationLin', self.lastForceExcitationLin, false, false );
            end

            if self.loggingSettings.forceExcitationNonLin
                self.logger.logVal ( 'ForceExcitationNonLin', self.lastForceExcitationNonLin, false, false );
            end

            if self.loggingSettings.forceRadiationDamping
                self.logger.logVal ( 'ForceRadiationDamping', self.lastForceRadiationDamping, false, false );
            end

            if self.loggingSettings.forceRestoring
                self.logger.logVal ( 'ForceRestoring', self.lastForceRestoring, false, false );
            end

            if self.loggingSettings.forceMorrison
                self.logger.logVal ( 'ForceMorrison', self.lastForceMorrison, false, false );
            end

            if self.loggingSettings.forceViscousDamping
                self.logger.logVal ( 'ForceViscousDamping', self.lastForceViscousDamping, false, false );
            end

            if self.loggingSettings.forceAddedMass || self.loggingSettings.forceAddedMassUncorrected
                self.logger.logVal ( 'ForceAddedMassUncorrected', self.lastForceAddedMassUncorrected, false, false );
            end

            % log the moments

            if self.loggingSettings.nodeMomentsUncorrected ...
                    || self.loggingSettings.nodeForces ...
                    || self.loggingSettings.forceAddedMass ...
                    || self.loggingSettings.nodeMoments ...
                    || self.loggingSettings.momentAddedMass

                self.logger.logVal ( 'NodeMomentsUncorrected', self.lastNodeMomentsUncorrected, false, false );

            end

            if self.loggingSettings.momentHydro ...
                || self.loggingSettings.nodeForces ...
                || self.loggingSettings.forceAddedMass ...
                || self.loggingSettings.nodeMoments ...
                || self.loggingSettings.momentAddedMass
                self.logger.logVal ( 'MomentHydro', self.lastMomentHydro, false, false );
            end

            if self.loggingSettings.momentExcitation
                self.logger.logVal ( 'MomentExcitation', self.lastMomentExcitation, false, false );
            end

            if self.loggingSettings.momentExcitationRamp
                self.logger.logVal ( 'MomentExcitationRamp', self.lastMomentExcitationRamp, false, false );
            end

            if self.loggingSettings.momentExcitationLin
                self.logger.logVal ( 'MomentExcitationLin', self.lastMomentExcitationLin, false, false );
            end

            if self.loggingSettings.momentExcitationNonLin
                self.logger.logVal ( 'MomentExcitationNonLin', self.lastMomentExcitationNonLin, false, false );
            end

            if self.loggingSettings.momentRadiationDamping
                self.logger.logVal ( 'MomentRadiationDamping', self.lastMomentRadiationDamping, false, false );
            end

            if self.loggingSettings.momentRestoring
                self.logger.logVal ( 'MomentRestoring', self.lastMomentRestoring, false, false );
            end

            if self.loggingSettings.momentMorrison
                self.logger.logVal ( 'MomentMorrison', self.lastMomentMorrison, false, false );
            end

            if self.loggingSettings.momentViscousDamping
                self.logger.logVal ( 'MomentViscousDamping', self.lastMomentViscousDamping, false, false );
            end

            if self.loggingSettings.momentAddedMass || self.loggingSettings.momentAddedMassUncorrected
                self.logger.logVal ( 'MomentAddedMassUncorrected', self.lastMomentAddedMassUncorrected, false, false );
            end

            if self.loggingSettings.forceIterations
                self.logger.logVal ( 'ForceIterations', self.lastForceIterations, false, false );
            end


        end

        function advanceStep (self)

            self.logData ();

            % accept the last data into the time history of solutions
            % for the hydrodynamic system and advance
            self.hydroSystem.advanceStep ( self.lastTime, ...
                                           [ self.lastVelocities(:,self.hydroNodeIndexMap(:,1)); ...
                                             self.lastAngularVelocities(:,self.hydroNodeIndexMap(:,1)) ], ...
                                           [ self.lastAccelerations(:,self.hydroNodeIndexMap(:,1)); ...
                                             self.lastAngularAccelerations(:,self.hydroNodeIndexMap(:,1)) ] );

            for ptoind = 1:numel (self.powerTakeOffs)
                % PTOs handle their own data logging which is triggered by
                % calling advanceStep on each PTO object
                self.powerTakeOffs{ptoind}.advanceStep (self.lastTime);

            end

            % advance the controller if there is one
            if ~isempty (self.wecController)
                self.wecController.advanceStep ();
            end

        end

        function finishPTOs (self)
            % tell all the PTOs that the simulation is complete

            for ptoind = 1:numel (self.powerTakeOffs)

                self.powerTakeOffs{ptoind}.finish (self.lastTime);

            end

        end

        function displayLastNLinesOfFile (self, filename, nlines)
            % display the last n lines of a text file on the command line

            lines = cell (nlines,1);

            fid = fopen (filename);

            if fid ~= -1

                CC = onCleanup (@() fclose (fid));

                line1ind = nlines + 1;

                while ~feof(fid)

                    line1ind = max (1, line1ind - 1);

                    lines = circshift (lines, -1, 1);

                    lines{end,1} = fgetl(fid);

                end

                % only print the lines we actually read in if there's less
                % than the max allowed
                for ind = line1ind:nlines
                    fprintf (1, '%s\n', lines{ind});
                end

            end

        end

        function makeAxes (self)
            % create axes and transform object

            figure;
            self.drawAxesH = axes;

        end

        function checkAxes (self, hax)
            % checks if there is a valid set of axes to plot to, and if
            % not, creates them

            % try to figure out if there is a valid set of axes to plot to
            if isempty (hax)

                % plot in the existing axes if possible as no new axes
                % supplied
                if mbdyn.pre.base.isAxesHandle (self.drawAxesH)
                    % use existing axes

                    if ~ishghandle (self.drawAxesH)
                        % axes no longer exists, or figure has been closed
                        self.drawAxesH = [];
%                         self.deleteAllDrawnObjects ();
                        % make a new set of axes to draw to
                        self.makeAxes ();
                    end

                elseif isempty (self.drawAxesH)
                    % make figure and axes
                    self.makeAxes ();
%                     self.needsRedraw = true;

                else
                    error ('drawAxesH property is not empty or an axes handle');
                end

            elseif mbdyn.pre.base.isAxesHandle (hax)
                % an axes has been supplied, so we plot to this new axes

                if isoctave
                    drawnow ();
                end

                if ~ishghandle (hax)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = hax;

%                 self.needsRedraw = true;

            else

                error ('hax must be an axes handle or empty');

            end

        end

        function cleanUpAfterError (self)

            if isempty (self.mBDynMBCNodal)
                mbdynpid = [];
                sockpath = [];
            else
                mbdynpid = self.mBDynMBCNodal.mBDynPID;
                sockpath = self.mBDynMBCNodal.path;
            end

            % clear mBDynMBCNodal to trigger the delete method on the
            % object and close sockets etc
            self.mBDynMBCNodal = [];

            if isunix

                % kill mbdyn process if we know the pid
                if ~isempty (mbdynpid)

                    cleansystem (sprintf ('kill %d', mbdynpid));

                end

                % remove local socket file if it exists
                if ~isempty (sockpath)
                    if exist (sockpath, 'file')
                        delete (sockpath);
                    end
                end

            end

        end

    end

end
