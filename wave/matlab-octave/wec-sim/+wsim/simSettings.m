classdef simSettings < handle
% Class which handles simulation parameters and settings for a WEC simulation
%
% Syntax
%
% ssobj = wsim.simSettings (casedir)
% ssobj = wsim.simSettings (..., 'Parameter', Value)
%
% Description
%
% wsim.simSettings handles simulation parameters and settings for a WEC
% simulation.
%
% wsim.simSettings Methods:
%
%   simSettings - wsim.simSettings constructor
%   checkInputs - Validate user input for simulation
%   listInfo - Prints simulation info to the command line
%   rhoDensitySetup - Assigns density and gravity values
%   setupSim - Sets simulation properties based on values specified in input file
%
%
% See Also: wsim.waveSettings
%


%
% Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation
% Copyright 2017-2018 The University of Edinburgh
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed und   er the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties (SetAccess = 'public', GetAccess = 'public')%input file
        
        multibodySolver = 'MBDyn'; % solver to use for mulitbody dynamics, can be SimMechanics or MBDyn
        startTime = 0; % Simulation start time (default = 0 s)
        endTime = 500; % Simulation end time (default = 500 s)
        dt = 0.1;  % Simulation time step (default = 0.1 s)
        dtMax = []; % Maximum simulation time step for variable step (default = 0.1 s) 
        dtOut = []; % Output sampling time (default = dt)
        dtFeNonlin = []; % Sample time to calculate nonlinear forces (default = dt)
        dtCITime = []; % Sample time to calculate Convolution Integral (default = dt)
        rampT = 100; % Ramp time for wave forcing (default = 100 s)
        domainSize = 200; % Size of free surface and seabed. This variable is only used for visualization (default = 200 m)
        CITime = 60; % Convolution integral time (default = 60 s)
        ssCalc = 0; % Option for convolution integral or state-space calculation: convolution integral->'0', state-space->'1', (default = 0)
        mode = 'normal'; %'normal','accelerator','rapid-accelerator' (default = 'normal')
        autoRateTranBlk = 'on'; % Automatically handle rate transition for data transfer
        zeroCrossCont = 'DisableAll';  % Disable zero cross control 
        rho = 1000; % Density of water (default = 1000 kg/m^3)
        g = 9.81; % Acceleration due to gravity (default = 9.81 m/s)
        nlHydro = 0; % Option for nonlinear hydrohanamics calculation: linear->'0', nonlinear->'1', (default = 0)
        b2b = false; % Option for body2body interactions: off->false, on->true, (default = false)
        paraview = 0; % Option for writing vtp files for paraview visualization.
        adjMassWeightFun = 2; % Weighting function for adjusting added mass term in the translational direction (default = 2)
        numIntMidTimeSteps = 5; % Number of intermidiate time steps (default = 5 for ode4 method)
        mcrCaseFile = []; % mat file that contain a list of the multiple conditions runs with given conditions  
        morrisonElement = 0; % Option for Morrison Element calculation: Off->'0', On->'1', (default = 0)
%         outputtxt = 0; % Option to save results as ASCII files.
        reloadH5Data = 0; % Option to re-load hydro data from hf5 file between runs: Off->'0', On->'1', (default = 0)
        numWecBodies = []; % Number of hydrodynamic bodies that comprise the WEC device (default = [])
        numPtos = []; % Number of power take-off elements in the model (default = [])
        numConstraints = []; % Number of contraints in the wec model (default = [])
        numMoorings = [];% Number of moorings in the wec model (default = [])
        
    end

    properties (SetAccess = 'protected', GetAccess = 'public')%internal
        
        verbose = false; % flag determining whether to print information to the command line
        version;  % WEC-Sim version
        simulationDate = datestr (now ()); % Simulation date and time
        time = 0; % Simulation time [s] (default = 0 s)
        caseDir = []; % WEC-Sim case directory
        CIkt = []; % Number of timesteps in the convolution integral length
        maxIt = []; % Total number of simulation time steps (default = dependent)        CIkt                                                               % Calculate the number of convolution integral timesteps (default = dependent)
        CTTime = []; % Convolution integral time series (default = dependent)
        
    end

    methods
        
        function obj = simSettings (casedir, varargin)
            % wsim.simSettings constructor
            %
            % Syntax
            %
            % ssobj = wsim.simSettings (casedir)
            % ssobj = wsim.simSettings (..., 'Parameter', Value)
            %
            % Description
            %
            % wsim.simSettings creates a new wsim.simSettings object. 
            %
            % Input
            %
            %  casedir - root directory containing the simulation data
            %   files. This directory is expected to have two
            %   subdirectories 'hydroData', where hydrodynamic data fiels
            %   will be searched for, and 'geometry', where geometry files
            %   (e.g. stl files will be searched for.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Verbose' - true/false flag indicating whether to print
            %    information (about the simSettings) to the command line.
            %    Does not effect the level of output from any other parts
            %    of the simulation system. Default is false.
            %
            % Output
            %
            %  ssobj - wsim.simSettings object
            %
            %
            %
            % See Also: wsim.waveSettings
            %
            
            options.Verbose = false;
            
            options = parse_pv_pairs (options, varargin);
            
            check.isLogicalScalar (options.Verbose, true, 'Verbose');
            
            obj.verbose = options.Verbose;
            
            obj.getVersion ();
            
            if obj.verbose
                fprintf ('WEC-Sim: An open-source code for simulating wave energy converters\n')
                fprintf ('Version: %s\n\n', obj.version)
                fprintf ('Initializing the Simulation Class...\n')
            end
            
            if nargin < 1
                obj.caseDir = pwd ();
            else
                obj.caseDir = casedir;
            end
            
            if obj.verbose, fprintf('\tCase Dir: %s \n', obj.caseDir); end
            
            obj.checkInputs ();
            
        end
        

        function setupSim (obj)
            % Sets simulation properties based on values specified in input file
            
            obj.time = obj.startTime:obj.dt:obj.endTime;
            obj.maxIt = floor((obj.endTime - obj.startTime) / obj.dt);
            
            % Set dtOut if it was not specificed in input file
            if isempty(obj.dtOut) || obj.dtOut < obj.dt
                obj.dtOut = obj.dt;
            end
            
            % Set dtFeNonlin if it was not specificed in input file
            if isempty(obj.dtFeNonlin) || obj.dtFeNonlin < obj.dt
                obj.dtFeNonlin = obj.dt;
            end
            
            % Set dtCITime if it was not specificed in input file
            if isempty(obj.dtCITime) || obj.dtCITime < obj.dt
                obj.dtCITime = obj.dt;
            end
            
            % Set dtMax if it was not specificed in input file
            if isempty(obj.dtMax) || obj.dtMax < obj.dt
                obj.dtMax = obj.dt;
            end
            
            obj.CTTime = 0:obj.dtCITime:obj.CITime;            
            obj.CIkt = length(obj.CTTime);
            
        end
        

        function checkInputs (obj)
            % Validate user input for simulation
            
            if exist (obj.caseDir, 'dir') ~= 7
                error ('The WEC case directory, %s, does not appear to exist', obj.caseDir);
            end
            
            if exist (fullfile (obj.caseDir, 'hydroData'), 'dir') ~= 7
                error ('The case directory does not contain a ''hydroData'' subfolder.');
            end
            
            if exist (fullfile (obj.caseDir, 'geometry'), 'dir') ~= 7
                error ('The case directory does not contain a ''geometry'' subfolder.');
            end
            
        end
        

        function rhoDensitySetup (obj, rho, g)
            % Assigns density and gravity values
            
            obj.rho = rho;
            obj.g   = g;
            
        end
        

        function listInfo (obj, waveTypeNum)
            % Prints simulation info to the command line
            
            fprintf ('\nWEC-Sim Simulation Settings:\n');
            %fprintf('\tTime Marching Solver                 = Fourth-Order Runge-Kutta Formula \n')
            fprintf ('\tStart Time                     (sec) = %G\n',obj.startTime);
            fprintf ('\tEnd Time                       (sec) = %G\n',obj.endTime);
            fprintf ('\tTime Step Size                 (sec) = %G\n',obj.dt);
            fprintf ('\tRamp Function Time             (sec) = %G\n',obj.rampT);
            if waveTypeNum > 10
                fprintf ('\tConvolution Integral Interval  (sec) = %G\n',obj.CITime);
            end
            fprintf ('\tTotal Number of Time Step            = %u \n',obj.maxIt);
            
        end


    end % methods
    
    
    methods (Access = protected)
        
        function getVersion (obj)
            % Determines software version being used
            %
            
            wecsim_version_file = fullfile (wecsim_rootdir (), 'version.txt');
            
            if exist (wecsim_version_file, 'file') == 2

                obj.version = strtrim (urlread (['file:///' wecsim_version_file]));
                
            else
                
                origdir = pwd ();
                
                CC = onCleanup (@() cd (origdir));
                
                cd (wecsim_rootdir ());
                
                [status,result] = system (sprintf ('hg id', wecsim_rootdir));
                
                if status == 0
                    obj.version = sprintf ('development sources, hg id: %s', strtrim (result));
                else
                    obj.version = 'unable to determine version';
                end
            
            end
            
        end
        
    end
    
end % classdef
