classdef powerTakeOff < handle
% base class for power take-off simulation classes
%
% Syntax
%
% pto = wsim.powerTakeOff (reference_node, other_node, logginginfo)
% pto = wsim.powerTakeOff (..., 'Parameter', value)
%
% Description
%
% wsim.powerTakeOff is a base class for the power take-off simulation
% classes. It has methods to help with setting up data logging and checking
% some standard inputs.
%
% wsim.powerTakeOff Methods:
%
%   powerTakeOff - wsim.powerTakeOff constructor
%   advanceStep - advance to the next simulation time step
%   forceAndMoment - 
%   logData - appends the internal variable data to the log
%   loggingSetup - sets up data logging for a wsim.linearPowerTakeOff object
%
%
% See Also: wsim.linearPowerTakeOff, wsim.rotaryPowerTakeOff
%

    properties
        
        id; % positive scalar integer for uniquely idenifying the object
        
    end
    
    properties (GetAccess=public, SetAccess=private)
        uniqueLoggingNames;
        referenceNode;
        otherNode;
        loggingInfo;
    end
    
    properties (GetAccess=protected, SetAccess=protected)
        logger;
        loggerReady;
        internalVariables;
        simulationInfo;
    end
    
    methods (Abstract)
        
        forceAndMoment (self);
        
    end
    
    methods
        
        function self = powerTakeOff (reference_node, other_node, logginginfo, varargin)
            % wsim.powerTakeOff constructor
            %
            % Syntax
            %
            % pto = wsim.powerTakeOff (reference_node, other_node, logginginfo)
            % pto = wsim.powerTakeOff (..., 'Parameter', value)
            %
            % Description
            %
            % wsim.powerTakeOff is a base class for the power take-off
            % simulation classes. It has methods to help with setting up
            % data logging and checking some standard inputs. 
            %
            % Input
            %
            %  reference_node - mbdyn.pre.structuralNode object
            %
            %  other_node - mbdyn.pre.structuralNode object
            %
            %  logginginfo - scalar structure containing the following
            %   fields:
            %
            %   AvailableNames : cell array of variable names that can be
            %    logged. To use the logData method of the powerTakeOff
            %    class, each name must correspond to a class property,
            %    names "last<var name>", where <var name> is replaced by
            %    the contents of each cell in AvailableNames. When calling
            %    logData, it is the contents of the class property which is
            %    appended to the log.
            %
            %   IndepVars : cell array the same size as the AvailableNames
            %    field containingthe name of the independant variable
            %    associated with the corresponding variable name in
            %    AvailableNames, e.g. 'Time'. If there is no independant
            %    variable, the cell should containin an empty character
            %    vector (or just be empty).
            %
            %   Sizes : cell array containing vector containing the size of
            %    each of the the logged variables.
            %
            %   Descriptions : cell array of character vectors containing
            %    descriptions of the variables. Each cell can be empty, but
            %    the cell array must be of the same length as the
            %    AvailableNames cell array.
            %
            %  From the AvailableNames, a unique variable name is created
            %  for use in the logger object based on the pto id. Data is
            %  actually logged to this unique name in the logger. This
            %  ensures data for mutiple PTO objects of the same type can be
            %  logged in the same logger object. The unique logging names
            %  are stored in the uniqueLoggingNames property of the class.
            %  The variables which are actually logged can be restricted by
            %  using the LoggedVars option (see below)
            %
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'LoggedVars' - cell array of character vectors containing
            %    the names of internal variables to be logged. If supplied,
            %    only those variables named in this cell array will
            %    actually have data logged. Alternatively, this can be a
            %    character vector, 'none', in which case no internal
            %    variables will be logged. Default is an empty cell array,
            %    which mans all available nternal variables will be logged.
            %
            % Output
            %
            %  pto - wsim.powerTakeOff object
            %
            %
            %
            % See Also: wsim.linearPowerTakeOff, wsim.rotaryPowerTakeOff
            %

            options.LoggedVars = {};
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isa (reference_node, 'mbdyn.pre.structuralNode')
                error ('reference_node must be an mbdyn.pre.structuralNode')
            end
            
            if ~isa (other_node, 'mbdyn.pre.structuralNode') 
                error ('other_node must be an mbdyn.pre.structuralNode')
            end
            
            check.structHasAllFields ( logginginfo, {'AvailableNames', 'IndepVars', 'Sizes', 'Descriptions'}, true, 'logginginfo');

            
            if isempty (options.LoggedVars)
                % log everything
                logginginfo.LoggedVarInds = 1:logginginfo.NAvailable;
                logginginfo.LoggedVariables = logginginfo.AvailableNames;
                
            else
                
                if ischar (options.LoggedVars)
                    
                    if ~strcmpi (options.LoggedVars, 'none')
                        error ('LoggedVars must be a character vector ''none'', or a cell array of character vectors')
                    end
                    logginginfo.LoggedVarInds = [];
                    
                elseif iscellstr (options.LoggedVars)
                    
                    logginginfo.LoggedVarInds = [];
                    
                    for ind = 1:numel (options.LoggedVars)
                        
                        check.allowedStringInputs ( options.LoggedVars{ind}, ...
                                                    logginginfo.AvailableNames, ...
                                                    true, ...
                                                    sprintf('LoggedVars{%d} is ''%s'' but', ind, options.LoggedVars{ind}) ...
                                                  );
                                              
                        for avnind = 1:logginginfo.NAvailable
                            if strcmp (options.LoggedVars{ind}, logginginfo.AvailableNames{avnind})
                                logginginfo.LoggedVarInds = [ logginginfo.LoggedVarInds, avnind ];
                                break;
                            end
                        end
                                              
                    end
                else
                    error ('LoggedVars must be a character vector ''none'', or a cell array of character vectors')
                end
                
                logginginfo.LoggedVariables = options.LoggedVars;
                
            end
            
            self.referenceNode = reference_node;
            self.otherNode = other_node;
            self.loggingInfo = logginginfo;
            
            self.loggerReady = false;
            self.logger = [];
            
        end
        
        function set.id (self, newid)
            % set new integer id for the power take-off
            
            check.isPositiveScalarInteger (newid, true, 'id');
            
            self.id = newid;
        end
        
        function start (self, siminfo)
            % initialise the pto simulation
            % replace h1 line
            %
            % Syntax
            %
            % start (pto, siminfo)
            %
            % Description
            %
            % wsim.powerTakeOff.start is called by wecSim at the beginning
            % of a simulation to allow any initialisation to take place.
            %
            % Input
            %
            %  pto - wsim.powerTakeOff object
            %
            %  siminfo - structure containing information about the
            %   simulation the PTO is part of. It should contain the
            %   following fields:
            %
            %   TStart : Start time of the simulation
            %   
            %   TEnd : End time of the simulation
            %
            %   TStep : Time step size
            %
            %   MBDynSystem : the MBDyn object used in the simulation
            %
            %   HydroSystem : the wsim.hydroSystem object used in the
            %    simulation
            %
            %
            % See Also: wsim.powerTakeOff.advanceStep
            %
            
            % TODO: put in restricted class access block when Octave supports this
            
            self.simulationInfo = siminfo;

        end
        
        function advanceStep (self, varargin)
            % advance to the next simulation time step
            %
            % Syntax
            %
            % advanceStep (pto)
            %
            % Description
            %
            % wsim.powerTakeOff.advanceStep is intended to be called when
            % the simulation is ready to advance to the next time step. It
            % is always called when the powerTakeOff is used with the
            % wsim.wecSim class to manage a simulation. The only tasks done
            % by wsim.powerTakeOff.advanceStep is to log data to the
            % wsim.logger object (internally it calls
            % wsim.powerTakeOff.logData). It is intended that subclasses of
            % wsim.powerTakeOff can extend this method to perform other
            % tasks.
            %
            % Input
            %
            %  pto - wsim.powerTakeOff object
            %
            %
            % See Also: wsim.powerTakeOff.logData
            %
            
            % TODO: put in restricted class access block when Octave supports this
            
            self.logData ();
            
        end
        
        function loggingSetup (self, logger)
            % sets up data logging for a wsim.linearPowerTakeOff object
            %
            % Syntax
            %
            % info = loggingSetup (lpto)
            % info = loggingSetup (lpto, logger)
            %
            % Description
            %
            % loggingSetup initialises
            %
            % Input
            %
            %  lpto - wsim.powerTakeOff object
            %
            %  logger - optional wsim.logger object
            %
            % Output
            %
            %  info - structure containing information on the data to be
            %   logged by the PTO. Will contain the following fields:
            %   AvailableNames, IndepVars, Sizes, Descriptions, and
            %   NAvailable. 
            %
            %
            % See Also: wsim.powerTakeOff.logData,
            %           wsim.powerTakeOff.advanceStep
            %

            if nargin < 2
                logger = [];
            end

            self.initLogging (self.loggingInfo, logger);
                              
        end
        
        function logData (self)
            % appends the internal variable data to the log
            %
            % Syntax
            %
            % logData (pto)
            %
            % Description
            %
            % logData appends the last recorded values of the internal
            % variables to a logger object.
            %
            % Input
            %
            %  pto - wsim.powerTakeOff object
            %
            % See Also: wsim.powerTakeOff.loggingSetup,
            %           wsim.powerTakeOff.advanceStep
            %
            
            if self.loggerReady
                for ind = 1:numel(self.loggingInfo.LoggedVarInds)
                    
                    fieldname = self.loggingInfo.AvailableNames{self.loggingInfo.LoggedVarInds(ind)};
                    
                    self.logger.logVal ( self.uniqueLoggingNames{self.loggingInfo.LoggedVarInds(ind)}, ...
                                         self.internalVariables.(fieldname) ...
                                       );
                end
            else
                error ('You have called logData, but logging has not been set up, have you called loggingSetup yet?');
            end
            
        end
        
    end
    
    methods (Access=protected)
        
        function initLogging (self, info, logger)
            
            for ind = 1:numel (info.AvailableNames)
                % initialise internal variable storage
            	self.internalVariables.(info.AvailableNames{ind}) = [];
                
                info.UniqueLoggingNames{ind} = sprintf ('PTO_%d_%s', self.id, info.AvailableNames{ind});
                
            end
            
            self.uniqueLoggingNames = info.UniqueLoggingNames;
            
            if ~isempty (logger)
                
                assert (isa (logger, 'wsim.logger'), ...
                    'logger must be a wsim.logger object (or empty)');
            
                initLoggerObj (self, info, logger)
                
            end

            
        end
        
        function initLoggerObj (self, info, logger)
            
            for ind = 1:info.NAvailable
                
                % putting IndepVar in will automatically allocate enough
                % space for variable, same length as indepvar, if there is
                % no Indep var, preallocated length will be 1
                
                if ischar (info.IndepVars{ind})
                    indepvar = info.IndepVars{ind};
                elseif check.isScalarInteger (info.IndepVars{ind}, false)
                    indepvar = info.UniqueLoggingNames{info.IndepVars{ind}};
                else
                    error ('Invalid PTO IndepVar');
                end
                
                logger.addVariable ( info.UniqueLoggingNames{ind}, ...
                                     info.Sizes{ind}, ...
                                     'Desc', info.Descriptions{ind}, ...
                                     'Indep', indepvar );
                
            end
            
            self.logger = logger;
            
            self.loggerReady = true;
            
        end
        
        
    end
    
end