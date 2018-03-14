classdef linearPMMachinePTO < wsim.powerTakeOff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = private, SetAccess = private)
        
        mbdynForceObj;
        
        % internal logging variables
        lastInternalForce;
        lastRelativeDisplacement;
        lastRelativeVelocity;
        
    end
    
    methods
        
        function self = linearPMMachinePTO (reference_node, other_node, axisNum, varargin)
            % construct a wsim.linearPMMachinePTO object
            %
            %
            % Syntax
            %
            % lpto = wsim.linearPMMachinePTO (reference_node, other_node, axisNum)
            % lpto = wsim.linearPMMachinePTO (..., 'Parameter', value)
            %
            % Description
            %
            % wsim.linearPMMachinePTO is a class representing a linear
            % power-take-off mechanism in a wave energy converter. It
            % facilitates sending the correct forces to an MBDyn multibody
            % simulation. wsim.linearPMMachinePTO applies forces between
            % two MBDyn nodes based on their relative displacement. Forces
            % are applied based on the relative displacement and velocity
            % of the two nodes along axis 3 in the reference frame of the
            % first node. It is assumed that the nodes motion is
            % constrained appropriately by other MBDyn elements (e.g. 
            %
            % Input
            %
            %  reference_node - mbdyn.pre.structuralNode6dof object
            %
            %  other_node - 
            %
            %  axisNum - 
            %
            % Additional options my be supplied as parameter-value pairs.
            % The avaialable options are:
            %
            %  'InitialDisplacementZero' - optional true/false flag
            %    indicating whether the intial relative displacement (along
            %    axis 3 of the reference node) in the global frame should
            %    be taken as the reference point for displacement during
            %    the simulation, i.e. the PTO starts with an initial
            %    displacement of zero for the purposes of force calulation,
            %    and future displacement is measured relative to this
            %    initial position. If false, the raw position is used
            %    instead. Default is true if not supplied.
            %
            % Output
            %
            %  lpto - a wsim.linearPMMachinePTO
            %
            %
            %
            % See Also: 
            %

            options.InitialDisplacementZero = true;
            options.ForceFcn = [];
            
            options = parse_pv_pairs (options, varargin);
            
            lginfo.AvailableNames = { 'InternalForce', ...
                                    'RelativeDisplacement', ...
                                    'RelativeVelocity' ...
                                  };
                              
            lginfo.IndepVars = { 'Time', ...
                               'Time', ...
                               'Time' };
                              
            lginfo.Sizes = { [1,1], [1,1], [1,1] };
            
            lginfo.Descriptions = {'', '', ''};
            
            lginfo.NAvailable = numel(lginfo.AvailableNames);
            
            self = self@wsim.powerTakeOff (reference_node, other_node, lginfo);
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero, ...
                                    'ForceFcn', [] );
            
        end
        
        function [FM, ptoforce, reldisp, relvel] = forceAndMoment (self)
            
            [FM, ptoforce, reldisp, relvel] = self.mbdynForceObj.force ();
            
            % need to add zero moments to forces
            FM = [FM; zeros(size (FM))];
            
            self.lastInternalForce = ptoforce;
            self.lastRelativeDisplacement = reldisp;
            self.lastRelativeVelocity = relvel;
            
        end
        
        function info = loggingSetup (self, logger)

            if nargin < 2
                logger = [];
            end
            

            
            self.initLogging (info, logger);
                              
        end
        
        function logData (self)
            % appends the internal variable data to the log
            
            if self.loggerReady
                self.logger.logVal (self.uniqueLoggingNames{1}, self.lastInternalForce);
                self.logger.logVal (self.uniqueLoggingNames{2}, self.lastRelativeDisplacement);
                self.logger.logVal (self.uniqueLoggingNames{3}, self.lastRelativeVelocity);
            else
                error ('You have called logData, but logging has not been set up, have you called loggingSetup yet?');
            end
            
        end
        
    end
    
end