classdef linearPowerTakeOff < wsim.powertakeoff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = private, SetAccess = private)
        
        mbdynForceObj;
        
        % internal logging variables
        lastInternalForce;
        lastRelativeDisplacement;
        lastRelativeVelocity;
        
    end
    
    methods
        
        function self = linearPowerTakeOff (reference_node, other_node, axisNum, varargin)
            
            options.InitialDisplacementZero = true;
            options.ForceFcn = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@wsim.powertakeoff (reference_node, other_node);
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero, ...
                                    'ForceFcn', options.ForceFcn );
            
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
            
            info.AvailableNames = { 'InternalForce', ...
                                    'RelativeDisplacement', ...
                                    'RelativeVelocity' ...
                                  };
                              
            info.IndepVars = { 'Time', ...
                               'Time', ...
                               'Time' };
                              
            info.Sizes = { [1,1], [1,1], [1,1] };
            
            info.Descriptions = {'', '', ''};
            
            info.NAvailable = numel(info.AvailableNames);
            
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