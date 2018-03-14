classdef rotaryPowerTakeOff < wsim.powerTakeOff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = private, SetAccess = private)
        mbdynMomentObj;
        joint;
        
        % logging variables
        lastInternalMoment;
        lastRelativeAngularDisplacement;
        lastRelativeAngularVelocity;
    end
    
    methods
        
        function self = rotaryPowerTakeOff (revolute_hinge, varargin)
            
            options.InitialThetaZero = true;
            options.TorqueFcn = [];
            options.ReferenceNode = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            momobj = mbdyn.mint.twoNodeTorque ( ...
                                    revolute_hinge, ...
                                    'ReferenceNode', options.ReferenceNode, ...
                                    'InitialThetaZero', options.InitialThetaZero, ...
                                    'TorqueFcn', options.TorqueFcn );
                                
            self = self@wsim.powerTakeOff ( momobj.nodes(1), ...
                                            momobj.nodes(2) );
            
            self.joint = revolute_hinge;
            self.mbdynMomentObj = momobj;
            
        end
        
        function [FM, ptotorque, reltheta, relomega] = forceAndMoment (self)
            
            [FM, ptotorque, reltheta, relomega] = self.mbdynForceObj.momentFromFcn ();
            
            % need to add zero forces to forces
            FM = [ zeros(size (FM));
                  T ];
              
            self.lastInternalMoment = ptotorque;
            self.lastRelativeAngularDisplacement = reltheta;
            self.lastRelativeAngularVelocity = relomega;
            
        end
        
        function n = forceSize (self)
            n = 1;
        end
        
        function info = loggingSetup (self, logger)
            
            if nargin < 2
                logger = [];
            end
            
            info.AvailableNames = { 'InternalMoment', ...
                                    'RelativeAngularDisplacement', ...
                                    'RelativeAngularVelocity' ...
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
                self.logger.logVal (self.uniqueLoggingNames{1}, self.lastInternalMoment);
                self.logger.logVal (self.uniqueLoggingNames{2}, self.lastRelativeAngularDisplacement);
                self.logger.logVal (self.uniqueLoggingNames{3}, self.lastRelativeAngularVelocity);
            else
                error ('You have called logData, but logging has not been set up, have you called loggingSetup yet?');
            end
            
        end
        
    end
    
end

