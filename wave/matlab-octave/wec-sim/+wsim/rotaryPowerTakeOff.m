classdef rotaryPowerTakeOff < wsim.powerTakeOff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = private, SetAccess = private)
        
        mbdynMomentObj;
        joint;

    end
    
    methods
        
        function self = rotaryPowerTakeOff (revolute_hinge, varargin)
            
            options.InitialThetaZero = true;
            options.TorqueFcn = [];
            options.ReferenceNode = 1;
            options.LoggedVars = {};
            
            options = parse_pv_pairs (options, varargin);
            
            momobj = mbdyn.mint.twoNodeTorque ( ...
                                    revolute_hinge, ...
                                    'ReferenceNode', options.ReferenceNode, ...
                                    'InitialThetaZero', options.InitialThetaZero, ...
                                    'TorqueFcn', options.TorqueFcn );
                                
            lginfo.AvailableNames = { 'InternalMoment', ...
                                      'RelativeAngularDisplacement', ...
                                      'RelativeAngularVelocity' ...
                                    };
                              
            lginfo.IndepVars = { 'Time', ...
                                 'Time', ...
                                 'Time' };
                              
            lginfo.Sizes = { [1,1], [1,1], [1,1] };
            
            lginfo.Descriptions = {'', '', ''};
            
            lginfo.NAvailable = numel(info.AvailableNames);
            
            
            self = self@wsim.powerTakeOff ( momobj.nodes(1), ...
                                            momobj.nodes(2), ...
                                            lginfo, ...
                                            'LoggedVars', options.LoggedVars );
            
            self.joint = revolute_hinge;
            self.mbdynMomentObj = momobj;
            
            % set up internal logging variables
            self.internalVariables.LastInternalMoment = [];
            self.internalVariables.LastRelativeAngularDisplacement = [];
            self.internalVariables.LastRelativeAngularVelocity = [];
            
        end
        
        function [FM, ptotorque, reltheta, relomega] = forceAndMoment (self)
            
            [FM, ptotorque, reltheta, relomega] = self.mbdynForceObj.momentFromFcn ();
            
            % need to add zero forces to forces
            FM = [ zeros(size (FM));
                  T ];
              
            self.internalVariables.LastInternalMoment = ptotorque;
            self.internalVariables.LastRelativeAngularDisplacement = reltheta;
            self.internalVariables.LastRelativeAngularVelocity = relomega;
            
        end
         
    end
    
end

