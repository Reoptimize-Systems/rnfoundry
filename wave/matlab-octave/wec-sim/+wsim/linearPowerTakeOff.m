classdef linearPowerTakeOff < wsim.powertakeoff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = private, SetAccess = private)
        mbdynForceObj;
    end
    
    methods
        
        function self = linearPowerTakeOff (reference_node, other_node, axisNum)
            
            options.InitialDisplacementZero = true;
            options.ForceFcn = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero, ...
                                    'ForceFcn', options.ForceFcn );
            
        end
        
        function [F, ptoforce, reldisp, relvel] = forceAndTorque (self)
            
            [F, ptoforce, reldisp, relvel] = self.mbdynForceObj.force ();
            
            % need to add zero moments to forces
            F = [F; zeros(size (F))];
            
        end
        
        function n = forceSize (self)
            n = 1;
        end
        
        function info = loggingInfo (self)
            
            info.NAvailable = 3;
            
            info.AvailableNames = { 'InternalForce', ...
                                    'RelativeDisplacement', ...
                                    'RelativeVelocity' ...
                                  };
                              
            info.Sizes = { [1,1], [1,1], [1,1] };
                              
        end
        
    end
    
end