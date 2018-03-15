classdef linearPowerTakeOff < wsim.powerTakeOff
   
    properties (GetAccess = private, SetAccess = private)
        
        mbdynForceObj;
        
    end
    
    methods
        
        function self = linearPowerTakeOff (reference_node, other_node, axisNum, force_fcn, varargin)
            % construct a wsim.linearPowerTakeOff object
            %
            % Syntax
            %
            % lpto = wsim.linearPowerTakeOff (reference_node, other_node, axisNum)
            % lpto = wsim.linearPowerTakeOff (..., 'Parameter', value)
            %
            % Description
            %
            % wsim.linearPowerTakeOff is a class representing a linear
            % power-take-off mechanism in a wave energy converter. It
            % facilitates sending the correct forces to an MBDyn multibody
            % simulation. wsim.linearPowerTakeOff applies forces between
            % two MBDyn nodes based on their relative displacement. Forces
            % are applied based on the relative displacement and velocity
            % of the two nodes along axis 3 in the reference frame of the
            % first node. It is assumed that the nodes motion is
            % constrained appropriately by other MBDyn elements (e.g. a
            % prismatic joint).
            %
            % wsim.linearPowerTakeOff Methods:
            %
            %  wsim.linearPowerTakeOff - constructor
            %  forceAndMoment - returns pto forces and moments on the 
            %    attached nodes
            %
            % Input
            %
            %  reference_node - mbdyn.pre.structuralNode6dof object
            %
            %  other_node - mbdyn.pre.structuralNode6dof object
            %
            %  axisNum - axis in the frame of the reference node. Forces
            %   will be applied to the node in a direction parallel to this
            %   axis.
            %
            %  force_fcn - function handle or string to be used to
            %   calculate the force to be applied between the two nodes
            %   making up the PTO. force_fcn is a function which takes two
            %   arguments with the following signature:
            %
            %        force_value = myfcn (reldisp, relvel)
            %
            %   where reldisp is the relative displacement of the two
            %   nodes along the specified axis in axisNum in the
            %   reference frame of the reference node, and relvel is the
            %   relative velocity of the two nodes in the same frame.
            %   force_value is expected to be a scalar value, the value of
            %   the force acting on the reference node parallel to the
            %   axis in forceAxis in the frame of the reference node.
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
            %  lpto - a wsim.linearPowerTakeOff
            %
            %
            %
            % See Also: wsim.rotaryPowerTakeOff
            %

            options.InitialDisplacementZero = true;
            options.LoggedVars = {};
            
            options = parse_pv_pairs (options, varargin);
            
            
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
            
            
            self = self@wsim.powerTakeOff ( reference_node, other_node, info, ...
                                            'LoggedVars', options.LoggedVars );
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero, ...
                                    'ForceFcn', force_fcn );
                                
            self.internalVariables.LastInternalForce = [];
            self.internalVariables.LastRelativeDisplacement = [];
            self.internalVariables.LastRelativeVelocity = [];
                                
        end
        
        function [FM, ptoforce, reldisp, relvel] = forceAndMoment (self)
            
            [FM, ptoforce, reldisp, relvel] = self.mbdynForceObj.forceFromFcn ();
            
            % need to add zero moments to forces
            FM = [FM; zeros(size (FM))];
            
            self.internalVariables.LastInternalForce = ptoforce;
            self.internalVariables.LastRelativeDisplacement = reldisp;
            self.internalVariables.LastRelativeVelocity = relvel;
            
        end
        
    end
    
end