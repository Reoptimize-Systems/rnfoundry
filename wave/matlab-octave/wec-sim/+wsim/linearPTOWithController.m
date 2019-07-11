classdef linearPTOWithController < wsim.powerTakeOff
% class representing a linear power take-off for a wave energy converter
%
% Syntax
%
% lpto = wsim.linearPTOWithController (reference_node, other_node, axisNum, controller)
% lpto = wsim.linearPTOWithController (..., 'Parameter', value)
%
% Description
%
% wsim.linearPTOWithController is a class representing a linear power-take-off
% mechanism in a wave energy converter. It facilitates sending the correct
% forces to an MBDyn multibody simulation. wsim.linearPTOWithController applies
% forces between two MBDyn nodes based on their relative displacement.
% Forces are applied based on the relative displacement and velocity of the
% two nodes along axis 3 in the reference frame of the first node. It is
% assumed that the nodes motion is constrained appropriately by other MBDyn
% elements (e.g. a prismatic joint).
%
% wsim.linearPTOWithController Methods:
%
%   linearPTOWithController - construct a wsim.linearPTOWithController object
%   forceAndMoment - returns the forces and moments applied by the PTO in 3D
%
%   Inheirited Methods
%
%   advanceStep - advance to the next simulation time step
%   logData - appends the internal variable data to the log
%   loggingSetup - sets up data logging for a wsim.linearPTOWithController object
%
%
% See Also: wsim.rotaryPowerTakeOff, wsim.powerTakeOff
%

    properties (GetAccess = protected, SetAccess = protected)
        
        mbdynForceObj;
        controller;
        
    end
    
    methods
        
        function self = linearPTOWithController (reference_node, other_node, axisNum, controller, varargin)
            % construct a wsim.linearPTOWithController object
            %
            % Syntax
            %
            % lpto = wsim.linearPTOWithController (reference_node, other_node, axisNum, controller)
            % lpto = wsim.linearPTOWithController (..., 'Parameter', value)
            %
            % Description
            %
            % wsim.linearPTOWithController is a class representing a linear
            % power-take-off mechanism in a wave energy converter. It
            % facilitates sending the correct forces to an MBDyn multibody
            % simulation. wsim.linearPTOWithController applies forces
            % between two MBDyn nodes based on a wsim.wecController's
            % output pto forces. Forces are applied based on the relative
            % displacement and velocity of the two nodes along axis 3 in
            % the reference frame of the first node. It is assumed that the
            % nodes motion is constrained appropriately by other MBDyn
            % elements (e.g. a prismatic joint).
            %
            % wsim.linearPTOWithController Methods:
            %
            %  wsim.linearPTOWithController - constructor
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
            %  controller - a wsim.wecController object
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
            %  'LoggedVars' - character vector, or cell array of character
            %    vectors indicating what internal variables are to be
            %    logged during a simulation. If a character vector, it must
            %    be 'none', meaning nothing will be logged. Otherwise it
            %    must be a cell array containing any combination of the
            %    following names:
            %
            %    'InternalForce' : This is the internal PTO force
            %      calculated using the controller
            %
            %    'RelativeDisplacement' : The relative displacement used to 
            %      calculate the PTO force
            %
            %    'RelativeVelocity' : The relative velocity used to 
            %      calculate the PTO force
            %
            %    The logged variables will be put in the wecSim wsim.logger
            %    object with unique names, created by adding a prefix
            %    'PTO_X_' where 'X' is replaced with an integer. The
            %    integer is incremented for each PTO in the system, e.g.
            %    PTO_1_, PTO_2_ etc. This allows multiple PTO objects with
            %    the same internal variable names to be used.
            %
            % Output
            %
            %  lpto - a wsim.linearPTOWithController object
            %
            %
            %
            % See Also: wsim.rotaryPowerTakeOff
            %

            options.InitialDisplacementZero = true;
            options.LoggedVars = {};
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (controller, 'wsim.wecController'), ...
                'controller must be a wsim.wecController or derived object');
            
            
            info.AvailableNames = { 'InternalForce', ...
                                    'RelativeDisplacement', ...
                                    'RelativeVelocity' ...
                                  };
                              
            info.IndepVars = { 'Time', ...
                               'Time', ...
                               'Time' };
                              
            info.Sizes = { [1,1], [1,1], [1,1] };
            
            info.Descriptions = { 'Force applied between the two PTO nodes', ...
                                  'Relative displacement of the two PTO nodes in a direction parallel to the reference node''s chosen axis', ...
                                  'Relative velocity of the two PTO nodes in a direction parallel to the reference node''s chosen axis' };
                              
            info.AxisLabels = { 'Force [N]', 'Displacement [m]', 'Velocity [ms^{-1}]' };
            
            
            % get the information of variables to be logged by the
            % controller
            ctrl_info = controller.loggingInfo ();
            % append it to the pto info
            f = fieldnames(info);
            for i = 1:length(f)
                info.(f{i}) = [info.(f{i}), ctrl_info.(f{i})];
            end
            
            info.NAvailable = numel(info.AvailableNames);
            
            self = self@wsim.powerTakeOff ( reference_node, other_node, info, ...
                                            'LoggedVars', options.LoggedVars );
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero );
                                
            self.controller = controller;
                                
        end
        
        function [FM, ptoforce, reldisp, relvel] = forceAndMoment (self, time)
            % returns the forces and moments applied by the PTO in 3D
            %
            % Syntax
            %
            % [FM, ptotorque, reltheta, relomega] = forceAndMoment (rot)
            %
            % Description
            %
            % wsim.linearPTOWithController.forceAndMoment calculates the forces
            % and moments applied in the global frame by the PTO to the two
            % nodes associated with the PTO.
            %
            % Input
            %
            %  rot - wsim.rotaryPowerTakeOff object
            %
            % Output
            %
            %  FM - (6 x 2) vector of forces and moments in the global
            %   frame which the PTO is applying to the two nodes. The first
            %   column is the forces applied to the reference node, the
            %   second column is the forces and moments applied to the
            %   other node.
            %
            %  ptoforce - scalar force value applied to the nodes parallel
            %   to axis three of the reference node
            %
            %  reldisp - relative displacement of the nodes in a direction
            %   parallel to axis three of the reference node
            %
            %  relvel - relative angular velocity of the nodes in a
            %   direction parallel to axis three of the reference node
            %
            %
            
            [reldisp, relvel] = self.mbdynForceObj.displacements ();

            self.internalVariables.RelativeDisplacement = reldisp;
            self.internalVariables.RelativeVelocity = relvel;
            
            % get the force, which is expected to be the control output
            % variable in this case
            [ptoforce, ctrl_log_vars] = self.controller.ptoControlOutput (self.id, time, self.internalVariables);
            
            FM = self.mbdynForceObj.force (ptoforce);
            
            % need to add zero moments to forces
            FM = [FM; zeros(size (FM))];
            
            % note that this internalVariables property is astructure which
            % is initialised by wsim.powerTakeOff
            self.internalVariables.InternalForce = ptoforce;
            
             % copy the controller logged variables int the PTO
             % internalVariables structure
             f = fieldnames(ctrl_log_vars);
             for i = 1:length(f)
                self.internalVariables.(f{i}) = ctrl_log_vars.(f{i});
             end
            
        end
        
    end
    
end