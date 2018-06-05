classdef rotaryPowerTakeOff < wsim.powerTakeOff
% power take-off from relative angular displacement and velocity of two nodes
%
% Syntax
%
% rot = wsim.rotaryPowerTakeOff (revolute_hinge)
% rot = wsim.rotaryPowerTakeOff (..., 'Parameter', Value)
%
% Description
%
% wsim.rotaryPowerTakeOff is a class representing a rotary power-take-off
% mechanism in a wave energy converter. It facilitates sending the correct
% forces to an MBDyn multibody simulation. wsim.rotaryPowerTakeOff applies
% opposite moments to two MBDyn nodes based attahed to a
% mbdyn.pre.revoluteHinge object from the MBDYn simulation toolbox. Moments
% are applied based on the relative angular displacement and velocity of
% the two nodes about the revolute hinge rotation axis.
%
% wsim.rotaryPowerTakeOff Methods:
%
%   rotaryPowerTakeOff - wsim.rotaryPowerTakeOff constructor
%   forceAndMoment - returns the forces and moments applied by the PTO in 3D
%
%   Inheirited Methods
%
%   advanceStep - advance to the next simulation time step
%   logData - appends the internal variable data to the log
%   loggingSetup - sets up data logging for a wsim.linearPowerTakeOff object
%
%
% See Also: wsim.linearPowerTakeOff, wsim.powerTakeOff
%

    properties (GetAccess = private, SetAccess = private)
        
        mbdynMomentObj;
        joint;

    end
    
    methods
        
        function self = rotaryPowerTakeOff (revolute_hinge, torque_fcn, varargin)
            % wsim.rotaryPowerTakeOff constructor
            %
            % Syntax
            %
            % rot = wsim.rotaryPowerTakeOff (revolute_hinge)
            % rot = wsim.rotaryPowerTakeOff (..., 'Parameter', Value)
            %
            % Description
            %
            % wsim.rotaryPowerTakeOff is a class representing a rotary
            % power-take-off mechanism in a wave energy converter. It
            % facilitates sending the correct forces to an MBDyn multibody
            % simulation. wsim.rotaryPowerTakeOff applies opposite moments
            % to two MBDyn nodes based attahed to a mbdyn.pre.revoluteHinge
            % object from the MBDYn simulation toolbox. Moments are applied
            % based on the relative displacement and velocity of the two
            % nodes about the revolute hinge rotation axis.
            %
            % Input
            %
            %  revolute_hinge - mbdyn.pre.revoluteHinge object
            %
            %  torque_fcn - function handle or string to be used to
            %   calculate the torque to be applied between the two nodes
            %   making up the PTO. torque_fcn is a function which takes two
            %   arguments with the following signature:
            %
            %        torque_value = myfcn (time, reltheta, relomega)
            %
            %   where reltheta is the relative angular displacement of the
            %   two nodes about the revolute hinge axis, and relomega is
            %   the relative angular velocity of the two nodes abot that
            %   axis. torque_value is expected to be a scalar value, the
            %   value of the torqe acting on the reference node. The
            %   opposite torque is applied to the other node. By default
            %   the reference node is the first node connected to the
            %   revolute hinge. See the 'ReferenceNode' options below to
            %   change this.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'InitialThetaZero' - optional true/false flag
            %    indicating whether the intial relative rotation angle
            %    (about the revolute hinge axis) should be taken as the
            %    reference point for relative rotation during the
            %    simulation, i.e. the PTO starts with an initial relative
            %    theta of zero for the purposes of torque calulation, and
            %    future displacement is measured relative to this initial
            %    position. If false, the raw relative angular displacement
            %    is used instead. Default is true if not supplied.
            %
            %  'ReferenceNode' - optional integer, 1 or 2. By default the
            %    torque returned by torque_fcn (see above) is applied to
            %    the first node attached to the revolute hinge. This can be
            %    changed with this option. The number supplied in
            %    ReferenceNode sets which node is the reference node.
            %    Default is 1 if not supplied.
            %
            %  'LoggedVars' - cell array of character vectors containing
            %    the names of internal variables to be logged. If supplied,
            %    only those variables named in this cell array will
            %    actually have data logged. The available variables are:
            %            
            %    'InternalMoment' : the scalar moment (torque) applied
            %      about the revolute hinge axis
            %
            %    'RelativeAngularDisplacement' : the relative angular
            %      displacement of the nodes attached to the revolute hinge
            %
            %    'RelativeAngularVelocity' : the relative angular
            %      velocity of the nodes attached to the revolute hinge
            %
            %    Alternatively, this can be a character vector, 'none', in
            %    which case no internal variables will be logged. Default
            %    is an empty cell array, which means all available nternal
            %    variables will be logged.
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
            %  rot - wsim.rotaryPowerTakeOff object
            %
            %
            %
            % See Also: wsim.linearPowerTakeOff, wsim.powerTakeOff
            %

            options.InitialTheta = 0;
            options.ReferenceNode = 1;
            options.LoggedVars = {};
            
            options = parse_pv_pairs (options, varargin);
            
            momobj = mbdyn.mint.twoNodeTorque ( ...
                                    revolute_hinge, ...
                                    'ReferenceNode', options.ReferenceNode, ...
                                    'InitialTheta', options.InitialTheta, ...
                                    'TorqueFcn', torque_fcn );
                                
            lginfo.AvailableNames = { 'InternalTorque', ...
                                      'RelativeAngularDisplacement', ...
                                      'RelativeAngularVelocity' ...
                                    };
                              
            lginfo.IndepVars = { 'Time', ...
                                 'Time', ...
                                 'Time' };
                              
            lginfo.Sizes = { [1,1], [1,1], [1,1] };
            
            lginfo.Descriptions = { 'Internal torque applied between the two PTo nodes', ...
                                    'Relative angular displacement of the other node relative to the reference node', ...
                                    'Relative angular velocity of the other node relative to the reference node'};
            
            lginfo.AxisLabels = {'Torque [Nm]', 'Theta [rad]', 'Omega [rads^{-1}]'};
            
            lginfo.NAvailable = numel(lginfo.AvailableNames);
            
            
            self = self@wsim.powerTakeOff ( momobj.nodes(1), ...
                                            momobj.nodes(2), ...
                                            lginfo, ...
                                            'LoggedVars', options.LoggedVars );
            
            self.joint = revolute_hinge;
            self.mbdynMomentObj = momobj;
            
        end
        
        function [FM, ptotorque, reltheta, relomega] = forceAndMoment (self, time)
            % returns the forces and moments applied by the PTO in 3D
            %
            % Syntax
            %
            % [FM, ptotorque, reltheta, relomega] = forceAndMoment (rot)
            %
            % Description
            %
            % wsim.rotaryPowerTakeOff.forceAndMoment calculates the forces
            % and moments applied in the global frame by the PTO to the two
            % nodes attached to the revolute hinge associated with the PTO.
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
            %   second node.
            %
            %  ptotorque - scalar torque value applied to the nodes about
            %   the revolute hinge axis
            %
            %  reltheta - relative angular displacement of the nodes about
            %   the revolute hinge axis
            %
            %  relomega - relative angular velocity of the nodes about
            %   the revolute hinge axis
            %
            %
            
            [M, ptotorque, reltheta, relomega] = self.mbdynMomentObj.momentFromFcn (time);
            
            % need to add zero forces to forces vector
            FM = [ zeros(size (M));
                   M ];
              
            self.internalVariables.InternalTorque = ptotorque;
            self.internalVariables.RelativeAngularDisplacement = reltheta;
            self.internalVariables.RelativeAngularVelocity = relomega;
            
        end
        
        function start (self, siminfo)
            
            self.mbdynMomentObj.initialise (siminfo.TStart);
            
            start@wsim.powerTakeOff (self, siminfo);

        end
        
        function advanceStep (self, time)
            
            self.mbdynMomentObj.advance (time);
            
            advanceStep@wsim.powerTakeOff (self);
            
        end
         
    end
    
end

