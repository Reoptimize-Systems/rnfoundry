classdef twoNodeTorque < mbdyn.mint.base
% class for helping with application of torques via external structural elements
%
% Syntax
%
% tnt = twoNodeTorque (revolute_hinge)
% tnt = twoNodeTorque (..., 'Parameter', Value)
%
% Description
%
% mbdyn.mint.twoNodeTorque is a helper class used to apply a torque
% calculated in Matlab to two nodes connected by a rotational joint
% (specifically, a revolute hinge) during an MBDyn multibody dynamics
% simulation.
%
%
% mbdyn.mint.twoNodeTorque Methods:
%
%   twoNodeTorque - mbdyn.mint.twoNodeTorque constructor
%   displacements - gets the relative displacement and velocity of the other node
%   moment - returns moments in the global frame on the two nodes
%   momentFromFcn - returns global forces on the two nodes evaluated from a function
%
%
% See Also: mbdyn.mint.twoNodeTranslationalForce
%

    properties (GetAccess = public, SetAccess = private)
        
        torqueFcn; % torque function
        referenceTheta; 
        joint;
        referenceNode;
        
    end
    
    properties 
        
        nodes;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        omegaIntegral;
        lastTime;
        lastTheta;
        lastOmega;
        
    end
    
    methods
        
        function self = twoNodeTorque (revolute_hinge, varargin)
            % mbdyn.mint.twoNodeTorque constructor
            %
            % Syntax
            %
            % tnt = twoNodeTorque (revolute_hinge)
            % tnt = twoNodeTorque (..., 'Parameter', Value)
            %
            % Description
            %
            % mbdyn.mint.twoNodeTorque is a helper class used to apply a
            % torque calculated in Matlab to two nodes connected by a
            % rotational joint (specifically, a revolute hinge) during an
            % MBDyn multibody dynamics simulation.
            %
            % Input
            %
            %  revolute_hinge - the mbdyn.pre.revoluteHinge object
            %   representing the revolute hinge joint connecting the two
            %   nodes
            %
            % Additional arguments may be supplied using parameter-value
            % pairs. The available  options are:
            %
            %  'InitialTheta' - optional scalar value which is the
            %    initial relative angular position in radians of the two
            %    nodes, all subsequent displacements considered relative to
            %    this initial position. Default is zero.
            %
            %  'TorqueFcn' - handle to matlab function which returns a
            %    scalar value, the torque on the reference node about axis
            %    3 of the revolute_hinge. The function must take as input
            %    two arguments, and it is passed the relative angle and the
            %    relative angular velocity between the nodes when the
            %    moment is to be calculated. Therefore the function handle
            %    must have a calling sequence like:
            %
            %       moment_value = myfcn (rel_theta, rel_omega)
            % 
            %  'ReferenceNode' - scalar integer, eithe 1 or 2. Indicates
            %    which node is to be considered the reference node for
            %    calculating the relative angular velocity.
            %
            % Output
            %
            %  tnt - mbdyn.mint.twoNodeTorque object
            %
            %
            %
            % See Also: mbdyn.mint.twoNodeTranslationalForce
            %           
            
            options.InitialTheta = 0;
            options.TorqueFcn = [];
            options.ReferenceNode = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkNumericScalar ( options.InitialTheta, ...
                                                true, 'InitialTheta' );
                                            
            mbdyn.pre.base.checkScalarInteger ( options.ReferenceNode, ...
                                                true, 'ReferenceNode' );
                                            
            assert ( options.ReferenceNode == 1 || options.ReferenceNode == 2, ...
                'ReferenceNode must be 1 or 2');
            
            self = self@mbdyn.mint.base ();
            
            if isa (revolute_hinge, 'mbdyn.pre.revoluteHinge')
                
                self.joint = revolute_hinge;
                
            else
                error ('revolute_hinge must be an mbdyn.pre.revoluteHinge object');
            end
            
            self.referenceNode = options.ReferenceNode;
            
            if self.referenceNode == 1
                self.nodes = [ self.joint.node1, self.joint.node2 ];
            else
                self.nodes = [ self.joint.node2, self.joint.node1 ];
            end
            
%             self.referenceTheta = 0;
%             if options.InitialThetaZero
%                 % make the initial displacement be set to zero so all
%                 % future displacements are reported as relative to this
%                 % initial position
%                 [reltheta, ~] = displacements (self);
%                 
%                 self.referenceTheta = reltheta;
%                 
%             end

            self.referenceTheta = options.InitialTheta;
            
            if ~isempty (options.TorqueFcn)
                assert (isa (options.TorqueFcn, 'function_handle'), ...
                    'TorqueFcn is not a function handle' ) ;
            end
            
            self.torqueFcn = options.TorqueFcn;
            
        end
        
        function [M, torqueval, reltheta, relomega] = momentFromFcn (self, time)
            % returns global forces on the two nodes evaluated from a function
            %
            % Syntax
            %
            % [M, torqueval, reltheta, relomega] = momentFromFcn (self)
            %
            % Description
            %
            % momentFromFcn evaluates the torque function supplied when
            % constructing the twoNodeTorque object, and which is stored in
            % the torqueFcn property. It then calculates the moments on the
            % reference and other nodes in the global frame after applying
            % this moment about the axis of the revolute hinge to which the
            % nodes are attached. The supplied moment is applied to the
            % NON-reference (OTHER) node, and the reverse of these moments
            % is also calculated to give the force on the reference node.
            %
            % torqueFcn is a function which takes three arguments with the
            % following signature:
            %
            %         torque_value = myfcn (time, reltheta, relomega)
            %
            % where time is the current simulation time, supplied as an
            % input to momentFromFcn, reltheta is the relative angular
            % displacement of the two nodes about the rotation axis of the
            % atached revolute hinge, and relomega is the relative angular
            % velocity of the two nodes in the same frame. torque_value is
            % expected to be a scalar value, the value of the torque acting
            % on the non-reference node.
            %
            % Input
            %
            %  tnf - mbdyn.mint.twoNodeTranslationalForce object
            %
            %  time - scalar value of the current simulation time
            %
            % Output
            %
            %  M - (3 X 2) matrix. The first column is the moments in the
            %   global frame on the reference node, the second is the
            %   moments on the other node. These are simply the inverse of
            %   the first column.
            %
            %  torqueval - scalar value of the torque calculated by
            %   evaluating the torqueFcn
            %
            %  reltheta - scalar value of the relative angular displacement
            %   of the two nodes about the attached revolute hinge rotation
            %   axis
            %
            %  relomega - scalar value of the relative angular velocity of
            %   the two nodes about the attached revolute hinge rotation
            %   axis
            %
            %
            % See Also: mbdyn.mint.twoNodeTranslationalForce.moment
            %
            
            [reltheta, relomega] = displacements (self, time);
            
            torqueval = feval ( self.torqueFcn, time, reltheta, relomega );
            
            M = moment (self, torqueval);
            
        end
        
        
        function M = moment (self, torqueval)
            % returns moments in the global frame on the two nodes 
            %
            % Syntax
            %
            % M = moment (tnf, torqueval)
            %
            % Description
            %
            % moment calculates the moments on the two nodes attached to
            % the joint node in the global frame after applying a torque
            % about the rotation axis of a revolute hinge. The reverse of
            % these forces is also calculated to give the moment on the
            % other node.
            %
            % Input
            %
            %  tnf - mbdyn.mint.twoNodeTranslationalForce object
            %
            %  torqueval - scalar value of the torque applied to the
            %   NON-reference node (the OTHER node), about the revolute
            %   hinge axis.
            %
            % Output
            %
            %  M - (3 X 2) matrix. The first column is the moments in the
            %   global frame on the reference node, the second is the
            %   moments on the other node. These are simply the inverse of
            %   the first column.
            %
            % See Also: mbdyn.mint.twoNodeTorque.momentFromFcn
            %
            
            M_pto_frame = [0; 0; 0];
            M_pto_frame(3) = torqueval;
            
            M = zeros (3,2);
            
            % apply the reverse of the absolute joint orientation to the
            % moments to put them in the global frame
            M(1:3,2) = (self.joint.absoluteJointOrientation.orientationMatrix' * M_pto_frame)';
            M(1:3,1) = -M(1:3,2);
            
        end
        
        function [reltheta, relomega] = displacements (self, time)
            % gets the relative displacement and velocity of the other node
            % relative to the reference node in the chosen axis of the
            % reference node coordinate frame
            %
            % Syntax
            %
            % [reltheta, relomega] = displacements (obj)
            %
            % Input
            %
            %  obj - mbdyn.mint.twoNodeTorque object
            %
            % Output
            %
            %  reltheta - displacement of the second node relative to the
            %    reference node in the chosen axis of the reference node's
            %    frame
            %
            %  relomega - velocity of the second node relative to the
            %    reference node in the chosen axis of the reference node's
            %    frame
            %
            %
            
    
            % relative angular velocity

% orig
%             R = self.nodes(2).absoluteOrientation.orientationMatrix ...
%                 * self.joint.node2FrameRelativeOrientation.orientationMatrix;

            R = self.nodes(1).absoluteOrientation.orientationMatrix * self.joint.relativeOrientation1.orientationMatrix;
            
            omegas = R * (self.nodes(2).absoluteAngularVelocity - self.nodes(1).absoluteAngularVelocity);
            
            relomega = omegas(3);
            
            % get theta by integrating omega
            reltheta = self.omegaIntegral + ( relomega * (time - self.lastTime) );
            
            % update the last calculated values of omega and theta
            self.lastOmega = relomega;
            self.lastTheta = reltheta;
            
        end
        
        function initialise (self, time)
            % intialise the twoNodeTorque object
            %
            % Syntax
            %
            % initialise (obj)
            %
            % Description
            %
            % intialise the twoNodeTorque object in preparation for a
            % simulation. This initialises the integral of omega, the
            % relative angular velocity to obtain theta, the relative
            % angular displacement.
            %
            % Input
            %
            %  obj - mbdyn.mint.twoNodeTorque object
            %
            %  time - the intial simulation time
            %
            %
            % See Also: 
            %

            mbdyn.pre.base.checkNumericScalar (time, true, 'time');
                       
            self.lastTime = time;
            self.lastTheta = self.referenceTheta;
            self.omegaIntegral = self.lastTheta;
            
            [~, relomega] = displacements (self, time);
            
            self.lastOmega = relomega;
            
        end
        
        function advance (self, time)
            % advance the twoNodeTorque calculation
            %
            % Syntax
            %
            % advance (obj, time)
            %
            % Description
            %
            % advance the twoNodeTorque calculation, accepting the last
            % clculated integral of omega, and the last simulation time
            % into the simulation history.
            %
            % Input
            %
            %  obj - mbdyn.mint.twoNodeTorque object
            %
            %  time - the current simulation time (used to calculate the
            %   last value of theta in a call to the displacements method
            %
            %
            % See Also: 
            %
            
            self.lastTime = time;
            self.omegaIntegral = self.lastTheta;
            
        end
        
    end
    
end