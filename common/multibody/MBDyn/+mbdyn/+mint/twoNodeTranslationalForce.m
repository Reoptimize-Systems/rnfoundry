classdef twoNodeTranslationalForce < mbdyn.mint.twoNodeForce
% replace h1 line
%
% Syntax
%
% tnf = mbdyn.mint.twoNodeTranslationalForce (reference_node, other_node, axisNum)
% tnf = mbdyn.mint.twoNodeTranslationalForce (..., 'Parameter', value)
%
% Description
%
% twoNodeTranslationalForce is a class used to help apply forces to two
% appropriately constrained structural nodes during an MBDyn simulation. It
% is intended to help apply forces from systems such as actuators, with two
% nodes constrained so that only motion along a common axis is possible.
% However, no check that the nodes are constrained in this way is
% performed.
%
% The forces are applied along the chosen axis of a reference node, in the
% reference frame of that node and the resulting forces calculated in the
% global frame. The reverse of these forces is calculated for the other
% node.
%
% mbdyn.mint.twoNodeTranslationalForce Methods:
%
%   twoNodeTranslationalForce - mbdyn.mint.twoNodeTranslationalForce 
%    constructor. See for full description of setup options.
%   displacements - gets the relative displacement and velocity of the 
%    other node
%   force - returns forces in the global frame on the two nodes
%   forceFromFcn - returns global forces on the two nodes evaluated from a 
%    function
%
%
% See Also: mbdyn.mint.twoNodeTorque
%
    
    properties (GetAccess=public, SetAccess = private)
        
        forceAxis;
        forceFcn;
        referenceDisp;
        
    end
    
    properties (GetAccess=private, SetAccess = private)
        
        hasForceFcn;
        
    end
    
    methods
        
        function self = twoNodeTranslationalForce (reference_node, other_node, axisNum, varargin)
            % mbdyn.mint.twoNodeTranslationalForce constructor
            %
            % Syntax
            %
            % tnf = mbdyn.mint.twoNodeTranslationalForce (reference_node, other_node, axisNum)
            % tnf = mbdyn.mint.twoNodeTranslationalForce (..., 'Parameter', value)
            %
            % Description
            %
            % twoNodeTranslationalForce is a class used to help apply
            % forces to two appropriately constrained structural nodes
            % during an MBDyn simulation. It is intended to help apply
            % forces from systems such as actuators, with two nodes
            % constrained so that only motion along a common axis is
            % possible. However, no check that the nodes are constrained in
            % this way is performed.
            %
            % The forces are applied along the chosen axis of a reference
            % node, in the reference frame of that node and the resulting
            % forces calculated in the global frame. The reverse of these
            % forces is calculated for the other node.
            %
            % Input
            %
            %  reference_node - mbdyn.pre.structuralNode6dof object
            %   representing the structural node in whose reference frame
            %   the forces will be applied along the chosen axis
            %
            %  other_node - mbdyn.pre.structuralNode6dof
            %
            %  axisNum - scalar integer indicating the axis of the
            %   reference node along which forces will be applied. Can be
            %   1, 2 or 3.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'InitialDisplacementZero' - optional scalar true/false flag
            %    indicating whether the initial position of the two nodes
            %    should be considered a displacement of zero, and all
            %    subsequent displacements considered relative to this
            %    initial position. If false, the actual initial and
            %    subsequent displacement in the direction of the reference
            %    node's chosen axis will be reported by the 'displacements'
            %    method (and used in the forceFromFcn method to evaluate
            %    the force function). Default is true.
            %    
            %  'ForceFcn' - optional function handle or string to be used
            %    to calculate the force to be applied between the two
            %    nodes. This is required to use the forceFromFcn method.
            %    forceFcn is a function which takes two arguments with the
            %    following signature:
            %
            %         force_value = myfcn (reldisp, relvel)
            %
            %    where reldisp is the relative displacement of the two
            %    nodes along the specified axis in forceAxis in the
            %    reference frame of the reference node, and relvel is the
            %    relative velocity of the two nodes in the same frame.
            %    force_value is expected to be a scalar value, the value of
            %    the force acting on the reference node parallel to the
            %    axis in forceAxis in the frame of the reference node.
            %
            % Output
            %
            %  tnf - mbdyn.mint.twoNodeTranslationalForce object
            %
            %
            %
            % See Also: mbdyn.mint.twoNodeTorque
            %
            
            options.InitialDisplacementZero = true;
            options.ForceFcn = [];
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkLogicalScalar ( options.InitialDisplacementZero, ...
                                                true, 'InitialDisplacementZero' );
            
            self = self@mbdyn.mint.twoNodeForce (reference_node, other_node);
            
            if mbdyn.pre.base.checkScalarInteger (axisNum, false) ...
                    && any (axisNum == [1,2,3])
                
                self.forceAxis = axisNum;
                
            else
                error ('axisNum contains invalid axis number, should be an integer 1, 2 or 3');
            end
            
            self.referenceDisp = 0;
            if options.InitialDisplacementZero
                % make the initial displacement be set to zero so all
                % future displacements are reported as relative to this
                % initial position
                [reldisp, ~] = displacements (self);
                
                self.referenceDisp = reldisp;
                
            end
            
            if isempty (options.ForceFcn)
                self.hasForceFcn = false;
            else
                assert (isa (options.ForceFcn, 'function_handle'), ...
                    'ForceFcn is not a function handle' );
                self.hasForceFcn = true;
            end
            
            self.forceFcn = options.ForceFcn;
            
        end
        
        function [F, ptoforce, reldisp, relvel] = forceFromFcn (self, time)
            % returns global forces on the two nodes evaluated from a function
            %
            % Syntax
            %
            % [F, ptoforce, reldisp, relvel] = forceFromFcn (self)
            %
            % Description
            %
            % forceFromFcn evaluates the force function supplied when
            % constructing the twoNodeTranslationalForce object, and which
            % is stored in the forceFcn property. It then calculates the
            % forces on the reference and other nodes in the global frame
            % after applying this force along the previously specified axis
            % (the value stored n the forceAxis property) in the frame of
            % the reference node. The supplied force is applied to the
            % NON-reference (OTHER) node, and the reverse of these forces
            % is also calculated to give the force on the reference node.
            %
            % forceFcn is a function which takes three arguments with the
            % following signature:
            %
            %         force_value = myfcn (time, reldisp, relvel)
            %
            % where  time is the current simulation time, supplied as an
            % input to forceFromFcn, reldisp is the relative displacement
            % of the two nodes along the specified axis in forceAxis in the
            % reference frame of the reference node, and relvel is the
            % relative velocity of the two nodes in the same frame.
            % force_value is expected to be a scalar value, the value of
            % the force acting on the non-reference node parallel to the
            % axis in forceAxis in the frame of the reference node.
            %
            % Input
            %
            %  tnf - mbdyn.mint.twoNodeTranslationalForce object
            %
            %  time - scalar value of the current simulation time
            %
            % Output
            %
            %  F - (3 X 2) matrix. The first column is the forces in the
            %   global frame on the reference node, the second is the
            %   forces on the other node. These are simply the inverse of
            %   the first column.
            %
            %  ptoforce - scalar value of the force calculated by
            %   evaluating the forceFcn
            %
            %  reldisp - scalar value of the relative displacement of the
            %   two nodes along the specified axis in forceAxis in the
            %   reference frame of the reference node
            %
            %  relvel - scalar value of the relative velocity of the
            %   two nodes along the specified axis in forceAxis in the
            %   reference frame of the reference node
            %
            %
            % See Also: mbdyn.mint.twoNodeTranslationalForce.force
            %
            
            assert (self.hasForceFcn, ...
                'There is no force function to evaluate, use the ''ForceFcn'' option when constructing');
            
            [reldisp, relvel] = displacements (self);
            
            ptoforce = feval ( self.forceFcn, time, reldisp, relvel );
            
            F = force (self, ptoforce);
            
        end
        
        function F = force (self, ptoforce)
            % returns forces in the global frame on the two nodes 
            %
            % Syntax
            %
            % F = force (tnf, ptoforce)
            %
            % Description
            %
            % force calculates the forces on the reference node in the
            % global frame after applying a force parallel to a previously
            % specified axis (the value stored n the forceAxis property) in
            % the frame of the reference node. The reverse of these forces
            % is also calculated to give the force on the other node.
            %
            % Input
            %
            %  tnf - mbdyn.mint.twoNodeTranslationalForce object
            %
            %  ptoforce - scalar value of the force applied to the
            %   NON-reference node (the OTHER node), parallel to the chosen
            %   axis of the reference node.
            %
            % Output
            %
            %  F - (3 X 2) matrix. The first column is the forces in the
            %   global frame on the reference node, the second is the
            %   forces on the other node. These are simply the inverse of
            %   the first column.
            %
            % See Also: mbdyn.mint.twoNodeTranslationalForce.forceFromFcn
            %

            F_pto_frame = [0; 0; 0];
            F_pto_frame(self.forceAxis) = ptoforce;
            
            F = zeros (3,2);
            F(1:3,2) = (self.referenceNode.absoluteOrientation.orientationMatrix * F_pto_frame)' ;
            F(1:3,1) = -F(1:3,2);
            
        end
        
        function [reldisp, relvel] = displacements (self)
            % gets the relative displacement and velocity of the other node
            % relative to the reference node in the chosen axis of the
            % reference node coordinate frame
            %
            % Syntax
            %
            % [reldisp, relvel] = displacements (obj)
            %
            % Input
            %
            %  obj - mbdyn.mint.twoNodeTranslationalForce object
            %
            % Output
            %
            %  reldisp - displacement of the second node relative to the
            %   reference node in the chosen axis of the reference node's
            %   frame
            %
            %  relvel - velocity of the second node relative to the
            %   reference node in the chosen axis of the reference node's
            %   frame
            %
            %
            
            % get the relative position and velocity of the other node
            % relative to the reference node
            xRvec = self.otherNode.absolutePosition - self.referenceNode.absolutePosition;
            vRvec = self.otherNode.absoluteVelocity - self.referenceNode.absoluteVelocity;
            
            % perform a rotation to bring the orientation of the reference
            % node into line with the global coordinate system.
            % 
            % we pre-multiply by the transpose of the reference node
            % orientation in the global frame to get the reverse rotation
            % in the global frame
            xRforceVec = self.referenceNode.absoluteOrientation.orientationMatrix.' * xRvec;
            vRforceVec = self.referenceNode.absoluteOrientation.orientationMatrix.' * vRvec;
            
            % velocities and displacements are the desired components of
            % the vectors in the translator coordinate system
            reldisp = xRforceVec(self.forceAxis) - self.referenceDisp;
            relvel = vRforceVec(self.forceAxis);
            
        end
        
    end
    
end