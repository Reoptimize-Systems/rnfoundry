classdef twoNodeTorque < mbdyn.mint.base
    
    properties (GetAccess = public, SetAccess = private)
        
        torqueFcn;
        referenceTheta;
        joint;
        referenceNode;
        
    end
    
    properties
        
        nodes;
        
    end
    
    methods
        
        function self = twoNodeTorque (revolute_hinge, varargin)
            % mbdyn.mint.twoNodeTorque constructor
            %
            %
            % Syntax
            %
            % tnt = twoNodeTorque (revolute_hinge, varargin)
            %
            % Description
            %
            % mbdyn.mint.twoNodeTorque is a helper class used to apply a
            % torque calculated in Matlab to a rotational joint conencting
            % two nodes during an MBDyn simulation.
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
            %  'InitialThetaZero' - 
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
            
            options.InitialThetaZero = true;
            options.TorqueFcn = [];
            options.ReferenceNode = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkLogicalScalar ( options.InitialThetaZero, ...
                                                true, 'InitialThetaZero' );
                                            
            mbdyn.pre.base.checkScalarInteger ( options.ReferenceNode, ...
                                                true, 'ReferenceNode' );
                                            
            assert ( options.ReferenceNode == 1 || options.ReferenceNode == 2, ...
                'ReferenceNode must be 1 or 2');
            
            
            
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
            
            self.referenceTheta = 0;
            if options.InitialThetaZero
                % make the initial displacement be set to zero so all
                % future displacements are reported as relative to this
                % initial position
                [reltheta, ~] = displacements (self);
                
                self.referenceTheta = reltheta;
                
            end
            
            if ~isempty (options.TorqueFcn)
                assert (isa (options.TorqueFcn, 'function_handle'), ...
                    'TorqueFcn is not a function handle' ) ;
            end
            
            self.torqueFcn = options.TorqueFcn;
            
        end
        
        function [M, momentval, reltheta, relomega] = momentFromFcn (self)
            
            [reltheta, relomega] = displacements (self);
            
            momentval = feval ( self.torqueFcn, reltheta, relomega );
            
            M = moment (self, momentval);
            
        end
        
        
        function M = moment (self, momentval)
            
            M = self.joint.node1FrameRelativeOrientation.orientationMatrix(:,3) * momentval;
            M(1:3,2) = -M;
            
        end
        
        function [reltheta, relomega] = displacements (self)
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
            
            % get the relative position and velocity of the other node
            % relative to the reference node in the 
%             [joint_ref_pos, joint_ref_orient] = self.joint.reference ();
            
            % 	Vec3 v(RotManip::VecRot((pNode1->GetRCurr()*R1h).MulTM(pNode2->GetRCurr()*R2h)));
            % 	doublereal dThetaTmp(v(3));
            % 
            % 

            v = self.vecRot ( ...
                ( self.nodes(1).absoluteOrientation.orientationMatrix * self.joint.node1FrameRelativeOrientation.orientationMatrix ) ...
                * ( self.nodes(2).absoluteOrientation.orientationMatrix * self.joint.node2FrameRelativeOrientation.orientationMatrix ) ...
                            );
                        
            reltheta = self.referenceTheta + v(3);
    
            % relative angular velocity
            %
            % perform a rotation to bring the orientation of the joint into
            % line with the global coordinate system.
            % 
            % we pre-multiply by the transpose of the joint orientation in
            % the global frame to get the reverse rotation in the global
            % frame
%             w_n_ref_joint_frame = joint_ref_orient.orient.orientationMatrix.' ...
%                 * self.nodes(1).absoluteAngularVelocity;
%             
%             w_n_other_joint_frame = joint_ref_orient.orient.orientationMatrix.' ...
%                 * self.nodes(2).absoluteAngularVelocity;
%             
%             relomega = w_n_other_joint_frame(3) - w_n_ref_joint_frame(3);

            R = self.nodes(2).absoluteOrientation.orientationMatrix ...
                * self.joint.node2FrameRelativeOrientation.orientationMatrix;
            
            relomega = R * (self.nodes(2).absoluteAngularVelocity - self.nodes(1).absoluteAngularVelocity);
            
        end
        
    end
    
    methods (Access=protected)
        
        function unit = vecRot (self, Phi)
            
%         Vec3 RotManip::VecRot(const Mat3x3 & Phi) 

            % Modified from Appendix 2.4 of
            %
            % author = {Marco Borri and Lorenzo Trainelli and Carlo L. Bottasso},
            % title = {On Representations and Parameterizations of Motion},
            % journal = {Multibody System Dynamics},
            % volume = {4},
            % pages = {129--193},
            % year = {2000}
            
            cosphi = (trace(Phi) - 1) / 2;
            
            if (cosphi > 0)
                
                unit = self.ax (Phi);
                sinphi = norm (unit);
                phi = atan2 (sinphi, cosphi);
                absphi = abs(phi);
                a = sin (absphi) / absphi;
                unit = unit ./ a;
                
            else
                % -1 <= cosphi <= 0
                eet = (Phi + Phi.')/2; % ensure matrix is symmetric
                eet(1, 1) = eet(1, 1) - cosphi;
                eet(2, 2) = eet(2, 2) - cosphi;
                eet(3, 3) = eet(3, 3) - cosphi;
                
                % largest (abs) component of unit vector phi/|phi|
                maxcol = 1;
                
                if (eet(2, 2) > eet(1, 1))
                    maxcol = 2;
                end
                
                if (eet(3, 3) > eet(maxcol, maxcol))
                    maxcol = 3;
                end
                
                unit = eet(:,maxcol) ./ sqrt (eet(maxcol, maxcol) * (1. - cosphi));
                
                % sinphi = -(Mat3x3(unit)*Phi).Trace()/2.;
                sinphi = -(trace (cross (unit,Phi)) ) ./ 2;
                
                unit = unit * atan2 (sinphi, cosphi);
                
            end
            
            function v = ax (self, M)
                
                v = 0.5 .* [ M(3,2)-M(2,3), ...
                             M(1,3)-M(3,1), ...
                             M(2,1)-M(1,2) ];
                
            end
            
        end
    end
    
end