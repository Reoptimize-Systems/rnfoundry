classdef twoNodeTranslationalForce < mbdyn.mint.twoNodeForce
    
    properties
        
        forceAxis;
        forceFcn;
        referenceDisp;
        
    end
    
    methods
        
        function self = twoNodeTranslationalForce (reference_node, other_node, axisNum, varargin)
            
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
            
            if ~isempty (options.ForceFcn)
                assert (isa (options.ForceFcn, 'function_handle'), ...
                    'ForceFcn is not a function handle' ) ;
            end
            
            self.forceFcn = options.ForceFcn;
            
        end
        
        function [F, ptoforce, reldisp, relvel] = force (self)
            
            [reldisp, relvel] = displacements (self);
            
            ptoforce = feval ( self.forceFcn, reldisp, relvel );
            
            F_pto_frame = [0, 0, 0];
            F_pto_frame(self.forceAxis) = ptoforce;
            
            F = (F_pto_frame * self.referenceNode.absoluteOrientation.orientationMatrix)' ;
            F(1:3,2) = -F;
            
%             F(1:3,1) = F(1:3,1) + FptoVec(1:3,1);
%             F(1:3,2,ind) = F(1:3,2) - FptoVec(1:3,1);
            
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
            %    reference node in the chosen axis of the reference node's
            %    frame
            %
            %  relvel - velocity of the second node relative to the
            %    reference node in the chosen axis of the reference node's
            %    frame
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