classdef twoNodeTranslationalForce < mbdyn.mint.twoNodeForce
    
    properties
        
        forceAxis;
        
    end
    
    methods
        
        function self = twoNodeTranslationalForce (reference_node, other_node, axisNum)
            
            self = self@mbdyn.mint.twoNodeForce (reference_node, other_node);
            
            if isint2eps (axisNum) && any (axisNum == [1,2,3])
                self.forceAxis = axisNum;
            else
                error ('Invalid axis number');
            end
            
        end
        
        function force (self)
            
            [reldisp, relvel] = displacements (self);
            
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
            xRforceVec = self.referenceNode.absoluteOrientation.' * xRvec;
            vRforceVec = self.referenceNode.absoluteOrientation.' * vRvec;
            
            % velocities and displacements are the desired components of
            % the vectors in the translator coordinate system
            reldisp = xRforceVec(self.forceAxis);
            relvel = vRforceVec(self.forceAxis);
            
        end
        
    end
    
end