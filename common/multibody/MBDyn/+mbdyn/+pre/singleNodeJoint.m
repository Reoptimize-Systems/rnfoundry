classdef singleNodeJoint < mbdyn.pre.joint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node;
        
    end
    
    methods
        function self = singleNodeJoint (node)
            % generic base class for joints which constrain a single node
        
            self.checkIsStructuralNode (node, true);
            
            self.node = node;
            
        end
        
        function str = generateOutputString (self)
            str = generateOutputString@mbdyn.pre.joint(self);
        end
        
    end
    
    methods (Access = protected)
        
        function ok = checkNodeReferenceType (self, ref, throw)
            % checks that the specified reference frame is valid
            %
            % Syntax
            %
            % ok = checkNodeReferenceType (jntobj, ref, throw)
            %
            % Input
            %
            %  jntobj - mbdyn.pre.singleNodeJoint object
            %
            %  ref - char array specifying the reference frame
            %    in which a position is defined realtive to a node in a
            %    single node joint. Valid strings are: 'node', 'local' and
            %    'global'.
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkNodeReferenceType if ref fails check
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            % See Also: 
            %
            
            allowedstrs = {'node', 'local', 'global'};
            
            ok = self.checkAllowedStringInputs (ref, allowedstrs, throw, 'node reference type');
            
        end
        
    end
    
end