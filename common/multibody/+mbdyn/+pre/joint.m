classdef joint < mbdyn.pre.element
    
    properties
        
    end
    
    methods
        
        function str = generateOutputString (self)
            str = sprintf ('    joint : %d, %s,', self.label, self.type);
        end
        
    end
    
    methods (Access = protected)
        
        function ok = checkNodeReferenceType (self, ref, throw)
            
            ok = false;
            if ischar (ref)
                switch ref

                    case 'node'
                        ok = true;
                    case 'local'
                        ok = true;
                    case 'global'
                        ok = true;

                end
            end
            
            if ~ok && throw
                error ('invalid reference type, should be ''node'' or ''global''');
            end
            
        end
        
    end
    
end