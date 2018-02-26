classdef force < mbdyn.pre.element
    
    
    properties (GetAccess = public, SetAccess = protected)
        subType;
    end
    
    properties (GetAccess = protected, SetAccess = protected)

    end
    
    methods
        
        function self = force ()
            
            self.type = 'force';
            
        end

        function str = generateOutputString (self)
            str = sprintf ('    %s : %d,', self.type, self.label);
        end
        
%         function setSize (self, sx, sy, sz)
%             self.sx = sx;
%             self.sy = sy;
%             self.sz = sz;
%         end
%         
%         function setColour (self, newcolour)
%             self.drawColour = newcolour;
%         end
        
    end
    
    methods (Access = protected)
        
    end
    
end