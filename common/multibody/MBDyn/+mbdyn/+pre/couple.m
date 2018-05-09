classdef couple < mbdyn.pre.element
    
    
    properties (GetAccess = public, SetAccess = protected)
       
    end
    
    properties (GetAccess = protected, SetAccess = protected)

    end
    
    methods
        
%         function self = force ()
%             
%             
%         end

        function str = generateMBDynInputString (self)
            str = sprintf ('    couple : %d,', self.label);
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