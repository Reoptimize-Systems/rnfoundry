classdef force < mbdyn.pre.element
    
    
    properties (GetAccess = public, SetAccess = protected)
        subType;
    end
    
    properties (GetAccess = protected, SetAccess = protected)

    end
    
    methods
        
        function self = force (varargin)
            
            [ options, nopass_list ] = mbdyn.pre.force.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.element (pvpairs{:} );
            
            self.type = 'force';
            self.netCDFName = 'force';
            
        end

        function str = generateMBDynInputString (self)
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
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.element.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end