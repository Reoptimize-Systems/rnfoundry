classdef lineSearchSolver < mbdyn.pre.newtonRaphsonSolver
    
    
    properties (GetAccess = public, SetAccess = protected)

    end
    
    methods
        
        function self = lineSearchSolver (varargin)
            
            options.Modified = {};
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.newtonRaphsonSolver ('Modified', options.Modified);
            
            self.type = 'line search';
            
        end

        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.newtonRaphsonSolver (self);
        end
         
    end
    
end