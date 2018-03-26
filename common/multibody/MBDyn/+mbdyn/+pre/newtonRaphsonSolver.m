classdef newtonRaphsonSolver < mbdyn.pre.nonlinearSolver
    
    
    properties (GetAccess = public, SetAccess = protected)
       modified;
    end
    
    methods
        
        function self = newtonRaphsonSolver (varargin)
            
            options.Modified = {};
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Modified)
               if (~iscell (options.Modified)) ...
                       || (numel (options.Modified) > 3)
                  error ('Modified must be a cell array of strings, of length less than 3');
               end
               
               self.checkNumericScalar (options.Modified{1}, true, 'iterations');
               
               allowedstrings = { 'keep jacobian matrix', ...
                                  'honor element requests' };
               for ind = 2:numel(options.Modified)
                   self.checkAllowedStringInputs (options.Modified{ind}, ...
                       allowedstrings, true, sprintf('Modified(%d)', ind));
                   
                   allowedstrings(1) = [];
               end
            end
            
            self.modified = options.Modified;
            
            self.type = 'newton raphson';
            
        end

        function str = generateMBDynInputString (self)
            
            args = {self.type};
            
            if ~isempty (self.modified)
                args = [args, {sprintf('modified, %d', self.modified{1})}];
                
                if numel (self.modified) > 1
                    args = [args, self.modified(2:end)];
                end
            end
            
            str = self.commaSepList (args{:});
            
        end
         
    end
    
end