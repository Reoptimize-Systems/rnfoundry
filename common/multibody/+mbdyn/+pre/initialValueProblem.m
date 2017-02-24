classdef initialValueProblem < mbdyn.pre.problem
    
    properties
        initialTime;
        finalTime;
        timeStep;
        maxIterations;
        tolerance;
    end
    
    methods
        
        function self = initialValueProblem (itime, ftime, tstep, varargin)
            
            options.MaxIterations = 10;
            options.Tolerance = 1e-6;
            
            options = parse_pv_pairs (options, varargin);
            
            check.multicheck (@(x) (isscalar(x) && isnumeric(x)), ...
                'initial time, final time, tolerance must all be scalar numeric values', '', ...
                itime, ftime, tstep, options.Tolerance );
            
            self.initialTime = itime;
            self.finalTime = ftime;
            self.timeStep = tstep;
            self.maxIterations = options.MaxIterations;
            self.tolerance = options.Tolerance;
            self.type = 'initial value';
            
        end
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 0, false, 'initial value problem');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str , 'begin: initial value;', 0, false);
            
            str = self.addOutputLine (str, sprintf('initial time: %.8f;', self.initialTime), 1, false);
            str = self.addOutputLine (str, sprintf('final time: %.8f;', self.finalTime), 1, false);
            str = self.addOutputLine (str, sprintf('time step: %.8f;', self.timeStep), 1, false);
            str = self.addOutputLine (str, sprintf('max iterations: %d;', self.maxIterations), 1, false);
            str = self.addOutputLine (str, sprintf('tolerance: %.8f;', self.tolerance), 1, false);
            
            str = self.addOutputLine (str, 'end: initial value;', 0, false);
            
        end
       
    end
    
end