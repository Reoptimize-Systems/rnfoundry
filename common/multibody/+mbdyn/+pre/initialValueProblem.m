classdef initialValueProblem < mbdyn.pre.problem
    
    properties
        initialTime;
        finalTime;
        timeStep;
        maxIterations;
        tolerance;
        derivativesTolerance;
    end
    
    methods
        
        function self = initialValueProblem (itime, ftime, tstep, varargin)
            
            options.MaxIterations = 20;
            options.Tolerance = 1e-9;
            options.SolutionTolerance = [];
            options.Method = [];
            options.DerivativesTolerance = [];
            options = parse_pv_pairs (options, varargin);
            
            check.multicheck (@(x) (isscalar(x) && isnumeric(x)), ...
                'initial time, final time, tolerance must all be scalar numeric values', '', ...
                itime, ftime, tstep, options.Tolerance );
            
            self.initialTime = itime;
            self.finalTime = ftime;
            self.timeStep = tstep;
            self.maxIterations = options.MaxIterations;
            self.tolerance = options.Tolerance;
            self.derivativesTolerance = options.DerivativesTolerance;
            self.type = 'initial value';
            
        end
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 0, false, 'initial value problem');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str , 'begin: initial value;', 0, false);
            
            str = self.addOutputLine (str, sprintf('initial time: %s;', self.formatNumber (self.initialTime)), 1, false);
            str = self.addOutputLine (str, sprintf('final time: %s;', self.formatNumber (self.finalTime)), 1, false);
            str = self.addOutputLine (str, sprintf('time step: %s;', self.formatNumber (self.timeStep)), 1, false);
            
            if ~isempty (self.maxIterations)
                str = self.addOutputLine (str, sprintf('max iterations: %d;', self.maxIterations), 1, false);
            end
            
            if ~isempty (self.tolerance)
                str = self.addOutputLine (str, sprintf('tolerance: %e;', self.tolerance), 1, false);
            end
            
            if ~isempty (self.derivativesTolerance)
                str = self.addOutputLine (str, sprintf('derivatives tolerance: %e;', self.derivativesTolerance), 1, false);
            end
            
            str = self.addOutputLine (str, 'end: initial value;', 0, false);
            
        end
       
    end
    
end