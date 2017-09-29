classdef initialValueProblem < mbdyn.pre.problem
    
    properties
        initialTime;
        finalTime;
        timeStep;
        maxIterations;
        tolerance;
        derivativesTolerance;
        output;
    end
    
    methods
        
        function self = initialValueProblem (itime, ftime, tstep, varargin)
            
            options.MaxIterations = 20;
            options.Tolerance = 1e-9;
            options.SolutionTolerance = [];
            options.Method = [];
            options.DerivativesTolerance = [];
            options.Output = {};
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (itime, true, 'itime');
            self.checkNumericScalar (ftime, true, 'ftime');
            self.checkNumericScalar (tstep, true, 'tstep');
            self.checkNumericScalar (options.Tolerance , true, 'Tolerance');
            
            if ~isempty (options.Output)
               if ~iscellstr (options.Output)
                  error ('Output must be a cell array of strings');
               end
               for ind = 1:numel(options.Output)
                   self.checkAllowedStringInputs (options.Output{ind}, ...
                       { 'iterations', ...
                         'residual', ...
                         'solution', ...
                         'jacobian matrix', ...
                         'messages', ...
                         'counter', ...
                         'bailout', ...
                         'matrix condition number', ...
                         'solver condition number', ...
                         'cpu time', ...
                         'none' }, ...
                       true, sprintf('Output(%d)', ind));
               end
            end
            
            self.initialTime = itime;
            self.finalTime = ftime;
            self.timeStep = tstep;
            self.maxIterations = options.MaxIterations;
            self.tolerance = options.Tolerance;
            self.derivativesTolerance = options.DerivativesTolerance;
            self.output = options.Output;
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
            
            if ~isempty (self.output)
                str = self.addOutputLine (str, sprintf('output: %s;', self.commaSepList(self.output{:})), 1, false);
            end
            
            str = self.addOutputLine (str, 'end: initial value;', 0, false);
            
        end
       
    end
    
end