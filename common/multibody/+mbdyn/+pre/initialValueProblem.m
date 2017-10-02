classdef initialValueProblem < mbdyn.pre.problem
    
    properties
        initialTime;
        finalTime;
        timeStep;
        maxIterations;
        tolerance;
        derivativesTolerance;
        output;
        nonlinearSolver;
    end
    
    methods
        
        function self = initialValueProblem (itime, ftime, tstep, varargin)
            
            options.MaxIterations = 20;
            options.ResidualTolerance = 1e-9;
            options.SolutionTolerance = {};
            options.Method = [];
            options.DerivativesTolerance = [];
            options.Output = {};
            options.NonlinearSolver = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (itime, true, 'itime');
            self.checkNumericScalar (ftime, true, 'ftime');
            self.checkNumericScalar (tstep, true, 'tstep');
            
            ResidualTolerance = self.checkResidualTolerance (options.ResidualTolerance);
            SolutionTolerance = self.checkSolutionTolerance (options.SolutionTolerance);
            
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
            
            if ~isempty (options.NonlinearSolver) ...
                    && ~isa (options.NonlinearSolver, 'mbdyn.pre.nonlinearSolver')
                
                error ('NonlinearSolver must be an object derived from the mbdyn.pre.nonlinearSolver class');
                
            end
            
            self.initialTime = itime;
            self.finalTime = ftime;
            self.timeStep = tstep;
            self.maxIterations = options.MaxIterations;
            self.tolerance = [ ResidualTolerance, SolutionTolerance ];
            self.derivativesTolerance = options.DerivativesTolerance;
            self.output = options.Output;
            self.nonlinearSolver = options.NonlinearSolver;
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
            
            if ~isempty (self.tolerance)
                str = self.addOutputLine (str, sprintf ('tolerance: %s;', self.commaSepList (self.tolerance{:})), 1, false);
            end
            
            if ~isempty (self.maxIterations)
                str = self.addOutputLine (str, sprintf('max iterations: %d;', self.maxIterations), 1, false);
            end
            
            if ~isempty (self.nonlinearSolver)
                str = self.addOutputLine (str, sprintf('nonlinear solver: %s ;', self.nonlinearSolver.generateOutputString ()), 1, false);
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
    
    methods (Access = protected)
        
        function ResidualTolerance = checkResidualTolerance (self, tolerance)
            
            if isempty (tolerance)
                
                % return empty cell array
                ResidualTolerance = {};
                
            else
                
                if ~iscell (tolerance)
                    tolerance = {tolerance};
                end

                if ~( (ischar (tolerance{1}) && strcmp (tolerance{1}, 'null')) ...
                        || self.checkNumericScalar (tolerance{1}, false, 'Tolerance') )

                    error ('First ResidualTolerance value must be a numeric scalar value (or the keyword ''null''');

                end

                if numel (tolerance) > 1
                    if strcmp (tolerance{2}, 'test')

                        self.checkAllowedStringInputs ( tolerance{3}, ...
                           { 'none', 'norm', 'minmax' }, true, 'residual tolerance test');

                       if numel (tolerance) > 3
                           if ~strcmp (tolerance{4}, 'scale')
                               error ('The last value in the residual tolerance specification is not the keyword ''scale''');
                           end
                       end

                    else
                        error ('Missing ''test'' keyword in residual tolerance cell array')
                    end
                end

                ResidualTolerance = tolerance;
            end
            
        end
        
        function SolutionTolerance = checkSolutionTolerance (self, tolerance)
            
            if isempty (tolerance)
                
                % return empty cell array
                SolutionTolerance = {};
                
            else
                
                if ~iscell (tolerance)
                    tolerance = {tolerance};
                end

                if ~( (ischar (tolerance{1}) && strcmp (tolerance{1}, 'null')) ...
                        || self.checkNumericScalar (tolerance{1}, false) )

                    error ('First SolutionTolerance value must be a numeric scalar value (or the keyword ''null''');

                end

                if numel (tolerance) > 1
                    if strcmp (tolerance{2}, 'test')

                        self.checkAllowedStringInputs ( tolerance{3}, ...
                           { 'none', 'norm', 'minmax' }, true, 'solution tolerance test' );

                       if numel (tolerance) > 3
                           error ('The SolutionTolerance has too many members, cell array length is greater than 3.');
                       end

                    else
                        error ('Missing ''test'' keyword in SolutionTolerance cell array')
                    end
                end

                SolutionTolerance = tolerance;
            end
            
        end
        
    end
    
end