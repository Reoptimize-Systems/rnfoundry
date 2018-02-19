classdef initialValueProblem < mbdyn.pre.problem
    
    properties
        initialTime;
        finalTime;
        timeStep;
        maxIterations;
        tolerance;
        method;
        derivativesTolerance;
        output;
        nonlinearSolver;
        linearSolver;
    end
    
    methods
        
        function self = initialValueProblem (itime, ftime, tstep, varargin)
            % constructs an mbdyn.pre.initialValueProblem object
            %
            % Syntax
            %
            % ivp = mbdyn.pre.initialValueProblem (itime, ftime, tstep)
            % ivp = mbdyn.pre.initialValueProblem (..., 'Parameter', value)
            %
            % Description
            %
            % mbdyn.pre.initialValueProblem creates an initial value
            % problem description object for use in setting up an MBDyn
            % multibody simulation. This solves initial value problems by
            % means of generic integration schemes that can be cast in a
            % broad family of multistep and, experimentally, Implicit Runge
            % Kutta-like schemes. For more details on any of the input
            % options described below, see the corresponding section of the
            % MBDyn manual.
            %
            % Input
            %
            %  itime - intial time of simulation
            %
            %  ftime - final time of simulation
            %
            %  tstep - The initial time step. This value is used throughout
            %   the simulation unless some variable time step strategy is
            %   defined (not currently implemented in this matlab
            %   interface).
            %
            % Additional optional inputs may be specified using
            % parameter-value pairs. The avaialable options are:
            %
            %  'MaxIterations' - Number of iterations to attempt without
            %    passing the convergence test before erroring out
            %
            %  'ResidualTolerance' - the tolerance used for the residual
            %    test.
            %
            %  'SolutionTolerance' - A tolerance to test the solution (the
            %    difference between the states at two iterations)
            %
            %  'Method' - cell array containing the time stepping method
            %    settings. See the MBDyn manual for more information on the
            %    various options. The following options are possible:
            %
            %    crank nicolson 
            %
            %    In this case Method must be a scalar cell array containing
            %    just the string 'crank nicolson'.
            %
            %    ms | hope
            %
            %    Method must be a cell array of length 2 or three. The
            %    first element must be either the string 'ms' or the string
            %    'hope' The second element must be an mbdyn.pre.drive
            %    object which returns the differential_radius of the
            %    algorithm. The optional third cell is also an
            %    mbdyn.pre.drive which returns the algebraic_radius of the
            %    algorithm
            %
            %    third order
            %
            %    Method must a cell array of length 2. The fisrt cell must
            %    contain the string 'third order'. The second must contain
            %    either the string 'ad hoc' or an mbdyn.pre.drive object
            %    which returns the differential_radius of the algorithm.
            %
            %    bdf
            %
            %    Method must be a cell array of length 1 or 3. The first
            %    cell must contain the sting 'bdf', If the second and thrid
            %    cell are supplied, they must be the string 'order' and a
            %    scalar integer giving the order of the algorithm.
            %    
            %    implicit euler
            %
            %    In this case Method must be a scalar cell array containing
            %    just the string 'implicit euler'
            %
            %  'DerivativesTolerance' - Right after the initial assembly
            %    and before the simulation starts, the so-called
            %    derivatives solution is performed. The system is solved
            %    with the kinematic unknowns constrained, in order to
            %    properly determine the dynamic unknowns, namely momenta
            %    and constraint reactions. For this purpose, the
            %    coefficient that relates the state perturbation to the
            %    derivative perturbation must be set to a value that is
            %    small enough to allow the determination of accurate
            %    derivatives with very small change in the states. This
            %    coefficient should be zero, but this leads to matrix
            %    singularity, so it must be chosen by the user, since it is
            %    highly problem dependent. A rule-of-thumb is: if the
            %    system has small stiffness and high inertia, the
            %    coefficient can be big, if the system has high stiffness
            %    and small inertia, the coefficient must be small.
            %
            %    The derivatives solution is always performed and cannot be
            %    disabled. If for any reason it does not converge, to
            %    practically disable it one can set a very large tolerance.
            %    Subsequent time steps may start with a less than ideal
            %    initialization, resulting in a rough transient.
            %
            %  'Output' - Cell array of strings requesting special output
            %    related to the solution phase. The following strings may
            %    be supplied: 'iterations', 'residual', 'solution',
            %    'jacobian matrix', 'messages', 'counter', 'bailout',
            %    'matrix condition number', 'solver condition number', 'cpu
            %    time', 'none'. The keyword 'counter' logs on the standard
            %    output a string that summarizes the statistics of the time
            %    step just completed. By default it is off; the same
            %    information, by default, is logged in the .out file. The
            %    keyword 'iterations' logs on stdout a detailed output of the
            %    error at each iteration inside a time step. The keywords
            %    'residual', 'solution' and 'jacobian matrix' log to stdout the
            %    residual, the solution and the Jacobian matrix; they are
            %    mainly intended for last resort debugging purposes, as the
            %    output can be extremely verbose. The item messages refers
            %    to all messaging on the .out file; by default, it is on.
            %    The special item 'bailout' instructs the solver to output
            %    the residual (as with residual) only in case of no
            %    convergence, before bailing out. The special item 'none'
            %    clears the output flags; the items are additive, so, for
            %    instance, to clear out the default and add the iterations
            %    output, use: {'none', 'iterations'}
            %
            %  'NonlinearSolver' - mbdyn.pre.nonlinearSolver object,
            %    defining the nonlinear solver method to use.
            %
            %  'LinearSolver' - mbdyn.pre.linearSolver object defining the
            %    nonlinear solver to use
            %
            % Output
            %
            %  ivp - mbdyn.pre.initialValueProblem object
            %
            %
            % See Also: mbdyn.pre.system
            %
            
            options.MaxIterations = 20;
            options.ResidualTolerance = 1e-9;
            options.SolutionTolerance = {};
            options.Method = {};
            options.DerivativesTolerance = [];
            options.Output = {};
            options.NonlinearSolver = [];
            options.LinearSolver = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (itime, true, 'itime');
            self.checkNumericScalar (ftime, true, 'ftime');
            self.checkNumericScalar (tstep, true, 'tstep');
            
            ResidualTolerance = self.checkResidualTolerance (options.ResidualTolerance);
            SolutionTolerance = self.checkSolutionTolerance (options.SolutionTolerance);
            
            if ~isempty (options.Method)
                self.checkAllowedStringInputs (options.Method{1}, ...
                    {'crank nicolson', 'ms', 'hope', 'third order', 'bdf', 'implicit euler'}, true, 'Method(1)');
                
                
                switch options.Method{1}
                    
                    case 'crank nicolson'
                        
                        if numel (options.Method) > 1
                            error ('Method length is > 1, but crank nicolson does not support any additional arguments');
                        end
                        
                    case {'ms', 'hope'}
                        
                        if ~any (numel (options.Method) == [2,3])
                            error ('For ms or hope method, Method cell array should be of length 2 or 3');
                        end
                        
                        self.checkDrive (options.Method{2}, true, 'Method(2)');
                        
                        if numel (options.Method) == 3
                            self.checkDrive (options.Method{3}, true, 'Method(2)');
                        end
                        
                    case 'third order'
                        
                        if numel (options.Method) ~= 2
                            error ('For third order method, Method cell array should be of length 2');
                        end
                        
                        if ischar (options.Method{2})
                            self.checkAllowedStringInputs (options.Method{2}, {'ad hoc'}, true, 'Method(2)');
                        elseif isa (options.Method{2}, 'mbdyn.pre.drive')
                            % do nothing
                        else
                            error ('For third order method, Method{2} must be the string ''ad hoc'' or an mbdyn.pre.drive object');
                        end
                        
                    case 'bdf'
                        
                        if numel (options.Method) > 3
                            error ('Method length is > 3, but bdf does not support any additional arguments');
                        end
                        
                        self.checkAllowedStringInputs (options.Method{2}, {'order'}, true, 'Method(2)');
                        
                        if ~ (options.Method{3} == 1 || options.Method{3} == 2)
                            error ('For bdf Method(3) if supplied, must be 1 or 2');
                        end
                        
                    case 'implicit euler'
                        
                        if numel (options.Method) > 1
                            error ('Method length is > 1, but implicit euler does not support any additional arguments');
                        end
                end
            end
            
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
            
            if ~isempty (options.LinearSolver) ...
                    && ~isa (options.LinearSolver, 'mbdyn.pre.linearSolver')
                
                error ('LinearSolver must be an object derived from the mbdyn.pre.linearSolver class');
                
            end
            
            self.initialTime = itime;
            self.finalTime = ftime;
            self.timeStep = tstep;
            self.maxIterations = options.MaxIterations;
            self.tolerance = [ ResidualTolerance, SolutionTolerance ];
            self.method = options.Method;
            self.derivativesTolerance = options.DerivativesTolerance;
            self.output = options.Output;
            self.nonlinearSolver = options.NonlinearSolver;
            self.linearSolver = options.LinearSolver;
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
            
            if ~isempty (self.method)
                
                switch self.method{1}
                    
                    case {'crank nicolson', 'bdf', 'implicit euler'}
                        
                        methodstr = self.commaSepList (self.method{:});
                        
                    case {'ms', 'hope'}
                        
                        args = {self.method{1}, self.method{2}.generateOutputString()};

                        if numel (self.method) == 3
                            args = [args, {self.method{3}.generateOutputString()}];
                        end
                        
                        methodstr = self.commaSepList (args{:});
                        
                    case 'third order'
                        
                        args = {self.method{1}};
                        
                        if ischar (self.method{2})
                            
                            args = [args, {self.method{2}}];
                            
                        elseif isa (options.Method{2}, 'mbdyn.pre.drive')
                            
                            args = [args, {self.method{2}.generateOutputString()}];
                            
                        end
                        
                        methodstr = self.commaSepList (args{:});
                        
                end
                
                str = self.addOutputLine (str, sprintf ('method: %s;', methodstr), 1, false);
                
            end
            
            if ~isempty (self.maxIterations)
                str = self.addOutputLine (str, sprintf('max iterations: %d;', self.maxIterations), 1, false);
            end
            
            if ~isempty (self.nonlinearSolver)
                str = self.addOutputLine (str, sprintf('nonlinear solver: %s ;', self.nonlinearSolver.generateOutputString ()), 1, false);
            end
            
            if ~isempty (self.linearSolver)
                str = self.addOutputLine (str, sprintf('linear solver: %s ;', self.linearSolver.generateOutputString ()), 1, false);
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