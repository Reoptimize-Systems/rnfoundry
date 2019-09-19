classdef timeIntegrator < handle
% class implementing the PID control algorithm
%
% Syntax
%
% self = timeIntegrator (Kp, Ki, Kd)
% self = timeIntegrator (..., 'Parameter', Value)
%
% Description
%
% timeIntegrator integrates a quantity in time
%      
% timeIntegrator Methods:
%
%  timeIntegrator - constructor
%  reset - reset the controller to it's initial state
%  calculate - calculate the pid output signal
%

    
    properties (GetAccess = public, SetAccess = private)
        
        maxOutput; % maximum possible output value from the controller
        minOutput; % minimum possible output value from the controller
        
        integral;  % current value of the input signal integral
        
        isLimited;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        prevIntegral; % previous value of the input signal integral
        prevT;
        initT;
        
    end
    
    methods
        
        function self = timeIntegrator (varargin)
            % timeIntegrator class constructor
            %
            % Syntax
            %
            % self = timeIntegrator ('Parameter', Value)
            %
            % Description
            %
            % timeIntegrator integrates a quantity in time, generally
            % intended for use with the ode solver routines
            %
            % Input
            %
            % All arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'InitialTime' - 
            %
            %  'IsLimited' - treu/false flag indicating whether to apply
            %    limits to the output of the integration in order to
            %    prevent windup. The limits to be applied are supplied
            %    using the
            %
            %  'MaxOut' - maximum possible value of the output. Default is
            %    inf if not supplied.
            %
            %  'MinOut' - mimimum possible value of the output. Default is
            %    -inf if not supplied.
            %
            %
            % Output
            %
            %  tiobj - timeIntegrator object
            %
            %
            % See Also:
            %
            
            options.MaxOut = inf;
            options.MinOut = -inf;
            options.InitialTime = [];
            options.IsLimited = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % input checking
            check.isNumericColVector (options.MaxOut, true, 'MaxOut');
            check.isNumericColVector (options.MinOut, true, 'MinOut');
            assert (all (options.MaxOut > options.MinOut), 'MaxOut must be greater than MinOut');
            check.isLogicalScalar (options.IsLimited, true, 'IsLimited');
            
            self.maxOutput = options.MaxOut;
            self.minOutput = options.MinOut;
            self.isLimited = options.IsLimited;
            
            if ~isempty (options.InitialTime)
                check.isNumericScalar (options.InitialTime, true, 'InitialTime');
                self.prevT = options.InitialTime;
                self.initT = options.InitialTime;
            end
            
            self.reset ();
            
        end
        
        function newval = toggleIsLimited (self)
            
            self.isLimited = ~self.isLimited;
            
            newval = self.isLimited;
            
        end
        
        function reset (self, t)
            % reset the controller to initial values
            %
            % Syntax
            %
            % reset (tiobj)
            %
            % Description
            %
            % Resets the timeIntegrator to it's initial state, this sets
            % the internally calculated integral to zero.
            %
            % Input
            %
            %  tiobj - timeIntegrator object
            %
            %  t - optional new value for the initial time of the 
            %   timeIntegrator. If not supplied, the intial time supplied
            %   (if any) when creating the timeIntegrator is used. If it is
            %   supplied, the value also replaces the stored initial time
            %   value (so future calls to reset without supplying a new
            %   value for t, will reset to the last value supplied here).
            %
            % See Also: timeIntegrator.calcDt
            %
            
            if nargin < 2
                t = self.initT;
            else
                self.initT = t;
            end
            
            self.prevT = t;
            self.integral = 0;
            self.prevIntegral = 0;
            
        end
        
        function integral = calculate (self, nextval, t)
            % calculate the value of the integral at the current time
            %
            % Syntax
            %
            % output = calculate (tiobj)
            % output = calculate (tiobj, dt)
            %
            % Description
            %
            % calculate calculates the value of the PID controller output
            % signal by conparing the actual and desired setpioint for a
            % variable under control and applying the PID algorithm. 
            %
            % Input
            %
            %  tiobj - timeIntegrator object
            %
            %  dt - optional value of the current time step which is used
            %   to calculate the integral and derivative terms. If dt is
            %   not supplied, the calculate method sill use the value in
            %   the class property fixed_dt. Therefore, if you wish to use
            %   a constant value time step you must set this value when
            %   constructing the timeIntegrator object. 
            %
            % Output
            %
            %  output - The calculated value of the integration
            %
            %
            %
            % See Also: 
            %

            dt = t - self.prevT;

            integral = self.prevIntegral + nextval .* dt;
            
            if self.isLimited
                integral = max (min (integral, self.maxOutput), self.minOutput);
            end
            
            self.integral = integral;
            
        end
        
        function advanceStep (self, t)
            % accept the current values into the time history of the integration
            %
            % Syntax
            %
            % advanceStep (tiobj, t)
            %
            % Description
            %
            % Accepts the current value of the integral into the history.
            % This must be called before advancing to the next time point.
            %
            % Input
            %
            %  tiobj - timeIntegrator object
            %
            %  t - current simulation time
            %
            %
            
            if nargin < 2
                t = nan;
            end
            
            self.prevT = t;

            % update the integral
            self.prevIntegral = self.integral;

        end
        
    end
    
end