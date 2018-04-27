classdef pidController < handle
% class implementing the PID control algorithm
%
% Syntax
%
% self = pidController (Kp, Ki, Kd)
% self = pidController (..., 'Parameter', Value)
%
% Description
%
% pidController implements a proportional-integral-derivative
% controller. Limits on the maximum and minimum values of the
% output signal can be set. To convert to a PI or PD controller
% one can simply set the values of Ki or Kd to zero.
%
% pidController Methods:
%
%  pidController - constructor
%  reset - reset the controller to it's initial state
%  calculate - calculate the pid output signal
%


    properties (GetAccess = public, SetAccess = public)
        
        Kp;        % proportional control gain 
        Kd;        % derivative control gain
        Ki;        % integral control gain
        taui;      % integral control factor
        saturationTimeConstant; % saturation time constant for antiwindup action
        
    end
    
    properties (GetAccess = public, SetAccess = private)
        
        fixed_dt;  % fixed time step value for when time step is not provided
        maxOutput; % maximum possible output value from the controller
        minOutput; % minimum possible output value from the controller
        antiWindup;
        
        error; % current value of the error 
        integral;  % current value of the input signal integral
        saturationError; 
        active;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        prevError; % previously calculated value of the error signal
        prevIntegral; % previous value of the input signal integral
        prevSaturationError;
        saturationIntegral;
        prevSaturationIntegral;
        prevT;
        initT;
        
    end
    
    methods
        
        function self = pidController (Kp, Ki, Kd, varargin)
            % PID controller class constructor
            %
            % Syntax
            %
            % self = pidController (Kp, Ki, Kd)
            % self = pidController (..., 'Parameter', Value)
            %
            % Description
            %
            % pidController implements a proportional-integral-derivative
            % controller. Limits on the maximum and minimum values of the
            % output signal can be set. To convert to a PI or PD controller
            % one can simply set the values of Ki or Kd to zero. 
            %
            % The PID controller class also has features to assist with
            % using it within a simulation being solved using the ode*
            % family of functions (e.g. ode45, ode15s etc). These funcitons
            % typically do not give access to the current time step, so it
            % must be calculated from the previous value to the simulation
            % time, which must be stored internally by the pidController
            % class.
            %
            % pidController Methods:
            %
            %  reset - reset the controller to it's initial state
            %  calculate - calculate the pid output signal
            %
            % Input
            %
            %  Kp - proportional control gain 
            %
            %  Ki - integral control gain
            %
            %  Kd - derivative control gain
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'MaxOut' - maximum possible value of the output. Default is
            %    inf if not supplied.
            %
            %  'MinOut' - mimimum possible value of the output. Default is
            %    -inf if not supplied.
            %
            %  'FixedDT' - for fixed time step applications one can set the
            %    fixed time step using this option. The value is used if
            %    calling the calculate method without providing a time
            %    step. The value provided here is put in the fixed_dt
            %    property of the object, and it's default value is NaN (Not
            %    a Number) if not supplied.
            %
            % Output
            %
            %  pidobj - pidController object
            %
            %
            % See Also:
            %
            
            options.MaxOut = inf;
            options.MinOut = -inf;
            options.FixedDT = nan;
            options.InitialTime = [];
            options.SaturationTimeConstant = [];
            options.AntiWindup = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % input checking
            check.isNumericScalar (options.MaxOut, true, 'MaxOut');
            check.isNumericScalar (options.MinOut, true, 'MinOut');
            check.isNumericScalar (Kp, true, 'Kp');
            check.isNumericScalar (Kd, true, 'Kd');
            check.isNumericScalar (Ki, true, 'Ki');
            assert (options.MaxOut > options.MinOut, 'MaxOut must be greater than MinOut');
            check.isLogicalScalar (options.AntiWindup, true, 'AntiWindup');
            
            self.maxOutput = options.MaxOut;
            self.minOutput = options.MinOut;
            self.Kp = Kp;
            self.Kd = Kd;
            self.Ki = Ki;
            
            self.taui = self.Kp / self.Ki;
            
            self.antiWindup = options.AntiWindup;
            
            if self.antiWindup == false
                self.saturationTimeConstant = inf;
            else
                
                if isempty (options.SaturationTimeConstant)
                    % calcuate the saturation time constant
                    options.SaturationTimeConstant = self.taui / 2;
                end
                
                check.isNumericScalar (options.SaturationTimeConstant, true, 'SaturationTimeConstant');
                
                self.saturationTimeConstant = options.SaturationTimeConstant;
                
            end
            
            if ~isempty (options.InitialTime)
                check.isNumericScalar (options.InitialTime, true, 'InitialTime');
                self.prevT = options.InitialTime;
                self.initT = options.InitialTime;
            end
            
            self.reset ();
            
        end
        
        function reset (self, t)
            % reset the controller to initial values
            %
            % Syntax
            %
            % reset (pidobj)
            %
            % Description
            %
            % Resets the PID controller to it's initial state, this sets
            % both the internally calculated integral and stored value of
            % the previously calculated error to zero.
            %
            % Input
            %
            %  pidobj - pidController object
            %
            %  t - optional new value for the initial time of the PID
            %   controller. If not supplied, the intial time supplied (if
            %   any) when creating the controller is used. If it is
            %   supplied, the value also replaces the stored initial time
            %   value (so future calls to reset without supplying a new
            %   value for t, will reset to the last value supplied here).
            %   The time is used in the calcDt method which is used when
            %   not using a fixed time step.
            %
            % See Also: pidController.calcDt
            %
            
            if nargin < 2
                t = self.initT;
            else
                self.initT = t;
            end
            
            self.prevT = t;
            self.error = 0;
            self.prevError = 0;
            self.integral = 0;
            self.prevIntegral = 0;
            self.prevSaturationError = 0;
            self.saturationIntegral = 0;
            self.prevSaturationIntegral = 0;
            self.active = true;
            
        end
        
        function toggleOnOff (self)
            self.active = ~self.active;
        end
        
        function output = calculate (self, actual, setpoint, dt)
            % calculate the value of the PID controller output signal
            %
            % Syntax
            %
            % output = calculate (pidobj)
            % output = calculate (pidobj, dt)
            %
            % Description
            %
            % calculate calculates the value of the PID controller output
            % signal by conparing the actual and desired setpioint for a
            % variable under control and applying the PID algorithm. 
            %
            % Input
            %
            %  pidobj - pidController object
            %
            %  actual - current real value of the variable under control
            %
            %  setpoint - target value of the variable under control
            %
            %  dt - optional value of the current time step which is used
            %   to calculate the integral and derivative terms. If dt is
            %   not supplied, the calculate method sill use the value in
            %   the class property fixed_dt. Therefore, if you wish to use
            %   a constant value time step you must set this value when
            %   constructing the pidController object. 
            %
            % Output
            %
            %  output - The output (or manipulated variable) from the PID
            %   control algorithm
            %
            %
            %
            % See Also: 
            %

            if self.active
                
                if nargin < 4
                    dt = self.fixed_dt;
                end

                % Calculate error
                self.error = setpoint - actual;

                % Proportional term
                Pout = self.Kp * self.error;

                % Integral term
                self.integral = self.prevIntegral + self.error * dt;
                Iout = self.Ki * self.integral;

                % Derivative term
                derivative = (self.error - self.prevError) / dt;
                Dout = self.Kd * derivative;

                % Calculate total output
                total_action = Pout + Iout + Dout;

                % Restrict to max/min
                if ( total_action > self.maxOutput )
                    output = self.maxOutput;
                elseif ( total_action < self.minOutput )
                    output = self.minOutput;
                else
                    output = total_action;
                end

                if self.antiWindup == true

                    self.saturationError = output - total_action;

                    self.saturationIntegral = self.prevSaturationIntegral + ...
                                  (1/self.saturationTimeConstant) * self.saturationError * dt;

                    output = output + self.saturationIntegral;

                end
            
            else
                output = 0;
            end
        end
        
        function advanceStep (self, t)
            % accept the current values into the time history of PID states
            %
            % Syntax
            %
            % advanceStep (pidobj, t)
            % Description
            %
            % Accepts the current value of the error and integral into the
            % PID history. This must be called before advancing to the next
            % time point. 
            %
            % Input
            %
            %  pidobj - pidController object
            %
            %  t - current simulation time
            %
            %
            
            if nargin < 2
                t = nan;
            end
            
            self.prevT = t;
            
            % Save current error to previous error
            self.prevError = self.error;
            
            % update the integral
            self.prevIntegral = self.integral;
            
            self.prevSaturationError = self.saturationError;
            
            self.prevSaturationIntegral = self.saturationIntegral;
            
        end
        
        function dt = calcDt (self, t)
            % calculates the step since the previous time
            
            dt = t - self.prevT;
            
        end
        
    end
    
end