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

            
    properties (GetAccess = public, SetAccess = private)
        
        fixed_dt;  % fixed time step value for when time step is not provided
        maxOutput; % maximum possible output value from the controller
        minOutput; % minimum possible output value from the controller
        Kp;        % proportional control gain 
        Kd;        % derivative control gain
        Ki;        % integral control gain
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        error; % current value of the error 
        prevError; % previously calculated value of the error signal
        integral;  % current value of the input signal integral
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
            
            options = parse_pv_pairs (options, varargin);
            
            % input checking
            check.isNumericScalar (options.MaxOut, true, 'MaxOut');
            check.isNumericScalar (options.MinOut, true, 'MinOut');
            check.isNumericScalar (Kp, true, 'Kp');
            check.isNumericScalar (Kd, true, 'Kd');
            check.isNumericScalar (Ki, true, 'Ki');
            assert (options.MaxOut > options.MinOut, 'MaxOut must be greater than MinOut');
            
            self.maxOutput = options.MaxOut;
            self.minOutput = options.MinOut;
            self.Kp = Kp;
            self.Kd = Kd;
            self.Ki = Ki;
            
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
            %
            % See Also:
            %
            
            if nargin < 2
                t = self.initT;
            end
            
            self.prevT = t;
            self.prevError = 0;
            self.integral = 0;
            
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

            if nargin < 4
                dt = self.fixed_dt;
            end
            
            % Calculate error
            self.error = setpoint - actual;
            
            % Proportional term
            Pout = self.Kp * self.error;
            
            % Integral term
            self.integral = self.integral + self.error * dt;
            Iout = self.Ki * self.integral;
            
            % Derivative term
            derivative = (self.error - self.prevError) / dt;
            Dout = self.Kd * derivative;
            
            % Calculate total output
            output = Pout + Iout + Dout;
            
            % Restrict to max/min
            if ( output > self.maxOutput )
                output = self.maxOutput;
            elseif ( output < self.minOutput )
                output = self.minOutput;
            end
            
            
            
        end
        
        function advanceStep (self, t)
            
            if nargin < 2
                t = nan;
            end
            
            self.prevT = t;
            
            % Save error to previous error
            self.prevError = self.error;
            
        end
        
        function dt = calcDt (self, t)
            % calculates the step since the previous time
            
            dt = t - self.prevT;
            
        end
        
    end
    
%     methods (Static)
%        
%         function status = outputFcn (t, y, flag)
%             
%             switch flag
%                 
%                 case []
%                     
%                     
%                     
%                 case 'init'
%                     
%                 case 'done'
%                     
%             end
%             
%         end
%         
%     end
    
end