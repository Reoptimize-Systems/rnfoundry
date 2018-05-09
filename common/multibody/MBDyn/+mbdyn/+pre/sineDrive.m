classdef sineDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        numberOfCycles;
        initialTime;
        omega;
        amplitude;
        initialValue;
    end
    
    methods
        
        function self = sineDrive (t_init, omega, amplitude, num_cycles, initial_val)
            % Constructor for sine drive
            %
            % Syntax
            %
            % sd = mbdyn.pre.sineDrive (t_init, omega, amplitude, num_cycles, initial_val)
            %
            % Description
            %
            % where angular_velocity is 2*pi/T. This drive actually computes
            %
            %  f(t) = initial_value + amplitude * sin (angular_velocity * (t − initial_time))
            %
            % The value of number_of_cycles determines the behavior of the
            % drive. If it is positive, num_cycles-1/2 oscillations
            % are performed. If it is negative, the oscillations end after
            % num_cycles-3/4 cycles at the top of the sine, with null
            % tangent. Special keywords can be used for num_cycles.
            %
            % Input
            %
            %  t_init - time when sine wave starts operating
            %
            %  omega - frequency in rad/s
            %
            %  amplitude - mean to peak amplitude of the sime wave
            %
            %  num_cycles - either a scalar value or a character vector
            %   which can be one of:
            %
            %   'forever' : the oscillation never stops
            %   
            %   'one' : exactly half period is performed (equivalent to number_of_cycles = 1);
            %
            %   'half' : exactly a quarter of period is performed
            %     (equivalent to number_of_cycles = -1), so the function
            %     stops at
            %
            %  initial_val - offset of sine wave (mean value)
            %
            % Output
            %
            %  sd - mbdyn.pre.sineDrive object
            %
            %
            %
            % See Also: 
            %
            
            self.checkNumericScalar (t_init, true, 't_init');

            self.checkNumericScalar (omega, true, 'omega');
            self.checkNumericScalar (amplitude, true, 'amplitude');
            
            if ischar (num_cycles)
                ok = self.checkAllowedStringInputs (num_cycles, {'forever', 'one', 'half'}, false);
            else
                ok = self.checkNumericScalar (t_final, false, 't_final');
            end
            
            assert (ok, 'num_cycles msut be a scalar numeric value of a character vector, ''forever'' | ''one'' | ''half''');
            
            self.checkNumericScalar (initial_val, true, 'initial_val');
            
            self.type = 'sine';
            self.numberOfCycles = num_cycles;
            self.initialTime = t_init;
            self.omega = omega;
            self.amplitude = amplitude;
            self.initialValue = initial_val;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = [ self.type, ',' ];
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList ( self.initialTime, ...
                                                           self.omega, ...
                                                           self.amplitude, ...
                                                           self.numberOfCycles, ...
                                                           self.initialValue ), ...
                                       1, ...
                                       false );
            
        end
        
    end
    
end