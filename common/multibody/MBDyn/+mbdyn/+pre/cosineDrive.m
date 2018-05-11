classdef cosineDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        numberOfCycles;
        initialTime;
        omega;
        amplitude;
        initialValue;
    end
    
    methods
        
        function self = cosineDrive (initial_time, omega, amplitude, num_cycles, varargin)
            % Constructor for cosine drive
            %
            % Syntax
            %
            % cd = mbdyn.pre.cosineDrive (t_init, omega, amplitude, num_cycles, initial_val)
            %
            % Description
            %
            % where angular_velocity is 2*pi/T. This drive actually computes
            %
            %  f(t) = initial_value + amplitude * (1 − cos (angular_velocity * (t − initial_time))).
            %
            % The value of number_of_cycles determines the behavior of the
            % drive. If it is positive, num_cycles-1/2 oscillations
            % are performed. If it is negative, the oscillations end after
            % num_cycles-3/4 cycles at the top of the cosine, with null
            % tangent. Special keywords can be used for num_cycles.
            %
            % Input
            %
            %  t_init - time when cosine wave starts operating
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
            %   'one' : exactly half period is performed (equivalent to
            %     number_of_cycles = 1);
            %
            %   'half' : exactly a quarter of period is performed
            %     (equivalent to number_of_cycles = -1), so the function
            %     stops at
            %
            %   If num_cylces is a numeric value and it is positive,
            %   num_cycles-1/2 oscillations are performed. If it is
            %   negative, the oscillations end after num_cycles-3/4 cycles
            %   at the top of the cosine, with null tangent.
            %
            % Additional arguments may be supplied as parameter-value
            % pairs. The available options are:
            %
            %  'InitialValue' - offset of sine wave from the x axis (i.e. 
            %    the mean value of the sine wave). Default is zero if not
            %    supplied.
            %
            % Output
            %
            %  cd - mbdyn.pre.cosineDrive object
            %
            %
            %
            % See Also: 
            %
            
            options.InitialValue = 0;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (initial_time, true, 't_init');

            self.checkNumericScalar (omega, true, 'omega');
            assert (omega > 0, 'omega must be greater than zero');
            self.checkNumericScalar (amplitude, true, 'amplitude');
            assert (amplitude > 0, 'amplitude must be greater than zero');
            
            if ischar (num_cycles)
                ok = self.checkAllowedStringInputs (num_cycles, {'forever', 'one', 'half'}, false);
            else
                ok = self.checkScalarInteger (num_cycles, false, 'num_cycles');
            end
            
            assert (ok, 'num_cycles must be a scalar integer or a character vector, ''forever'' | ''one'' | ''half''');
            
            self.checkNumericScalar (options.InitialValue, true, 'InitialValue');
            
            self.type = 'cosine';
            self.numberOfCycles = num_cycles;
            self.initialTime = initial_time;
            self.omega = omega;
            self.amplitude = amplitude;
            self.initialValue = options.InitialValue;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for cosineDrive
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (cd)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  cd - mbdyn.pre.cosineDrive object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = [ self.type, ',' ];
            
            if ischar (self.numberOfCycles)
                num_cycles = self.numberOfCycles;
            else
                num_cycles = self.formatInteger (self.numberOfCycles);
            end
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList ( self.initialTime, ...
                                                           self.omega, ...
                                                           self.amplitude, ...
                                                           num_cycles, ...
                                                           self.initialValue ), ...
                                       1, ...
                                       false );
            
        end
        
    end
    
end