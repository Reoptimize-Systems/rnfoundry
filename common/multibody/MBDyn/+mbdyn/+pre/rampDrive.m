classdef rampDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        slope;
        initialTime;
        finalTime;
        initialValue;
    end
    
    methods
        
        function self = rampDrive (t_init, t_final, slope, init_val)
            % Constructor for ramp drive
            %
            % Syntax
            %
            % rd = mbdyn.pre.rampDrive (t_init, t_final, slope, init_val)
            %
            % Description
            %
            % Drive which increases from an initial value with a given
            % slope over a given time period.
            %
            % Input
            %
            %  t_init - time at which ramp begins to be applied
            %
            %  t_final - time at which ramp is no longer applied, and the
            %    output remains at whatever value is reached. Can also be a
            %    string 'forever', in which case it never stops changing.
            %
            %  slope - slope with which the ramp changes.
            %
            %  init_val - initial value of the drive from which it
            %   increases  or decreases with slope in the specified times.
            %
            % Output
            %
            %  rd - mbdyn.pre.rampDrive object
            %
            %
            %
            % See Also: 
            %
            
            self.checkNumericScalar (t_init, true, 't_init');
            if ischar (t_final)
                assert (self.checkAllowedStringInputs (t_final, {'forever'}, false), ...
                    'If t_final is a char array, it must be ''forever''');
            else
                self.checkNumericScalar (t_final, true, 't_final');
            end
            self.checkNumericScalar (slope, true, 'slope');
            self.checkNumericScalar (init_val, true, 'init_val');
            
            self.slope = slope;
            self.initialTime = t_init;
            self.finalTime = t_final;
            self.initialValue = init_val;
            self.type = 'ramp';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = [ self.type, ',' ];
            
            str = self.addOutputLine ( str, ...
                self.commaSepList ( self.slope, self.initialTime, self.finalTime,  self.initialValue), ...
                                       1, ...
                                       false );
            
        end
        
    end
    
end