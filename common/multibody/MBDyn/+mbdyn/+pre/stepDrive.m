classdef stepDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        initialTime;
        stepValue;
        initialValue;
    end
    
    methods
        
        function self = stepDrive (t_init, step_val, init_val)
            % Constructor for step drive
            %
            % Syntax
            %
            % sd = mbdyn.pre.stepDrive (t_init, step_val, init_val)
            %
            % Description
            %
            % Drive which increases from an initial value with a given
            % slope over a given time period.
            %
            % Input
            %
            %  t_init - time at which step begins to be applied
            %
            %  step_val - value of the drive after time t_init
            %
            %  init_val - initial value of the drive until t_init
            %
            % Output
            %
            %  sd - mbdyn.pre.stepDrive object
            %
            %
            %
            % See Also: mbdyn.pre.rampDrive
            %
            
            self.checkNumericScalar (t_init, true, 't_init');
            self.checkNumericScalar (step_val, true, 'step_val');
            self.checkNumericScalar (init_val, true, 'init_val');
            
            self.stepValue = step_val;
            self.initialTime = t_init;
            self.initialValue = init_val;
            self.type = 'step';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = [ self.type, ',' ];
            
            str = self.addOutputLine ( str, ...
                self.commaSepList ( self.initialTime, self.stepValue,  self.initialValue), ...
                                       1, ...
                                       false );
            
        end
        
    end
    
end