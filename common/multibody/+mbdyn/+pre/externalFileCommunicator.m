classdef externalFileCommunicator < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        sleepTime;
        coupling;
        precision;
        sendAfterPredict;
    end
    
    methods
        
        function self = externalFileCommunicator (varargin)
            
            options.SleepTime = [];
            options.Precision = [];
            options.Coupling = 'loose';
            options.SendAfterPredict = 'yes';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Coupling)
                if ischar (options.Coupling)
                    
                    self.checkAllowedStringInputs ( options.Coupling, ...
                        {'staggared', 'loose', 'tight'}, ...
                        true, 'Coupling' );
                    
                elseif ~(isnumeric (options.Coupling) && isscalar (options.Coupling) && isint2eps (options.Coupling))
                    
                    error ('Coupling must be a string or integer number of steps');
                    
                end
            end
            
            if ~( ( isnumeric (options.SleepTime) ...
                        && isscalar (options.SleepTime) ...
                        && options.SleepTime >= 0 ) ...
                    || isempty (options.SleepTime) )
                
                error ('SleepTime must be a scalar numeric value >= 0');
                
            end
            
            if ~( ( isnumeric (options.Precision) && isscalar (options.Precision) && isint2eps (options.Precision) ) ...
                    || isempty (options.Precision) )
                
                error ('Precision must be an integer');
                
            end
            
            if ~isempty (options.SendAfterPredict)
                self.checkAllowedStringInputs (options.SendAfterPredict, {'yes', 'no'}, true, 'SendAfterPredict');
            end
            
            self.sleepTime = options.SleepTime;
            self.coupling = options.Coupling;
            self.precision = options.Precision;
            self.sendAfterPredict = options.SendAfterPredict;
            
        end
        
        function str = generateOutputString (self)
            
            str = sprintf ('%s,', self.type);
            
        end
        
    end
    
    
end