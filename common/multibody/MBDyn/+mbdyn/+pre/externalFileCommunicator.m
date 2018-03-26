classdef externalFileCommunicator < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        sleepTime;
        coupling;
        precision;
        sendAfterPredict;
        commMethod;
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
                    
                elseif ~self.checkScalarInteger (options.Coupling, false)
                    
                    error ('Coupling must be a string or integer number of steps');
                    
                end
            end
            
            if ~( ( isnumeric (options.SleepTime) ...
                        && isscalar (options.SleepTime) ...
                        && options.SleepTime >= 0 ) ...
                    || isempty (options.SleepTime) )
                
                error ('SleepTime must be a scalar numeric value >= 0');
                
            end
            
            if ~( ( self.checkScalarInteger (options.Precision, false) ) ...
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
        
        function str = generateMBDynInputString (self)
            
            str = sprintf ('%s,', self.type);
            
        end
        
%         function comminfo = commInfo (self)
%             % gets communication info for the socket communicator.
%             %
%             %
%             
%             comminfo.commMethod = self.commMethod;
%             
%             if strcmp (self.commMethod, 'local socket')
%                 comminfo.path = self.path;
%             elseif strcmp (self.commMethod, 'inet socket')
%                 comminfo.host = self.host;
%                 comminfo.port = self.port;
%             else
%                 error ('unrecognised communication type');
%             end
% 
%         end
        
    end
    
    
end