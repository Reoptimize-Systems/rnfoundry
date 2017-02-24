classdef tplDriveCaller < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = protected)
        
        type;
        driveCallers;
        direction;
    end
    
    methods
        
        function self = tplDriveCaller (type, drivecallers, varargin)
            
            switch type
                
                case 'null'
                    
                case 'single'
                    options.Direction = [];
                case 'component'
                    
                case 'array'
                    
            end
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkCartesianVector (options.Direction);
            
            switch type
                
                case 'null'
                    
                case 'single'
                    if numel (drivecallers) ~= 1
                        error ('single type tplDriveCaller should have only drivecallers of length 1');
                    end
                case 'component'
                    if numel (drivecallers) < 1
                        error ('component type tplDriveCaller should have drivecallers of at least length 1');
                    end
                case 'array'
                    
                otherwise
                    error ('tplDriveCaller type must be null | single | component');
            end
            
            self.type = type;
            self.driveCallers = drivecallers;
            self.direction = options.Direction;
            
        end
        
        
        function str = generateOutputString (self)

            % delete newline character and space from start
            
            str = self.type;
            
            switch self.type
                
                case 'null'
                    % do nothing further
                    
                case 'single'
                    if ~isempty (self.direction)
                        str = sprintf ('%s, %s', str, self.commaSepList (self.direction));
                    end
                    str = sprintf ('%s, %s', str, self.driveCallers.generateOutputString ());
                    
                case 'component'
                    
                case 'array'
                    
            end
            
        end
        
    end
    
end