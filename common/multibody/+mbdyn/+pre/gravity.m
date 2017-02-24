classdef gravity < mbdyn.pre.element
    
    properties
        gravityAcceleration;
        gravityModel;
%         origin;
    end
    
    methods
        
        function self = gravity (varargin)
            
            options.GravityModel = 'uniform';
            options.GravityAcceleration = mbdyn.pre.tplDriveCaller ('single', mbdyn.pre.const (9.81), 'Direction', [0;0;-1]);
            options.Origin = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if isnumeric (options.GravityAcceleration) ...
                if numel (options.GravityAcceleration) == 4
                    dirvec = options.GravityAcceleration(3:end);
                    dirvec = dirvec(:);
                elseif numel (options.GravityAcceleration) == 1
                    dirvec = [ 0; 0; -1];
                else
                    error ('GravityAcceleration must be a scalar value, a 4 element vector or a mbdyn.pre.drive');
                end
                gaccel = options.GravityAcceleration(end);
                
                options.GravityAcceleration = ...
                    mbdyn.pre.tplDriveCaller ('single', mbdyn.pre.const (gaccel), 'Direction', dirvec);
            end
            
            if isa (options.GravityAcceleration, 'mbdyn.pre.tplDriveCaller')
                self.gravityAcceleration = options.GravityAcceleration;
            else
                error ('GravityAcceleration must be a template drive caller')
            end

        end
        
        function str = generateOutputString (self)
            
            str = '';
            
            str = self.addOutputLine (str , sprintf ('gravity : %s;', self.gravityAcceleration.generateOutputString()), 1, false);
            
            str = str (2:end);
            
        end
        
        function draw (self, varargin)
            
           % do nothing
           
        end
    end
    
end