classdef gravity < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = private)
        gravityAcceleration;
        gravityModel;
        origin;
        mass
        G;
%         origin;
    end
    
    methods
        
        function self = gravity (varargin)
            
            options.GravityModel = 'uniform';
            options.GravityAcceleration = mbdyn.pre.tplDriveCaller ('single', mbdyn.pre.const (9.81), 'Direction', [0;0;-1]);
            options.Origin = [0;0;0];
            options.Mass = nan;
            options.G = 'si';
            
            options = parse_pv_pairs (options, varargin);
            
            switch options.GravityModel
                
                case 'uniform'
                    
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
                    
                case 'central'
                    
                    self.checkCartesianVector (options.Origin, true);
                     
                    if ~( (isscalar (options.G) && isnumeric (options.G)) ...
                            || (ischar (options.G) && ~strcmp(options.G, 'si')) ...
                        )
                    
                        error ('G must be the keyword ''si'' or a scalar numeric value')
                    end
                    
                    if ~(isscalar (options.Mass) && isnumeric (options.Mass))
                    
                        error ('Mass a scalar numeric value')
                    end
                    
                    self.origin = options.Origin;
                    self.mass = options.Mass;
                    self.G = options.G;
                    
                otherwise
                    
                    error ('Gravity model must be ''uniform'', or ''central''');
            end
            
            self.gravityModel = options.GravityModel;

        end
        
        function str = generateOutputString (self)
            
            str = '';
            
            switch self.gravityModel
                
                case 'uniform'
                    str = self.addOutputLine (str , sprintf ('gravity : %s, %s;', ...
                                                     self.gravityModel, ...
                                                     self.gravityAcceleration.generateOutputString() ...
                                                            ), ...
                                              1, false);
                    
                case 'central'
                    str = self.addOutputLine (str , sprintf ('gravity : %s, %s;', ...
                                                     self.gravityModel, ...
                                                     self.commaSepList ('origin', self.origin, 'mass', self.mass, 'G', self.G) ...
                                                            ), ...
                                              1, false);
                    
            end
            
            str = str (2:end);
            
        end
        
        function draw (self, varargin)
            
           % do nothing
           
        end
    end
    
end