classdef gravity < mbdyn.pre.element
% class to add gravity simulation to an MBDyn problem
%
% Syntax
%
% [g] = mbdyn.pre.gravity ()
% [g] = mbdyn.pre.gravity ('Parameter', Value)
%
% Description
%
% mbdyn.pre.gravity adds a gravitational field to an MBDyn
% simulation. Gravity can be a uniform fild or a point source.
% The direction and magnitude can also be controlled.
%
% mbdyn.pre.gravity Methods:
%
%  gravity - constructor, see for details of class options
%  generateOutputString - generate string for MBDyn file
%  draw - called when plotting currently does nothing
%
%
            
            
    properties (GetAccess = public, SetAccess = private)
        gravityAcceleration;
        gravityModel;
        origin;
        mass
        G;
    end
    
    methods
        
        function self = gravity (varargin)
            % mbdyn.pre.gravity constructor
            %
            % Syntax
            %
            % [g] = mbdyn.pre.gravity ()
            % [g] = mbdyn.pre.gravity ('Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.gravity adds a gravitational field to an MBDyn
            % simulation. Gravity can be a uniform fild or a point source.
            % The direction and magnitude can also be controlled.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
            %
            %  'GravityModel' - character vector which sets the tye of
            %    gravitational field. Can be either 'uniform' or 'central'.
            %    Default is 'uniform' if not supplied.
            %
            %  'GravityAcceleration' - The acceleration due to graity. If
            %    not suplied, the default is a contstant value of 9.81.
            %    Otherwise This can be one of:
            %
            %    1. A scalar value, representing a constant field value. In
            %    this case the direction vector of the field is [0;0;-1],
            %    i.e. in the negative z direction.
            %
            %    2. A four element vector, the first three elements are the
            %    direction vector of the field, and the last value is the
            %    acceleration value.
            %
            %    3. A mbdyn.pre.tplDriveCaller of size 3, can also be
            %    mbdyn.pre.singleTplDriveCaller
            %
            %  'Origin' - only relevant if using 'central' as the
            %    GravityModel. In this case Origin is a three element
            %    vector, the location of the origin of the field.
            %
            %  'Mass' - scalar value of the mass producing the field.
            %
            %  'G' - either a scalar numeric value of a string 'si'. If a
            %    numeric value, it is the universal gravity constantto be
            %    used. If the keyword 'si' is used, the value 
            %    6.67384x10^−11 m^3 kg^-1 s^−2 is used
            %
            % Output
            %
            %  g - mbdyn.pre.gravity object
            %
            %
            %
            % See Also: 
            %

            options.GravityModel = 'uniform';
            options.GravityAcceleration = mbdyn.pre.singleTplDriveCaller ([0;0;-1], mbdyn.pre.const (9.81));
            options.Origin = [0;0;0];
            options.Mass = nan;
            options.G = 'si';
            
            options = parse_pv_pairs (options, varargin);
            
            switch options.GravityModel
                
                case 'uniform'
                    
                    if isnumeric (options.GravityAcceleration) ...
                        if numel (options.GravityAcceleration) == 4
                            dirvec = options.GravityAcceleration(1:end-1);
                            dirvec = dirvec(:);
                        elseif numel (options.GravityAcceleration) == 1
                            dirvec = [ 0; 0; -1];
                        else
                            error ('GravityAcceleration must be a scalar value, a 4 element vector or a mbdyn.pre.drive');
                        end
                        gaccel = options.GravityAcceleration(end);

                        options.GravityAcceleration = ...
                            mbdyn.pre.singleTplDriveCaller (dirvec, mbdyn.pre.const (gaccel));
                    end
                    
                    if isa (options.GravityAcceleration, 'mbdyn.pre.driveCaller')
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
                    
                    self.checkNumericScalar (options.Mass, true, 'Mass');
                    
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