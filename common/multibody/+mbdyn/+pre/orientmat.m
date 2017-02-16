classdef orientmat
    
    properties (GetAccess = public, SetAccess = protected)
        
        orientationMatrix;
    end
    
    methods
        
        function this = orientmat (spectype, spec)
        % orentmat constructor
        %
        % Syntax
        %
        % om = orientmat (spectype, spec)
        %
        % Input
        %
        % 
        
            switch spectype
                
                case 'orientation'
                    
                    this.orientationMatrix = spec;
                    
                case 'euler'
                    
                    this.orientationMatrix = SpinCalc('EA123toDCM', rad2deg (spec), eps (), 1);
                    
                case 'euler123'
                    
                    this.orientationMatrix = SpinCalc('EA123toDCM', rad2deg (spec), eps (), 1);
                    
                case 'euler321'
                    
                    this.orientationMatrix = SpinCalc('EA321toDCM', rad2deg (spec), eps (), 1);
                    
                case 'vector'
                    % axis and angle (angle in rad = norm of matrix)
                    wcrs = [ 0         spec(3) -spec(2)
                            -spec(3)        0   spec(1)
                             spec(2)  -spec(1)       0] ;   
             
                    this.orientationMatrix = expm (wcrs);
                    
                case '2vectors'
                    
                    % normalise the fisr vector
                    spec.vec1 = this.unit (spec.vec1);
                    spec.vec2 = this.unit (spec.vec2);
                    
                    spec.vec3 = cross (spec.vec1, spec.vec2);
                    
                    spec.vec2 = this.unit (cross (this.unit (spec.vec3), spec.vec1));
                    
                    switch spec.vec1axis
                        
                        case 1
                            X = spec.vec1;
                            if spec.vec2axis == 2
                                Y = spec.vec2;
                                Z = spec.vec3;
                            elseif spec.vec2axis == 3
                                Y = spec.vec3;
                                Z = spec.vec2;
                            end
                            
                        case 2
                            Y = spec.vec1;
                            if spec.vec2axis == 1
                                X = spec.vec2;
                                Z = spec.vec3;
                            elseif spec.vec2axis == 3
                                X = spec.vec3;
                                Z = spec.vec2;
                            end
                            
                        case 3
                            Z = spec.vec1;
                            if spec.vec2axis == 2
                                X = spec.vec2;
                                Y = spec.vec3;
                            elseif spec.vec2axis == 3
                                X = spec.vec3;
                                Y = spec.vec2;
                            end
                        
                    end
                    
                    this.orientationMatrix = [ X, Y, Z ];
                    
            end 
            
        end
        
    end
    
    % operator overloading
    methods 
        
        function om = plus (om1, om2)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix + om2.orientationMatrix);
            
        end
        
        function om = minus (om1, om2)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix - om2.orientationMatrix);
            
        end
        
        function om = times (om1, om2)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix .* om2.orientationMatrix);
            
        end
        
        function om = mtimes (om1, om2)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix * om2.orientationMatrix);
            
        end
        
        function om = double (om1)
            
            om = om1.orientationMatrix;
            
        end
        
        function om = uminus (om1)
            
            om = mbdyn.pre.orientmat ('orientation', -om1.orientationMatrix);
            
        end
        
        function om = uplus (om1)
            
            om = mbdyn.pre.orientmat ('orientation', +om1.orientationMatrix);
            
        end
        
        function om = transpose (om1)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix.');
            
        end
        
        function om = ctranspose (om1)
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix');
            
        end
        
    end
    
    methods (Access = private)
       
        function out = unit (self, vec)
            
            out = vec ./ norm (vec);
            
        end
        
    end

    
end