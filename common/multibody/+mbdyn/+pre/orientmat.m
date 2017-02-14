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

    
end