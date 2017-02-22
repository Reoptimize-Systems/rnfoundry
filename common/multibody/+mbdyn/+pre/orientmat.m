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
                    
                    om = mbdyn.pre.orientmat ('euler123', spec);
                    this.orientationMatrix = om.orientationMatrix;
                    clear om;
                    
                case 'euler123'
                    
                    d = spec(1);
                    dCosAlpha = cos(d);
                    dSinAlpha = sin(d);
                    d = spec(2);
                    dCosBeta = cos(d);
                    dSinBeta = sin(d);
                    d = spec(3);
                    dCosGamma = cos(d);
                    dSinGamma = sin(d);
                    
                    this.orientationMatrix = ...
                        [ dCosBeta*dCosGamma, dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma, dSinAlpha*dSinGamma - dCosAlpha*dSinBeta*dCosGamma;
                          -dCosBeta*dSinGamma, dCosAlpha*dCosGamma - dSinAlpha*dSinBeta*dSinGamma, dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma;
                          dSinBeta, -dSinAlpha*dCosBeta, dCosAlpha*dCosBeta ];
                    
                case 'euler321'
               
                    d = spec(1);
                    dCosAlpha = cos(d);
                    dSinAlpha = sin(d);
                    d = spec(2);
                    dCosBeta = cos(d);
                    dSinBeta = sin(d);
                    d = spec(3);
                    dCosGamma = cos(d);
                    dSinGamma = sin(d);

                    this.orientationMatrix = ...
                        [ dCosAlpha*dCosBeta,  dSinAlpha*dCosBeta, -dSinBeta;
                         -dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma, dCosAlpha*dCosGamma + dSinAlpha*dSinBeta*dSinGamma, dCosBeta*dSinGamma;
                         dSinAlpha*dSinGamma + dCosAlpha*dSinBeta*dCosGamma, -dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma, dCosBeta*dCosGamma ];
                    
                case 'euler313'
                    
                    d = spec(1);
                    dCosAlpha = cos(d);
                    dSinAlpha = sin(d);
                    d = spec(2);
                    dCosBeta = cos(d);
                    dSinBeta = sin(d);
                    d = spec(3);
                    dCosGamma = cos(d);
                    dSinGamma = sin(d);

                    this.orientationMatrix = ...
                        [ dCosAlpha*dCosGamma - dSinAlpha*dCosBeta*dSinGamma, dSinAlpha*dCosGamma + dCosAlpha*dCosBeta*dSinGamma, dSinBeta*dSinGamma; 
                          -dCosAlpha*dSinGamma - dSinAlpha*dCosBeta*dCosGamma, -dSinAlpha*dSinGamma + dCosAlpha*dCosBeta*dCosGamma, dSinBeta*dCosGamma;
                          dSinAlpha*dSinBeta, -dCosAlpha*dSinBeta, dCosBeta ];
                      
                case 'vector'
                    % axis and angle (angle in rad = norm of matrix)
                    wcrs = [ 0         spec(3) -spec(2)
                            -spec(3)        0   spec(1)
                             spec(2)  -spec(1)       0] ;   
             
                    this.orientationMatrix = expm (wcrs);
                    
                case '2vectors'
                    
                    if (spec.ia < 1 || spec.ia > 3)
                        error ('Axis index A must be 1, 2, or 3')
                    end
                    
                    r = cell([1,3]);
                    
                    i1 = spec.ia - 1;
                    i2 = mod (spec.ia, 3);
                    i3 = mod ((spec.ia+1), 3);
                    
                    i1 = i1 + 1;
                    i2 = i2 + 1;
                    i3 = i3 + 1;

                    if (spec.ib == mod (spec.ia,3)+1)
                        d = norm(spec.vecA);
                        if (d <= eps())
                        	error ('first vector must be non-null');
                        end
                        r{i1} = spec.vecA / d;
                        d = norm(spec.vecB);
                        if (d <= eps())
                        	error ('second vector must be non-null');
                        end
                        r{i3} = cross(r{i1}, spec.vecB);
                        d = dot(r{i3}, r{i3});
                        if (d <= eps())
                        	error ('vectors must be distinct');
                        end
                        d = sqrt(d);
                        r{i3} = r{i3} / d;
                        r{i2} = cross(r{i3}, r{i1});

                    elseif (spec.ib == (mod (spec.ia+1,3)+1))
                        d = norm(spec.vecA);
                        if (d <= eps())
                        	error ('first vector must be non-null');
                        end
                        
                        r{i1} = spec.vecA / d;
                        d = norm(spec.vecB);
                        if (d <= eps())
                        	error ('second vector must be non-null');
                        end

                        r{i2} = cross(spec.vecB, r{i1});
                        d = dot(r{i2}, r{i2});
                        if (d <= eps())
                        	error ('vectors must be distinct');
                        end
                        
                        d = sqrt(d);
                        r{i2} = r{i2} / d;
                        r{i3} = cross(r{i1}, r{i2});
                        
                    else 
                        error ('second index is illegal');
                    end
               
                    this.orientationMatrix = [ r{1}.'; r{2}.'; r{3}.' ];
                    
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