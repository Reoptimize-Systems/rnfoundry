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
        %  spectype - string stating how the orientation matrix is to be
        %    specified, possible options are: 'orientation', 'euler',
        %    'euler123', 'euler321', 'euler313', 'vector' and '2vectors'
        %
        %  spec - orientation, the value of spec is dementdent on the value
        %    of spectype with the following possibilities:
        %
        %    'orientation' - spec should be a (3 x 3) full orientation
        %      matrix aka direction-cosine aka rotation matrix
        %
        %    'euler' - same as for 'euler123'
        %
        %    'euler123' - spec should be a 3 element vector representing
        %      the euler angles using the 123 convention
        %
        %    'euler321' - spec should be a 3 element vector representing
        %      the euler angles using the 321 convention
        %
        %    'euler313' - spec should be a 3 element vector representing
        %      the euler angles using the 313 convention
        %
        %    'vector' - spec should be a 3 element vector. The rotation is
        %      the magnitude of the vector and is performed around the axis
        %      formed by the vector direction.
        %
        %    '2vectors' - spec should be a structure containing the fields
        %      'ia', 'vecA', 'ib' and 'vecB'.
        %
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
                    
                    % make sure vectors etc are good
                    if numel (spec.vecA) ~= 3 ...
                            || numel (spec.vecB) ~= 3 ...
                            || ~(isscalar (spec.ia) && isnumeric(spec.ia) && isint2eps (spec.ia)) ...
                            || ~(isscalar (spec.ib) && isnumeric(spec.ib) && isint2eps (spec.ib))
                        error ('ia and ib must be integers and vecA and vecB must be 3 element vectors');
                    end
                    
                    spec.vecA = spec.vecA(:);
                    spec.vecB = spec.vecB(:);
                    
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
                otherwise
                    
                    error ('unrecognised specification');
            end 
            
        end
        
    end
    
    methods
        
       function eul = euler123 (this)
           % returns the extrinsic euler123 angles corresponding to the
           % orientation matrix

           alpha = -atan2 ( this.orientationMatrix(2,3) , this.orientationMatrix(3,3) );
           
           beta = atan2 ( this.orientationMatrix(1,3) ...
                            , ( cos (alpha)*this.orientationMatrix(3,3) - sin (alpha)*this.orientationMatrix(2,3) ) );
           
           gamma = atan2 ( ( cos (alpha)*this.orientationMatrix(2,1) - sin (alpha)*this.orientationMatrix(3,1) ) ...
                           , ( cos (alpha)*this.orientationMatrix(2,2) - sin (alpha)*this.orientationMatrix(3,2) ) );
           
           eul = [ alpha;
                   beta;
                   gamma ];
       end
       
       function [hquiv, hax] = draw (this, varargin)
           
           options.PlotAxes = [];
           options.Title = true;
           options.Parent = [];
           options.Offset = [];
           options.DrawGlobal = true;
           options.Scale = 1;
           
           options = parse_pv_pairs (options, varargin);
           
           if isempty (options.PlotAxes)
               figure;
               hax = axes;
           else
               hax = options.PlotAxes;
           end
           
           x = options.Scale * [1;0;0];
           y = options.Scale * [0;1;0];
           z = options.Scale * [0;0;1];
           
           ox = this.orientationMatrix * (0.5*x);
           oy = this.orientationMatrix * (0.5*y);
           oz = this.orientationMatrix * (0.5*z);
           
           % orientation frame
           hquiv(1) = vect.plotvec3 (ox, options.Offset, 'Properties', {'Color', 'r', 'LineWidth', 4, 'LineStyle', ':'}, 'PlotAxes', hax);           
           hquiv(2) = vect.plotvec3 (oy, options.Offset, 'Properties', {'Color', 'g', 'LineWidth', 4, 'LineStyle', ':'}, 'PlotAxes', hax);           
           hquiv(3) = vect.plotvec3 (oz, options.Offset, 'Properties', {'Color', 'b', 'LineWidth', 4, 'LineStyle', ':'}, 'PlotAxes', hax);
           
           if options.DrawGlobal
               % global frame
               hquiv(4) = vect.plotvec3 (x, options.Offset, 'Properties', {'Color', 'r'}, 'PlotAxes', hax);
               hquiv(5) = vect.plotvec3 (y, options.Offset, 'Properties', {'Color', 'g'}, 'PlotAxes', hax);
               hquiv(6) = vect.plotvec3 (z, options.Offset, 'Properties', {'Color', 'b'}, 'PlotAxes', hax);
           end
           
           xlabel (hax, 'x');
           ylabel (hax, 'y');
           zlabel (hax, 'z');
           
           set (hax, 'XLim', [-1.1, 1.1]);
           set (hax, 'YLim', [-1.1, 1.1]);
           set (hax, 'ZLim', [-1.1, 1.1]);
           
           if options.Title
               title ('Orientation Matrix Plot')
           end
           
           axis equal;
           
           view (3);
           
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