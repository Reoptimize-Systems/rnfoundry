classdef orientmat
    % class representing orientations in 3D space
    %
    % Syntax
    %
    % om = orientmat (spectype, spec)
    %
    % Description
    %
    % The orientmat class is used to represent orientations in 3D space
    % and construct orientation matrices from a variety of methods of
    % describing orientations, e.g. from euler angles, orientation vectors
    % etc. See the hlep for the constructor for full details and examples.
    %
    % Input
    %
    %  spectype - string stating how the orientation matrix is to be
    %    specified, possible options are: 'eye', 'orientation matrix',
    %    'orientation', 'matrix', 'euler', 'euler123', 'euler321',
    %    'euler313', 'orientation vector, 'vector' and '2vectors'.
    %
    %  spec - orientation, the value of spec is dementdent on the value
    %    of spectype with the following possibilities:
    %
    %    'eye' - spec is ignores, an object representing the (3 x 3)
    %      identity matrix is returned. Equivalent to:
    %
    %      om = mbdyn.pre.orientmat ('orientation matrix', eye(3))
    %
    %    'orientation matrix' - spec should be a (3 x 3) full
    %      orientation matrix aka direction-cosine aka rotation matrix
    %
    %    'orientation' - same as 'orientation matrix'
    %
    %    'matrix' - same as 'orientation matrix'
    %
    %    'euler123' - spec should be a 3 element vector representing
    %      the euler angles using the 123 convention
    %
    %    'euler' - same as for 'euler123'
    %
    %    'euler321' - spec should be a 3 element vector representing
    %      the euler angles using the 321 convention
    %
    %    'euler313' - spec should be a 3 element vector representing
    %      the euler angles using the 313 convention
    %
    %    'orientation vector' - spec should be a 3 element vector. The
    %      rotation is the magnitude of the vector and is performed
    %      around the axis formed by the vector direction.
    %
    %    'vector' - same as 'orientation vector'
    %
    %    '2vectors' - spec should be a structure containing the
    %      fields 'ia', 'vecA', 'ib' and 'vecB'. The vecA and vecB fields
    %      are three element vectors, ia and ib are integers of value 1, 2
    %      or 3. vecA in this case is a vector representing the axis of a
    %      three-axis coordinate system with the index given in 'ia'.
    %      'vecB' is another vector which is not parallel to vecA which is
    %      used to define a plane in 3D space. Once the plane is defined,
    %      another axis is defined as the vector orthogonal to vecA on this
    %      3D plane. This axis is assigned to be the index given in 'ib' in
    %      the 3D coordinate system. The remaining axis is the vector
    %      orthogonal to the other two axes found previously, and given the
    %      index not yet assigned to either of those axes.
    %
    %      If only one direction matters, vecB may be empty, or the string
    %      'guess'. In this case a random vector is generated that is
    %      orthogonal to vecA, which is used as the direction indicated by
    %      the index ib.
    %
    %      Some examples may be found below which clarify the use
    %      of the '2vectors' syntax.
    %
    %
    % orientmat Methods:
    %
    %  orientmat - constructor 
    %  euler123 - constructs euler123 angles from orientation matrix
    %  draw - plots the orientation matrix as 3 axes in the global frame
    %  plus - add two orientation matrices
    %  minus - subtract two orientation matrices
    %
    %
    %
    % Examples
    %
    % The following example represents the identity matrix, i.e. no
    % rotation occurs with respect to the global reference frame. direction
    % 1 in the local frame is parallel to (1,0,0), which represents
    % direction 1 in the global frame, while direction 2 in the local frame
    % is parallel to (0,1,0), which represents direction 2 in the global
    % frame.
    %
    % om = mbdyn.pre.orientmat ( '2vectors', ...
    %                            struct ('ia', 1, ...
    %                                    'vecA', [1.,0.,0.], ...
    %                                    'ib', 2, ...
    %                                    'vecB', [0.,1.,0.] ) );
    % om.draw ()
    %
    %
    % The exact same result as the previous example can be had choosing
    % vector B to be (0.5,0.5,0.0).
    %
    % om = mbdyn.pre.orientmat ( '2vectors', ...
    %                            struct ('ia', 1, ...
    %                                    'vecA', [1.,0.,0.], ...
    %                                    'ib', 2, ...
    %                                    'vecB', [0.5,0.5,0.] ) );
    % om.draw ()
    %
    %
    % The second example describes a rotation of π/6 radian about global
    % direction 3: direction 1 in the local frame results from composing
    % cos(pi/6.) in global direction 1 and sin(pi/6.) in global direction
    % 2, while direction 3 in the local frame remains parallel to 0.,0.,1.,
    % which represents direction 3 in the global frame.
    %
    % alpha = pi/6;
    % om = mbdyn.pre.orientmat ( '2vectors', ...
    %                            struct ('ia', 1, ...
    %                                    'vecA', [cos(alpha), sin(alpha), 0], ...
    %                                    'ib', 3, ...
    %                                    'vecB', [0,0,1] ) );
    % om.draw ()
    %
    % Guess the second vector when only one direction matters
    %
    % om = mbdyn.pre.orientmat ( '2vectors', ...
    %                            struct ('ia', 1, ...
    %                                    'vecA', [1.,0.,0.], ...
    %                                    'ib', 2, ...
    %                                    'vecB', 'guess' ) );
    % om.draw ()
    %
    %
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        orientationMatrix; % full (3 x 3) orientation matrix
        
    end
    
    methods
        
        function this = orientmat (spectype, spec)
            % orientmat constructor
            %
            % Syntax
            %
            % om = orientmat ('eye')
            % om = orientmat (spectype, spec)
            %
            % Input
            %
            %  spectype - string stating how the orientation matrix is to be
            %    specified, possible options are: 'eye', 'orientation
            %    matrix', 'orientation', 'matrix', 'euler', 'euler123',
            %    'euler321', 'euler313', 'orientation vector, 'vector' and
            %    '2vectors'.
            %
            %  spec - orientation, the value of spec is dependent on the 
            %    value of spectype with the following possibilities:
            %
            %    'eye' - spec is ignores, an object representing the 
            %      (3 x 3) identity matrix is returned. Equivalent to 
            %      om = mbdyn.pre.orientmat ('orientation matrix', eye(3))
            %
            %    'orientation matrix' - spec should be a (3 x 3) full
            %      orientation matrix aka direction-cosine aka rotation
            %      matrix.
            %
            %    'orientation' - same as 'orientation matrix'
            %
            %    'matrix' - same as 'orientation matrix'
            %
            %    'euler123' - spec should be a 3 element vector
            %      representing the euler angles using the 123 convention
            %
            %    'euler' - same as for 'euler123'
            %
            %    'euler321' - spec should be a 3 element vector
            %      representing the euler angles using the 321 convention
            %
            %    'euler313' - spec should be a 3 element vector
            %      representing the euler angles using the 313 convention
            %
            %    'orientation vector' - spec should be a 3 element vector.
            %      The rotation is the magnitude of the vector and is
            %      performed around the axis formed by the vector
            %      direction.
            %
            %    'vector' - same as 'orientation vector'
            %
            %    '2vectors' - spec should be a structure containing the
            %      fields 'ia', 'vecA', 'ib' and 'vecB'. The vecA and vecB
            %      fields are three element vectors, ia and ib are integers
            %      of value 1, 2 or 3. vecA in this case is a vector
            %      representing the axis of a three-axis coordinate system
            %      with the index given in 'ia'. 'vecB' is another vector
            %      which is not parallel to vecA which is used to define a
            %      plane in 3D space. Once the plane is defined, another
            %      axis is defined as the vector orthogonal to vecA on this
            %      3D plane. This axis is assigned to be the index given in
            %      'ib' in the 3D coordinate system. The remaining axis is
            %      the vector orthogonal to the other two axes found
            %      previously, and given the index not yet assigned to
            %      either of those axes. 
            %
            %      If only one direction matters, vecB may be empty, or the
            %      string 'guess'. In this case a random vector is
            %      generated that is orthogonal to vecA, which is
            %      used as the direction indicated by the index ib.
            %
            %      Some examples may be found below which clarify the use
            %      of the '2vectors' syntax.
            %
            % Output
            %
            %  om - mbdyn.pre.orientmat object
            %
            %
            % Examples
            %
            % The following example represents the identity matrix, i.e. no
            % rotation occurs with respect to the global reference frame.
            % direction 1 in the local frame is parallel to (1,0,0), which
            % represents direction 1 in the global frame, while direction 2
            % in the local frame is parallel to (0,1,0), which represents
            % direction 2 in the global frame.
            %
            % om = mbdyn.pre.orientmat ( '2vectors', ...
            %                            struct ('ia', 1, ...
            %                                    'vecA', [1.,0.,0.], ...
            %                                    'ib', 2, ...
            %                                    'vecB', [0.,1.,0.] ) );
            % om.draw ()
            %
            %
            % The exact same result as the previous example can be had
            % choosing vector B to be (0.5,0.5,0.0).
            %
            % om = mbdyn.pre.orientmat ( '2vectors', ...
            %                            struct ('ia', 1, ...
            %                                    'vecA', [1.,0.,0.], ...
            %                                    'ib', 2, ...
            %                                    'vecB', [0.5,0.5,0.] ) );
            % om.draw ()
            %
            %
            % The second example describes a rotation of π/6 radian about
            % global direction 3: direction 1 in the local frame results
            % from composing cos(pi/6.) in global direction 1 and
            % sin(pi/6.) in global direction 2, while direction 3 in the
            % local frame remains parallel to 0.,0.,1., which represents
            % direction 3 in the global frame.
            %
            % alpha = pi/6;
            % om = mbdyn.pre.orientmat ( '2vectors', ...
            %                            struct ('ia', 1, ...
            %                                    'vecA', [cos(alpha), sin(alpha), 0], ...
            %                                    'ib', 3, ...
            %                                    'vecB', [0,0,1] ) );
            % om.draw ()
            %
            % Guess the second vector when only one direction matters
            %
            % om = mbdyn.pre.orientmat ( '2vectors', ...
            %                            struct ('ia', 1, ...
            %                                    'vecA', [1.,0.,0.], ...
            %                                    'ib', 2, ...
            %                                    'vecB', 'guess' ) );
            % om.draw ()
            %
            %
        
            switch spectype
                
                case 'eye'
                    
                    this.orientationMatrix = eye (3);
                
                case {'orientation', 'orientation matrix', 'matrix'}
                    
                    assert (~isempty (spec), 'Matrix cannot be empty.');
                        
                    mbdyn.pre.base.check3X3Matrix (spec, true, ...
                        'when using ''orientation'', ''orientation matrix'' or ''matrix'' keyword, spec' );
                    
                    this.orientationMatrix = spec;
                    
                case 'euler'
                    
                    assert (~isempty (spec), 'Input vector cannot be empty.');
                    
                    mbdyn.pre.base.check3ElementNumericVector (spec, true, 'when using ''euler'', spec');
                    
                    om = mbdyn.pre.orientmat ('euler123', spec);
                    this.orientationMatrix = om.orientationMatrix;
                    clear om;
                    
                case 'euler123'
                    
                    assert (~isempty (spec), 'Input vector cannot be empty.');
                    
                    mbdyn.pre.base.check3ElementNumericVector (spec, true, 'when using ''euler123'', spec');
                    
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
                        [ dCosBeta*dCosGamma,                                -dCosBeta*dSinGamma,                                 dSinBeta;
                          dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma, dCosAlpha*dCosGamma - dSinAlpha*dSinBeta*dSinGamma, -dSinAlpha*dCosBeta;
                          dSinAlpha*dSinGamma - dCosAlpha*dSinBeta*dCosGamma, dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma, dCosAlpha*dCosBeta; ];
                    
                case 'euler321'
                    
                    assert (~isempty (spec), 'Input vector cannot be empty.');
                    
                    mbdyn.pre.base.check3ElementNumericVector (spec, true, 'when using ''euler321'', spec');
               
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
                        [ dCosAlpha*dCosBeta,  -dSinAlpha*dCosGamma + dCosAlpha*dSinBeta*dSinGamma, dSinAlpha*dSinGamma + dCosAlpha*dSinBeta*dCosGamma
                          dSinAlpha*dCosBeta,  dCosAlpha*dCosGamma + dSinAlpha*dSinBeta*dSinGamma, -dCosAlpha*dSinGamma + dSinAlpha*dSinBeta*dCosGamma
                          -dSinBeta,           dCosBeta*dSinGamma,                                 dCosBeta*dCosGamma ];
                    
                case 'euler313'
                    
                    assert (~isempty (spec), 'Input vector cannot be empty.');
                    
                    mbdyn.pre.base.check3ElementNumericVector (spec, true, 'when using ''euler313'', spec');
                    
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
                        [ dCosAlpha*dCosGamma - dSinAlpha*dCosBeta*dSinGamma,  -dCosAlpha*dSinGamma - dSinAlpha*dCosBeta*dCosGamma, dSinAlpha*dSinBeta; 
                          dSinAlpha*dCosGamma + dCosAlpha*dCosBeta*dSinGamma,  -dSinAlpha*dSinGamma + dCosAlpha*dCosBeta*dCosGamma, -dCosAlpha*dSinBeta;
                          dSinBeta*dSinGamma,                                  dSinBeta*dCosGamma,                                  dCosBeta; ];
                      
                case {'vector', 'orientation vector'}
                    % axis and angle (angle in rad = norm of matrix)
                    
                    assert (~isempty (spec), 'Input vector cannot be empty.');
                    
                    mbdyn.pre.base.check3ElementNumericVector (spec, true, 'when using ''vector'' or ''orientation vector'', spec');
                    
                    wcrs = [ 0         spec(3) -spec(2)
                            -spec(3)        0   spec(1)
                             spec(2)  -spec(1)       0] ;   
             
                    this.orientationMatrix = expm (wcrs);
                    
                case '2vectors'
                    
                    % chech input
                    mbdyn.pre.base.checkScalarInteger (spec.ia, true, 'ia');
                    assert (~isempty (spec.vecA), 'vecA vector cannot be empty.');
                    mbdyn.pre.base.check3ElementNumericVector (spec.vecA, true, 'vecA');
                    mbdyn.pre.base.checkScalarInteger (spec.ib, true, 'ib');
                    
                    if isempty (spec.vecB) ...
                            || (ischar (spec.vecB) && strcmp (spec.vecB, 'guess'))
                        
                        spec.vecB = [0,0,0];
                        
                        i_max = 1;
                        i_min = 1;
                        
                        if abs(spec.vecA(2)) > abs(spec.vecA(1))
                            i_max = 2;
                        else
                            i_min = 2;
                        end
                        
                        if abs(spec.vecA(3)) > abs(spec.vecA(i_max))
                            
                            i_max = 3;
                            
                        elseif abs(spec.vecA(3)) < abs(spec.vecA(i_min))
                            
                            i_min = 3;
                            
                        end
                        
                        spec.vecB(i_min) = 1.0;
                        spec.vecB(i_max) = -spec.vecA(i_min) / spec.vecA(i_max);
                    
                    else
                        mbdyn.pre.base.check3ElementNumericVector (spec.vecB, true, 'vecB');
                    end
                    
                    % make sure they're column vectors
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
               
                    this.orientationMatrix = [ r{1}, r{2}, r{3} ];
                    
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
       
       function [hquiv, hax, hhiddenquiv] = draw (this, varargin)
           % plot an orientation matrix representing it as three axes
           %
           % Syntax
           %
           % [hquiv, hax] = draw (om)
           % [hquiv, hax] = draw (om, 'Parameter', value)
           %
           % Description
           %
           % plots the orientation matrix, representing it as a 3-axis
           % coordinate system rotated in the global coordinate system.
           % Axis 1 is represented as a dashed red arrow, axis 2 as blue
           % and axis 3 as green.
           %
           % Input
           %
           %  om - mbdyn.pre.orientmat object
           %
           % Additional optional arguments may be supplied as
           % parameter-value pairs. The avaialable options are:
           %
           %  'PlotAxes' - handle to matlab figure axes in shich to create
           %    the plot. If not supplied, a new figure and axes are
           %    created.
           %
           %  'Title' - logical (true/false) flag indicating whether to add
           %    a title to the plot. Title will be 'Orientation Matrix
           %    Plot'. Default is true, the title will be added.
           %
           %  'Offset' - 3 element vector. The origin of the coordinate
           %    axes representing the orientation will be offset from the
           %    global coordinate system origin by this vector in the plot.
           %    However, note that is 'DrawGlobal' is true, the global
           %    orientation coordinate axes will also be offset by this
           %    value (this is generally more useful for examining the
           %    orientation visually).
           %
           %  'DrawGlobal' - logical (true/false) flag indicating whether
           %    the global coordinate axes should also be plotted. Defautl
           %    is true.
           %
           %  'Scale' - scalar value, the length of the orientation matrix
           %    axes are 1 units by default. The length of the axes will
           %    instead be scaled by this factor in the plot.
           %
           %  'GlobalScale' - scalar value, the length of the global
           %    axes in the plot are by default half the length of the
           %    orientation matrix axes. The length of all axes will
           %    instead be scaled by this factor in the plot. If this
           %    option is supplied, the scaling is no longer related to the
           %    orientation matrix axes lengths, but is a scaling relative
           %    to a length of 1.
           %
           %  'DrawHidden' - logical (true/false) flag. When the axes are
           %    drawn, by default, 3 hidden axes are also drawn which point in
           %    the opposite direction of the drawn axes. The purpose of
           %    this is to force Matlab to choose sensible axis limits for
           %    the plot. This behaviour can be contorlled using this
           %    options. If DrawHidden is true, hidden axes is drawn, if
           %    false they are not. Default is true. Handles to the hidden
           %    quiver objects are returned in hhiddenquiv.
           %
           % Output
           %
           %  hquiv - 3 or 6 element array of handles to quiver objects. If
           %    DrawGlobal is false, it will be 3 elements, the 3 quiver
           %    plots for the orientation matrix axes plot. If DrawGlobal
           %    is true, it will be 6 elements, the additional 3 elements
           %    are the handles to the global axes quiver plots.
           %
           %  hax - handle to axes in which the plot is made
           %
           %  hhiddenquiv - 3 or 6 element array of handles to hidden quiver
           %    objects (i.e. with 'LineStyle' set to 'none'). See the
           %    'DrawHidden' option above for explanation. If DrawHidden is
           %    true, and DrawGlobal is false, it will be 3 elements, the 3
           %    hidden quiver plots for the orientation matrix axes plot.
           %    If DrawGlobal is true, it will be 6 elements, the
           %    additional 3 elements are the handles to the global axes
           %    hidden quiver plots. If DrawHidden is false, it will be
           %    empty.
           %
           %
           
           options.PlotAxes = [];
           options.Title = true;
%            options.Parent = [];
           options.Offset = [];
           options.DrawGlobal = true;
           options.Scale = 1;
           options.GlobalScale = [];
           options.DrawHidden = true;
           
           options = parse_pv_pairs (options, varargin);
           
           mbdyn.pre.base.checkLogicalScalar (options.Title, true, 'Title');
           if ~isempty (options.Offset)
               mbdyn.pre.base.check3ElementNumericVector (options.Offset, true, 'Offset');
           end
           mbdyn.pre.base.checkLogicalScalar (options.DrawGlobal, true, 'DrawGlobal');
           mbdyn.pre.base.checkNumericScalar (options.Scale, true, 'Scale');
           if isempty (options.GlobalScale)
               options.GlobalScale = options.Scale * 0.5;
           else
               mbdyn.pre.base.checkNumericScalar (options.GlobalScale, true, 'GlobalScale');
           end
           mbdyn.pre.base.checkLogicalScalar (options.DrawHidden, true, 'DrawHidden');
           
           if isempty (options.PlotAxes)
               figure;
               hax = axes;
               view (hax, 3);
               axis (hax, 'equal');
           else
               hax = options.PlotAxes;
           end
           
           x = [1;0;0];
           y = [0;1;0];
           z = [0;0;1];
           
           ox = this.orientationMatrix * x;
           oy = this.orientationMatrix * y;
           oz = this.orientationMatrix * z;
           
           ax1colour = [0.635,0.078,0.184]; % red
           ax2colour = [0,0.447,0.741]; % blue
           ax3colour = [0.466,0.674,0.188]; % green
           
           % orientation frame
           hquiv(1) = vect.plotvec3 (options.Scale .* ox, options.Offset(:), 'Properties', {'Color', ax1colour, 'LineWidth', 2, 'LineStyle', ':'}, 'PlotAxes', hax);           
           hquiv(2) = vect.plotvec3 (options.Scale .* oy, options.Offset(:), 'Properties', {'Color', ax2colour, 'LineWidth', 2, 'LineStyle', ':'}, 'PlotAxes', hax);   
           hquiv(3) = vect.plotvec3 (options.Scale .* oz, options.Offset(:), 'Properties', {'Color', ax3colour, 'LineWidth', 2, 'LineStyle', ':'}, 'PlotAxes', hax);
           
           if options.DrawHidden
               % plot reverse axes, but invisible so plot axes limits are set
               % nicely
               hhiddenquiv(1) = vect.plotvec3 (options.Scale .* -ox, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
               hhiddenquiv(2) = vect.plotvec3 (options.Scale .* -oy, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
               hhiddenquiv(3) = vect.plotvec3 (options.Scale .* -oz, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
           else
               hhiddenquiv = [];
           end
           
           if options.DrawGlobal
               % global frame
               hquiv(4) = vect.plotvec3 (options.GlobalScale .* x, options.Offset(:), 'Properties', {'Color', ax1colour}, 'PlotAxes', hax);
               hquiv(5) = vect.plotvec3 (options.GlobalScale .* y, options.Offset(:), 'Properties', {'Color', ax2colour}, 'PlotAxes', hax);
               hquiv(6) = vect.plotvec3 (options.GlobalScale .* z, options.Offset(:), 'Properties', {'Color', ax3colour}, 'PlotAxes', hax);
               
               if options.DrawHidden
                   % plot reverse axes, but invisible so plot axes limits are set
                   % nicely
                   hhiddenquiv(4) = vect.plotvec3 (options.GlobalScale .* -x, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);           
                   hhiddenquiv(5) = vect.plotvec3 (options.GlobalScale .* -y, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);           
                   hhiddenquiv(6) = vect.plotvec3 (options.GlobalScale .* -z, options.Offset(:), 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
               end
           end
           
           xlabel (hax, 'x');
           ylabel (hax, 'y');
           zlabel (hax, 'z');
                      
           if options.Title
               title ('Orientation Matrix Plot')
           end
           
       end
        
    end
    
    % operator overloading
    methods 
        
        function om = plus (om1, om2)
            % adds the orientation matrices om1 and om2
            %
            % Syntax
            %
            % om = plus (om1, om2)
            % om = om1 + om2
            %
            % Input
            %
            %  om1 - mbdyn.pre.orientmat object
            %
            %  om2 - mbdyn.pre.orientmat object
            %
            % Output
            %
            %  om - mbdyn.pre.orientmat object with (3 x 3) orientation
            %   matrix = om1.orientationMatrix + om2.orientationMatrix
            %
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix + om2.orientationMatrix);
            
        end
        
        function om = minus (om1, om2)
            % subtracts the orientation matrix om2 from om1
            %
            % Syntax
            %
            % om = minus (om1, om2)
            % om = om1 - om2
            %
            % Input
            %
            %  om1 - mbdyn.pre.orientmat object
            %
            %  om2 - mbdyn.pre.orientmat object
            %
            % Output
            %
            %  om - mbdyn.pre.orientmat object with (3 x 3) orientation
            %   matrix = om1.orientationMatrix - om2.orientationMatrix
            %
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix - om2.orientationMatrix);
            
        end
        
        function om = times (om1, om2)
            % element-by-element multiplication of orientation matrices
            %
            % Syntax
            %
            % om = times (om1, om2)
            % om = om1 .* om2
            %
            % Description
            %
            % creates a new orientmat object with orientationMatrix equal
            % to the elementwise multiplication of the input orientmat
            % orinetation matrices, i.e. the new orientmat object has
            % orientation matrix given by:
            %
            % om1.orientationMatrix .* om2.orientationMatrix
            %
            % Input
            %
            %  om1 - mbdyn.pre.orientmat object
            %
            %  om2 - mbdyn.pre.orientmat object
            %
            % Output
            %
            %  om - mbdyn.pre.orientmat object with (3 x 3) orientation
            %   matrix = om1.orientationMatrix .* om2.orientationMatrix
            %
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix .* om2.orientationMatrix);
            
        end
        
        function out = mtimes (arg1, arg2)
            % matrix multiplication of orientation matrices
            %
            % Syntax
            %
            % om = mtimes (om1, om2)
            % om = om1 * om2
            % om = mtimes (vec, om)
            % om = vec * om
            % om = mtimes (om, vec)
            % om = om * vec
            %
            % Description
            %
            % creates a new orientmat object with orientationMatrix equal
            % to the multiplication of the input orientmat orientation
            % matrices in the matrix sense, i.e. the new orientmat object
            % has orientation matrix given by:
            %
            % om1.orientationMatrix * om2.orientationMatrix
            %
            % Input
            %
            %  om1 - mbdyn.pre.orientmat object
            %
            %  om2 - mbdyn.pre.orientmat object
            %
            % Output
            %
            %  om - mbdyn.pre.orientmat object with (3 x 3) orientation
            %   matrix = om1.orientationMatrix * om2.orientationMatrix
            %
        
            if isa (arg1, 'mbdyn.pre.orientmat') && isa (arg2, 'mbdyn.pre.orientmat')
                
                out = mbdyn.pre.orientmat ('orientation', arg1.orientationMatrix * arg2.orientationMatrix);
                
            elseif isa (arg1, 'mbdyn.pre.orientmat') && isvector (arg2) && numel (arg2) == 3
                
                out = arg1.orientationMatrix * arg2;
                
            elseif isvector (arg1) && isa (arg2, 'mbdyn.pre.orientmat') && numel (arg1) == 3
                
                out = arg1 * arg2.orientationMatrix;
                
            else
                error ('matrix multiplication only defined for other orientation matrices and 3 element vectors')
            end
            
        end
        
        function om = double (om1)
            % returns underlying (3 x 3) orientation matrix

            om = om1.orientationMatrix;
            
        end
        
        function om = uminus (om1)
            
            om = mbdyn.pre.orientmat ('orientation', -om1.orientationMatrix);
            
        end
        
        function om = uplus (om1)
            
            om = mbdyn.pre.orientmat ('orientation', +om1.orientationMatrix);
            
        end
        
        function om = transpose (om1)
            % non-conjugate transpose of orientmat object
            %
            % Syntax
            %
            % om = transpose (om1)
            % om = om1.'
            %
            
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix.');
            
        end
        
        function om = ctranspose (om1)
            % complex conjugate transpose of orientmat object
            %
            % Syntax
            %
            % om = transpose (om1)
            % om = om1'
            %
        
            om = mbdyn.pre.orientmat ('orientation', om1.orientationMatrix');
            
        end
        
    end
    
    methods (Access = private)
       
        function out = unit (self, vec)
            
            out = vec ./ norm (vec);
            
        end
        
    end

    
end