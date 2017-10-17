classdef reference < handle
    
    properties (GetAccess = public, SetAccess = private)
        
        name;
        positionParent;
        orientParent;
        velocityParent;
        omegaParent;
        
        dpos;
        dorientm;
        dv;
        domega;
        
    end
    
    properties (Dependent)
        
        pos;
        orientm;
        v;
        omega;
        
    end
    
    methods
        
        function this = reference (dpos, dorientm, dv, domega, varargin)
            % constructor for mbdyn.pre.reference class
            %
            % Syntax
            %
            % ref = reference (dpos, dorientm, dv, domega)
            % ref = reference (..., 'Parameter', value)
            %
            % Description
            %
            % The reference class is intended to ease the finding of
            % positions and orientations in the global coordinate system.
            % The reference class allows you to define locations and
            % orientation in the coordinate space relative to another
            % position and orientation. References can be defined relative
            % to other references, so the describe positions in frame of
            % that reference. The purpose of this is to allow the
            % modification of the orientation or position of a thing,
            % keeping orientations and positions defined relative to it
            % consistent.
            %
            % Different aspects of a referenc can also be defined relative
            % to different references. For example, the position can be
            % defined relative to reference 1, while the velocity s defined
            % relative to another reference, reference 2.
            %
            %
            % Input
            %
            %  dpos - displacement from the origin in the coordinate frame
            %    of the parent reference (see the 'Parent' and
            %    'PositionParent' option below). 
            %
            %  dorientm - orientation in the coordinate frame of the parent
            %    reference (see the 'Parent' and 'PositionParent' option
            %    below).
            %
            %  dv - velocity relative to the frame of the parent reference
            %
            %  domega - angular velocity
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            % 'Parent' - The parent reference. The position and orientation
            %   of a reference are defined in the coordinate frame of its
            %   parent reference. A reference's absolute position and
            %   orientation are found by evaluating the position and
            %   orientation of the parent reference, displacing by dpos in
            %   this coordinate frame, and then applying the rotation
            %   defined by dorientm in this coordinate frame. If not
            %   supplied, the parent is the global reference frame, given
            %   by the mbdyn.pre.globalref object. 
            %
            % 'PositionParent' - To use a different parent reference for
            %   the position. PositionParent is another mbdyn.pre.reference
            %   object.
            %
            % 'OrientParent' - To use a different parent reference for
            %   the orientation. OrientParent is another mbdyn.pre.reference
            %   object.
            %
            % 'VelParent' - To use a different parent reference for
            %   the velocity. VelParent is another mbdyn.pre.reference
            %   object.
            %
            % 'OmegaParent' - To use a different parent reference for
            %   the angular velocity. OmegaParent is another
            %   mbdyn.pre.reference object.
            %
            % 'Name' - string. Name for this reference, used in plotting.
            %
            % Output
            %
            %  ref - an mbdyn.pre.reference object
            %
            %
            % See Also: mbdyn.pre.orientmat, 
            %           mbdyn.pre.base.drawReferences
            %
            
            options.Parent = mbdyn.pre.globalref ();
            options.PositionParent = [];
            options.OrientParent = [];
            options.VelParent = [];
            options.OmegaParent = [];
            options.Name = '';
            
            options = parse_pv_pairs (options, varargin);
            
            options.Parent = processParentInput (options.Parent);
            options.PositionParent = processParentInput (options.PositionParent, options.Parent);
            options.OrientParent = processParentInput (options.OrientParent, options.Parent);
            options.VelParent = processParentInput (options.VelParent, options.Parent);
            options.OmegaParent = processParentInput (options.OmegaParent, options.Parent);
            
            this.positionParent = options.PositionParent;
            this.orientParent = options.OrientParent;
            this.velocityParent = options.VelParent;
            this.omegaParent = options.OmegaParent;
            
            if isempty (dpos) || strcmp (dpos, 'null')
                dpos = [ 0; 0; 0 ];
            end
            
            if isempty (dorientm) || strcmp (dorientm, 'null')
                this.dorientm = mbdyn.pre.orientmat ('orientation', eye (3));
            elseif ~isa (dorientm, 'mbdyn.pre.orientmat')
                error ('RENEWNET:mbdyn:badreforientation', ...
                    'dorientm should be a mbdyn.pre.orientmat  object or empty' );
            else
                this.dorientm = dorientm;
            end
            
            if isempty (dv) || strcmp (dv, 'null')
                dv = [ 0; 0; 0 ];
            end
            
            if isempty (domega) || strcmp (domega, 'null')
                domega = [ 0; 0; 0 ];
            end
            
            check.multicheck ( @(x) (isnumeric(x) && (size (x,1) == 3) && (size (x,2) == 1) ), ...
                'dpos, dv, domega must all be numeric column vectors 3 elements (or empty)', ...
                'RENEWNET:mbdyn:badrefvalues', ...
                dpos, dv, domega );
            
            if ischar (options.Name)
                this.name = options.Name;
            else
                error ('Name must be a char array.')
            end
            
            this.dpos = dpos;
            this.dv = dv;
            this.domega = domega;
            
        end
        
        function value = get.pos (this)
            value =  this.get_pos ();
        end
        
        function value = get.orientm (this)
            value = this.get_orientm ();
        end
        
        function value = get.v (this)
            value = this.get_v ();
        end
        
        function value = get.omega (this)
            value = this.get_omega ();
        end
        
        function [hax, hquiv, hhiddenquiv] = draw (this, varargin)
            % draw the reference in a figure
            %
            % Syntax
            %
            % draw (ref, 'Parameter', value)
            %
            % Description
            %
            % Draws a reference objects. The reference is represented as a
            % three-axis coordinate system.
            %
            % Input
            %
            %  ref - mbdyn.pre.reference object
            %
            % Additional optional arguments can be provided using
            % parameter-value pairs. The available options are:
            %
            %  'PlotAxes' - axes in which to create the plot. If not
            %    supplied a new figure and axes will be created.
            %
            %  'Title' - flag determining whether to add a title to the
            %    plot. Default is true.
            %
            %  'DrawGlobal' - flag determining whether the global axes will
            %    be drawn in the plot for reference. Default is true. The
            %    global axes is always drawn at the origin.
            %
            %  'Scale' - scalar value. The reference is drawn with a
            %    default size, this option may be used to adjust this by
            %    scaling the length of its axes up or down. Default is 1,
            %    no scaling.
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
            options.DrawGlobal = true;
            options.Scale = 1;
            options.DrawHidden = true;
            options.GlobalScale = [];
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkLogicalScalar (options.Title, true, 'Title');
            mbdyn.pre.base.checkLogicalScalar (options.DrawGlobal, true, 'DrawGlobal');
            mbdyn.pre.base.checkNumericScalar (options.Scale, true, 'Scale');
            if isempty (options.GlobalScale)
                options.GlobalScale = options.Scale * 0.5;
            else
                mbdyn.pre.base.checkNumericScalar (options.GlobalScale, true, 'GlobalScale');
            end
            mbdyn.pre.base.checkLogicalScalar (options.DrawHidden, true, 'DrawHidden');
            
            [hquiv, hax, hhiddenquiv] = draw ( this.orientm, ...
                                               'PlotAxes', options.PlotAxes, ...
                                               'Title', false, ...
                                               'Offset', this.pos, ...
                                               'DrawGlobal', false, ...
                                               'Scale', options.Scale, ...
                                               'DrawHidden', options.DrawHidden );
            
            x = options.Scale * [1;0;0];
            y = options.Scale * [0;1;0];
            z = options.Scale * [0;0;1];
            
            if options.DrawGlobal
                ax1colour = [0.635,0.078,0.184]; % red
                ax2colour = [0,0.447,0.741]; % blue
                ax3colour = [0.466,0.674,0.188]; % green
                
                % global frame
                hquiv(4) = vect.plotvec3 (x, [], 'Properties', {'Color', ax1colour}, 'PlotAxes', hax);
                hquiv(5) = vect.plotvec3 (y, [], 'Properties', {'Color', ax2colour}, 'PlotAxes', hax);
                hquiv(6) = vect.plotvec3 (z, [], 'Properties', {'Color', ax3colour}, 'PlotAxes', hax);
                
                if options.DrawHidden
                    % plot reverse axes, but invisible so plot axes limits are set
                    % nicely
                    hhiddenquiv(4) = vect.plotvec3 (options.GlobalScale .* -x, [], 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
                    hhiddenquiv(5) = vect.plotvec3 (options.GlobalScale .* -y, [], 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
                    hhiddenquiv(6) = vect.plotvec3 (options.GlobalScale .* -z, [], 'Properties', {'LineStyle', 'none'}, 'PlotAxes', hax);
                end
            end
            
            if options.Title
                title ('Reference Plot')
            end
            
            axis (hax, 'equal');
            
        end
        
    end
    
    methods (Access = private)
        % need this complication as you cannot access dependent properties
        % of a class property in Matlab, i.e. the x, orientm, v and omega
        % properties of the parent reference from this class
        
        function value = get_pos (this)
            % return absolute cartesian position in global frame
            %
            
            value = this.positionParent.pos ...
                + this.orientParent.orientm.orientationMatrix * this.dpos;
            
        end
        
        function value = get_orientm (this)
            % get the absolute orientation of this reference in the global
            % frame
            %
            
            value =  this.orientParent.orientm * this.dorientm;
            
        end
        
        function value = get_v (this)
            
%           rfOut.GetR()*GetVec3()
%            +rfOut.GetV()
% 			+rfOut.GetW().Cross(x-rfOut.GetX());
            
            value = this.velocityParent.orientm.orientationMatrix * this.dv ...
                + this.velocityParent.v ...
                + cross (this.velocityParent.omega, this.dpos);
        end
        
        function value = get_omega (this)
            % 
            
            value = this.orientParent.orientm.orientationMatrix * this.domega + this.omegaParent.omega;
        end
        
        function parentref = processParentInput (parentinput, defaultparent)
            
            if nargin > 1
                if isempty (parentinput)
                    parentinput = defaultparent;
                end
            end
            
            if isa (parentinput, 'mbdyn.pre.structuralNode3dof')
                
                % replace node with reference
                
                parentref = mbdyn.pre.reference ( ...
                    parentinput.absolutePosition, ...
                    parentinput.absoluteOrientation, ...
                    parentinput.absoluteVelocity, ...
                    parentinput.absoluteAngularVelocity );
                
            elseif isa (parentinput, 'mbdyn.pre.structuralNode6dof')
                
                % create reference from node absolute positions etc
                parentref = mbdyn.pre.reference ( ...
                    parentinput.absolutePosition, ...
                    parentinput.absoluteOrientation, ...
                    parentinput.absoluteVelocity, ...
                    parentinput.absoluteAngularVelocity );
                
            elseif ~ ( isa (parentinput, 'mbdyn.pre.reference') ...
                    || isa (parentinput, 'mbdyn.pre.globalref') )
                
                error ('RENEWNET:mbdyn:badrefparent', ...
                    'The parent must be another reference object, the global object or a node');
                
            end
            
        end
        
    end
    
end