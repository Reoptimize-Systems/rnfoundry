classdef element < mbdyn.pre.base
% base class for all MBDyn elements
%
% Syntax
%
% el = mbdyn.pre.element ('Parameter', Value)
%
% Description
%
% mbdyn.pre.element is the base class for all other element types in the
% toolbox. It contains methods and properties common to all elements such
% as bodies, joints and forces etc. It is not intended to be used directly
% by ordinary users.
%
% mbdyn.pre.element Methods:
%
%   element - mbdyn.pre.element constructor
%   draw - plot the element in a figure
%   loadSTL - load an STL file which will be used when visualising the element
%   setColour - set the colour of the element in plots
%   setSize - set the size of the element in plots
%
%
    
    properties (GetAccess = public, SetAccess = protected)
       
        name; % name of the element
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       
        stlLoaded;
        defaultShape;
        defaultShapeOrientation;
        stlNormals;
        
    end
    
    methods
        
        function self = element (varargin)
            % mbdyn.pre.element constructor
            %
            % Syntax
            %
            % el = mbdyn.pre.element ('Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.element is the base class for all other element
            % types in the toolbox. It contains methods and properties
            % common to all elements such as bodies, joints and forces etc.
            % It is not intended to be used directly by ordinary users.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
            %
            %  'STLFile' - character vector containing the full path to an
            %    STL file which will be used when visualising/plotting the
            %    element.
            %
            %  'UseSTLName' - true/false flag indicating whether to use the
            %    name embedded in the STL file as the name for this
            %    element, in which case it will be used in the 'name'
            %    property of the element.
            %
            % Output
            %
            %  el - mbdyn.pre.element
            %
            
            [options, ~] = mbdyn.pre.element.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkAllowedStringInputs ( options.DefaultShape, ...
                                            { 'none', 'cuboid', 'box', 'cylinder', 'sphere', 'ellipsoid', 'tube', 'pipe', 'annularcylinder' }, ...
                                            true, ...
                                            'DefaultShape' );
                                        
            self.checkOrientationMatrix ( options.DefaultShapeOrientation, ...
                                          true, ...
                                          'DefaultShapeOrientation' );
            
            self.stlLoaded = false;
            self.drawColour = [0.8, 0.8, 1.0];
            self.shapeData = struct([]);
            self.shapeObjects = {};
            self.drawAxesH = [];
            self.defaultShape = options.DefaultShape;
            self.defaultShapeOrientation = options.DefaultShapeOrientation;
            
            switch self.defaultShape
                
                case 'none'
                    % do nothing
                    
                case {'cuboid', 'box'}
                    
                    % cuboid, 3 arguments expected, x, y and z dimensions
                    self.setSize (1, 1, 1);

                case 'cylinder'
                    
                    % cylinder, two arguments expected, radius and axial
                    % length
                    self.setSize (1, 2);

                case {'tube', 'pipe', 'annularcylinder'}

                    % tube, 3 arguments expected, router, rinner and
                    % axiallength dimensions
                    self.setSize (1, 0.5, 2);

                case 'sphere'


                case 'ellipsoid'


                otherwise
                    error ('Bad defaultShape string');
                        
            end
            
            if ~isempty (options.STLFile)
                if exist (options.STLFile, 'file')
                    self.loadSTL (options.STLFile, options.UseSTLName);
                else
                    error ('Supplied STL file %s does not appear to exist', options.STLFile);
                end
            end
            
        end
        
        function loadSTL (self, filename, usename)
            % load an STL file which will be used when visualising the element
            %
            % Syntax
            %
            % mbdyn.pre.element (self, filename, usename)
            %
            % Description
            %
            % load an STL file which will be used when visualising the
            % element instead of the default shapes.
            %
            % Input
            %
            %  el - mbdyn.pre.element object
            %
            %  filename - full path to the STL file to be loaded
            %
            %  usename - true/false 
            %
            %
            % See Also: mbdyn.pre.element.node
            %

            self.checkLogicalScalar (usename, true, 'usename');
            
            self.shapeData{1} = struct ();
            
            [self.shapeData{1}.Vertices, self.shapeData{1}.Faces, self.stlNormals, stlname] = stl.read(filename);
            
            if usename
                self.name = stlname;
            end
            
            self.stlLoaded = true;
            self.needsRedraw = true;
            
            setSize (self, ...
                max(self.shapeData{1}.Vertices(:,1)) - min(self.shapeData{1}.Vertices(:,1)), ...
                max(self.shapeData{1}.Vertices(:,2)) - min(self.shapeData{1}.Vertices(:,2)), ...
                max(self.shapeData{1}.Vertices(:,3)) - min(self.shapeData{1}.Vertices(:,3)) )
            
        end
        
        function setSize (self, varargin)
            % set the size of the element in plots
            %
            % Syntax
            %
            % setSize (el, sx, sy, sz)
            % setSize (el, radius, axiallength)
            % setSize (el, router, rinner, axiallength)
            %
            % Description
            %
            % setSize is used to set the size of the default element shape
            % for plotting the element in a figure. This is used when no
            % STL file is avaialable, of the subclassed elemnt does not
            % provide it's own drawing of the element. The inputs to
            % setSize depend on what the element's chosen shape is.
            %
            % Input
            %
            %  el - mbdyn.pre.element object
            %
            %  sx - used when the shape is a box/cuboid, this is the
            %   length along the x axis
            %
            %  sy - used when the shape is a box/cuboid, this is the
            %   length along the y axis
            %
            %  sz - used when the shape is a box/cuboid, this is the
            %   length along the z axis
            %
            %  radius - used when the shape is a cylinder, this is the
            %   radius of the cylinder
            %
            %  axiallength - used when the shape is a cylinder, this is the
            %   axial length of the cylinder
            %
            %  router - used when the shape is a tube/pipe/annularcylinder,
            %   this is the outer radius of the tube.
            %
            %  rinne - used when the shape is a tube/pipe/annularcylinder,
            %   this is the inner radius of the tube.
            %
            %  axiallength - used when the shape is a tube/pipe/annularcylinder,
            %   this is the axial length of the tube.
            %
            % Output
            %
            %
            %
            % See Also: 
            %

            if self.stlLoaded
                
                % cuboid, 3 arguments expected, x, y and z dimensions
                assert (numel (varargin) == 3, ...
                        'setSize requires 3 size input arguments when the shape is from an STL file, sx, sy and sz, which represent the bounding box of the shape');

                self.checkNumericScalar (varargin{1}, true, 'sx');
                self.checkNumericScalar (varargin{2}, true, 'sy');
                self.checkNumericScalar (varargin{3}, true, 'sz');

                assert (varargin{1} > 0, 'sx must be greater than zero');
                assert (varargin{2} > 0, 'sy must be greater than zero');
                assert (varargin{3} > 0, 'sz must be greater than zero');

                self.shapeParameters(1) = varargin{1};
                self.shapeParameters(2) = varargin{2};
                self.shapeParameters(3) = varargin{3};
                        
            else
                
                switch self.defaultShape

                    case 'none'

                        warning ('Default shape is set to ''none'', setting the size has no effect');

                    case {'cuboid', 'box'}

                        % cuboid, 3 arguments expected, x, y and z dimensions
                        assert (numel (varargin) == 3, ...
                                'setSize requires 3 size input arguments when the shape is a box/cuboid, sx, sy and sz');

                        self.checkNumericScalar (varargin{1}, true, 'sx');
                        self.checkNumericScalar (varargin{2}, true, 'sy');
                        self.checkNumericScalar (varargin{3}, true, 'sz');

                        assert (varargin{1} > 0, 'sx must be greater than zero');
                        assert (varargin{2} > 0, 'sy must be greater than zero');
                        assert (varargin{3} > 0, 'sz must be greater than zero');

                        self.shapeParameters(1) = varargin{1};
                        self.shapeParameters(2) = varargin{2};
                        self.shapeParameters(3) = varargin{3};

                    case 'cylinder'

                        % cylinder, two arguments expected, radius and axial
                        % length
                        assert (numel (varargin) == 2, ...
                                'setSize requires 2 size input arguments when the shape is a cylinder, radius, axiallength');

                        self.checkNumericScalar (varargin{1}, true, 'radius');
                        self.checkNumericScalar (varargin{2}, true, 'axiallength');

                        assert (varargin{1} > 0, 'radius must be greater than zero');
                        assert (varargin{2} > 0, 'axiallength must be greater than zero');

                        self.shapeParameters(1) = varargin{1};
                        self.shapeParameters(2) = varargin{2};


                    case {'tube', 'pipe', 'annularcylinder'}

                        % tube, 3 arguments expected, router, rinner and
                        % axiallength dimensions
                        assert (numel (varargin) == 3, ...
                                'setSize requires 3 size input arguments when the shape is a tube/pipe/annularcylinder, router, rinner and axiallength');

                        self.checkNumericScalar (varargin{1}, true, 'router');
                        self.checkNumericScalar (varargin{2}, true, 'rinner');
                        self.checkNumericScalar (varargin{3}, true, 'axiallength');

                        assert (varargin{1} > 0, 'router must be greater than zero');
                        assert (varargin{2} > 0, 'rinner must be greater than zero');
                        assert (varargin{1} > varargin{2}, 'router must be greater than rinner');
                        assert (varargin{3} > 0, 'axiallength must be greater than zero');

                        self.shapeParameters(1) = varargin{1};
                        self.shapeParameters(2) = varargin{2};
                        self.shapeParameters(3) = varargin{3};

                    case 'sphere'


                    case 'ellipsoid'


                    otherwise
                        error ('Bad defaultShape string');

                end

                % set the shapedata to empty so it is recreated with the new
                % sizes when draw is next called
                self.shapeData = [];
            
            end
            
        end
        
        function setColour (self, newcolour)
            % set the colour of the element in plots
            
            self.drawColour = newcolour;
        end
        
        function hax = draw (self, varargin)
            % plot the element in a figure
            %
            % Syntax
            %
            % hax = draw (el)
            % hax = draw (..., 'Parameter', Value)
            %
            % Description
            %
            % The draw method creates a visualisation of the element in a
            % figure. If an STL file has previously be added to the
            % element, this will be plotted. Otherwise a standard shape
            % will be used.
            %
            % The draw method tries to be efficient. Each element is
            % associated with a hgtransform object. If the shape of the
            % object has not changed, the transform matrix is simply
            % updated to adjust the location and orientation on the plot. 
            %
            % Input
            %
            %  el - mbdyn.pre.element object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'AxesHandle' - optional handle to axes in which to plot the
            %    element. If not supplied, a new figure and axes will be
            %    created. The first time the element is drawn the handle to
            %    the axes will be stored internally and future calls to
            %    draw will plot to the same axes, unless the axes are
            %    destroyed, or this option is used to override it.
            %
            %  'ForceRedraw' - true/false flag indicating whether to force
            %    a full redraw of the object (rather than just update the
            %    transform matrix), even if the element does not think it
            %    needs it.
            %
            %  'Mode' - character vector determining the style in which the
            %    element will be plotted. Can be one of 'solid',
            %    'wiresolid', 'ghost', 'wireframe', 'wireghost'. Default is
            %    'solid'.
            %
            %  'Light' - deterined whether the scene should have light
            %    source
            %
            % Output
            %
            %
            %
            % See Also: 
            %

            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar ( options.ForceRedraw, true, 'ForceRedraw' );
            self.checkAllowedStringInputs ( options.Mode, ...
                                            { 'solid', ...
                                              'wiresolid', ...
                                              'ghost', ...
                                              'wireframe', ...
                                              'wireghost' }, ...
                                            true, ...
                                            'Mode' );
            self.checkLogicalScalar ( options.Light, true, 'Light' );
            
            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            if isempty (self.shapeData)
                
                switch self.defaultShape
                    
                    case 'none'
                        
                        self.shapeData = {};
                    
                    case {'cuboid', 'box'}
                        
                        self.shapeData{1} = struct ();
                        
                        % make a unit box by default for drawing
                        self.shapeData{1}.Vertices = [ -self.shapeParameters(1)/2, -self.shapeParameters(2)/2, -self.shapeParameters(3)/2;
                                                        self.shapeParameters(1)/2, -self.shapeParameters(2)/2, -self.shapeParameters(3)/2;
                                                        self.shapeParameters(1)/2,  self.shapeParameters(2)/2, -self.shapeParameters(3)/2;
                                                       -self.shapeParameters(1)/2,  self.shapeParameters(2)/2, -self.shapeParameters(3)/2;
                                                       -self.shapeParameters(1)/2, -self.shapeParameters(2)/2,  self.shapeParameters(3)/2;
                                                        self.shapeParameters(1)/2, -self.shapeParameters(2)/2,  self.shapeParameters(3)/2;
                                                        self.shapeParameters(1)/2,  self.shapeParameters(2)/2,  self.shapeParameters(3)/2;
                                                       -self.shapeParameters(1)/2,  self.shapeParameters(2)/2,  self.shapeParameters(3)/2; ];

                        self.shapeData{1}.Faces = [ 1, 4, 3, 2;
                                                    1, 5, 6, 2;
                                                    2, 6, 7, 3;
                                                    7, 8, 4, 3;
                                                    8, 5, 1, 4;
                                                    8, 7, 6, 5 ];
                                     
                    case 'cylinder'
                        
                        npnts = 30;
                        [X,Y,Z] = cylinder (self.shapeParameters(1), npnts-1);
                        Z = Z .* self.shapeParameters(2);
                        Z = Z - self.shapeParameters(3)/2;
                        
                        % rotate
                        XYZtemp = [ X(1,:);
                                    Y(1,:)
                                    Z(1,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        X(1,:) = XYZtemp(1,:);
                        Y(1,:) = XYZtemp(2,:);
                        Z(1,:) = XYZtemp(3,:);
                        
                        XYZtemp = [ Xo(2,:);
                                    Yo(2,:)
                                    Zo(2,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        X(2,:) = XYZtemp(1,:);
                        Y(2,:) = XYZtemp(2,:);
                        Z(2,:) = XYZtemp(3,:);
                        
                        
                        self.shapeData{1} = struct ();
                        self.shapeData{1}.Vertices = [];
                        self.shapeData{1}.Faces = [];
                        
                        self.shapeData{1}.Vertices = [ X(1,:)', Y(1,:)', Z(1,:)';
                                                       X(2,:)', Y(2,:)', Z(2,:)'; ];
                                                   
                        self.shapeData{1}.Faces = [ (1:npnts)', (1:npnts)' + 1, (1:npnts)' + 1 + npnts, (1:npnts)' + npnts ];
                        self.shapeData{1}.Faces (end, 2) = 1;
                        self.shapeData{1}.Faces (end, 3) = 1 + npnts;
                        
                        self.shapeData{2} = struct ();
                        self.shapeData{2}.Vertices = [ X(1,:)', Y(1,:)', Z(1,:)' ];
                        self.shapeData{2}.Faces = 1:npnts;
                        
                        self.shapeData{3} = struct ();
                        self.shapeData{3}.Vertices = [ X(2,:)', Y(2,:)', Z(2,:)' ];
                        self.shapeData{3}.Faces = 1:npnts;
                        
                    case {'tube', 'pipe', 'annularcylinder'}
                        
                        npnts = 20;
                        [Xo,Yo,Zo] = cylinder (self.shapeParameters(1), npnts-1);
                        Zo = Zo .* self.shapeParameters(3);
                        
                        [Xi,Yi,Zi] = cylinder (self.shapeParameters(2), npnts-1);
                        Zi = Zi .* self.shapeParameters(3);
                        
                        Zo = Zo - self.shapeParameters(3)/2;
                        Zi = Zi - self.shapeParameters(3)/2;
                        
                        % rotate
                        XYZtemp = [ Xo(1,:);
                                    Yo(1,:)
                                    Zo(1,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        Xo(1,:) = XYZtemp(1,:);
                        Yo(1,:) = XYZtemp(2,:);
                        Zo(1,:) = XYZtemp(3,:);
                        
                        XYZtemp = [ Xo(2,:);
                                    Yo(2,:)
                                    Zo(2,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        Xo(2,:) = XYZtemp(1,:);
                        Yo(2,:) = XYZtemp(2,:);
                        Zo(2,:) = XYZtemp(3,:);
                        
                        XYZtemp = [ Xi(1,:);
                                    Yi(1,:)
                                    Zi(1,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        Xi(1,:) = XYZtemp(1,:);
                        Yi(1,:) = XYZtemp(2,:);
                        Zi(1,:) = XYZtemp(3,:);
                        
                        XYZtemp = [ Xi(2,:);
                                    Yi(2,:)
                                    Zi(2,:) ];
                                
                        XYZtemp = self.defaultShapeOrientation.orientationMatrix * XYZtemp;
                        
                        Xi(2,:) = XYZtemp(1,:);
                        Yi(2,:) = XYZtemp(2,:);
                        Zi(2,:) = XYZtemp(3,:);
                        
                        
                        self.shapeData{1} = struct ();
                        self.shapeData{1}.Vertices = [];
                        self.shapeData{1}.Faces = [];
                        
                        self.shapeData{1}.Vertices = [ Xo(1,:)', Yo(1,:)', Zo(1,:)';
                                                       Xo(2,:)', Yo(2,:)', Zo(2,:)'; ];
                                                   
                        self.shapeData{1}.Faces = [ (1:npnts)', (1:npnts)' + 1, (1:npnts)' + 1 + npnts, (1:npnts)' + npnts ];
                        self.shapeData{1}.Faces (end, 2) = 1;
                        self.shapeData{1}.Faces (end, 3) = 1 + npnts;
                        
                        self.shapeData{2}.Vertices = [ Xi(1,:)', Yi(1,:)', Zi(1,:)';
                                                       Xi(2,:)', Yi(2,:)', Zi(2,:)'; ];
                                                   
                        self.shapeData{2}.Faces = [ (1:npnts)', (1:npnts)' + 1, (1:npnts)' + 1 + npnts, (1:npnts)' + npnts ];
                        self.shapeData{2}.Faces (end, 2) = 1;
                        self.shapeData{2}.Faces (end, 3) = 1 + npnts;
                        
                        self.shapeData{3} = struct ();
                        self.shapeData{3}.Vertices = [ Xo(1,:)', Yo(1,:)', Zo(1,:)';
                                                       Xi(1,:)', Yi(1,:)', Zi(1,:)'; ];
                        self.shapeData{3}.Faces = [ 1:npnts, 1, (1:npnts) + npnts, 1 + npnts ];
                        
                        self.shapeData{4} = struct ();
                        self.shapeData{4}.Vertices = [ Xo(2,:)', Yo(2,:)', Zo(2,:)';
                                                       Xi(2,:)', Yi(2,:)', Zi(2,:)'; ];
                        self.shapeData{4}.Faces = [ 1:npnts, 1, (1:npnts) + npnts, 1 + npnts ];
                        
                        
                    case 'sphere'
                        
                        
                    case 'ellipsoid'
                        
                        
                    otherwise
                        error ('Bad defaultShape string');
                        
                end
                                     
                self.needsRedraw = true;
                
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                self.shapeObjects = {};
                
                for ind = 1:numel (self.shapeData)
                    if all ( isfield (self.shapeData{ind}, {'Faces', 'Vertices'})) ...
                        || all (isfield (self.shapeData{ind}, {'XData', 'YData', 'ZData'}))

                        self.shapeData{ind}.FaceLighting = 'gouraud';
                        self.shapeData{ind}.AmbientStrength = 0.15;
                        self.shapeData{ind}.Parent = self.transformObject;
                        
                        self.shapeObjects = [ self.shapeObjects, ...
                                              { patch( self.drawAxesH, ...
                                                       self.shapeData{ind} ) } ...
                                             ];
                                             
%                        self.shapeObjects = [ self.shapeObjects, ...
%                                               { patch( self.drawAxesH, ...
%                                                        self.shapeData{ind}, ...
%                                                        'FaceLighting', 'gouraud', ...
%                                                        'AmbientStrength', 0.15, ...
%                                                        'Parent', self.transformObject ) ...
%                                               } ...
%                                              ];

                    else
                        error ('Invalid shape data');
                    end
                end
                
                self.needsRedraw = false;
                
                if options.Light
                    light (self.drawAxesH);
                end
                
            end
            
            for ind = 1:numel (self.shapeObjects)
                
                switch options.Mode

                    case 'solid'
                        set (self.shapeObjects{ind}, 'FaceAlpha', 1.0);
                        set (self.shapeObjects{ind}, 'FaceColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'EdgeColor', 'none');
                    case 'wiresolid'
                        set (self.shapeObjects{ind}, 'FaceAlpha', 1.0);
                        set (self.shapeObjects{ind}, 'FaceColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'EdgeColor', self.drawColour);
                    case 'ghost'
                        set (self.shapeObjects{ind}, 'FaceAlpha', 0.25);
                        set (self.shapeObjects{ind}, 'FaceColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'EdgeColor', 'none');
                    case 'wireframe'
                        set (self.shapeObjects{ind}, 'EdgeColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'FaceColor', 'none');
                    case 'wireghost'
                        set (self.shapeObjects{ind}, 'EdgeColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'FaceColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'FaceAlpha', 0.25);

                    otherwise
                        set (self.shapeObjects{ind}, 'FaceAlpha', 1.0);
                        set (self.shapeObjects{ind}, 'FaceColor', self.drawColour);
                        set (self.shapeObjects{ind}, 'EdgeColor', 'none');

                end
            
            end
            
            if nargout > 0
                hax = self.drawAxesH;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function ok = checkInertiaMatrix (self, mat, throw)
            
            ok = self.check3X3Matrix (mat, false);
            
            if ~ok && throw
                error ('Inertia matrix must be a 3 x 3 numeric matrix');
            end
            
        end
        
        function ok = checkCOGVector (self, cog, throw)
            
            if isempty (cog) || (ischar (cog) && strcmp (cog, 'null'))
                ok = true;
            else
                ok = self.checkCartesianVector (cog, false);
            end
            
            if ~ok && throw
                error ('Centre of gravity offset must 3 element numeric column vector, or keyword ''null'' or empty');
            end
            
            
        end

        
    end
    
    methods (Static)
                
        function comment = nodeLabelComment (node)
            
            if isempty (node.name)
                comment = sprintf ('node label');
            else 
                comment = sprintf ('node label: %s', node.name);
            end
            
        end
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options.STLFile = '';
            options.UseSTLName = false;
            options.DefaultShape = 'cuboid';
            options.DefaultShapeOrientation = mbdyn.pre.orientmat ('eye');
            
            nopass_list = {};
            
        end

    end
    
end