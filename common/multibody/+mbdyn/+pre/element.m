classdef element < mbdyn.pre.base
    
    
    properties (GetAccess = public, SetAccess = protected)
       
        drawColour;
        name;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       
        shapeData;
        sx;
        sy;
        sz;
        shapeDataChanged;
        stlLoaded;
        patchObject;   % patch object for the body shape
        drawAxesH; % handle to figure for plotting
        transformObject;
        defaultShape;
        
    end
    
    methods
        
        function self = element (varargin)
            
            options.STLFile = '';
            options.UseSTLName = false;
            
            options = parse_pv_pairs (options, varargin);
            
            self.stlLoaded = false;
            self.drawColour = [0.8, 0.8, 1.0];
            self.shapeData = struct([]);
            self.drawAxesH = [];
            self.sx = 1;
            self.sy = 1;
            self.sz = 1;
            
            if ~isempty (options.STLFile)
                if exist (options.STLFile, 'file')
                    self.loadSTL (options.STLFile, options.UseSTLName);
                else
                    error ('Supplied STL file %s does not appear to exist', options.STLFile);
                end
            end
            
        end
        
        function loadSTL (self, filename, usename)
            
            [self.shapeData(1).Vertices, self.shapeData(1).Faces, self.shapeData(1).N, stlname] = stl.read(filename);
            
            if usename
                self.name = stlname;
            end
            
            self.stlLoaded = true;
            self.shapeDataChanged = true;
            
        end
        
        function setSize (self, sx, sy, sz)
            self.sx = sx;
            self.sy = sy;
            self.sz = sz;
        end
        
        function setColour (self, newcolour)
            self.drawColour = newcolour;
        end
        
        function hax = draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % try to figure out if there is a valid axees to plot to
            if isa (options.AxesHandle, 'matlab.graphics.axis.Axes')
                if ~isvalid (options.AxesHandle)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = options.AxesHandle;
            elseif isempty (options.AxesHandle)
                % plot in the existing axes if possible and no new axes
                % supplied
                if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
                    if ~isvalid (self.drawAxesH)
                        self.drawAxesH = [];
                    end
                end
            end
            
            if isempty (self.shapeData)
                % make a unit box by default for drawing
                self.shapeData(1).Vertices = [ -self.sx/2, -self.sy/2, -self.sz/2;
                                                self.sx/2, -self.sy/2, -self.sz/2;
                                                self.sx/2,  self.sy/2, -self.sz/2;
                                               -self.sx/2,  self.sy/2, -self.sz/2;
                                               -self.sx/2, -self.sy/2,  self.sz/2;
                                                self.sx/2, -self.sy/2,  self.sz/2;
                                                self.sx/2,  self.sy/2,  self.sz/2;
                                               -self.sx/2,  self.sy/2,  self.sz/2; ];
                                        
                self.shapeData.Faces = [ 1, 4, 3, 2;
                                         1, 5, 6, 2;
                                         2, 6, 7, 3;
                                         7, 8, 4, 3;
                                         8, 5, 1, 4;
                                         8, 7, 6, 5 ];
                                     
                self.shapeDataChanged = true;
                
            end
            
            if isempty (self.patchObject) ...
                    || ~isvalid (self.patchObject) ...
                    || options.ForceRedraw ...
                    || self.shapeDataChanged
                
                % delete the current patch object
                if ~isempty (self.patchObject) && isvalid (self.patchObject) 
                    delete (self.patchObject);
                end
                self.patchObject = [];
                
                % make figure and axes if necessary
                if isempty (self.drawAxesH)
                    figure;
                    self.drawAxesH = axes;
                    if ~isempty (self.transformObject) && isvalid (self.transformObject)
                        delete (self.transformObject);
                    end
                    self.transformObject = [];
                end

                if isempty (self.transformObject) || ~isvalid (self.transformObject)
                    self.transformObject = hgtransform (self.drawAxesH);
                end
                
                if all ( isfield (self.shapeData, {'Faces', 'Vertices'}))

                    self.patchObject = patch (self.drawAxesH, ...
                                        'Faces', self.shapeData.Faces, ...
                                        'Vertices', self.shapeData.Vertices, ...
                                        'FaceLighting', 'gouraud',     ...
                                        'AmbientStrength', 0.15, ...
                                        'Parent', self.transformObject);
                                
                elseif all (isfield (self.shapeData, {'XData', 'YData', 'ZData'}))
                    
                    self.patchObject = patch (self.drawAxesH, ...
                                        'XData', self.shapeData.XData, ...
                                        'YData', self.shapeData.YData, ...
                                        'ZData', self.shapeData.ZData, ...
                                        'FaceLighting', 'gouraud',     ...
                                        'AmbientStrength', 0.15, ...
                                        'Parent', self.transformObject);
                                    
                else
                    error ('Invalid shape data');
                end
                
                self.shapeDataChanged = false;
                
                if options.Light
                    light (self.drawAxesH);
                end
                
                
                
            end
            
            switch options.Mode
                
                case 'solid'
                    set (self.patchObject, 'FaceAlpha', 1.0);
                    set (self.patchObject, 'FaceColor', self.drawColour);
                    set (self.patchObject, 'EdgeColor', 'none');
                case 'wiresolid'
                    set (self.patchObject, 'FaceAlpha', 1.0);
                    set (self.patchObject, 'FaceColor', self.drawColour);
                    set (self.patchObject, 'EdgeColor', self.drawColour);
                case 'ghost'
                    set (self.patchObject, 'FaceAlpha', 0.25);
                    set (self.patchObject, 'FaceColor', self.drawColour);
                    set (self.patchObject, 'EdgeColor', 'none');
                case 'wireframe'
                    set (self.patchObject, 'EdgeColor', self.drawColour);
                    set (self.patchObject, 'FaceColor', 'none');
                case 'wireghost'
                    set (self.patchObject, 'EdgeColor', self.drawColour);
                    set (self.patchObject, 'FaceColor', self.drawColour);
                    set (self.patchObject, 'FaceAlpha', 0.25);
                    
                otherwise
                    set (self.patchObject, 'FaceAlpha', 1.0);
                    set (self.patchObject, 'FaceColor', self.drawColour);
                    set (self.patchObject, 'EdgeColor', 'none');
                    
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
    
end