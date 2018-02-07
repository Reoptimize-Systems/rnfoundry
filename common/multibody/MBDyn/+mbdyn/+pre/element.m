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
        stlLoaded;
        shapeObjects;   % object for the element shape drawing
        needsRedraw; % track whether abject needs to be redrawn
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
            self.shapeObjects = {};
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
            self.needsRedraw = true;
            
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
            
            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
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
                                     
                self.needsRedraw = true;
                
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                
                if all ( isfield (self.shapeData, {'Faces', 'Vertices'}))

                    self.shapeObjects = { patch( self.drawAxesH, ...
                                                 'Faces', self.shapeData.Faces, ...
                                                 'Vertices', self.shapeData.Vertices, ...
                                                 'FaceLighting', 'gouraud', ...
                                                 'AmbientStrength', 0.15, ...
                                                 'Parent', self.transformObject ) ...
                                        };
                                
                elseif all (isfield (self.shapeData, {'XData', 'YData', 'ZData'}))
                    
                    self.shapeObjects = { patch( self.drawAxesH, ...
                                                 'XData', self.shapeData.XData, ...
                                                 'YData', self.shapeData.YData, ...
                                                 'ZData', self.shapeData.ZData, ...
                                                 'FaceLighting', 'gouraud',     ...
                                                 'AmbientStrength', 0.15, ...
                                                 'Parent', self.transformObject) ...
                                        };
                                    
                else
                    error ('Invalid shape data');
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
        
        function checkAxes (self, hax)
            % checks if there is a valid set of axes to plot to, and if
            % not, creates them
            
            % try to figure out if there is a valid set of axes to plot to
            if isempty (hax)
                
                % plot in the existing axes if possible as no new axes
                % supplied
                if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
                    % use existing axes
                    
                    if ~isvalid (self.drawAxesH)
                        % axes no longer exists, or figure has been closed
                        self.drawAxesH = [];
                        self.deleteAllDrawnObjects ();
                        self.transformObject = [];
                        self.needsRedraw = true;
                        % make a new set of axes to draw to
                        self.makeAxes ();
                    end
                    
                elseif isempty (self.drawAxesH)
                    % make figure and axes
                    self.makeAxes ();
                    self.needsRedraw = true;
                    
                else
                    error ('drawAxesH property is not empty or an axes handle');
                end
            
            elseif isa (hax, 'matlab.graphics.axis.Axes')
                % an axes has been supplied, so we plot to this new axes
                
                if ~isvalid (hax)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = hax;
                % abandon existing shape objects and transform object
                % TODO: should we delet objects in old axes?
                self.shapeObjects = {};
                self.transformObject = [];
                % we need to redraw, as we're plotting in a different set
                % of axes
                self.needsRedraw = true;
                
            else
                
                error ('hax must be an axes handle or empty');
                
            end
            
            if isempty (self.transformObject) || ~isvalid (self.transformObject)
                self.transformObject = hgtransform (self.drawAxesH);
            end
            
        end
        
        function makeAxes (self)
            % create axes and transform object
            
            figure;
            self.drawAxesH = axes;
            if ~isempty (self.transformObject) && isvalid (self.transformObject)
                delete (self.transformObject);
            end
            self.transformObject = [];
            self.needsRedraw = true;
            
        end
        
        function deleteAllDrawnObjects (self)
            
            for ind = 1:numel (self.shapeObjects)
                if ~isempty (self.shapeObjects{ind}) ...
                        && isvalid (self.shapeObjects{ind})

                    delete (self.shapeObjects{ind});

                end
            end
            self.shapeObjects = {};
            
        end
        
    end
    
end