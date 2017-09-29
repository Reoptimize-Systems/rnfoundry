classdef couple < mbdyn.pre.element
    
    
    properties (GetAccess = public, SetAccess = protected)
       
    end
    
    properties (GetAccess = protected, SetAccess = protected)

    end
    
    methods
        
%         function self = force ()
%             
%             
%         end

        function str = generateOutputString (self)
            str = sprintf ('    couple : %d, %s,', self.label, self.type);
        end
        
        function setSize (self, sx, sy, sz)
            self.sx = sx;
            self.sy = sy;
            self.sz = sz;
        end
        
        function setColour (self, newcolour)
            self.drawColour = newcolour;
        end
        
%         function draw (self, varargin)
%             
%             options.AxesHandle = [];
%             options.ForceRedraw = false;
%             options.Mode = 'solid';
%             
%             options = parse_pv_pairs (options, varargin);
%             
%             % try to figure out if there is a valid axees to plot to
%             if isa (options.AxesHandle, 'matlab.graphics.axis.Axes')
%                 if ~isvalid (options.AxesHandle)
%                     error ('provided axes object is not valid');
%                 end
%                 self.drawAxesH = options.AxesHandle;
%             elseif isempty (options.AxesHandle)
%                 % plot in the existing axes if possible and no new axes
%                 % supplied
%                 if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
%                     if ~isvalid (self.drawAxesH)
%                         self.drawAxesH = [];
%                     end
%                 end
%             end
%             
%             if isempty (self.shapeData)
%                 % make a unit box by default for drawing
%                 self.shapeData(1).Vertices = [ -self.sx/2, -self.sy/2, -self.sz/2;
%                                                 self.sx/2, -self.sy/2, -self.sz/2;
%                                                 self.sx/2,  self.sy/2, -self.sz/2;
%                                                -self.sx/2,  self.sy/2, -self.sz/2;
%                                                -self.sx/2, -self.sy/2,  self.sz/2;
%                                                 self.sx/2, -self.sy/2,  self.sz/2;
%                                                 self.sx/2,  self.sy/2,  self.sz/2;
%                                                -self.sx/2,  self.sy/2,  self.sz/2; ];
%                                         
%                 self.shapeData.Faces = [ 1, 4, 3, 2;
%                                          1, 5, 6, 2;
%                                          2, 6, 7, 3;
%                                          7, 8, 4, 3;
%                                          8, 5, 1, 4;
%                                          8, 7, 6, 5 ];
%                                      
%                 self.shapeDataChanged = true;
%                 
%             end
%             
%             if isempty (self.patchObject) ...
%                     || ~isvalid (self.patchObject) ...
%                     || options.ForceRedraw ...
%                     || self.shapeDataChanged
%                 
%                 % delete the current patch object
%                 if ~isempty (self.patchObject) && isvalid (self.patchObject) 
%                     delete (self.patchObject);
%                 end
%                 self.patchObject = [];
%                 
%                 % make figure and axes if necessary
%                 if isempty (self.drawAxesH)
%                     figure;
%                     self.drawAxesH = axes;
%                     if ~isempty (self.transformObject) && isvalid (self.transformObject)
%                         delete (self.transformObject);
%                     end
%                     self.transformObject = [];
%                 end
% 
%                 if isempty (self.transformObject) || ~isvalid (self.transformObject)
%                     self.transformObject = hgtransform (self.drawAxesH);
%                 end
%                 
%                 if all ( isfield (self.shapeData, {'Faces', 'Vertices'}))
% 
%                     self.patchObject = patch (self.drawAxesH, ...
%                                         'Faces', self.shapeData.Faces, ...
%                                         'Vertices', self.shapeData.Vertices, ...
%                                         'FaceLighting', 'gouraud',     ...
%                                         'AmbientStrength', 0.15, ...
%                                         'Parent', self.transformObject);
%                                 
%                 elseif all (isfield (self.shapeData, {'XData', 'YData', 'ZData'}))
%                     
%                     self.patchObject = patch (self.drawAxesH, ...
%                                         'XData', self.shapeData.XData, ...
%                                         'YData', self.shapeData.YData, ...
%                                         'ZData', self.shapeData.ZData, ...
%                                         'FaceLighting', 'gouraud',     ...
%                                         'AmbientStrength', 0.15, ...
%                                         'Parent', self.transformObject);
%                                     
%                 else
%                     error ('Invalid shape data');
%                 end
%                 
%                 self.shapeDataChanged = false;
%                 
%             end
%             
%             switch options.Mode
%                 
%                 case 'solid'
%                     set (self.patchObject, 'FaceAlpha', 1.0);
%                     set (self.patchObject, 'FaceColor', self.drawColour);
%                     set (self.patchObject, 'EdgeColor', 'none');
%                 case 'wiresolid'
%                     set (self.patchObject, 'FaceAlpha', 1.0);
%                     set (self.patchObject, 'FaceColor', self.drawColour);
%                     set (self.patchObject, 'EdgeColor', self.drawColour);
%                 case 'ghost'
%                     set (self.patchObject, 'FaceAlpha', 0.25);
%                     set (self.patchObject, 'FaceColor', self.drawColour);
%                     set (self.patchObject, 'EdgeColor', 'none');
%                 case 'wireframe'
%                     set (self.patchObject, 'EdgeColor', self.drawColour);
%                     set (self.patchObject, 'FaceColor', 'none');
%                 case 'wireghost'
%                     set (self.patchObject, 'EdgeColor', self.drawColour);
%                     set (self.patchObject, 'FaceColor', self.drawColour);
%                     set (self.patchObject, 'FaceAlpha', 0.25);
%                     
%                 otherwise
%                     set (self.patchObject, 'FaceAlpha', 1.0);
%                     set (self.patchObject, 'FaceColor', self.drawColour);
%                     set (self.patchObject, 'EdgeColor', 'none');
%                     
%             end
%             
%         end
        
    end
    
    methods (Access = protected)
        
    end
    
end