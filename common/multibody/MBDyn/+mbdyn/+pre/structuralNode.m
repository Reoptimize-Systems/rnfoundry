classdef structuralNode < mbdyn.pre.node
    % generic structural node class, ancestor of all structural type nodes 
    
    
    properties (GetAccess = public, SetAccess = public)
        
        % we can set this later
        absolutePosition;
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        
        absoluteVelocity;
        
        % assembly
        positionInitialStiffness;
        velocityInitialStiffness;
        
        accelerations;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       
        sx;
        sy;
        sz;
        drawObjects; % line objects for the node drawing
        drawAxesH; % handle to figure for plotting
        transformObject;
        drawColour;
        
    end
    
    methods
        
        function self = structuralNode (varargin)
            
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteVelocity = [0;0;0];
            options.Accelerations = [];
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.node ( ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output );
            
            if ~isempty (options.Accelerations)
                if islogical (options.Accelerations) && isscalar (options.Accelerations)
                    self.accelerations = options.Accelerations;
                else
                    error ('Accelerations should be a scalar boolean true/false');
                end
            end
%             self.checkCartesianVector (options.AbsolutePosition, true);
            self.checkCartesianVector (options.AbsoluteVelocity, true);
            self.absolutePosition = options.AbsolutePosition;
            self.absoluteVelocity = options.AbsoluteVelocity;
            
            % TODO: find out what values of initial stiffness can be
            % supplied
            self.positionInitialStiffness = [];
            self.velocityInitialStiffness = [];
            
            self.drawColour = [0.8, 0.1, 0.1];
            self.drawAxesH = [];
            self.sx = 1;
            self.sy = 1;
            self.sz = 1;
            
        end
        
        function str = generateOutputString (self)
            str = generateOutputString@mbdyn.pre.node (self);
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            
            options = parse_pv_pairs (options, varargin);
            
            if options.ForceRedraw
                % delete the current drawing object
                if ~isempty (self.drawObjects)
                    for ind = 1:numel (self.drawObjects)
                        if isvalid (self.drawObjects(ind)) 
                            delete (self.drawObjects(ind));
                        end
                    end
                end
                self.drawObjects = [];
            end
            
            % try to figure out if there is a valid set of axes to plot to
            if isa (options.AxesHandle, 'matlab.graphics.axis.Axes')
                
                if ~isvalid (options.AxesHandle)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = options.AxesHandle;
                % we need to redraw, as we're plotting in a new set of axes
                options.ForceRedraw = true;
                % abandon existing patch objects 
                self.drawObjects = [];
                % abandon existing transform object
                self.transformObject = [];
                
            elseif isempty (options.AxesHandle)
                
                % plot in the existing axes if possible as no new axes
                % supplied
                if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
                    if ~isvalid (self.drawAxesH)
                        self.drawAxesH = [];
                        options.ForceRedraw = true;
                        self.transformObject = [];
                    end
                end
                
            end
            
            if isempty (self.drawObjects) ... % we need to create new plot as there's nothing in exsiting one
                    || ~all(isvalid (self.drawObjects)) ... % we need to create new plot as existing pathces have been destroyed
                    || options.ForceRedraw ... % new plot is forced by caller
                    || isempty (self.drawAxesH) % if self.drawAxesH is empty we need to create new axes and new plot
                
                % delete the current patch object
                if ~isempty (self.drawObjects)
                    for ind = 1:numel (self.drawObjects)
                        if isvalid (self.drawObjects(ind)) 
                            delete (self.drawObjects(ind));
                        end
                    end
                end
                self.drawObjects = [];
                
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
                
                self.drawObjects = line ([-self.sx/2,     0,           0; 
                                           self.sx/2,     0,           0; ], ...
                                         [  0 ,       -self.sy/2,      0;
                                            0,         self.sy/2,      0  ], ...
                                         [  0 ,           0,       -self.sz/2;
                                            0,            0,        self.sz/2 ], ...
                                          'Parent', self.transformObject, ...
                                          'Color', self.drawColour );
            end
            
            self.setTransform ();
            
        end
        
        function setSize (self, sx, sy, sz)
            self.sx = sx;
            self.sy = sy;
            self.sz = sz;
        end
        
        function setColour (self, newcolour)
            self.drawColour = newcolour;
        end
        
    end
    
    % getters/setters
    methods
        function set.absolutePosition (self, newpos)
            % set the absolute position of the structural node
            
            self.checkCartesianVector (newpos, true);
            
            self.absolutePosition = newpos;
            
        end
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            
        end
        
    end
    
end

