classdef structuralNode < mbdyn.pre.node
% generic structural node class, ancestor of all structural type nodes
%
% Syntax
%
% [str] = mbdyn.pre.structuralNode ()
% [str] = mbdyn.pre.structuralNode (..., 'Parameter', Value)
%
% Description
%
% Base class for the structural node types. Not generally
% intended to be used directly by normal users of the toolbox.
% Instead, the mbdyn.pre.structuralNode6dof and
% mbdyn.pre.structuralNode3dof classes should be used.  
%
% mbdyn.pre.structuralNode Methods:
%
%  structuralNode - constructor
%  generateOutputString - make MBDyn output string for node
%  draw - draw the node in a figure
%  setSize - set the size of the node when plotted in a figure
%  setColour - set the colour of the node when plotted in a figure
%
%
% See Also: mbdyn.pre.structuralNode6dof, 
%           mbdyn.pre.structuralNode3dof
%

    properties (GetAccess = public, SetAccess = public)
        
        absolutePosition; % absolute position of the node in the global frame
        absoluteVelocity; % absolute velocity of the node in the global frame
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        
        positionInitialStiffness;
        velocityInitialStiffness;
        
        accelerations; % flag indicating whether to output acceleration data
        
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
            % mbdyn.pre.structuralNode constructor
            %
            % Syntax
            %
            % [str] = mbdyn.pre.structuralNode ()
            % [str] = mbdyn.pre.structuralNode (..., 'Parameter', Value)
            %
            % Description
            %
            % Base class for the structural node types. Not generally
            % intended to be used directly by normal users of the toolbox.
            % Instead, the mbdyn.pre.structuralNode6dof and
            % mbdyn.pre.structuralNode3dof classes should be used.
            %
            % Input
            %
            % Addtional arguments may be supplied as parameter-value pairs. The available options are:
            %
            %  'AbsolutePosition' - optional (3 x 1) vector containing the
            %    intial position of the node in the global frame. Default
            %    is [0;0;0] if not supplied.
            %
            %  'AbsoluteVelocity' - optional (3 x 1) vector containing the
            %    intial velocity of the node in the global frame. Default
            %    is [0;0;0] if not supplied.
            %
            %  'Accelerations' - true/false flag, or a character vector
            %    which must be 'yes' of 'no'. Determines whether this node
            %    will output acceleration data.
            %
            % 'HumanReadableLabel' - a text string intended to provide a
            %   meaningful label. For some node types this may optionally
            %   be displayed when they they are drawn.
            %
            %  'Scale' - optional. Used to control the scaling of the
            %    residual for the node equations before performing testing
            %    whether the required tolerance has been met. For more
            %    information see the help for mbdyn.pre.initialValueProblem
            %    (see the 'ScaleResidual' option in the constructor), and
            %    mbdyn.pre.system (see the 'DefaultScales' option in the
            %    constructor).
            %
            %  'Output' - true/false flag, or a character vector which must
            %    be 'yes' of 'no'. Determines wheter this node will produce
            %    output. By default output will be produced.
            %
            %
            % Output
            %
            %  str - bdyn.pre.structuralNode object
            %
            %
            %
            % See Also: mbdyn.pre.structuralNode6dof, 
            %           mbdyn.pre.structuralNode3dof
            %

            options.AbsolutePosition = [0;0;0];
            options.AbsoluteVelocity = [0;0;0];
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            options.Accelerations = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.node ( ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output );
            
            self.checkCartesianVector (options.AbsolutePosition, true, 'AbsolutePosition');
            self.checkCartesianVector (options.AbsoluteVelocity, true, 'AbsoluteVelocity');
            
            self.absolutePosition = options.AbsolutePosition;
            self.absoluteVelocity = options.AbsoluteVelocity;
            
            if ~isempty (options.Accelerations)
                
                if self.checkLogicalScalar ( options.Accelerations, false )
                    if options.Accelerations
                        options.Accelerations = 'yes';
                    else
                        options.Accelerations = 'no';
                    end
                elseif self.checkAllowedStringInputs ( options.Accelerations, {'yes', 'no'}, false)
                    % do nothing
                else
                    error ('Accelerations should be a lgical true/false flag, or a string ''yes'', or ''no''') 
                end
                
            end
            
            self.accelerations = options.Accelerations;
            
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
            
            % report name as absolutePosition as this is what user will see
            self.checkCartesianVector (newpos, true, 'absolutePosition');
            
            self.absolutePosition = newpos;
            
        end
        
        function set.absoluteVelocity (self, newvel)
            % set the absolute position of the structural node
            
            % report name as absoluteVelocity as this is what user will see
            self.checkCartesianVector (newvel, true, 'absoluteVelocity');
            
            self.absoluteVelocity = newvel;
            
        end
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            
        end
        
    end
    
end

