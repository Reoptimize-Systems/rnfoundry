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
%  generateMBDynInputString - make MBDyn output string for node
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
        sx;
        sy;
        sz;
        
        accelerations; % flag indicating whether to output acceleration data
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)

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
            
            self.setColour ( [ 0.635, 0.078, 0.184; % red
                               0, 0.447, 0.741; % blue
                               0.466, 0.674, 0.188 ] ); % green
           
            self.drawAxesH = [];
            self.sx = 1;
            self.sy = 1;
            self.sz = 1;
            
            self.netCDFName = 'struct';
            
        end
        
        function str = generateMBDynInputString (self)
            str = generateMBDynInputString@mbdyn.pre.node (self);
        end
        
        function draw (self, varargin)
            % plot the node in a figure
            %
            % Syntax
            %
            % hax = draw (nd)
            % hax = draw (..., 'Parameter', Value)
            %
            % Description
            %
            % The draw method creates a visualisation of the node in a
            % figure. The node is plotted as three lines. These are the 3
            % axes in the node frame of reference and rotate with the node.
            % The lines are red, blue and green by default.
            %
            % The draw method tries to be efficient. Each node is
            % associated with a hgtransform object. If the shape of the
            % node has not changed, the transform matrix is simply
            % updated to adjust the location and orientation on the plot. 
            %
            % Input
            %
            %  nd - mbdyn.pre.structuralNode object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'AxesHandle' - optional handle to axes in which to plot the
            %    node. If not supplied, a new figure and axes will be
            %    created. The first time the node is drawn the handle to
            %    the axes will be stored internally and future calls to
            %    draw will plot to the same axes, unless the axes are
            %    destroyed, or this option is used to override it.
            %
            %  'ForceRedraw' - true/false flag indicating whether to force
            %    a full redraw of the node (rather than just update the
            %    transform matrix), even if the node does not think it
            %    needs it.
            %
            %
            % See Also: mbdyn.pre.element.draw, 
            %           mbdyn.pre.system.draw
            %
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            
            options = parse_pv_pairs (options, varargin);
            
            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects). We will recreate all
                % the objects
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                self.shapeObjects = {};
                
                self.shapeObjects = { line([-self.sx/2; 
                                             self.sx/2; ], ...
                                          [  0;
                                             0 ], ...
                                          [  0;
                                             0 ], ...
                                          'Parent', self.transformObject, ...
                                          'Color', self.drawColour(1,1:3) ), ...
                                      line([ 0; 
                                             0; ], ...
                                          [  -self.sy/2;
                                             self.sy/2  ], ...
                                          [  0;
                                             0 ], ...
                                          'Parent', self.transformObject, ...
                                          'Color', self.drawColour(2,1:3) ), ...
                                      line([0; 
                                            0; ], ...
                                          [ 0;
                                            0  ], ...
                                          [ -self.sz/2;
                                            self.sz/2 ], ...
                                          'Parent', self.transformObject, ...
                                          'Color', self.drawColour(3,1:3) ) ...
                                     };
            end
            
            self.setTransform ();
            
        end
        
        function setSize (self, sx, sy, sz)
            % sets the size of an mbdyn.pre.structuralNode when plotted in a figure
            %
            % Syntax
            %
            % setSize (nd, sx, sy, sz)
            %
            % Description
            %
            % setSize is used to set the size of the node when plotting the
            % element in a figure. This sets the length of the three lines
            % representing the node in the plot.
            %
            % Input
            %
            %  nd - mbdyn.pre.structuralNode object
            %
            %  sx - length along the node x axis
            %
            %  sy - length along the node y axis
            %
            %  sz - length along the node z axis
            %
            %
            
            self.checkNumericScalar (sx, true, 'sx');
            self.checkNumericScalar (sy, true, 'sy');
            self.checkNumericScalar (sz, true, 'sz');

            assert (sx > 0, 'sx must be greater than zero');
            assert (sy > 0, 'sy must be greater than zero');
            assert (sz > 0, 'sz must be greater than zero');
            
            self.sx = sx;
            self.sy = sy;
            self.sz = sz;
            
            % size has changes so objects must be recreated
            self.needsRedraw = true;
            
        end
        
        function setColour (self, newcolour)
            
            ok = self.check3X3Matrix (newcolour, false);
            
            if ~ok
                error ('newcolour should be a 3x3 where each row is the colour of the x, y and z node axes respectively');
            end
            
            self.drawColour = newcolour;
        end
        
    end
    
    % getters/setters
    methods
%         function set.absolutePosition (self, newpos)
%             % set the absolute position of the structural node
%             
%             % report name as absolutePosition as this is what user will see
%             self.checkCartesianVector (newpos, true, 'absolutePosition');
%             
%             self.absolutePosition = newpos;
%             
%         ended
%         
%         function set.absoluteVelocity (self, newvel)
%             % set the absolute position of the structural node
%             
%             % report name as absoluteVelocity as this is what user will see
%             self.checkCartesianVector (newvel, true, 'absoluteVelocity');
%             
%             self.absoluteVelocity = newvel;
%             
%         end
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            % setTransform stub function
            %
            % Description
            %
            % Does nothing, child classes are expected to overload this
            % function to set the transform matrix as appropriate for
            % plotting the structural node
            %
            
            
        end
        
    end
    
end

