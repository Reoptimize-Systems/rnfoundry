classdef viscousBody < mbdyn.pre.singleNodeOffsetJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
 
        constitutiveLaw;
        
    end
    
    methods
        
        function self = viscousBody (node, const_law, relative_offset, varargin)
            % viscousBody constructor
            %
            % Syntax
            %
            %  rh = viscousBody (node1, const_law)
            %  rh = viscousBody (..., 'Parameter', value)
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode (or derived class) object
            %    representing the node the joint connects
            %
            %  const_law - mbdyn.pre.constituativeLaw object representing a
            %   6D constituative law.
            %
            %  relative_offset - (3 x 1) vector containing the relative
            %   offset of the force from the node, or the character vector
            %   'null', which is equivalent to [0;0;0]. By default this
            %   offset is expressed in the reference frame of the node,
            %   alternative reference frames can be selected using the
            %   'PositionReference' option described below.
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %  'RelativeOrientation' - mbdyn.pre.orientation object
            %    defining the an orientation offset of the force from the
            %    node. By default this offset is expressed in the reference
            %    frame of the node, alternative reference frames can be
            %    selected using the 'OrientationReference' option described
            %    below.
            %
            %  'PositionReference' - by default the position provided in
            %    RelativeOffset is relative to the node reference frame. An
            %    alternative reference frame can be provided using this
            %    argument. Possible value for this are:
            %      'node'          : the default behaviour
            %      'global'        : the global reference frame
            %      'local'         : same as 'node'  
            %
            %  'OrientationReference' - by default the position provided in
            %    RelativeOrientation is relative to the node reference
            %    frame. An alternative reference frame can be provided
            %    using this argument. Possible value for this are:
            %      'node'          : the default behaviour
            %      'global'        : the global reference frame
            %      'local'         : same as 'node'  
            %
            % Output
            %
            %  rh - mbdyn.pre.viscousBody object
            %
            %
            % See also: mbdyn.pre.linearViscousGenericConstituativeLaw, 
            %           mbdyn.pre.linearViscousIsotropicConstituativeLaw
            %
            
            [ options, nopass_list ] = mbdyn.pre.viscousBody.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeOffsetJoint (node, pvpairs{:}, 'RelativeOffset', relative_offset);
                    
            assert (isa (const_law, 'mbdyn.pre.constituativeLaw'), ...
                'law must be an mbdyn.pre.constituativeLaw' );
            
            self.constitutiveLaw = const_law;
            self.type = 'viscous body';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for viscousBody joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (sp)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  sp - mbdyn.pre.viscousBody object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %

            str = generateMBDynInputString@mbdyn.pre.singleNodeOffsetJoint (self, true);
                    
            str = self.addOutputLine (str, self.constitutiveLaw.generateMBDynInputString (), 2, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
            str = self.addRegularization (str);
            
        end
        
%         function draw (self, varargin)
%             
%             options.AxesHandle = [];
%             options.ForceRedraw = false;
%             options.Mode = 'solid';
%             
%             options = parse_pv_pairs (options, varargin);
%             
% %             draw@mbdyn.pre.twoNodeOffsetJoint ( self, ...
% %                 'AxesHandle', options.AxesHandle, ...
% %                 'ForceRedraw', options.ForceRedraw, ...
% %                 'Mode', options.Mode );
% 
%             if options.ForceRedraw
%                 self.needsRedraw = true;
%             end
%             
%             self.checkAxes (options.AxesHandle);
%             
%             node1pos = self.node.absolutePosition;
%             jref = self.reference ();
%             jpos = jref.pos ();
%                 
%             if ~self.needsRedraw
%                 % always have to redraw the line connecting the two points.
%                 % This changes shape, so we can't just transform the line
%                 % object
%                 delete (self.shapeObjects{1})
%                 self.shapeObjects{1} =  line ( self.drawAxesH, ...
%                                                [ node1pos(1), jpos(1) ], ...
%                                                [ node1pos(2), jpos(2) ], ...
%                                                [ node1pos(3), jpos(3) ], ...
%                                                'Color', self.drawColour );
%                                        
%             end
%             
%             if isempty (self.shapeObjects) ...
%                     || self.needsRedraw
%                 % a full redraw is needed (and not just a modification of
%                 % transform matrices for the objects).
%                 
%                 % delete the current patch object
%                 self.deleteAllDrawnObjects ();
%                 
%                 self.shapeObjects = { line( self.drawAxesH, ...
%                                             [ node1pos(1), jpos(1) ], ...
%                                             [ node1pos(2), jpos(2) ], ...
%                                             [ node1pos(3), jpos(3) ], ...
%                                             'Color', self.drawColour ), ...
%                                       line( self.drawAxesH, ...
%                                             [ 0, 0 ], ...
%                                             [ 0, 0 ], ...
%                                             [ -self.shapeParameters(1)/2, self.shapeParameters(3)/2 ], ...
%                                             'Parent', self.transformObject, ...
%                                             'Color', self.drawColour, ...
%                                             'LineStyle', '--' )
%                                     };
%                 
%                 self.needsRedraw = false;
% 
% %                 if options.Light
% %                     light (self.drawAxesH);
% %                 end
%                 
%             end
%             
%             self.setTransform ();
% 
% %             self.setTransform ();
%             
%         end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
%             
%             om = self.absoluteJointOrientation;
%             
%             M = [ om.orientationMatrix, self.absoluteJointPosition; ...
%                   0, 0, 0, 1 ];
%                   
%             set ( self.transformObject, 'Matrix', M );
%             
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.singleNodeOffsetJoint.defaultConstructorOptions ();
            
            nopass_list = {'RelativeOffset'};
            
        end
        
    end
    
end