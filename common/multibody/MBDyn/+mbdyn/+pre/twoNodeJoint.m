classdef twoNodeJoint < mbdyn.pre.joint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node1;
        node2;
        
    end
    
    methods
        function self = twoNodeJoint (node1, node2, varargin)
            % generic base class for joints which constrains two nodes
            
            [options, nopass_list] = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.joint ( pvpairs{:} );
            
            self.checkIsStructuralNode (node1, true);
            self.checkIsStructuralNode (node2, true);
            
            self.node1 = node1;
            self.node2 = node2;
            
        end
    end
    
    methods
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for twoNodeJoint object
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (tnj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  tnj - mbdyn.pre.twoNodeJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.joint (self);

        end
        
        function abspos = offset2AbsolutePosition (self, offset, refstr, localnodenum)
            % gets an offset in global frame
            
            switch refstr

                case 'global'

                    abspos = offset;

                case {'node', 'local'}

                    if localnodenum == 1
                        abspos = self.node1.relativeToAbsolutePosition (offset);
                    else
                        abspos = self.node2.relativeToAbsolutePosition (offset);
                    end
                    
                case 'other node'
                    
                    if localnodenum == 2
                        abspos = self.node1.relativeToAbsolutePosition (offset);
                    else
                        abspos = self.node2.relativeToAbsolutePosition (offset);
                    end

            end
            
        end
        
        function absorientm = orient2AbsoluteOrientation (self, dorient, refstr, localnodenum)
            % gets an orientation in global frame
            
            switch refstr

                case 'global'

                    absorientm = dorient;

                case {'node', 'local'}

                    if localnodenum == 1
                        absorientm = self.node1.relativeToAbsoluteOrientation (dorient);
                    else
                        absorientm = self.node2.relativeToAbsoluteOrientation (dorient);
                    end
                    
                case 'other node'
                    
                    if localnodenum == 2
                        absorientm = self.node1.relativeToAbsoluteOrientation (dorient);
                    else
                        absorientm = self.node2.relativeToAbsoluteOrientation (dorient);
                    end

            end
            
        end
        
    end
    
    methods (Access = protected)
        
        function processed = checkJointPositionOffset (self, offset)
            % checks and processes the joint position reference frame
            % string
            %
            % Syntax
            %
            % ok = checkJointPositionOffset (jntobj, offset)
            %
            % Description
            %
            %
            %
            % Input
            %
            %  jntobj - mbdyn.pre.twoNodeJoint object
            %
            %  ref - char array specifying the reference frame
            %    in which a position is defined realtive to a node in a
            %    single node joint. Valid strings are: 'node', 'local' and
            %    'global'.
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkNodeReferenceType if ref fails check
            %
            % Output
            %
            %  processed - 
            %
            % See Also: 
            %
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            
                            self.checkAllowedStringInputs (offset{1}, {'global', 'node', 'local', 'other position', 'other node'}, true, 'Position Offset');
                            
                            if ischar (offset{2})
                                if ~strcmp (offset{2}, 'null')
                                    error ('unrecognised offset string (not ''null'')');
                                end
                            else
                                self.checkCartesianVector (offset{2});
                            end
                            
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                    processed = [{'reference'}, offset];
                else
                    self.checkCartesianVector (offset);
                    processed = offset;
                end
            end
            
        end
        
        function processed = checkJointOrientationOffset (self, offset)
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            
                            self.checkAllowedStringInputs (offset{1}, {'global', 'node', 'local', 'other orientation', 'other node'}, true, 'Orientation Offset');
                            
                            if ischar (offset{2})
                                if ~strcmp (offset{2}, 'null')
                                    error ('unrecognised offset string (not ''null'')');
                                end
                            else
                                self.checkOrientationMatrix (offset{2});
                                offset{2} = self.getOrientationMatrix (offset{2});
                            end
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                    processed = [{'reference'}, offset];
                else
                    self.checkCartesianVector (offset);
                    processed = offset;
                end
            end
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.joint.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end