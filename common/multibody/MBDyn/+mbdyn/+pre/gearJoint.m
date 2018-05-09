classdef gearJoint < mbdyn.pre.userDefined
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node1;
        node2;
        referenceNode;
        relativeOrientation1;
        relativeOrientation2;
        orientation1Reference;
        orientation2Reference;
        ratio;
        r1;
        r2;
        output;
        
    end
    
    methods
        function self = gearJoint (node1, node2, refnode, varargin)
            % mbdyn.pre.gearJoint constructor
            
            [options, nopass_list] = mbdyn.pre.gearJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            loadable_mod = mbdyn.pre.loadableModule ('libmodule-fabricate.so');
            
            self = self@mbdyn.pre.userDefined ( loadable_mod, 'GearJoint', pvpairs{:} );
            
            self.checkIsStructuralNode (node1, true);
            self.checkIsStructuralNode (node2, true);
            self.checkIsStructuralNode (refnode, true);
            if ~isempty (options.Orientation1)
                self.checkOrientationMatrix (options.Orientation1, true, 'Orientation1');
            end
            if ~isempty (options.Orientation2)
                self.checkOrientationMatrix (options.Orientation2, true, 'Orientation2');
            end
            allowedorientrefstrs = {'global', 'node', 'local'};
            self.checkAllowedStringInputs ( options.Orientation1Reference, allowedorientrefstrs, true, 'Orientation1Reference');
            self.checkAllowedStringInputs ( options.Orientation2Reference, allowedorientrefstrs, true, 'Orientation2Reference');
            allowedstrs = {'yes', 'no'};
            self.checkAllowedStringInputs ( options.Output, allowedstrs, true, 'Output');
            
            if ~isempty (options.Ratio)
                if ~(isempty (options.R1) && isempty (options.R2))
                    error ('You have supplied a ratio for the gear joint, but either R1 or R2 or both have also been supplied');
                end
                
                self.checkNumericScalar (options.Ratio, true, 'Ratio');
                
                self.ratio = options.Ratio;
                self.r1 = [];
                self.r2 = [];
                
            elseif isempty (options.R1) 
                
                error ('Only R1 has been supplied');
                
            elseif isempty (options.R2)
                
                error ('Only R2 has been supplied');
                
            else
                
                self.checkNumericScalar (options.R1, true, 'R1');
                self.checkNumericScalar (options.R1, true, 'R1');
                
                self.ratio = [];
                self.r1 = options.R1;
                self.r2 = options.R2;
                
            end
            
            self.node1 = node1;
            self.node2 = node2;
            self.referenceNode = refnode;
            self.relativeOrientation1 = options.Orientation1;
            self.relativeOrientation2 = options.Orientation2;
            self.orientation1Reference = options.Orientation1Reference;
            self.orientation2Reference = options.Orientation2Reference;
            self.output = options.Output;
            
        end
    end
    
    methods
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.userDefined(self);
            
            str = self.addOutputLine ( str, ...
                                       sprintf ('%d', self.node1.label), ...
                                       2, ...
                                       true, ...
                                       'gear 1 node' );
                                   
            if ~isempty (self.relativeOrientation1)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'orientation', ...
                                                               'reference', ...
                                                               self.orientation1Reference, ...
                                                               self.relativeOrientation1 ), ...
                                           3, ...
                                           true );
            end
            
            str = self.addOutputLine ( str, ...
                                       sprintf ('%d', self.node2.label), ...
                                       2, ...
                                       true, ...
                                       'gear 2 node' );
                                      
            if ~isempty (self.relativeOrientation2)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'orientation', ...
                                                               'reference', ...
                                                               self.orientation2Reference, ...
                                                               self.relativeOrientation2 ), ...
                                           3, ...
                                           true );
            end
            
            str = self.addOutputLine ( str, ...
                                       sprintf ('%d', self.referenceNode.label), ...
                                       2, ...
                                       true, ...
                                       'reference node' );
                                   
            if isempty (self.ratio)
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( self.r1, ...
                                                               self.r2 ), ...
                                           2, ...
                                           true );

            else
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'ratio', ...
                                                               self.ratio ), ...
                                           2, ...
                                           true );
                                       
            end
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList ( 'output', ...
                                                           self.output), ...
                                       2, ...
                                       false );
                                   
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.userTypeName));

        end
        
    end
    
%     methods (Access = protected)
%         
%         function processed = checkJointOrientationOffset (self, offset)
%             
%             if~isempty (offset)
%                 if iscell (offset)
%                     if numel (offset) == 2
% 
%                         if ischar (offset{1})
%                             
%                             self.checkAllowedStringInputs (offset{1}, {'global', 'node', 'local', 'other orientation', 'other node'}, true, 'Orientation Offset');
%                             
%                             if ischar (offset{2})
%                                 if ~strcmp (offset{2}, 'null')
%                                     error ('unrecognised offset string (not ''null'')');
%                                 end
%                             else
%                                 self.checkOrientationMatrix (offset{2});
%                                 offset{2} = self.getOrientationMatrix (offset{2});
%                             end
%                         else
%                             error ('First offset value must be a char array.')
%                         end
% 
%                     else
%                         error ('If offset is supplied as a cell array it must have only 2 elements')
%                     end
%                     processed = [{'reference'}, offset];
%                 else
%                     self.checkCartesianVector (offset);
%                     processed = offset;
%                 end
%             end
%         end
%         
%     end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.userDefined.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.Ratio = [];
            options.R1 = [];
            options.R2 = [];
            options.Orientation1 = [];
            options.Orientation1Reference = 'node';
            options.Orientation2 = [];
            options.Orientation2Reference = 'node';
            options.Output = 'no';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end