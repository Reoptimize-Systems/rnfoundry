classdef stateSpaceMIMO < mbdyn.pre.stateSpaceFilter
    
    properties (GetAccess = public, SetAccess = protected)
        
        numberOfOutputs;
        outputNodeDOFs;
        
        numberOfInputs;
        inputList;
        
    end
    
    methods
        
        function self = stateSpaceMIMO (state_order, A, B, C, output_node_list, input_list, varargin)
            
            options.E = [];
            options.D = [];
            options.Gain = []; % If a gain is supplied, all the coefficients of B are multiplied by the gain.
            options.Balance = '';
            options.Value = [];
            options.Derivative = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                            'E', options.E, ...
                            'D', options.D, ...
                            'Gain', options.Gain, ...
                            'Balance', options.Balance, ...
                            'Value', options.Value, ...
                            'Derivative', options.Derivative);
                        
                        
             if ~isa (output_node_list, 'mbdyn.pre.nodeDOF') 
                error ('output_node_list must be an array of one or more mbdyn.pre.nodeDOF objects')
             end
            
             if ~iscell (input_list)
                error ('input_list must be a cell array')
             end
             
             for ind = 1:numel (input_list)
                 if ~ (isa (input_list{ind}, 'mbdyn.pre.nodeDOF') ...
                       || isa (input_list{ind}, 'mbdyn.pre.drive') )
                    error ('input_list must be a cell array of mbdyn.pre.nodeDOF objects and/or mbdyn.pre.driveCaller objects')
                 end
             end
             
             self.numberOfOutputs = numel (output_node_list);
             self.outputNodeDOFs = output_node_list;
             
             self.numberOfInputs = numel (input_list);
             self.inputList = input_list;
             
             self.type = 'state space MIMO';
            
        end
        
        function str = generateOutputString (self)
            
            % base indent level is 0ne
            str = generateOutputString@mbdyn.pre.genel(self);
            
            str = self.addOutputLine (str, self.commaSepList (self.numberOfOutputs), 2, true, 'number of outputs');
            
            for ind = 1:self.numberOfOutputs
                
                str = self.addOutputLine (str, self.commaSepList (self.outputNodeDOFs(ind).generateOutputString ()), 3, true);
                
            end
            
            str = self.addOutputLine (str, self.commaSepList (self.numberOfInputs), 2, true, 'number of inputs');
            
            for ind = 1:self.numberOfInputs
                
                if isa (self.inputList{ind}, 'mbdyn.pre.nodeDOF')
                    
                    str = self.addOutputLine (str, self.commaSepList ('node dof', self.outputNodeDOFs(ind).generateOutputString ()), 3, true);
                    
                elseif isa (self.inputList{ind}, 'mbdyn.pre.drive') 
                    
                    str = self.addOutputLine (str, self.commaSepList ('drive', self.outputNodeDOFs(ind).generateOutputString ()), 3, true);
                    
                end
                
            end
            
            % generate common part of stste-space filter
            str = self.addOutputLine (str, generateOutputString@mbdyn.pre.stateSpaceFilter(self), 2, false); 
            
            str = self.addOutputLine (str, ';', 1, false, 'end state space MIMO');
            
        end
        
%         function hax = draw (self, varargin)
%             
%             options.AxesHandle = [];
%             options.ForceRedraw = false;
%             options.Mode = 'solid';
%             options.Light = false;
%             
%             options = parse_pv_pairs (options, varargin);
%             
%             hax = draw@mbdyn.pre.element ( self, ...
%                     'AxesHandle', options.AxesHandle, ...
%                     'ForceRedraw', options.ForceRedraw, ...
%                     'Mode', options.Mode, ...
%                     'Light', options.Light );
% 
%             self.setTransform ();
%             
%         end
        
    end
    
    methods (Access = protected)
        
        
        
    end
    
end