classdef joint < mbdyn.pre.element
    
    properties
        regularizeCoeff;
        regularizationType;
    end
    
    methods
        
        function self = joint (varargin)
            % generic base class for joints which constrains two nodes
            
            [options, nopass_list] = mbdyn.pre.joint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.element ( pvpairs{:} );
            
            mbdyn.pre.base.checkAllowedStringInputs (options.RegularizationType, {'tikhonov'}, true, 'RegularizationType');
            if ~isempty (options.RegularizeCoeff)
                self.regularizeCoeff = options.RegularizeCoeff;
            end
            
            self.netCDFName = 'joint';
            self.regularizationType = options.RegularizationType;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for common joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (jnt)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  jnt - mbdyn.pre.joint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.element (self);
            
%             if ~isempty (self.regularizeCoeff)
%                 str = sprintf ('    joint regularization : %d, tikhonov, %s;\n\n', self.label, self.formatNumber (self.regularizeCoeff));
%             end
            str = self.addOutputLine (str, sprintf ('joint : %d, %s', self.label, self.type), 1);
            
%             str = sprintf ('%s    joint : %d, %s,', str, self.label, self.type);
        end
        
    end
    
    methods (Access = protected)
        
        function str = addRegularization (self, str)
            
            if ~isempty (self.regularizeCoeff)
                
                if numel (self.regularizeCoeff) == 1
                    
                    str = self.addOutputLine ( str, ...
                                               sprintf ( 'joint regularization : %d, %s, %s;', ...
                                                         self.label, ...
                                                         self.regularizationType, ...
                                                         self.formatNumber (self.regularizeCoeff) ), ...
                                               1, ...
                                               false );
                else
                    
                    str = self.addOutputLine ( str, ...
                                               sprintf ( 'joint regularization : %d, %s, list, %s;', ...
                                                         self.label, ...
                                                         self.regularizationType, ...
                                                         self.commaSepList (self.regularizeCoeff) ), ...
                                               1, ...
                                               false );
                end
                
            end
            
        end
        

        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.element.defaultConstructorOptions ();
            
            options.RegularizeCoeff = [];
            options.RegularizationType = 'tikhonov';
            
            nopass_list = {'RegularizeCoeff', 'RegularizationType'};
            
        end
        
    end
    
end