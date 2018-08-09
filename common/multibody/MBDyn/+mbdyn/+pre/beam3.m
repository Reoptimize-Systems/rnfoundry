classdef beam3 < mbdyn.pre.beam
    
    properties
        
        relativeOffset1;
        relativeOrientation1;
        offset1Reference;
        orientation1Reference;
            
        relativeOffset2;
        relativeOrientation2;
        offset2Reference;
        orientation2Reference;
        
        relativeOffset3;
        relativeOrientation3;
        offset3Reference;
        orientation3Reference;
        
        node1;
        node2;
        node3;
        
        constituativeLawSectionI;
        constituativeLawSectionII;
        
        orientationMatrixSection_I;
        orientationMatrixSection_II;
        
        customOutput;
        
    end
    
    methods
        
        function self = beam3 (node1, node2, node3, orientm_section_I, beam_const_law_I, orientm_section_II, beam_const_law_II, varargin)
            % mbdyn.pre.beam3 constructor
            
            % <elem_type> ::= beam3
            % <normal_arglist> ::=
            % <node_1_label> ,
            % [ position , ] (Vec3) <relative_offset_1> ,
            % [ orientation , (OrientationMatrix) <relative_orientation_1> , ]
            % <node_2_label> ,
            % [ position , ] (Vec3) <relative_offset_2> ,
            % [ orientation , (OrientationMatrix) <relative_orientation_2> , ]
            % <node_3> ,
            % [ position , ] (Vec3) <relative_offset_3> ,
            % [ orientation , (OrientationMatrix) <relative_orientation_3> , ]
            % { (OrientationMatrix) <orientation_matrix_section_I>
            % | from nodes } ,
            % (ConstitutiveLaw<6D>) <constitutive_law_section_I> ,
            % { same
            % | (OrientationMatrix) <orientation_matrix_section_II>
            % | from nodes } ,
            % { same
            % | (ConstitutiveLaw<6D>) <constitutive_law_section_II> }
            % [ , <custom_output> ]


            [options, nopass_list] = mbdyn.pre.beam3.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.beam ( pvpairs{:} );
            
            self.checkIsStructuralNode (node1, true);
            self.checkIsStructuralNode (node2, true);
            self.checkIsStructuralNode (node3, true);
            
            if ischar (orientm_section_I)
                allowed_orientm_section_I_strs = {'from nodes'};
                self.checkAllowedStringInputs ( orientm_section_I, allowed_orientm_section_I_strs, true, 'orientm_section_I');
            else
                self.checkOrientationMatrix (orientm_section_I, true, 'orientm_section_I');
            end
            
            assert (isa (beam_const_law_I, 'mbdyn.pre.constituativeLaw'), ...
                'beam_const_law_I must be an mbdyn.pre.constituativeLaw object or derived' );
            
            if ischar (orientm_section_II)
                allowed_orientm_section_II_strs = {'from nodes', 'same'};
                self.checkAllowedStringInputs ( orientm_section_II, allowed_orientm_section_II_strs, true, 'orientm_section_II');
            else
                self.checkOrientationMatrix (orientm_section_II, true, 'orientm_section_II');
            end
            
            if ischar (beam_const_law_II)
                allowed_beam_const_law_II_strs = {'same'};
                self.checkAllowedStringInputs ( beam_const_law_II, allowed_beam_const_law_II_strs, true, 'beam_const_law_II');
            else
                assert (isa (beam_const_law_II, 'mbdyn.pre.constituativeLaw'), ...
                    'beam_const_law_II must be an mbdyn.pre.constituativeLaw object or derived' );
            end
            
            allowedposrefstrs = {'global', 'node', 'local'};
            allowedorientrefstrs = {'global', 'node', 'local'};
            self.checkAllowedStringInputs ( options.Offset1Reference, allowedposrefstrs, true, 'Offset1Reference');
            self.checkAllowedStringInputs ( options.Offset2Reference, allowedposrefstrs, true, 'Offset2Reference');
            self.checkAllowedStringInputs ( options.Offset3Reference, allowedposrefstrs, true, 'Offset3Reference');
            self.checkAllowedStringInputs ( options.Orientation1Reference, allowedorientrefstrs, true, 'Orientation1Reference');
            self.checkAllowedStringInputs ( options.Orientation2Reference, allowedorientrefstrs, true, 'Orientation2Reference');
            self.checkAllowedStringInputs ( options.Orientation3Reference, allowedorientrefstrs, true, 'Orientation3Reference');
            self.checkCartesianVector (options.RelativeOffset1, true, 'RelativeOffset1');
            self.checkCartesianVector (options.RelativeOffset2, true, 'RelativeOffset2');
            self.checkCartesianVector (options.RelativeOffset3, true, 'RelativeOffset3');
            self.checkOrientationMatrix (options.RelativeOrientation1, true, 'RelativeOrientation1');
            self.checkOrientationMatrix (options.RelativeOrientation2, true, 'RelativeOrientation2');
            self.checkOrientationMatrix (options.RelativeOrientation3, true, 'RelativeOrientation3');
            
            if ~isempty (options.CustomOutput)
                
                allowedoutputstrs = { ...
                                      'position', ...
                                      'orientation', ... [ , <orientation_description> ]
                                      'configuration', ... % same as: position, orientation
                                      'force', ...
                                      'moment', ...
                                      'forces', ... % same as: force, moment
                                      'linear strain', ...
                                      'angular strain', ...
                                      'strains', ... % same as: linear strain, angular strain
                                      'linear strain rate', ...
                                      'angular strain rate', ...
                                      'strain rates', ... % same as: linear strain rate, angular strain rate
                                      'all' };
                                  
                if ischar (options.CustomOutput)
                    
                    self.checkAllowedStringInputs ( options.CustomOutput, allowedoutputstrs, true, 'CustomOutput');
                    
                elseif iscellstr (options.CustomOutput)
                    
                    k = 1;
                    while k < numel (options.CustomOutput)
                        
                        if strcmp (options.CustomOutput{k}, 'orientation')
                            
                            % check if the next thing is a valid
                            % orientation description
                            isorientdesc = self.checkOrientationDescription (options.CustomOutput{k+1}, false);
                            
                            if isorientdesc
                                % skip forward 1 and move to the next item
                                k = k + 2;

                                continue
                            
                            end
                            
                        end
                        
                        self.checkAllowedStringInputs ( options.CustomOutput{k}, allowedoutputstrs, true, sprintf ('CustomOutput{%d}', k));
                        
                        k = k + 1;
                        
                    end
                    
                else
                    
                end
                                  
            end
            
            self.relativeOffset1 = options.RelativeOffset1;
            self.relativeOrientation1 = options.RelativeOrientation1;
            self.relativeOffset2 = options.RelativeOffset2;
            self.relativeOrientation2 = options.RelativeOrientation2;
            self.relativeOffset3 = options.RelativeOffset3;
            self.relativeOrientation3 = options.RelativeOrientation3;
            self.offset1Reference = options.Offset1Reference;
            self.orientation1Reference = options.Orientation1Reference;
            self.offset2Reference = options.Offset2Reference;
            self.orientation2Reference = options.Orientation2Reference;
            self.offset3Reference = options.Offset3Reference;
            self.orientation3Reference = options.Orientation3Reference;
            
            self.node1 = node1;
            self.node2 = node2;
            self.node3 = node3;

            self.constituativeLawSectionI = beam_const_law_I;
            self.constituativeLawSectionII = beam_const_law_II;
        
            self.orientationMatrixSection_I = orientm_section_I;
            self.orientationMatrixSection_II = orientm_section_II;
            
            self.customOutput = options.CustomOutput;
            
            self.type = 'beam3';
            
%             self.netCDFName = 'joint';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for beam3
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (beam3)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  beam3 - mbdyn.pre.beam3 object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = sprintf ('    %s : %d,', self.type, self.label);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeOffset1)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset1Reference, self.relativeOffset1), 3, true);
            end
            
            if ~isempty (self.relativeOrientation1)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation1Reference, self.relativeOrientation1), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
            if ~isempty (self.relativeOffset2)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset2Reference, self.relativeOffset2), 3, true);
            end
            
            if ~isempty (self.relativeOrientation2)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation2Reference, self.relativeOrientation2), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node3.label), 2, true, 'node 3 label');
            
            if ~isempty (self.relativeOffset3)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset3Reference, self.relativeOffset3), 3, true);
            end
            
            if ~isempty (self.relativeOrientation3)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation3Reference, self.relativeOrientation3), 3, true);
            end
            
            str = self.addOutputLine (str, self.orientationMatrixSection_I, 2, true);
            
            str = self.addOutputLine (str, self.constituativeLawSectionI.generateMBDynInputString (), 2, true);
            
            str = self.addOutputLine (str, self.orientationMatrixSection_II, 2, true);
            
            addcomma = ~isempty (self.customOutput);
            
            if ischar (self.constituativeLawSectionII)
                str = self.addOutputLine (str, self.constituativeLawSectionII, 2, addcomma);
            else
                str = self.addOutputLine (str, self.constituativeLawSectionII.generateMBDynInputString (), 2, addcomma);
            end
            
            if ~isempty (self.customOutput)
                str = self.addOutputLine (str, self.commaSepList ('custom output', self.customOutput{:}), 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));

        end
        
    end
    
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.beam.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.RelativeOffset1 = [];
            options.RelativeOffset2 = [];
            options.RelativeOffset3 = [];
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.RelativeOrientation3 =  [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Offset3Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            options.Orientation3Reference = 'node';
            options.CustomOutput = {};
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end