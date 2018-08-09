classdef beamSlider < mbdyn.pre.singleNodeJoint
    
    properties
        
        beams;
        
        sliderType;
        relativeOffset;
        relativeOrientation;
        offsetReference;
        orientationReference;
        initialBeam;
        initialNode;
        smearingFactor;
        
    end
    
    methods
        
        function self = beamSlider (node, rel_offset, beams, varargin)
            % mbdyn.pre.beamSlider constructor

            % <joint_type> ::= kinematic
            % <joint_arglist> ::=
            % <slider_node_label> ,
            % (Vec3) <relative_offset> ,
            % [ hinge , (OrientationMatrix) <relative_orientation> ] ,
            % [ type , { spherical | classic | spline } , ]
            % <beam_number> ,
            % <3_node_beam> ,
            % { same | (Vec3) <first_node_offset> } ,
            % [ hinge , { same | (OrientationMatrix) <first_node_orientation> ] , }
            % (Vec3) <mid_node_offset> ,
            % [ hinge , (OrientationMatrix) <mid_node_orientation> ] ,
            % 226
            % (Vec3) <end_node_offset> ,
            % [ hinge , (OrientationMatrix) <end_node_orientation> ] ,
            % [ ... ]
            % [ , initial beam , <initial_beam> ]
            % [ , initial node , <initial_node> ]
            % [ , smearing , <smearing_factor> ]


            [options, nopass_list] = mbdyn.pre.beamSlider.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.singleNodeJoint ( node, pvpairs{:} );
            
            self.checkCartesianVector (rel_offset, true, 'rel_offset');
            
            assert (isstruct (isstruct), 'beams should be a structure');
            
            self.checkStructHasAllFields ( beams, ...
                                           fieldnames (self.makeSliderBeamStruct ()), ...
                                           true, ...
                                           'beams');
            
            for ind = 1:numel (beams)
                
                beamspec = beams(ind);
                
                assert (isa (beamspec.Beam, 'mbdyn.pre.beam'), ...
                    'beams(%d).Beam is not an mbdyn.pre.beam object', ind );

                assert (isstruct (beamspec), 'beams is not a structure array', ind);
                assert (isfield (beamspec, 'FirstNodeOffset'), 'FirstNodeOffset field is missing for beam %d', ind);
                assert (isfield (beamspec, 'MidNodeOffset'), 'MidNodeOffset field is missing for beam %d', ind);
                assert (isfield (beamspec, 'EndNodeOffset'), 'EndNodeOffset field is missing for beam %d', ind);
                
                if ischar (beamspec.FirstNodeOffset)
                    
                    if strcmp (beamspec.FirstNodeOffset, 'same')
                   
                        if ind == 1
                            error ('First beam cannot have a FirstNodeOffset field containing a keyword');
                        end
                        
                    else
                        % could be 'null' keyword
                        self.checkCartesianVector (beamspec.FirstNodeOffset, true, sprintf ('beams(%d).FirstNodeOffset', ind));
                    end
                    
                else
                    self.checkCartesianVector (beamspec.FirstNodeOffset, true, sprintf ('beams(%d).FirstNodeOffset', ind));
                end
                
                self.checkCartesianVector (beamspec.MidNodeOffset, true, sprintf ('beams(%d).MidNodeOffset', ind));
                self.checkCartesianVector (beamspec.EndNodeOffset, true, sprintf ('beams(%d).EndNodeOffset', ind));
                
                allowedrefstrs = {'global', 'node', 'local'};
                
                if isfield (beamspec.FirstNodeOffsetReference) && ~isempty (beamspec.FirstNodeOffsetReference)
                    
                    self.checkAllowedStringInputs ( beamspec.FirstNodeOffsetReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).FirstNodeOffsetReference', ind) );
                                                
                end

                if isfield (beamspec.MidNodeOffsetReference) && ~isempty (beamspec.MidNodeOffsetReference)
                    self.checkAllowedStringInputs ( beamspec.MidNodeOffsetReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).MidNodeOffsetReference', ind) );
                end
                
                if isfield (beamspec.EndNodeOffsetReference) && ~isempty (beamspec.EndNodeOffsetReference)
                    self.checkAllowedStringInputs ( beamspec.EndNodeOffsetReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).EndNodeOffsetReference', ind) );
                end
                
               
                if isfield (beamspec, 'FirstNodeOrientation') && ~isempty (beamspec.FirstNodeOrientation)
                    
                    self.checkOrientationMatrix ( beamspec.FirstNodeOrientation, ...
                                                  true, ...
                                                  sprintf ('beams(%d).FirstNodeOrientation', ind) );
                    
                end
                
                if isfield (beamspec, 'MidNodeOrientation') && ~isempty (beamspec.MidNodeOrientation)
                    
                    self.checkOrientationMatrix ( beamspec.MidNodeOrientation, ...
                                                  true, ...
                                                  sprintf ('beams(%d).MidNodeOrientation', ind) );
                    
                end
                
                if isfield (beamspec, 'EndNodeOrientation') && ~isempty (beamspec.EndNodeOrientation)
                    
                    self.checkOrientationMatrix ( beamspec.EndNodeOrientation, ...
                                                  true, ...
                                                  sprintf ('beams(%d).EndNodeOrientation', ind) );
                    
                end
                
                if isfield (beamspec.FirstNodeOrientationReference) && ~isempty (beamspec.FirstNodeOrientationReference)
                    
                    self.checkAllowedStringInputs ( beamspec.FirstNodeOrientationReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).FirstNodeOrientationReference', ind) );
                                                
                end

                if isfield (beamspec.MidNodeOrientationReference) && ~isempty (beamspec.MidNodeOrientationReference)
                    self.checkAllowedStringInputs ( beamspec.MidNodeOrientationReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).MidNodeOrientationReference', ind) );
                end
                
                if isfield (beamspec.EndNodeOrientationReference) && ~isempty (beamspec.EndNodeOrientationReference)
                    self.checkAllowedStringInputs ( beamspec.EndNodeOrientationReference, ...
                                                    allowedrefstrs, ...
                                                    true, ...
                                                    sprintf ('beams(%d).EndNodeOrientationReference', ind) );
                end
                
            end
            
            
            
            if ~isempty (options.SliderType)
                self.checkAllowedStringInputs ( options.SliderType, {'spherical', 'classic', 'spline'}, true, 'type');
            end

            if ~isempty (options.RelativeOrientation)
                self.checkOrientationMatrix (options.RelativeOrientation, true, 'RelativeOrientation');
            end

            self.checkAllowedStringInputs ( options.OffsetReference, {'global', 'node', 'local'}, true, 'OffsetReference');
            self.checkAllowedStringInputs ( options.OffsetReference, {'global', 'node', 'local'}, true, 'OffsetReference');

            if ~isempty (options.InitialBeam)
                self.checkScalarInteger (options.InitialBeam, true, 'InitialBeam');
            end

            if ~isempty (options.InitialNode)
                self.checkScalarInteger (options.InitialNode, true, 'InitialNode');
            end

            if ~isempty (options.SmearingFactor)
                self.checkNumericScalar (options.SmearingFactor, true, 'SmearingFactor');
            end
            
            self.relativeOffset = rel_offset;
            self.beams = beams;
            self.sliderType = options.SliderType;
            self.relativeOrientation = options.RelativeOrientation;
            self.offsetReference = options.OffsetReference;
            self.orientationReference = options.OrientationReference;
            self.initialBeam = options.InitialBeam;
            self.initialNode = options.InitialNode;
            self.smearingFactor = options.SmearingFactor;
            
            
            self.type = 'beam slider';
            
%             self.netCDFName = 'joint';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for beamSlider
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (bs)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  bs - mbdyn.pre.beamSlider object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint ();
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'slider node label');
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList ('reference', self.offsetReference, self.relativeOffset), ...
                                       2, ...
                                       true, ...
                                       'relative offset' );
                                   
            if ~isempty (self.relativeOrientation)
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('hinge', 'reference', self.orientationReference, self.relativeOrientation), ...
                                           2, ...
                                           true, ...
                                           'relative orientation' );
            end
            
            if ~isempty (self.type)
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('type', self.sliderType), ...
                                           2, ...
                                           true, ...
                                           'type' );
                                       
            end
            
            str = self.addOutputLine ( str, ...
                                       sprintf ('%d', numel (self.beams)), ...
                                       2, ...
                                       true, ...
                                       'number of beams' );
                                       
            for ind = 1:numel (self.beams)
                
                str = self.addOutputLine ( str, ...
                                           sprintf ('%d', self.beams(ind).Beam.label), ...
                                           3, ...
                                           true, ...
                                           sprintf ('beam %d label', ind) );
                                       
                if ind == 1 || ~ischar (self.beams(ind).FirstNodeOffset)
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList ('reference', self.beams(ind).FirstNodeOffsetReference, self.beams(ind).FirstNodeOffset), ...
                                               4, ...
                                               true, ...
                                               sprintf ('beam %d first_node_offset', ind) );
                else
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList (self.beams(ind).FirstNodeOffset), ...
                                               4, ...
                                               true, ...
                                               sprintf ('beam %d first_node_offset', ind) );
                    
                end
                
                if ~isempty (self.beams.FirstNodeOrientation)
                    
                    if ind == 1 || ~ischar (self.beams(ind).FirstNodeOrientation)
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', 'reference', self.beams(ind).FirstNodeOrientationReference, self.beams(ind).FirstNodeOrientation), ...
                                                   4, ...
                                                   true, ...
                                                   sprintf ('beam %d first_node_orientation', ind) );
                    else
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', self.beams(ind).FirstNodeOrientation), ...
                                                   4, ...
                                                   true, ...
                                                   sprintf ('beam %d first_node_orientation', ind) );

                    end
                    
                end
                
                if ischar (self.beams(ind).MidNodeOffset)
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList (self.beams(ind).MidNodeOffset), ...
                                               4, ...
                                               true, ...
                                               sprintf ('beam %d mid_node_offset', ind) );
                else
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList ('reference', self.beams(ind).MidNodeOffsetReference, self.beams(ind).MidNodeOffset), ...
                                               4, ...
                                               true, ...
                                               sprintf ('beam %d mid_node_offset', ind) );
                                           

                    
                end
                
                if ~isempty (self.beams.MidNodeOrientation)
                    
                    if ind == 1 || ~ischar (self.beams(ind).MidNodeOrientation)
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', 'reference', self.beams(ind).MidNodeOrientationReference, self.beams(ind).MidNodeOrientation), ...
                                                   4, ...
                                                   true, ...
                                                   sprintf ('beam %d mid_node_orientation', ind) );
                    else
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', self.beams(ind).MidNodeOrientation), ...
                                                   4, ...
                                                   true, ...
                                                   sprintf ('beam %d mid_node_orientation', ind) );

                    end
                    
                end
                
                addcomma = ~isempty (self.beams(ind).EndNodeOrientation) ...
                            || ind < numel (self.beams) ...
                            || ~isempty (self.initialBeam) ...
                            || ~isempty (self.initialNode) ...
                            || ~isempty (self.smearingFactor) ;
                        
                if ischar (self.beams(ind).EndNodeOffset)
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList (self.beams(ind).EndNodeOffset), ...
                                               4, ...
                                               addcomma, ...
                                               sprintf ('beam %d end_node_offset', ind) );
                else
                    str = self.addOutputLine ( str, ...
                                               self.commaSepList ('reference', self.beams(ind).EndNodeOffsetReference, self.beams(ind).EndNodeOffset), ...
                                               4, ...
                                               addcomma, ...
                                               sprintf ('beam %d end_node_offset', ind) );
                end
                
                
                
                if ~isempty (self.beams.EndNodeOrientation)
                    
                    addcomma = ind < numel (self.beams) ...
                            || ~isempty (self.initialBeam) ...
                            || ~isempty (self.initialNode) ...
                            || ~isempty (self.smearingFactor) ;
                        
                    
                    if ind == 1 || ~ischar (self.beams(ind).EndNodeOrientation)
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', 'reference', self.beams(ind).EndNodeOrientationReference, self.beams(ind).EndNodeOrientation), ...
                                                   4, ...
                                                   addcomma, ...
                                                   sprintf ('beam %d end_node_orientation', ind) );
                    else
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ('hinge', self.beams(ind).EndNodeOrientation), ...
                                                   4, ...
                                                   addcomma, ...
                                                   sprintf ('beam %d end_node_orientation', ind) );

                    end
                    
                end
                

            end

            if ~isempty (self.initialBeam)
                
                addcomma =  ~isempty (self.initialNode) ...
                            || ~isempty (self.smearingFactor) ;

                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('initial beam', self.initialBeam), ...
                                           2, ...
                                           addcomma, ...
                                           sprintf ('initial beam') );
            end

            if ~isempty (self.initialNode)
                
                addcomma = ~isempty (self.smearingFactor) ;
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('initial node', self.initialNode), ...
                                           2, ...
                                           addcomma, ...
                                           sprintf ('initial node') );

            end

            if ~isempty (self.smearingFactor)
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('smearing', self.smearingFactor), ...
                                           2, ...
                                           false, ...
                                           sprintf ('smearing factor') );

            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));

        end
        
    end
    
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.singleNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.SliderType = [];
            options.RelativeOrientation = [];
            options.OffsetReference = 'node';
            options.OrientationReference = 'node';
            options.InitialBeam = [];
            options.InitialNode = [];
            options.SmearingFactor = [];
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
        function beam = makeSliderBeamStruct ()
            
            beam = struct ( 'Beam', [], ...
                            'FirstNodeOffset', [], ...
                            'MidNodeOffset', [], ...
                            'EndNodeOffset', [], ...
                            'FirstNodeOffsetReference', [], ...
                            'MidNodeOffsetReference', [], ...
                            'EndNodeOffsetReference', [], ...
                            'FirstNodeOrientation', [], ...
                            'MidNodeOrientation', [], ...
                            'EndNodeOrientation', [], ...
                            'FirstNodeOrientationReference', [], ...
                            'MidNodeOrientationReference', [], ...
                            'EndNodeOrientationReference', [] );
            
        end
        
    end
    
end