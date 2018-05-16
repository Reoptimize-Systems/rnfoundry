classdef baseSystem < mbdyn.pre.base
    % class representing an mbdyn multibody dynamics system
    %
    % Syntax
    %
    % 
    %
    % Description
    %
    % This class represents a multibody system and associated problem data
    % and can be used to generate an input file for MBDyn. This class
    % depends on the many other classes in the mbdyn.pre namespace which
    % represent the various things which make up an mbdyn problem
    % description, such as nodes, forces, drives, bodies etc.
    %
    % The class can also be used to create a visual representaiton of the
    % system, and node positons and orientations can be set, allowing
    % simulations to be visualised at any time point.
    %
    % system Methods:
    %   addDrivers - adds drives to the system
    %   addElements - adds elements to the system (bodies, joints, forces, 
    %     drives etc.)
    %   addNodes - adds nodes to the system
    %   addProblems - adds problems to the system
    %   draw - draws the system in a figure window
    %   externalStructuralCommInfo - returns information about the
    %     communication method used in any external structural force element
    %   externalStructuralInfo - returns information about any external
    %     structural force element in the system
    %   generateMBDynInputFile - gerates an mbdyn input file for the system
    %   generateMBDynInputStr - gerates a string representing an mbdyn
    %     input file contents
    %   setNodeOrientation - sets a node's orientation (useful for drawing)
    %   setNodePosition - sets a node's position (useful for drawing)
    %   setStructuralNodeSize - sets the size of all nodes for drawing
    %   system - constructor
    %
    % About MBDyn
    % -----------
    %
    % MBDyn is the first and possibly the only free general purpose
    % Multibody Dynamics analysis software, released under GNU's GPL 2.1
    % (get a cached copy here).
    % 
    % It has been developed at the Dipartimento di Scienze e Tecnologie
    % Aerospaziali (formerly Dipartimento di Ingegneria Aerospaziale) of
    % the University "Politecnico di Milano", Italy.
    % 
    % MBDyn features the integrated multidisciplinary simulation of
    % multibody, multiphysics systems, including nonlinear mechanics of
    % rigid and flexible bodies (geometrically exact & composite-ready beam
    % and shell finite elements, component mode synthesis elements, lumped
    % elements) subjected to kinematic constraints, along with smart
    % materials, electric networks, active control, hydraulic networks, and
    % essential fixed-wing and rotorcraft aerodynamics.
    %
    
    
    properties
        
       references;
       data;
       nodes;
       drivers;
       elements;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        drawAxesH; % handle to figure for plotting
        hasSubsystems;
    end
    
    methods
        
        function self = baseSystem (problems, varargin)
            
            options.Nodes = {};
            options.Elements = {};
            options.Drivers = {};
            options.References = {};
            
            options = parse_pv_pairs (options, varargin);
            
            self.hasSubsystems = false;
            self.nodes = {};
            self.elements = {};
            self.drivers = {};
            self.references = {};
            
            self.addProblems (problems);
            
            if ~isempty (options.Nodes)
                self.addNodes (options.Nodes)
            end
            
            if ~isempty (options.Elements)
                self.addElements (options.Elements)
            end
            
            if ~isempty (options.Drivers)
                self.addDrivers (options.Drivers)
            end
            
            if ~isempty (options.References)
                self.addReferences (options.References)
            end           
            
        end
        
        function addProblems (self, problems)
            
            problems = self.makeCellIfNot (problems);
            
            self.checkCellArrayClass (problems, 'mbdyn.pre.problem');
            
            self.problems = [self.problems, problems];
            
        end
        
        function addNodes (self, nodes)
            
            nodes = self.makeCellIfNot (nodes);
            
            self.checkCellArrayClass (nodes, 'mbdyn.pre.node');
            
            self.nodes = [self.nodes, nodes];
            
        end
        
        function addElements (self, elements)
            
            elements = self.makeCellIfNot (elements);
            
%             for ind = 1:numel (elements)
%                 if isa (elements{ind}, 'mbdyn.pre.externalStructuralForce')
%                     error ();
%                 end
%             end
            
            self.checkCellArrayClass (elements, 'mbdyn.pre.element');
            
            self.elements = [self.elements, elements];
            
            for ind = 1:numel (elements)
                if isa (elements{ind}, 'mbdyn.pre.baseSystem')
                    self.hasSubsystems = true;
                end
            end
            
        end
        
        function addDrivers (self, drivers)
            
            drivers = self.makeCellIfNot (drivers);
            
            self.checkCellArrayClass (drivers, 'mbdyn.pre.driver');
            
            self.drivers = [self.drivers, drivers];
            
        end
        
        function addReferences (self, refs)
            
            refs = self.makeCellIfNot (refs);
            
            self.checkCellArrayClass (refs, 'mbdyn.pre.reference');
            
            self.references = [self.references, refs];
            
        end
        
        function setNodePosition (self, label, newposition)
            
            for ind = 1:numel (self.nodes)
               
                if self.nodes{ind}.label == label
                    self.nodes{ind}.absolutePosition = newposition;
                    return;
                end
                
            end
            
        end
        
        function setNodeOrientation (self, label, neworientation)
            
            for ind = 1:numel (self.nodes)
               
                if self.nodes{ind}.label == label
                    self.nodes{ind}.absoluteOrientation = neworientation;
                    return;
                end
                
            end
            
        end
        
        
        function draw (self, varargin)
            % draw the mbdyn system
            %
            % Syntax
            %
            %
            % Description
            %
            %
            %
            % Input
            %
            %  mbsys - mbdyn.pre.system object
            %
            % Additional optional arguments may be supplied using parameter
            %-value pairs. The available options are:
            %
            % 'AxesHandle' - self.drawAxesH
            %
            % 'ForceRedraw' - Default is false
            %
            % 'Mode' - Default is 'solid'
            %
            % 'Bodies' - Default is true
            %
            % 'StructuralNodes' - Default is true
            %
            % 'Joints' - Default is true
            %
            % 'Light' - Default is false
            %
            % 'AxLims' - (3 x 2) matrix containing new axis limites for the
            %   X,Y and Z axes. Each row has the new limits for the X,Y and
            %   Z axes respectively. Default is empty matrix meaning axis
            %   limites are not modified from the automatic values.
            %
            % 'References' - logical flag determining whether to draw the
            %   references (if there are any present in the system object).
            %
            % Output
            %
            %
            %
            % See Also: 
            %
            
            if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
                if ~ishghandle (self.drawAxesH)
                    self.drawAxesH = [];
                end
            end
            
            options.AxesHandle = self.drawAxesH;
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Bodies = true;
            options.StructuralNodes = true;
            options.Joints = true;
            options.Light = false;
            options.AxLims = [];
            options.References = false;
            options.ReferenceScale = 1;
            options.GlobalReference = true;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.ForceRedraw, true, 'ForceRedraw');
            self.checkLogicalScalar (options.Bodies, true, 'Bodies');
            self.checkLogicalScalar (options.Joints, true, 'Joints');
            self.checkLogicalScalar (options.Light, true, 'Light');
            self.checkLogicalScalar (options.References, true, 'References');
            self.checkNumericScalar (options.ReferenceScale, true, 'ReferenceScale');
            self.checkLogicalScalar (options.GlobalReference, true, 'GlobalReference');
            
            % make figure and axes if necessary
            if isempty (options.AxesHandle)
                figure;
                self.drawAxesH = axes;
            else
                self.drawAxesH = options.AxesHandle;
            end
            
            if options.References
                if ~isempty (self.references)
                    self.drawReferences ( self.references, ...
                                          'Scale', options.ReferenceScale, ...
                                          'Title', false, ...
                                          'PlotAxes', self.drawAxesH, ...
                                          'DrawGlobal', options.GlobalReference );
                end
            end
            
            if islogical (options.StructuralNodes)
                if options.StructuralNodes == false
                    options.StructuralNodes = [];
                else
                    options.StructuralNodes = 1:numel(self.nodes);
                end
            else
                assert (isnumeric (options.StructuralNodes) && isvector (options.StructuralNodes), ...
                    'error StructuralNodes must be a logical true/false flag, or a vector of integers indicating which nodes to plot');
            end

            for ind = 1:numel (options.StructuralNodes)
                if isa (self.nodes{options.StructuralNodes(ind)}, 'mbdyn.pre.structuralNode')
                    draw (self.nodes{options.StructuralNodes(ind)}, ...
                        'AxesHandle', self.drawAxesH, ...
                        'ForceRedraw', options.ForceRedraw);
                end
            end
            
            for ind = 1:numel (self.elements)
                if isa (self.elements{ind}, 'mbdyn.pre.body') && options.Bodies
                    draw (self.elements{ind}, ...
                        'AxesHandle', self.drawAxesH, ...
                        'ForceRedraw', options.ForceRedraw, ...
                        'Mode', options.Mode );
                elseif isa (self.elements{ind}, 'mbdyn.pre.joint') && options.Joints
                    draw (self.elements{ind}, ...
                        'AxesHandle', self.drawAxesH, ...
                        'ForceRedraw', options.ForceRedraw, ...
                        'Mode', options.Mode );
                elseif isa (self.elements{ind}, 'mbdyn.pre.baseSystem')
                    % recursively draw subsystems
                    draw (self.elements{ind}, ...
                        'AxesHandle', self.drawAxesH, ...
                        'ForceRedraw', options.ForceRedraw, ...
                        'Mode', options.Mode, ...
                        'Bodies', options.Bodies, ...
                        'StructuralNodes', options.StructuralNodes , ...
                        'Joints', options.Joints , ...
                        'Light', options.Light , ...
                        'AxLims', options.AxLims , ...
                        'References', options.References , ...
                        'ReferenceScale', options.ReferenceScale, ...
                        'GlobalReference', false );
                end
            end
            
            if options.Light
                camHandle = findobj(gcf,'Type','light');
                if isempty (camHandle)
                    light (self.drawAxesH);
                end
            end
            
%             axis equal;
            xlabel ('x'); ylabel ('y'); zlabel('z'); 
            
            if ~isempty (options.AxLims)
                set (self.drawAxesH, 'XLim', options.AxLims (1,1:2));
                set (self.drawAxesH, 'YLim', options.AxLims (2,1:2));
                set (self.drawAxesH, 'Zlim', options.AxLims (3,1:2));
                axis (self.drawAxesH, 'equal')
            else
                axis (self.drawAxesH, 'equal')
            end
            
            view (3);
            
        end
        
        function setStructuralNodeSize (self, sx, sy, sz)
            % set the size of all nodes for drawing
            
            for ind = 1:numel (self.nodes)
                if isa (self.nodes{ind}, 'mbdyn.pre.structuralNode')
                    setSize (self.nodes{ind}, sx, sy, sz);
                end
            end
        end
        
        function str = generateMBDynInputStr (self)
            % creates a string representing the contents of an MBDyn input file
            %
            % Syntax
            %
            % str = generateMBDynInputStr (mbs)
            %
            % Input
            %
            %  mbs - mbdyn.pre.system object
            %
            % Output
            %
            %  str - the contents of the MBDyn input representing the
            %    system file as a string
            %
            
            % make sure labels are set
            self.setLabels ();
            
            str = '';
            
            %% data section
            str = self.addOutputLine (str , 'begin: data;', 0, false);
            for ind = 1:numel (self.problems)
                str = self.addOutputLine (str , sprintf('problem: %s;', self.problems{ind}.type), 1, false);
            end
            str = self.addOutputLine (str , 'end: data;', 0, false);
            str = sprintf ('%s\n', str);
            
            %% problems section
            % write out each problem section
            for ind = 1:numel (self.problems)
                str = sprintf ('%s\n%s\n', str, self.problems{ind}.generateOutputString ());
            end
            str = sprintf ('%s\n', str);
            
            elcount = self.countControlElements ();
            
            %% control data section
            str = self.addOutputLine (str , 'begin: control data;', 0, false);
            
            if elcount.StructuralNodes > 0
                str = self.addOutputLine (str , sprintf('structural nodes: %d;', elcount.StructuralNodes), 1, false);
            end
            
            if elcount.AbstractNodes > 0
                str = self.addOutputLine (str , sprintf('abstract nodes: %d;', elcount.AbstractNodes), 1, false);
            end
            
            if elcount.RigidBodies > 0
                str = self.addOutputLine (str , sprintf('rigid bodies: %d;', elcount.RigidBodies), 1, false);
            end
            
            if elcount.Joints > 0
                str = self.addOutputLine (str , sprintf('joints: %d;', elcount.Joints), 1, false);
            end
            
            if elcount.Forces > 0
                str = self.addOutputLine (str , sprintf('forces: %d;', elcount.Forces), 1, false);
            end
            
            if elcount.Genels > 0
                str = self.addOutputLine (str , sprintf('genels: %d;', elcount.Genels), 1, false);
            end
            
            if elcount.Gravity
                str = self.addOutputLine (str , 'gravity;', 1, false);
            end
            
            if ~isempty (self.controlData.DefaultOutput)
                str = self.addOutputLine (str, sprintf ('default output: %s;', self.commaSepList (self.controlData.DefaultOutput{:})), 1, false);
            end
            
            if ~isempty (self.controlData.DefaultOrientation)
                str = self.addOutputLine (str, sprintf ('default orientation: %s;', self.controlData.DefaultOrientation), 1, false);
            end

            str = self.addOutputLine (str , 'end: control data;', 0, false);
            str = sprintf ('%s\n', str);
            
            %% drivers section
            if numel (self.drivers) > 0
                str = self.addOutputLine (str , 'begin: drivers;', 0, false);
                for ind = 1:numel (self.drivers)
                    str = sprintf ('%s\n', self.drivers{ind}.generateOutputString ());
                end
                str = self.addOutputLine (str , 'end: drivers;', 0, false);
                str = sprintf ('%s\n', str);
            end
            
            %% nodes section
            str = self.addOutputLine (str , 'begin: nodes;', 0, false);
            str = sprintf ('%s\n', str);
            
            for ind = 1:numel (self.nodes)
                str = sprintf ('%s\n%s\n', str, self.nodes{ind}.generateOutputString ());
            end
            str = self.addOutputLine (str , 'end: nodes;', 0, false);
            
            str = sprintf ('%s\n', str);
            
            %% elements section
            
            str = self.addOutputLine (str , 'begin: elements;', 0, false);
            
            str = sprintf ('%s\n', str);
            
            for ind = 1:numel (self.elements)
                str = sprintf ('%s\n%s\n', str, self.elements{ind}.generateOutputString ());
            end

            str = self.addOutputLine (str , 'end: elements;', 0, false);
            
        end
        
        
        function [filename, str] = generateMBDynInputFile (self, filename)
            % creates an input file for mbdyn based on the system
            %
            % Syntax
            %
            % [filename, str] = generateMBDynInputFile (mbs, filename)
            %
            % Input
            %
            %  mbs - mbdyn.pre.system object
            %
            %  filename - (optional) string containing path where mbdyn input file
            %    will be written. If not provided, the file will be created
            %    in an automatically generated temporary location.
            %
            % Output
            %
            %  filename - path to generated MBDyn input file
            %
            %  str - the contents of the file as a string
            %
            
            str = self.generateMBDynInputStr ();
            
            % write out the string to a file
            if nargin < 2
                filename = [tempname, '.mbd'];
            end
            
            [fid, errmsg] = fopen (filename, 'wt');
            
            CC = onCleanup (@() fclose(fid));
            
            fprintf (fid, '%s', str);
            
        end
        
    end
    
    methods (Access = protected)
        
        function elcount = countControlElements (self)
            
            elcount.StructuralNodes = 0;
            elcount.RigidBodies = 0;
            elcount.Joints = 0;
            elcount.Forces = 0;
            elcount.FileDrivers = 0;
            elcount.Gravity = false;
            elcount.Genels = 0;
            elcount.AbstractNodes = 0;
            
            % not yet implemented
            elcount.ElectricNodes = 0;
            elcount.HydraulicNodes = 0;
            elcount.ParameterNodes = 0;
            elcount.ThermalNodes = 0;
            elcount.AerodynamicElements = 0;
            elcount.Aeromodals = 0;
            elcount.AirProperties = 0;
            elcount.AutomaticStructuralElements = 0;
            elcount.Beams = 0;
            elcount.BulkElements = 0;
            elcount.ElectricBulkElements = 0;
            elcount.ElectricElements = 0;
            elcount.ExternalElements = 0;
            
            elcount.HydraulicElements = 0;
            elcount.LoadableElements = 0;
            elcount.OutputElements = 0;
            elcount.InducedVelocityElements = 0;
            
            for ind = 1:numel (self.nodes)
                
                if isa (self.nodes{ind}, 'mbdyn.pre.structuralNode')
                    elcount.StructuralNodes = elcount.StructuralNodes + 1;
                end
                
                if isa (self.nodes{ind}, 'mbdyn.pre.abstractNode')
                    elcount.AbstractNodes = elcount.AbstractNodes + 1;
                end
                
            end
            
            for ind = 1:numel (self.elements)
                
                if isa (self.elements{ind}, 'mbdyn.pre.joint')
                    elcount.Joints = elcount.Joints + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.force') ...
                        || isa (self.elements{ind}, 'mbdyn.pre.couple')
                    elcount.Forces = elcount.Forces + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.genel')
                    elcount.Genels = elcount.Genels + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.body')
                    elcount.RigidBodies = elcount.RigidBodies + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.gravity')
                    if elcount.Gravity
                        error ('Gravity present more than once');
                    end
                    elcount.Gravity = true;
                end
                
            end
            
        end
        
        function ok = checkCellArrayClass (self, CC, classname, throw)
            
            ok = true;
            if iscell (CC)
                for ind = 1:numel (CC)
                    if isa (CC{ind}, classname)
                        
                    else
                        ok = false;
                        if throw
                            error ('input cell array member %d was not of the class type %s', ind, classname);
                        end
                    end
                end
            elseif isa (CC, classname)
                % do nothing
                ok = true;
            else
                ok = false;
                if throw
                    error ('input was not of the class type %s', classname);
                end
            end
            
        end
        
        function setLabels (self)
            
            label = 1;
            
            for ind = 1:numel (self.nodes)
                self.nodes{ind}.label = label;
                label = label + 1;
            end
            
            for ind = 1:numel (self.drivers)
                self.drivers{ind}.label = label;
                label = label + 1;
            end
            
            for ind = 1:numel (self.elements)
                self.elements{ind}.label = label;
                label = label + 1;
            end
            
        end
        
    end
    
end