classdef system < mbdyn.pre.base
    % class representing an mbdyn system
    
    properties
        
%        references;
       data;
       problems;
       controlData;
       nodes;
       drivers;
       elements;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        drawAxesH; % handle to figure for plotting
    end
    
    methods
        
        function self = system (problems, varargin)
            
            options.Nodes = {};
            options.Elements = {};
            options.Drivers = {};
            
            options = parse_pv_pairs (options, varargin);
            
            self.problems = {};
            self.nodes = {};
            self.elements = {};
            self.drivers = {};
            
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
            
            self.checkCellArrayClass (elements, 'mbdyn.pre.element');
            
            self.elements = [self.elements, elements];
            
        end
        
        function addDrivers (self, drivers)
            
            drivers = self.makeCellIfNot (drivers);
            
            self.checkCellArrayClass (drivers, 'mbdyn.pre.driver');
            
            self.drivers = [self.drivers, drivers];
            
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
            
            if isa (self.drawAxesH, 'matlab.graphics.axis.Axes')
                if ~isvalid (self.drawAxesH)
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
            
            options = parse_pv_pairs (options, varargin);
            
            % make figure and axes if necessary
            if isempty (options.AxesHandle)
                figure;
                self.drawAxesH = axes;
            else
                self.drawAxesH = options.AxesHandle;
            end
            
%             hold all
            
            for ind = 1:numel (self.nodes)
                if isa (self.nodes{ind}, 'mbdyn.pre.structuralNode') && options.StructuralNodes
                    draw (self.nodes{ind}, ...
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
                end
            end
%             hold off
            
            if options.Light
                camHandle = findobj(gcf,'Type','light');
                if isempty (camHandle)
                    light (self.drawAxesH);
                end
            end
            
%             axis equal;
            xlabel ('x'); ylabel ('y'); zlabel('z'); 
            view (3);
            axis equal
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
            
             % make sure labels are set
            self.setLabels ();
            
%             begin: data;
%                 problem: initial value; # the default
%             end: data;
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
            
            if elcount.RigidBodies > 0
                str = self.addOutputLine (str , sprintf('rigid bodies: %d;', elcount.RigidBodies), 1, false);
            end
            
            if elcount.Joints > 0
                str = self.addOutputLine (str , sprintf('joints: %d;', elcount.Joints), 1, false);
            end
            
            if elcount.Forces > 0
                str = self.addOutputLine (str , sprintf('forces: %d;', elcount.Forces), 1, false);
            end
            
            if elcount.Gravity
                str = self.addOutputLine (str , 'gravity;', 1, false);
            end
            
            str = self.addOutputLine (str, 'default orientation: orientation matrix;', 1, false);

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
            
            str = self.generateMBDynInputStr ();
            
            % write out the string to a file
            if nargin < 2
                filename = [tempname, '.mbd'];
            end
            
            [fid, errmsg] = fopen (filename, 'wt');
            
            CC = onCleanup (@() fclose(fid));
            
            fprintf (fid, '%s', str);
            
        end
        
        function elcount = countControlElements (self)
            
            elcount.StructuralNodes = 0;
            elcount.RigidBodies = 0;
            elcount.Joints = 0;
            elcount.Forces = 0;
            elcount.FileDrivers = 0;
            elcount.Gravity = false;
           
            % not yet implemented
             elcount.AbstractNodes = 0;
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
             elcount.Genels = 0;        
             elcount.HydraulicElements = 0;
             elcount.LoadableElements = 0;
             elcount.OutputElements = 0;
             elcount.InducedVelocityElements = 0;
            
            for ind = 1:numel (self.nodes)
                if isa (self.nodes{ind}, 'mbdyn.pre.structuralNode')
                    elcount.StructuralNodes = elcount.StructuralNodes + 1;
                end
            end
            
            for ind = 1:numel (self.elements)
                
                if isa (self.elements{ind}, 'mbdyn.pre.joint')
                    elcount.Joints = elcount.Joints + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.force')
                    elcount.Forces = elcount.Forces + 1;
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
        
    end
    
    methods (Access = protected)
        
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