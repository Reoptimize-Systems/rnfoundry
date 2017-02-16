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
    
    methods
        
        function self = system ()
            
            self.problems = { struct('tyep', 'initial value') };
            
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
        
        function [filename, str] = generateMBDynInputFile (self, filename)
            
            
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
            
            %% problems section
            % write out each problem section
            for ind = 1:numel (self.problems)
                
                str = self.addOutputLine (str , sprintf('begin: %s;', self.problems{ind}.type), 0, false);
                
                str = self.addOutputLine (str , sprintf('end: %s;', self.problems{ind}.type), 0, false);
                
            end
            
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

            str = self.addOutputLine (str , 'end: control data;', 0, false);
            
            %% drivers section
            if numel (self.drivers) > 0
                str = self.addOutputLine (str , 'begin: drivers;', 0, false);
                for ind = 1:numel (self.drivers)
                    str = sprintf ('%s\n', self.drivers{ind}.generateOutputString ());
                end
                str = self.addOutputLine (str , 'end: drivers;', 0, false);
            end
            
            %% nodes section
            str = self.addOutputLine (str , 'begin: nodes;', 0, false);
            for ind = 1:numel (self.nodes)
                str = sprintf ('%s\n', self.nodes{ind}.generateOutputString ());
            end
            str = self.addOutputLine (str , 'end: nodes;', 0, false);
            
            %% elements section
            str = self.addOutputLine (str , 'begin: elements;', 0, false);
            
            str = sprintf ('%s\n', str);
            
            for ind = 1:numel (self.elements)
                str = sprintf ('%s\n', self.elements{ind}.generateOutputString ());
            end

            str = self.addOutputLine (str , 'end: elements;', 0, false);
            
            %% write out the string to a file
            if nargin < 2
                filename = [tempname, '.mbd'];
            end
            
            fid = fopen (filename);
            
            CC = onCleanup (@() fclose(fid));
            
            fprintf (fid, str);
            
        end
        
        function elcount = countControlElements (self)
            
            elcount.StructuralNodes = 0;
            elcount.RigidBodies = 0;
            elcount.Joints = 0;
            elcount.Forces = 0;
            elcount.FileDrivers = 0;
            
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