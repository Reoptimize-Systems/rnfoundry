classdef system < mbdyn.pre.base
% class representing an mbdyn multibody dynamics system
%
% Syntax
%
% sys = mbdyn.pre.system (problems)
% sys = mbdyn.pre.system (..., 'Parameter', Value)
% 
%
% Description
%
% This class represents a multibody system and associated problem data.
% It's primary use is to generate an input file for MBDyn. This class
% depends on the many other classes in the mbdyn.pre namespace which
% represent the various things which make up an mbdyn problem description,
% such as nodes, forces, drives, bodies etc.
%
% The class can also be used to create a visual representaiton of the
% system, and node positons and orientations can be set, allowing
% simulations to be visualised at any time point.
%
% The mbdyn.pre.system object can also assist with cosimulation (updating
% node postitions etc. in the system) and with postprocessing.
%
% mbdyn.pre.system Methods:
%
%   system - constructor, see constructor help for more detail on the main
%     class options
%   addDriveCallers - add drives to the system 
%   addDrivers - adds drivers to the system (like file drivers)
%   addElements - adds elements to the system (bodies, joints, forces, 
%     drives etc.)
%   addNodes - adds nodes to the system
%   addProblems - adds problems to the system
%   addVariables - adds variables to the system
%   addPluginVariables - 
%   addScalarFunctions - 
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
%
%
% About MBDyn
% -----------
%
% MBDyn is the first and possibly the only free general purpose Multibody
% Dynamics analysis software, released under GNU's GPL 2.1.
% 
% It has been developed at the Dipartimento di Scienze e Tecnologie
% Aerospaziali (formerly Dipartimento di Ingegneria Aerospaziale) of the
% University "Politecnico di Milano", Italy.
% 
% MBDyn features the integrated multidisciplinary simulation of multibody,
% multiphysics systems, including nonlinear mechanics of rigid and flexible
% bodies (geometrically exact & composite-ready beam and shell finite
% elements, component mode synthesis elements, lumped elements) subjected
% to kinematic constraints, along with smart materials, electric networks,
% active control, hydraulic networks, and essential fixed-wing and
% rotorcraft aerodynamics.
%
%
% See Also: mbdyn.mint.MBCNodal, mbdyn.postproc
%
    
    properties
        
       references;
       data;
       problems;
       controlData;
       nodes;
       driveCallers;
       drivers;
       elements;
       variables;
       pluginVariables;
       scalarFunctions;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
%         drawAxesH; % handle to figure for plotting
    end
    
    methods
        
        function self = system (problems, varargin)
            % mbdyn.pre.system constructor
            %
            % Syntax
            %
            % sys = mbdyn.pre.system (problems)
            % sys = mbdyn.pre.system (..., 'Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.system act as a container for all the components of
            % an MBDyn mutibody dynamics and multiphysics problem. It is
            % used to visualise the complete problem, generate MBDyn input
            % files and assist with cosimulation and pstprocessing.
            %
            % Input
            %
            %  problems - mbdyn.pre.problem object (or derived object such
            %   as mbdyn.pre.initialValueProblem), or a cell array of these
            %   objects. 
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Nodes' - cell array of mbdyn.pre.node objects (or derived
            %    objects, such as mbdyn.pre.structuralNode6dof).
            %
            %  'Elements' - cell array of mbdyn.pre.element objects (or
            %    derived objects, such as mbdyn.pre.body).
            %
            %  'DriveCallers' - cell array of mbdyn.pre.drive objects to add 
            %    to the MBDyn input file. These define drives which may be
            %    reused throughout the input file. See Section 2.4.3 'Drive
            %    Caller' of the MBDyn manual for further details.
            %
            %  'Drivers' - cell array mbdyn.pre.fileDriver objects to
            %    add to the drivers section of an MBDyn input file, see the
            %    section on Drivers in Chapter 7 of the MBDyn manual for
            %    further details.
            %
            %  'DefaultOutput' - cell array of strings indicating the
            %    default output flag for a type of node or element. It can
            %    be overridden for each entity either when it is created or
            %    later, for entity aggregates, in each entity module, by
            %    means of the output directive for nodes and elements.
            %
            %    Special values are: 
            %
            %      'all' : enables output of all entities
            %
            %      'none' : disables output of all entities; 
            %
            %      'reference frames' : by default, reference frames are
            %        not output. When enabled, a special file .rfm is
            %        generated, which contains all the reference frames
            %        defined, using the syntax of the .mov file.
            %
            %      'accelerations' enables output of linear and angular
            %        accelerations for dynamic structural nodes, which are
            %        disabled by default. Accelerations output can be
            %        enabled on a node by node basis; see
            %        mbdyn.pre.structuralNode (or derived classes) for
            %        details.
            %
            %    Other possible values are:
            %
            %      'abstract nodes'
            %      'electric nodes'
            %      'hydraulic nodes'
            %      'structural nodes'
            %      'aerodynamic elements'
            %      'air properties'
            %      'beams'
            %      'electric elements'
            %      'forces'
            %      'genels'
            %      'gravity'
            %      'hydraulic elements'
            %      'joints'
            %      'rigid bodies'
            %      'induced velocity elements'
            %
            %    Values are additive, except for 'none', so to select only
            %    specific entities use 'none' first, followed by a list of
            %    the entities whose output should be activated.
            %
            %  'DefaultOrientation' - character vector which is the default
            %    orientation output format. Can be one of 'euler123',
            %    'euler313', 'euler321', 'orientation vector', or
            %    'orientation matrix'
            %
            %  'DefaultScales' - optional structure or array of structures
            %    used to set default scale values to apply to the residuals
            %    of different dof owners before applying the tolerance
            %    test. The structure should contain only the fields
            %    'ScaleItem' and 'ScaleFactor'. Scaleitem is a a character
            %    vector, which can be one of the following: 'all', 'none',
            %    'abstract nodes', 'electric nodes', 'hydraulic nodes',
            %    'structural nodes' or 'thermal nodes'. The ScaleFactor
            %    field must contain a scalar value with is the scale factor
            %    to apply to the corresponding item in the ScaleItem field.
            %
            %  'OutputResults' - this option can be used to enable the
            %    output of data in the NetCDF format. It can be a character
            %    vector, or a cell array of character vectors. If a
            %    character vector, it must  be the keyword 'netcdf', which
            %    enable netcdf output in addition to the normal text
            %    output. If a cell array, it must be of length 1, 2 or
            %    three elements. The first cell must contain the keyword
            %    'netcdf'. The second cell can contain the keywords 'sync'
            %    or 'no text'. The 'no text' means only the netcdf file
            %    will be produced, otherwise both the standard text format
            %    output files are produced as well as the netcdf format
            %    files. The sync option is currently undocumented. If the
            %    second cell contains 'no text' no third cell may be
            %    present, if the second cell contains 'sync', the third
            %    cell may contain 'no text'. So in summary, the possible
            %    values for OutputResults are:
            %
            %    'netcdf'
            %    { 'netcdf' }
            %    { 'netcdf', 'no text' }
            %    { 'netcdf', 'sync' }
            %    { 'netcdf', 'sync', 'no text' }
            %
            %    By default OutputResults is { 'netcdf', 'no text' },
            %    meaning only the netcdf format output is produced by
            %    MBDyn.
            %
            %  'References' - optional cell array of mbdyn.pre.reference
            %    objects. If supplied, these are used purely for
            %    visualisation and have no effect on the simulation in any
            %    way.
            %
            % Output
            %
            %  sys - mbdyn.pre.system object
            %
            %
            %
            % See Also: mbdyn.min.MBCNodal, mbdyn.postproc
            %

            options.Nodes = {};
            options.Elements = {};
            options.Variables = {};
            options.PluginVariables = {};
            options.DriveCallers = {};
            options.ScalarFunctions = {};
            options.DefaultOutput = {};
            options.DefaultOrientation = '';
            options.References = {};
            options.DefaultScales = [];
            options.Drivers = {};
            options.OutputResults = { 'netcdf', 'no text' };
            options.UseInAssembly = {};
            options.AssemblyTolerance = [];
            options.AssemblyMaxIterations = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.problems = {};
            self.nodes = {};
            self.elements = {};
            self.driveCallers = {};
            self.references = {};
            self.variables = {};
            self.drivers = {};
            
            self.addProblems (problems);
            
            if ~isempty (options.Nodes)
                self.addNodes (options.Nodes);
            end
            
            if ~isempty (options.Elements)
                self.addElements (options.Elements);
            end
            
            if ~isempty (options.Variables)
                self.addVariables (options.Variables);
            end
            
            if ~isempty (options.PluginVariables)
                self.addPluginVariables (options.PluginVariables);
            end
            
            if ~isempty (options.DriveCallers)
                self.addDriveCallers (options.DriveCallers);
            end
            
            if ~isempty (options.ScalarFunctions)
                self.addScalarFunctions (options.ScalarFunctions);
            end
            
            if ~isempty (options.References)
                self.addReferences (options.References);
            end
            
            if ~isempty (options.Drivers)
                self.addDrivers (options.Drivers);
            end
            
            if ~isempty (options.DefaultOutput)
                if ~iscellstr (options.DefaultOutput)
                    error ('DefaultOutput must be a cell array of strings if supplied');
                else
                    okstrings = { 'all', ...
                                  'none', ...
                                  'reference frames', ...
                                  'abstract nodes', ...
                                  'electric nodes', ...
                                  'hydraulic nodes', ...
                                  'structural nodes', ...
                                  'accelerations', ...
                                  'aerodynamic elements', ...
                                  'air properties', ...
                                  'beams', ...
                                  'electric elements', ...
                                  'forces', ...
                                  'genels', ...
                                  'gravity', ...
                                  'hydraulic elements', ...
                                  'joints', ...
                                  'rigid bodies', ...
                                  'induced velocity elements', ...
                                };
                            
                    for ind = 1:numel (options.DefaultOutput)
                        self.checkAllowedStringInputs ( options.DefaultOutput{ind}, ...
                                                        okstrings, ...
                                                        true, ...
                                                        sprintf ('DefaultOutput{%d}', ind) );
                    end
                end
            end
            
            if ~isempty (options.DefaultOrientation)
                self.checkAllowedStringInputs (options.DefaultOrientation, ...
                    {'euler123', 'euler313', 'euler321', 'orientation vector', 'orientation matrix'}, ...
                    true, 'DefaultOrientation');
            end
            
            if ~isempty (options.DefaultScales)
                
                assert (isstruct (options.DefaultScales), ...
                    'DefaultScales should be an aray of one or more structures.');
                
                if ~all (isfield (options.DefaultScales, {'ScaleItem', 'ScaleFactor'})) ...
                        || ( numel(fieldnames (options.DefaultScales)) > 2)
                    error ('DefaultScales should be an aray of one or more structures with only the fields ''ScaleItem'' and ''ScaleFactor''.');
                end
                
                for ind = 1:numel (options.DefaultScales)
                    
                    self.checkAllowedStringInputs ( ...
                        options.DefaultScales(ind).ScaleItem, ...
                        { 'all', ...
                          'none', ...
                          'abstract nodes', ...
                          'electric nodes', ...
                          'hydraulic nodes', ...
                          'structural nodes', ...
                          'thermal nodes' ...
                        }, ...
                        true, ...
                        sprintf ('DefaultScales(%d).ScaleItem', ind) ...
                        );
                        
                    self.checkNumericScalar ( ...
                        options.DefaultScales(ind).ScaleFactor, ...
                        true, ...
                        sprintf ('DefaultScales(%d).ScaleFactor', ind) );

                end
                
            end
            
            if ~isempty (options.OutputResults)
                
                if ischar (options.OutputResults)
                    
                    ok = self.checkAllowedStringInputs (options.OutputResults, {'netcdf'}, false);
                    
                    assert (ok, 'If OutputResults is a character vector, it must be ''netcdf''');
                    
                    % convert to cell array
                    options.OutputResults = {options.OutputResults};
                    
                elseif ~iscellstr (options.OutputResults)
                    
                    error ('OutputResults must be a character vector (''netcdf'') or a cell array of character vectors of length 1, 2 or 3');
                    
                end
                
                self.checkAllowedStringInputs (options.OutputResults{1}, {'netcdf'}, true, 'OutputResults{1}');
                if numel (options.OutputResults) > 1
                    self.checkAllowedStringInputs (options.OutputResults{2}, {'sync', 'no text'}, true, 'OutputResults{2}');
                end
                if numel (options.OutputResults) > 2
                    if strcmp (options.OutputResults{2}, 'no text')
                        error ('If the second cell in OutputResults is ''no text'', there can be no third cell (if using ''sync'' it must come before ''no text'')');
                    end
                    self.checkAllowedStringInputs (options.OutputResults{3}, {'no text'}, true, 'OutputResults{3}');
                end

            end
            
            if ~isempty (options.UseInAssembly)
                
                assert (iscellstr (options.UseInAssembly), ...
                        'UseInAssembly must be a cell array of character vectors' );
                
                valid_use_in_assembly = { 'rigid bodies', ...
                                          'gravity', ...
                                          'forces', ...
                                          'beams', ...
                                          'aerodynamic elements', ...
                                          'loadable elements' };
                
                for ind = 1:numel (options.UseInAssembly)
                    self.checkAllowedStringInputs ( options.UseInAssembly{ind}, ...
                              valid_use_in_assembly, ...
                              true, ...
                              sprintf ('UseInAssembly{%d}', ind) );
                end

            end
            
            if ~isempty (options.AssemblyTolerance)
                self.checkNumericScalar (options.AssemblyTolerance, true, 'AssemblyTolerance');
            end
            
            if ~isempty (options.AssemblyMaxIterations)
                self.checkScalarInteger (options.AssemblyMaxIterations, true, 'AssemblyMaxIterations');
            end

            self.controlData.DefaultOutput = options.DefaultOutput;
            self.controlData.DefaultOrientation = options.DefaultOrientation;
            self.controlData.DefaultScales = options.DefaultScales;
            self.controlData.OutputResults = options.OutputResults;
            self.controlData.UseInAssembly = options.UseInAssembly;
            self.controlData.AssemblyTolerance = options.AssemblyTolerance;
            self.controlData.AssemblyMaxIterations = options.AssemblyMaxIterations;
            
        end
        
        function addProblems (self, problems)
            
            problems = self.makeCellIfNot (problems);
            
            % ensure it's a row vector
            problems = reshape (problems, 1, []);
            
            self.checkCellArrayClass (problems, 'mbdyn.pre.problem');
            
            self.problems = [self.problems, problems];
            
            self.problems = self.uniqueCells (self.problems);
            
        end
        
        function addNodes (self, nodes)
            
            nodes = self.makeCellIfNot (nodes);
            
            % ensure it's a row vector
            nodes = reshape (nodes, 1, []);
            
            % remove empty
            nodes(cellfun('isempty',nodes)) = [];
            
            self.checkCellArrayClass (nodes, 'mbdyn.pre.node');
            
            self.nodes = [self.nodes, nodes];
            
            self.nodes = self.uniqueCells (self.nodes);
            
        end
        
        function addElements (self, elements)
            
            elements = self.makeCellIfNot (elements);
            
            % ensure it's a row vector
            elements = reshape (elements, 1, []);
            
            % remove empty
            elements(cellfun('isempty',elements)) = [];
            
            self.checkCellArrayClass (elements, 'mbdyn.pre.element');
            
            self.elements = [self.elements, elements];
            
            self.elements = self.uniqueCells (self.elements);
            
        end
        
        function addVariables (self, variables)
            
            variables = self.makeCellIfNot (variables);
            
            % ensure it's a row vector
            variables = reshape (variables, 1, []);
            
            % remove empty
            variables(cellfun('isempty',variables)) = [];
            
            self.checkCellArrayClass (variables, 'mbdyn.pre.variable');
            
            self.variables = [self.variables, variables];
            
        end
        
        function addPluginVariables (self, plugin_variables)
            
            plugin_variables = self.makeCellIfNot (plugin_variables);
            
            % ensure it's a row vector
            plugin_variables = reshape (plugin_variables, 1, []);
            
            % remove empty
            plugin_variables(cellfun('isempty',plugin_variables)) = [];
            
            self.checkCellArrayClass (plugin_variables, 'mbdyn.pre.pluginVariable');
            
            self.pluginVariables = [self.pluginVariables, plugin_variables];
            
        end
        
        function addDriveCallers (self, drives)
            
            drives = self.makeCellIfNot (drives);
            
            % ensure it's a row vector
            drives = reshape (drives, 1, []);
            
            % remove empty
            drives(cellfun('isempty', drives)) = [];
            
            self.checkCellArrayClass (drives, 'mbdyn.pre.drive');
            
            self.driveCallers = [self.driveCallers, drives];
            
            self.driveCallers = self.uniqueCells (self.driveCallers);
            
        end
        
        function addScalarFunctions (self, sf)
            
            sf = self.makeCellIfNot (sf);
            
            % ensure it's a row vector
            sf = reshape (sf, 1, []);
            
            % remove empty
            sf(cellfun('isempty', sf)) = [];

            self.checkCellArrayClass (sf, 'mbdyn.pre.scalarFunction');
            
            self.scalarFunctions = [self.scalarFunctions, sf];
            
            self.scalarFunctions = self.uniqueCells (self.scalarFunctions);
            
        end
        
        function addReferences (self, refs)
            
            refs = self.makeCellIfNot (refs);
            
            % ensure it's a row vector
            refs = reshape (refs, 1, []);
            
            % remove empty
            refs(cellfun('isempty',refs)) = [];
            
            self.checkCellArrayClass (refs, 'mbdyn.pre.reference');
            
            self.references = [self.references, refs];
            
            self.references = self.uniqueCells (self.references);
            
        end
        
        function setNodePosition (self, label, newposition)
            % set the absolute position of node with the given label
            %
            % Syntax
            %
            % setNodePosition (sys, label, newposition)
            %
            % Description
            %
            % setNodePosition is used to set the absolute position of one
            % of the nodes in the system in the gloabl frame. Typically
            % this is used for post-processing data from a simulation to
            % create a visualisation of the system from the positions
            % output by MBDyn.
            %
            % Input
            %
            %  sys - mbdyn.pre.system object
            %
            %  label - the label of the node for which the new absolute
            %   position is to be set
            %
            %  newposition - (3 x 1) vector containing the new position of
            %   the node in the global coordinate system.
            %
            %
            % See Also: mbdyn.pre.system.setNodeOrientation
            %

            self.checkScalarInteger (label, true, 'label');
            
            % the nodes check that newposition is valid, so no need to
            % check here
            for ind = 1:numel (self.nodes)
               
                if self.nodes{ind}.label == label
                    self.nodes{ind}.absolutePosition = newposition;
                    return;
                end
                
            end
            
        end
        
        function setNodeOrientation (self, label, neworientation)
            % set the absolute orientation of node with the given label
            %
            % Syntax
            %
            % setNodeOrientation (sys, label, newposition)
            %
            % Description
            %
            % setNodePosition is used to set the absolute orientation of
            % one of the nodes in the system in the gloabl frame. Typically
            % this is used for post-processing data from a simulation to
            % create a visualisation of the system from the orientations
            % output by MBDyn.
            %
            % Input
            %
            %  sys - mbdyn.pre.system object
            %
            %  label - the label of the node for which the new absolute
            %   orientation is to be set
            %
            %  neworientation - (3 x 1) vector containing the new 
            %   orientation of the node in the global coordinate system.
            %
            %
            % See Also: mbdyn.pre.system.setNodePosition
            %
            
            self.checkScalarInteger (label, true, 'label');
            
            % the nodes check that neworientation is valid, so no need to
            % check here
            for ind = 1:numel (self.nodes)
               
                if self.nodes{ind}.label == label
                    self.nodes{ind}.setAbsoluteOrientation ( neworientation );
                    return;
                end
                
            end
            
        end
        
        function comminfo = externalStructuralCommInfo (self)
            % gets communication info for an external structural force.
            %
            % Syntax
            %
            % comminfo = externalStructuralCommInfo (sys)
            %
            % Description
            %
            % externalStructuralCommInfo gets information about the
            % communication method used by an external structural force
            % (mbdyn.pre.structuralExternal object). Throws an error if
            % there is no external structural force in the system.
            %
            % Input
            %
            %  sys - mbdyn.pre.system object
            %
            % Output
            %
            %  comminfo - structure containing information about the
            %   communication method used by the external structural force.
            %   It always contains the field 'commMethod' which contains a
            %   string describing the type of communation method used by
            %   the external structural force. Example strings include:
            %   'inet socket', 'local socket' or 'shared memory'. The other
            %   fields in the structure (if any) are dependent on the type
            %   of communication. To generate this structure,
            %   externalStructuralCommInfo simply calles the commInfo
            %   method of the external structural force object's
            %   communicator. So to see the expected contents of the
            %   sturcture for the different communicator types, see the
            %   help for the commInfo method for that class.
            %
            %
            % See Also: mbdyn.pre.system.externalStructuralInfo
            %
            
            for ind = 1:numel (self.elements)
                if isa (self.elements{ind}, 'mbdyn.pre.externalStructuralForce')
                    
                    el = self.elements{ind};
                    
                    comminfo = el.communicator.commInfo ();
                    
                    return;
                end
            end
            
            error ('No external structural force was in the system.');
            
        end
        
        function extforceinfo = externalStructuralInfo (self)
            % gets information on any external structural force in the system
            %
            % Syntax
            %
            % comminfo = externalStructuralInfo (sys)
            %
            % Description
            %
            % externalStructuralInfo gets information about the an external
            % structural force (mbdyn.pre.structuralExternal object)
            % present in the system. Throws an error if there is no
            % external structural force in the system.
            %
            % Input
            %
            %  sys - mbdyn.pre.system object
            %
            % Output
            %
            %  extforceinfo - structure containing information about the
            %   external structural force. Contains the following fields:
            %
            %   UseLabels : true/false flag indicating whether MBDyn will
            %    use the node labels when communicating.
            %
            %   UseAccelerations : true/false flag indicating whether MBDyn
            %    will return acceleration data when communicating.
            %
            %   NodeOrientationType : string containing the node
            %    orientation description, i.e. one of 'none',
            %    'orientation matrix', 'orientation vector', or 'euler
            %    123'
            %
            %   NNodes : The number of nodes referenced by the external
            %    structural force
            %
            %   Nodes : An array of mbdyn.pre.structuralNode objects
            %    representing the referenced nodes in the external
            %    structural force
            %
            %
            % See Also: mbdyn.pre.system.externalStructuralInfo
            %
            
            testyesno = @(x) strcmp (x, 'yes');
            
            for ind = 1:numel (self.elements)
                if isa (self.elements{ind}, 'mbdyn.pre.externalStructuralForce')
                    
                    extf = self.elements{ind};
                    
                    if isempty (extf.labels) 
                        extforceinfo.UseLabels = false;
                    else
                        extforceinfo.UseLabels = testyesno (extf.labels);
                    end
                    
                    if isempty (extf.accelerations) 
                        extforceinfo.UseAccelerations = false;
                    else
                        extforceinfo.UseAccelerations = testyesno (extf.accelerations);
                    end
                    
                    if isempty (extf.useReferenceNodeForces) 
                        extforceinfo.UseRefNode = false;
                    else
                        extforceinfo.UseRefNode = testyesno (extf.useReferenceNodeForces);
                    end
                    
                    if isempty (extf.orientation) 
                        extforceinfo.NodeOrientationType = false;
                    else
                        extforceinfo.NodeOrientationType = extf.orientation;
                    end
                    
                    extforceinfo.NNodes = numel (extf.nodes);
                    
                    extforceinfo.Nodes = extf.nodes;
                    
                    return;
                end
            end
            
            error ('No external structural force was in the system.');
            
        end
        
        function [hax, hfig] = draw (self, varargin)
            % draw the mbdyn system
            %
            % Syntax
            %
            % [hax, hfig] = draw (mbsys)
            %
            % Description
            %
            % mbdyn.pre.system.draw draws the complete mbdyn system in a
            % figure.
            %
            % Input
            %
            %  mbsys - mbdyn.pre.system object
            %
            % Additional optional arguments may be supplied using parameter
            %-value pairs. The available options are:
            %
            %
            %  'AxesHandle' - optional handle to axes in which to plot the
            %    system. If not supplied, a new figure and axes will be
            %    created. The first time the system is drawn the handle to
            %    the axes will be stored internally and future calls to
            %    draw will plot to the same axes, unless the axes are
            %    destroyed, or this option is used to override it.
            %
            %  'ForceRedraw' - true/false flag indicating whether to force
            %    a full redraw of the system (rather than just update the
            %    transform matrices of the elements), even if the system
            %    itself does not think a redraw is required.
            %
            %  'Bodies' - provides methods of indicating which bodies to 
            %    draw, can be one of:
            %
            %    1. a logical true/false flag determining whether to draw
            %       the system bodies. If true, all bodies are drawn. If
            %       false, no bodies are drawn.
            %    
            %    2. A vector of body indices, i.e. the number of the body
            %       in the order in which all bodies were added to the
            %       system (starting from 1). Only the bodies with these
            %       indices will be drawn.
            %
            %    3. A character vector containing the name of the body to
            %       be drawn. The first body with the same name in its
            %       'name' property will be drawn
            %
            %    4. A cell array of character vectors, each with a name of
            %       a body, as described in 3 above. only the bodies with
            %       these names will be drawn. For each cell the first
            %       body with that name will be drawn.
            %
            %    Default value of the 'Bodies' option is logical true,
            %    so all bodies will be drawn.
            %
            %  'StructuralNodes' - provides methods of indicating which 
            %    structural nodes to draw, can be one of:
            %
            %    1. a logical true/false flag determining whether to draw
            %       the system structural nodes. If true, all structural
            %       nodes are drawn. If false, no structural nodes are
            %       drawn.
            %    
            %    2. A vector of structural nodes indices, i.e. the number
            %       of the structural nodes in the order in which all
            %       bodies were added to the system (starting from 1). Only
            %       the structural nodes with these indices will be drawn.
            %
            %  'Joints' - provides methods of indicating which joints to 
            %    draw, can be one of:
            %
            %    1. a logical true/false flag determining whether to draw
            %       the system joints. If true, all joints are drawn. If
            %       false, no joints are drawn.
            %    
            %    2. A vector of joint indices, i.e. the number of the joint
            %       in the order in which all joints were added to the
            %       system (starting from 1). Only the joints with these
            %       indices will be drawn.
            %
            %    3. A character vector containing the name of the joint to
            %       be drawn. The first joint with the same name in its
            %       'name' property will be drawn
            %
            %    4. A cell array of character vectors, each with a name of
            %       a joint, as described in 3 above. only the joints with
            %       these names will be drawn. For each cell the first
            %       joint with that name will be drawn.
            %
            %    Default value of the 'Joints' option is logical true,
            %    so all joints will be drawn.
            %
            %  'Mode' - character vector determining the style in which the
            %    elements will be plotted. Can be one of 'solid',
            %    'wiresolid', 'ghost', 'wireframe', 'wireghost'. Default is
            %    'solid'. This option mainly effects only bodies, and
            %    joints with only a default shape.
            %
            %  'Light' - deterined whether the scene should have light
            %    source
            %
            %  'AxLims' - (3 x 2) matrix containing new axis limites for
            %    the X,Y and Z axes. Each row has the new limits for the
            %    X,Y and Z axes respectively. Default is empty matrix
            %    meaning axis limites are not modified from the automatic
            %    values.
            %
            %  'References' - logical flag determining whether to draw the
            %    references (if there are any present in the system object).
            %
            %  'ReferenceScale' - scalar value. The length of the arrows
            %    representing the axes of the plotted reference objects
            %    will be scaled by this value. Default is 1 if not
            %    supplied.
            %
            %  'GlobalReference' - logical flag determining whether to draw 
            %    the global reference, when the 'References' option is set
            %    to true. Default is true.
            %
            % Output
            %
            %  hax - axes handle for axes in which the system is drawn
            %
            %  hfig - figure handle for figure in which the system is drawn
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
            options.NodeLabels = false;
            options.Joints = true;
            options.Light = false;
            options.AxLims = [];
            options.References = false;
            options.ReferenceScale = 1;
            options.GlobalReference = true;
            options.View = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.ForceRedraw, true, 'ForceRedraw');
            self.checkLogicalScalar (options.Light, true, 'Light');
            self.checkLogicalScalar (options.References, true, 'References');
            self.checkNumericScalar (options.ReferenceScale, true, 'ReferenceScale');
            self.checkLogicalScalar (options.GlobalReference, true, 'GlobalReference');
            
            elcount = self.countControlElements ();

            options.Bodies = processDrawInds ( self, ...
                                               options.Bodies, ...
                                               'body', ...
                                               'Bodies', ...
                                               1:elcount.RigidBodies );
            
            options.Joints = processDrawInds ( self, ...
                                               options.Joints, ...
                                               'joint', ...
                                               'Joints', ...
                                               1:elcount.Joints );
            
            % make figure and axes if necessary
            self.checkAxes (options.AxesHandle);
            hax = self.drawAxesH;
            hfig = self.drawFigureH;
            
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

%             elseif iscell (options.StructuralNodes)
%                 
%                 drawinds = zeros (1,numel (options.StructuralNodes));
%                 
%                 for ind = 1:numel (options.StructuralNodes)
%                     
%                     if isa (options.StructuralNodes{ind}, 'mbdyn.pre.structuralNode')
%                         
%                         
%                         
%                     elseif self.checkValidScalarIndex (options.StructuralNodes{ind}, false)
%                         
%                     elseif ischar (options.StructuralNodes{ind})
%                         
%                     end
%                     
%                 end
                
            else
                assert (isnumeric (options.StructuralNodes) && isvector (options.StructuralNodes), ...
                    'error StructuralNodes must be a logical true/false flag, or a vector of integers indicating which nodes to plot');
            end

            for ind = 1:numel (options.StructuralNodes)
                if isa (self.nodes{options.StructuralNodes(ind)}, 'mbdyn.pre.structuralNode')
                    draw (self.nodes{options.StructuralNodes(ind)}, ...
                        'AxesHandle', self.drawAxesH, ...
                        'ForceRedraw', options.ForceRedraw, ...
                        'Label', options.NodeLabels);
                end
            end
            
            joint_ind = 1;
            joint_draw_ind = 1;
            body_ind = 1;
            body_draw_ind = 1;
            for ind = 1:numel (self.elements)
                
                if isa (self.elements{ind}, 'mbdyn.pre.body') && ~isempty (options.Bodies)
                    
                    if body_draw_ind <= numel (options.Bodies) ...
                            && body_ind == options.Bodies(body_draw_ind)
                        
                        draw (self.elements{ind}, ...
                            'AxesHandle', self.drawAxesH, ...
                            'ForceRedraw', options.ForceRedraw, ...
                            'Mode', options.Mode );
                        
                        body_draw_ind = body_draw_ind + 1;
                        
                    end
                    
                    body_ind = body_ind + 1;
                    
                elseif isa (self.elements{ind}, 'mbdyn.pre.joint') && ~isempty (options.Joints)
                    
                    if joint_draw_ind <= numel (options.Joints) ...
                            && joint_ind == options.Joints(joint_draw_ind)
                        
                        draw (self.elements{ind}, ...
                            'AxesHandle', self.drawAxesH, ...
                            'ForceRedraw', options.ForceRedraw, ...
                            'Mode', options.Mode );
                        
                        joint_draw_ind = joint_draw_ind + 1;
                        
                    end
                    
                    joint_ind = joint_ind + 1;
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
            
            if isempty (options.View)
                view (self.drawAxesH, 3);
            else
                view (self.drawAxesH, options.View);
            end
            
        end
        
        function setStructuralNodeSize (self, sx, sy, sz)
            % set the size of all nodes for drawing
            %
            % Syntax
            %
            % mbdyn.pre.system.setStructuralNodeSize (s)
            % mbdyn.pre.system.setStructuralNodeSize (sx, sy, sz)
            %
            % Description
            %
            % mbdyn.pre.system.setStructuralNodeSize sets the length of the
            % lines used to plot structural nodes. Nodes are drawn as three
            % lines representing the local coordinate axes of the nodes.
            % The length of all three axis lines can be set to the same
            % value, or they can be set individually.
            %
            % Input
            %
            %  s - if a single input is supplied, 's', this is the size of
            %   the nodes representations in all dimentions 
            %
            %  sx - if three values are supplied, sx is the length of the 
            %   node line representing the node's axis 1 (the X axis)
            %
            %  sy - if three values are supplied, sy is the length of the 
            %   node line representing the node's axis 2 (the Y axis)
            %
            %  sz - if three values are supplied, sz is the length of the 
            %   node line representing the node's axis 3 (the Z axis)
            %
            
            if nargin == 2
                sy = sx;
                sz = sx;
            elseif nargin < 4
                error ('Input should be either a single value, ''s'' or three values, ''sx'', ''sy'' and ''sz''');
            end
            
            for ind = 1:numel (self.nodes)
                if isa (self.nodes{ind}, 'mbdyn.pre.structuralNode')
                    setSize (self.nodes{ind}, sx, sy, sz);
                end
            end
        end
        
        function el = getElementByLabel (self, label)
            
            el = [];
            for ind = 1:numel (self.elements)
                if self.elements{ind}.label == label
                    el = self.elements{ind};
                end
            end
            
        end
        
        function node = getNodeByLabel (self, label)
            
            node = [];
            for ind = 1:numel (self.nodes)
                if self.nodes{ind}.label == label
                    node = self.nodes{ind};
                end
            end
            
        end
        
        function out = getComponentByLabel (self, label)
            
            out = getElementByLabel (self, label);
            
            if isempty (out)
                out = getNodeByLabel (self, label);
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
                str = sprintf ('%s\n%s\n', str, self.problems{ind}.generateMBDynInputString ());
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
            
            if elcount.AddedMasses > 0
                str = self.addOutputLine (str , sprintf('added masses: %d;', elcount.AddedMasses), 1, false);
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
            
            if elcount.LoadableElements > 0
                str = self.addOutputLine (str , sprintf('loadable elements: %d;', elcount.LoadableElements), 1, false);
            end
            
            if elcount.Drivers > 0
                str = self.addOutputLine (str , sprintf('file drivers: %d;', elcount.Drivers), 1, false);
            end
            
            if elcount.Gravity
                str = self.addOutputLine (str , 'gravity;', 1, false);
            end
            
            if ~isempty (self.controlData.DefaultOutput)
                str = self.addOutputLine ( str, ...
                                           sprintf ( 'default output: %s;', ...
                                                     self.commaSepList (self.controlData.DefaultOutput{:}) ), ...
                                           1, ...
                                           false );
            end
            
            if ~isempty (self.controlData.DefaultOrientation)
                str = self.addOutputLine ( str, sprintf ('default orientation: %s;', self.controlData.DefaultOrientation), 1, false);
            end
            
            if ~isempty (self.controlData.DefaultScales)
                
                str = self.addOutputLine (str, 'default scale: ', 1, true);
                
                for ind = 1:numel (self.controlData.DefaultScales)
                    
                    addcomma = ind < numel (self.controlData.DefaultScales);
                    
                    str = self.addOutputLine ( str, ...
                                    self.commaSepList ( self.controlData.DefaultScales.ScaleItem, ...
                                                        self.controlData.DefaultScales.ScaleFactor ), ...
                                               2, ...
                                               addcomma );
                    
                end
                
                str = self.addOutputLine ( str, ';', 1, false);
                
            end
            
            if ~isempty (self.controlData.OutputResults)
                
                str = self.addOutputLine ( str, ...
                                           [ 'output results: ', ...
                                             self.commaSepList(self.controlData.OutputResults{:}), ...
                                             ';' ], ...
                                           1, ...
                                           false );
                                       
            end
            
            if ~isempty (self.controlData.UseInAssembly)
                str = self.addOutputLine ( str, ...
                                           [ 'use : ', ...
                                             self.commaSepList(self.controlData.UseInAssembly{:}), ...
                                             ', in assembly ;' ], ...
                                           1, ...
                                           false );
            end
            
            if ~isempty (self.controlData.AssemblyTolerance)
                str = self.addOutputLine (str , sprintf('tolerance: %s;',self.formatNumber (self.controlData.AssemblyTolerance)), 1, false);
            end
            
            if ~isempty (self.controlData.AssemblyMaxIterations)
                str = self.addOutputLine (str , sprintf('max iterations: %s;', self.formatInteger (self.controlData.AssemblyMaxIterations)), 1, false);
            end

            str = self.addOutputLine (str , 'end: control data;', 0, false);
            str = sprintf ('%s\n', str);
            
            %% Loadable modules
            if elcount.LoadableElements > 0
                % gather up all the loadable modules from the user defined
                % elements and generate the module-load commands
                loadable_mods = mbdyn.pre.loadableModule.empty();
                for ind = 1:numel (self.elements)

                    if isa (self.elements{ind}, 'mbdyn.pre.userDefined')
                        loadable_mods(end+1) = self.elements{ind}.loadableModule;
                    end
                    
                end
                
                loadable_mods = unique (loadable_mods);
                
                for ind = 1:numel (loadable_mods)
                    
                    str = sprintf ('%s\n%s\n', str, loadable_mods(ind).generateMBDynInputString ());
                    
                end
                
                str = sprintf ('%s\n', str);
                
            end
           

            %% nodes section
            str = self.addOutputLine (str , 'begin: nodes;', 0, false);
            str = sprintf ('%s\n', str);
            
            for ind = 1:numel (self.nodes)
                str = sprintf ('%s\n%s\n', str, self.nodes{ind}.generateMBDynInputString ());
            end
            str = self.addOutputLine (str , 'end: nodes;', 0, false);
            
            str = sprintf ('%s\n', str);
            
            %% plugin variables
            if numel (self.pluginVariables) > 0

                for ind = 1:numel (self.pluginVariables)
                    str = sprintf ('%s\n%s\n', str, self.pluginVariables{ind}.generateMBDynInputString ());
                end
                
                str = sprintf ('%s\n', str);
            end
            
            %% scalar functions
            if numel (self.scalarFunctions ) > 0

                for ind = 1:numel (self.scalarFunctions)
                    str = sprintf ('%s\nscalar function : %s ;\n', str, self.scalarFunctions{ind}.generateMBDynInputString ());
                end
                
                str = sprintf ('%s\n', str);
            end

            %% variables
            if numel (self.variables) > 0

                for ind = 1:numel (self.variables)
                    str = sprintf ('%s\n%s\n', str, self.variables{ind}.generateMBDynInputString ());
                end
                
                str = sprintf ('%s\n', str);
            end
            
            %% file drivers section
            if numel (self.drivers) > 0
                str = self.addOutputLine (str , 'begin: drivers;', 0, false);
                for ind = 1:numel (self.drivers)
                    str = sprintf ('%s\n%s\n', str, self.drivers{ind}.generateMBDynInputString ());
                end
                str = self.addOutputLine (str , 'end: drivers;', 0, false);
                str = sprintf ('%s\n', str);
            end
            
            % drives (not file drivers)
            if numel (self.driveCallers) > 0
                str = sprintf ('%s\n#DriveCallers\n\n', str);
                for ind = 1:numel (self.driveCallers)
                    str = sprintf ('%s\ndrive caller : %d, %s;\n', str, self.driveCallers{ind}.label, self.driveCallers{ind}.generateMBDynInputString ());
                end
                str = sprintf ('%s\n', str);
            end
            
            %% elements section
            
            str = self.addOutputLine (str , 'begin: elements;', 0, false);
            
            str = sprintf ('%s\n', str);
            
            for ind = 1:numel (self.elements)
                str = sprintf ('%s\n%s\n', str, self.elements{ind}.generateMBDynInputString ());
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
            elcount.Drivers = numel (self.drivers);
            elcount.Gravity = false;
            elcount.Genels = 0;
            elcount.AbstractNodes = 0;
            elcount.AddedMasses = 0;
            
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
                
                if isa (self.elements{ind}, 'mbdyn.pre.addedMassAndInertia')
                    elcount.AddedMasses = elcount.AddedMasses + 1;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.gravity')
                    if elcount.Gravity
                        error ('Gravity present more than once');
                    end
                    elcount.Gravity = true;
                end
                
                if isa (self.elements{ind}, 'mbdyn.pre.userDefined')
                    elcount.LoadableElements = elcount.LoadableElements + 1;
                end
                
            end
            
        end
        
        function ok = checkCellArrayClass (self, CC, classname, throw)
            
            if nargin < 4
                throw = true;
            end
            
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
            
            for ind = 1:numel (self.driveCallers)
                self.driveCallers{ind}.label = label;
                label = label + 1;
            end
            
            for ind = 1:numel (self.elements)
                self.elements{ind}.label = label;
                label = label + 1;
            end
            
        end
        
        function elinds = processDrawInds (self, elspec, eltype, optname, allelinds)
            
            if ischar (elspec)
                elspec = {elspec};
            end
            
            if islogical (elspec)
                
                self.checkLogicalScalar (elspec, true, optname);
                if elspec
                    elinds = allelinds;
                else
                    elinds = [];
                end
                
            elseif isnumeric (elspec) && isvector (elspec)
                
                elinds = elspec;
                
            elseif iscellstr (elspec)
                
                elinds = [];
                
                for joptind = 1:numel (elspec)
                    % loop through the cell array of joint names. For each
                    % name, check all joints to see if the name matches
                    foundjointname = false;
                    
                    thisjointind = 0;
                    
                    for elind = 1:numel (self.elements)
                        
                        if isa (self.elements{elind}, sprintf ('mbdyn.pre.%s', eltype))
                            
                            thisjointind = thisjointind + 1;
                            fprintf (1, 'joint name: %s\n', self.elements{elind}.name)
                            if strcmp (self.elements{elind}.name, elspec{joptind})
                                
                                % note that we matched up the joint name
                                % this time
                                foundjointname = true;
                                
                                % add the joint index to the list of joints
                                % to plot
                                elinds = [elinds, allelinds(thisjointind)];
                                
                                % exit the inner element loop
                                break;
                                
                            end  
                            
                        end
                        
                    end
                    
                    if foundjointname == false
                        error ('Joint name %s was not found', elspec{joptind});
                    end
                    
                end
            else
                error ('Joints must be a logical scalar value or a vector of joint indices to plot');
            end
            
            if isvector (elinds)
                
                for ind = 1:numel (elinds)
                    self.checkScalarInteger (elinds(ind), true, sprintf ('%s(%d)', optname, ind));
                    assert (elinds(ind) >=1, sprintf ('%s(%d) is < 1', optname, ind));
                    assert (elinds(ind) <= numel (allelinds), sprintf ('%s(%d) is < 1', optname, ind));
                end
                
                % get the unique entries
                elinds = unique(elinds);
                
            end
                
            elinds = sort (elinds);
            
        end
        
    end
    
end