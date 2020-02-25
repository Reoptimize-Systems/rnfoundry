classdef pluginVariable < mbdyn.pre.base
% represents an MBDyn input file variable
%
% Syntax
%
% variable (varType, varName)
% variable (..., 'Parameter', Value)
%
% Description
%
% represents a variable in an mbdyn input file. Variable are
% values or expressions which are evaluated before a simulation
% starts and may be used in other places, such as string drive
% expressions during the simulation. Their values are not
% updated or changed as the simulation progresses, only once at
% startup. See the mbdyn manual for more information on
% variables.
%
%
% mbdyn.pre.variable Methods:
%
%   variable - mbdyn.pre.variable constructor
%   generateMBDynInputString - generate an mbdyn input file string for the element
%
%
% See Also: mbdyn.pre.system
%

    properties (GetAccess = public, SetAccess = protected)
        
        component; % element or node from which data is to be extracted into a variable
        componentType; % character vector containing the component type 'element' or 'node'
        varName; % name of the variable
        dataName; % name of the data to be 'plugged in' to the variable
        
    end

    methods
        
        function self = pluginVariable (component, varName, dataName)
            % mbdyn.pre.variable constructor
            %
            % Syntax
            %
            % pv = pluginVariable (component, varName, dataName)
            %
            % Description
            %
            % Plugin variables are a way of binding a variable to some
            % means to dynamically generate its value. Typically it is a
            % way to access some data about a node or element dynamically
            % during simulation. There are built-in plugins
            % that allow to link variables to the value of degrees of
            % freedom and to special parameters computed by nodes and
            % elements.
            %
            % Input
            %
            %  component - mbdyn.pre.node object or mbdyn.pre.element
            %   object
            %
            %  varName - character vector containing the name of the 
            %   variable that will be available in the input file which
            %   will be return the value of the component data specified in
            %   dataName. This is the variable name which is reused
            %   elsewhere in the MBDyn input file.
            %
            %  dataName - character vector containing the name of the 
            %   component data to be accessed through the variable. This is
            %   the identifier for the data from the element of node to be
            %   accessed, e.g. for a node it might be 'X[3]' for the Z
            %   position of the node in the global frame.
            %
            % Output
            %
            %  pv - mbdyn.pre.pluginVariable object
            %
            % See Also: mbdyn.pre.system, mbdyn.pre.variable
            %
            
            self = self@mbdyn.pre.base ();
                   

            assert (ischar (varName), 'varName must be a character vector');
            assert (ischar (dataName), 'dataName must be a character vector');

            if isa (component, 'mbdyn.pre.element')
                self.componentType = 'element';
            elseif isa (component, 'mbdyn.pre.node')
                self.componentType = 'node';
            else
                error ('component must be an mbdyn.pre.element or mbdyn.pre.node object');
            end
            
            self.varName = varName;
            self.dataName = dataName;
            self.component = component;

        end
        
        function str = generateMBDynInputString (self)
            % generate an mbdyn input file string for a pluginVariable
            %
            % Syntax
            %
            % str = generateMBDynInputString (var)
            %
            % Description
            %
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %
            % Input
            %
            %  var - mbdyn.pre.pluginVariable
            %
            % Output
            %
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            if isa (self.component, 'mbdyn.pre.element')
                type = self.component.netCDFName;
            elseif isa (self.component, 'mbdyn.pre.node')
                type = self.component.type;
            end
            
            str = sprintf('set: [ %s, %s, %d, %s, string="%s" ];', ...
                           self.componentType, ...
                           self.varName, ...
                           self.component.label, ...
                           type, ...
                           self.dataName );
            
        end
        
    end

end
