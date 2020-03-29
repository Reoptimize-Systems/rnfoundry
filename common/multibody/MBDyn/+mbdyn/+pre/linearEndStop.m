classdef linearEndStop < mbdyn.pre.base
% implements a an end stop for linear motion
%
% Syntax
%
%
% Description
%
% mbdyn.pre.linearEndStop creates a linearEndStop object. The linearEndStop applies
% forces which 
%      _
%     |  0               ,  (dl + x_p) <= 0
% F = |
%     |_ k * (dl + x_p)  ,  (dl + x_p) > 0
%
% Where k is a spring constant and x_p is a prestrain which can
% be be adjusted using the 'Prestrain' option (or implicitly
% via the 'Prestress' option).
%
% In other words, the tether acts like an ideal spring when
% extended beyond its initial length, and provides zero force
% when shorter, i.e. it goes slack. The default spring constant
% is 1e3, and can be adjusted using the 'SpringConstant'
% option.
%
% This class does not directly correspond to an MBDyn input
% file component, but instead is used to generate all the
% required components to create the tether using the
% generateMBDynSystemInputs method.
%
% mbdyn.pre.linearEndStop Methods:
%
%   linearEndStop - constructor for mbdyn.pre.linearEndStop object
%   generateMBDynSystemInputs - returns structure containing elements,
%     drives and variables for tether
%
%
% See Also: 
%
    
    properties
        
        node1;
        node2;
        startPoint;
        forceAtExt;
        ext;
        name;
        
    end
    
    
    methods
        
        function self = linearEndStop (node1, node2, start_point, forceatext, ext, varargin)
            % constructor for mbdyn.pre.tether object
            %
            % Syntax
            %
            % to = mbdyn.pre.tether (node1, node2)
            % to = mbdyn.pre.tether (..., 'Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.tether creates a tether object. The tether applies
            % forces which act like a simple weightless rope. The initial
            % length of the tether is determined from the intial position
            % of the two nodes. At each time step the difference between
            % this initial length and the current length, dl, is
            % determined. Then, a force is applied with the following
            % function:
            %      _
            %     |  0               ,  (dl + x_p) <= 0
            % F = |
            %     |_ k * (dl + x_p)  ,  (dl + x_p) > 0
            %
            % Where k is a spring constant and x_p is a prestrain which can
            % be be adjusted using the 'Prestrain' option (or implicitly
            % via the 'Prestress' option). The direction of the force in a
            % line pointing between the two nodes to which the tether is
            % attached, with equal and opposite forces applied to each node
            % (this is actually implemented using an internal structural
            % force, using the mbdyn.pre.structuralInternalForce class).
            %
            % In other words, the tether acts like an ideal spring when
            % extended beyond its initial length, and provides zero force
            % when shorter, i.e. it goes slack. The default spring constant
            % is 1e3, and can be adjusted using the 'SpringConstant'
            % option.
            %
            % This class does not directly correspond to an MBDyn input
            % file component, but instead is used to generate all the
            % required components to create the tether using the
            % generateMBDynSystemInputs method.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode object
            %
            %  node2 - mbdyn.pre.structuralNode object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Name' - 
            %
            % Output
            %
            %  to - mbdyn.pre.tether object
            %
            %
            %
            % See Also: 
            %
            
            options.Name = '';
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.base ();
            
            if isempty (options.Name)
                options.Name = sprintf ('linearEndStop_%d', self.uid);
            end
            
            mbdyn.pre.base.checkNumericScalar (start_point, true, 'start_point');
            assert (start_point > 0, 'start_point must be > 0');
            
            mbdyn.pre.base.checkNumericScalar (forceatext, true, 'forceatext');
            assert (forceatext > 0, 'forceatext must be > 0');
            
            mbdyn.pre.base.checkNumericScalar (ext, true, 'ext');
            assert (ext > 0, 'ext must be > 0');        
            
            self.node1 = node1;
            self.node2 = node2;
            self.startPoint = start_point;
            self.forceAtExt = forceatext;
            self.ext = ext;
            self.name = options.Name;
            
        end
        
        function sysinputs = generateMBDynSystemInputs (self)
            % returns structure containing elements, drives and variables for linearEndStop
            %
            % Syntax
            %
            % sysinputs = generateMBDynSystemInputs (leo)
            %
            % Description
            %
            % mndyn.pre.linearEndStop.generateMBDynSystemInputs returns a
            % structure containing the elements, drives and variables
            % required to create the linearEndStop in MBDyn. These objects
            % can be added to the mbdyn.pre.system object
            % which will construct the system.
            %
            % Input
            %
            %  leo - mbdyn.pre.linearEndStop object
            %
            % Output
            %
            %  sysinputs - structure containing the following fields:
            %   
            %   Elements : mbdyn.pre.element objects required for the 
            %    linearEndStop
            %
            %   Drives : mbdyn.pre.drive objects required for the 
            %    linearEndStop
            %
            %   Variables : mbdyn.pre.variable objects required for the 
            %    linearEndStop
            %
            %  
            %
            %
            
            
            coeff = self.forceAtExt / (self.ext^2);
            
            npoints = 15;
            x2 = self.startPoint + self.ext * ( exp ( linspace (log (1), log (10), npoints) ) - 1 );
%             x2 = linspace (start_point, 10*start_point, 10);
            
            force2 = coeff * (x2 - self.startPoint).^2;
            
            nzeros = 5;
            
            x = [linspace(0, self.startPoint, nzeros), x2(2:end)];
            
            x = [fliplr(-x(2:end)), x];
            
            force = [zeros(1, nzeros-1), force2 ];
            
            force = [fliplr(-force(2:end)), force];
            
%             figure; plot (x, force);

            sf = mbdyn.pre.cubicSplineScalarFunction ( sprintf ('%s_%d', self.name, self.uid), x, force );

            law = mbdyn.pre.scalarFunctionElasticIsotropicConstituativeLaw (sf);
            
            F_end_stop = mbdyn.pre.deformableDisplacementJoint (self.node1, self.node2, law, 'null', 'null', 'Offset1Reference', 'node', 'Offset2Reference', 'other node' );
            
            sysinputs.Elements = { F_end_stop };
            sysinputs.Drives = {};
            sysinputs.Variables = {};
            
        end
        
    end    
    
end