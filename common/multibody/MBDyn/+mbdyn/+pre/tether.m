classdef tether < mbdyn.pre.base
% implements a simple tether via forces
%
% Syntax
%
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
% mbdyn.pre.tether Methods:
%
%   tether - constructor for mbdyn.pre.tether object
%   generateMBDynSystemInputs - returns structure containing elements,
%     drives and variables for tether
%
%
% See Also: 
%
    
    properties
        
        preStrain;
        node1;
        node2;
        springConstant;
        zeroForceDistance;
        name;
        
    end
    
    
    methods
        
        function self = tether (node1, node2, varargin)
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
            %  'SpringConstant' - 
            %
            %  'PreStrain' - optional value of the initial prestrain of the
            %    tether. This will result in there being an initial force
            %    on the tether which is already 'stretched' by the
            %    prestrain value at the start of the simulation. By default
            %    the prestrain is empty, meaning it will be set to zero.
            %    The PreStrain option is mutually exclusive with the
            %    PreStress and ZeroForceDistance options. If supplied, the
            %    PreStrain must be greater than zero.
            %
            %  'PreStress' - optional value of the initial prestress on the
            %    tether. This will result in there being an initial force
            %    on the tether which is already 'stretched' at the start of
            %    the simulation. The Prestress is internall implemented by
            %    calculating the prestrain required to achieve this force
            %    and setting the prestrain to this value. By default the
            %    prestress is empty, meaning it will be set to zero. The
            %    PreStress option is mutually exclusive with the PreStrain
            %    and ZeroForceDistance options. If supplied, the PreStress
            %    must be greater than zero.
            %
            %  'ZeroForceDistance' - an alternative to the prestrain (or
            %    prestress), this option allows setting the distance
            %    between the nodes at which the tether force becomes zero.
            %    This is otherwise calculated from the intial node
            %    separation distance, and allows this to be overridden.
            %
            % Output
            %
            %  to - mbdyn.pre.tether object
            %
            %
            %
            % See Also: 
            %
            
            options.PreStrain = [];
            options.PreStress = [];
            options.SpringConstant = 1e3;
            options.ZeroForceDistance = [];
            options.Name = '';
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.base ();
            
            if isempty (options.Name)
                options.Name = sprintf ('tether_%d', self.uid);
            end
            
            mbdyn.pre.base.checkNumericScalar (options.SpringConstant, true, 'SpringConstant');
            assert (options.SpringConstant > 0, 'SpringConstant must be > 0');
            
            if ~isempty (options.ZeroForceDistance)
                if ~(isempty (options.PreStrain) && isempty (options.PreStress))
                    error ('Setting ZeroForceDistance is not compatible with also setting PreStrain or PreStress');
                end
                mbdyn.pre.base.checkNumericScalar (options.ZeroForceDistance, true, 'ZeroForceDistance');
                assert (options.ZeroForceDistance > 0, 'ZeroForceDistance must be > 0');
            end
                  
            if isempty (options.PreStrain) && isempty (options.PreStress)
                
                self.preStrain = 0;
                
            elseif isempty (options.PreStrain)
                
                mbdyn.pre.base.checkNumericScalar (options.PreStress, true, 'PreStress');
                assert (options.PreStress > 0, 'PreStress must be > 0');
                
                % calculate the prestrain (initial extension of the tether)
                % from the stress and spring constant
                self.preStrain = options.PreStress / options.SpringConstant;
                
            elseif isempty (options.PreStress)
                
                mbdyn.pre.base.checkNumericScalar (options.PreStrain, true, 'PreStrain');
                assert (options.PreStrain > 0, 'PreStrain must be > 0');
                
                self.preStrain = options.PreStrain;
                
            else
                
                mbdyn.pre.base.checkNumericScalar (options.Pred_tension_calcStress, true, 'PreStress');
                assert (options.PreStress > 0, 'PreStress must be > 0');
                mbdyn.pre.base.checkNumericScalar (options.PreStrain, true, 'PreStrain');
                assert (options.PreStrain > 0, 'PreStrain must be > 0');
                
                options.SpringConstant = options.PreStress / options.PreStrain;
                self.preStrain = options.PreStrain;
                
%                 error ('You have supplied values for both the prestrain and prestress, these are mutually exclusive options');
            end
            
            self.node1 = node1;
            self.node2 = node2;
            self.springConstant = options.SpringConstant;
            self.zeroForceDistance = options.ZeroForceDistance;
            
        end
        
        function sysinputs = generateMBDynSystemInputs (self)
            % returns structure containing elements, drives and variables for tether
            %
            % Syntax
            %
            % sysinputs = generateMBDynSystemInputs (to)
            %
            % Description
            %
            % mndyn.pre.tether.generateMBDynSystemInputs returns a
            % structure containing the elements, drives and variables
            % required to create the tether in MBDyn.
            %
            % Input
            %
            %  to - mbdyn.pre.tether object
            %
            % Output
            %
            %  sysinputs - structure containing the following fields:
            %   
            %   Elements : mbdyn.pre.element objects required for the 
            %    tether
            %
            %   Drives : mbdyn.pre.drive objects required for the tether
            %
            %   Variables : mbdyn.pre.variable objects required for the 
            %    tether
            %
            %   These objects can be added to the mbdyn.pre.system object
            %   which will construct the system.
            %
            %
            
            tether_dl_var_name = sprintf('dl_%d', self.uid);
            
            v_dl = mbdyn.pre.variable ( 'real', tether_dl_var_name, 'Value', 0);
    
            dist_calc_str = sprintf ('model::distance (UID:%d,UID:%d)', self.node1.uid, self.node2.uid);
            
            tether_length_var_name = sprintf('tether_length_%d', self.uid);
            
            if isempty (self.zeroForceDistance)
                % initial distance calculation to get the intial length of the
                % tether to help with calculating dl during the simulation
                v_init_dist = mbdyn.pre.variable ( 'real', tether_length_var_name, ...
                                                   'Value', dist_calc_str, ...
                                                   'LabelRepObjects', { self.node1, ...
                                                                        self.node2 } ...
                                                 );
                                          
            else
                
                % initial distance calculation to get the intial length of the
                % tether to help with calculating dl during the simulation
                v_init_dist = mbdyn.pre.variable ( 'real', tether_length_var_name, ...
                                                   'Value', mbdyn.pre.base.formatNumber (self.zeroForceDistance) ...
                                                 );
            end
            
            d_instant_dist = mbdyn.pre.stringDrive ( dist_calc_str, ...
                                           'LabelRepObjects', { self.node1, self.node2 } ...
                                                    );
                                     
            % distance - initial_length + prestrain
            dlength_calc_str = sprintf ( 'model::drive(UID:%d,Time) - %s + %.17g', ...
                                         d_instant_dist.uid, ...
                                         tether_length_var_name, ...
                                         self.preStrain);

            d_dlength_calc = mbdyn.pre.stringDrive ( dlength_calc_str, ...
                                                     'LabelRepObjects', { d_instant_dist }, ...
                                                     'Name', [self.name, '_d_dlength_calc'] );
            
%             tension_calc_str = sprintf ('%s = model::drive(UID:%d,Time); (%s > 0) * %.17g * %s', ...
%                                         tether_dl_var_name, ...
%                                         d_dlength_calc.uid, ...
%                                         tether_dl_var_name, ...
%                                         self.springConstant, ...
%                                         tether_dl_var_name );

            tension_calc_str = sprintf ('(model::drive(UID:%d,Time) > 0) * %.17g * model::drive(UID:%d,Time)', ...
                                        d_dlength_calc.uid, ...
                                        self.springConstant, ...
                                        d_dlength_calc.uid );
                                    
            d_tension_calc = mbdyn.pre.stringDrive ( tension_calc_str, ...
                                                     'LabelRepObjects', { d_dlength_calc }, ...
                                                     'Name', [self.name, '_d_tension_calc'] );
                                                 
                                                 
            xunitvec_str = sprintf ('model::xunitvec (UID:%d,UID:%d)', self.node1.uid, self.node2.uid);
            
            d_xunitvec = mbdyn.pre.stringDrive ( xunitvec_str, ...
                                                 'LabelRepObjects', { self.node1, self.node2 }, ...
                                                 'Name', [self.name, '_d_xunitvec'] ...
                                                    );
                                                
            yunitvec_str = sprintf ('model::yunitvec (UID:%d,UID:%d)', self.node1.uid, self.node2.uid);
            
            d_yunitvec = mbdyn.pre.stringDrive ( yunitvec_str, ...
                                                 'LabelRepObjects', { self.node1, self.node2 }, ...
                                                 'Name', [self.name, '_d_yunitvec'] ...
                                                    );
                                                
            zunitvec_str = sprintf ('model::zunitvec (UID:%d,UID:%d)', self.node1.uid, self.node2.uid);
            
            d_zunitvec = mbdyn.pre.stringDrive ( zunitvec_str, ...
                                                 'LabelRepObjects', { self.node1, self.node2 }, ...
                                                 'Name', [self.name, '_d_zunitvec'] ...
                                                    );
                                                
                                                
            F_tension = tensionForces (self, d_tension_calc, d_xunitvec, d_yunitvec, d_zunitvec);
            
            sysinputs.Elements = { F_tension };
            sysinputs.Drives = { d_instant_dist, d_dlength_calc, d_tension_calc, d_xunitvec, d_yunitvec, d_zunitvec };
            sysinputs.Variables = { v_dl, v_init_dist };
            
        end
        
    end
    
    methods (Access=private)
        
        function F_tension = tensionForces (self, d_tension_calc, d_xunitvec, d_yunitvec, d_zunitvec)

            Fx_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::drive(UID:%d,Time)', ...
                                    d_tension_calc.uid, ...
                                    d_xunitvec.uid );

            Fy_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::drive(UID:%d,Time)', ...
                                    d_tension_calc.uid, ...
                                    d_yunitvec.uid );

            Fz_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::drive(UID:%d,Time)', ...
                                    d_tension_calc.uid, ...
                                    d_zunitvec.uid );

            d_fx_string_drive = mbdyn.pre.stringDrive ( Fx_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, d_xunitvec} );

            d_fy_string_drive = mbdyn.pre.stringDrive ( Fy_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, d_yunitvec} );

            d_fz_string_drive = mbdyn.pre.stringDrive ( Fz_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, d_zunitvec} );

            d_force_drive = mbdyn.pre.componentTplDriveCaller ( { d_fx_string_drive, ...
                                                                  d_fy_string_drive, ...
                                                                  d_fz_string_drive } );

            F_tension = mbdyn.pre.structuralInternalForce ( self.node1, self.node2, 'absolute', d_force_drive, ...
                                      'Position1', 'null', ...
                                      'Position1Ref', 'node', ...
                                      'Position2', 'null', ...
                                      'Position2Ref', 'other node', ...
                                      'Name', [self.name, '_F_tension'] );

        end 
        
    end
    
    
end