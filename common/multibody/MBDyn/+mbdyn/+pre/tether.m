classdef tether < mbdyn.pre.base
    
    
    properties
        
        preStrain;
        node1;
        node2;
        springConstant;
        
    end
    
    
    methods
        
        function self = tether (node1, node2, varargin)
            
            options.PreStrain = [];
            options.PreStress = [];
            options.SpringConstant = 1e3;
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.base ();
            
            mbdyn.pre.base.checkNumericScalar (options.SpringConstant, true, 'SpringConstant');
            assert (options.SpringConstant > 0, 'SpringConstant must be > 0');
            
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
                
                mbdyn.pre.base.checkNumericScalar (options.PreStress, true, 'PreStress');
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
            
        end
        
        function sysinputs = generateMBDynSystemInputs (self)
            
            tether_dl_var_name = sprintf('dl_%d', self.uid);
            
            v_dl = mbdyn.pre.variable ( 'real', tether_dl_var_name, 'Value', 0);
    
            dist_calc_str = sprintf ('model::distance (UID:%d,UID:%d)', self.node1.uid, self.node2.uid);
            
            tether_length_var_name = sprintf('tether_length_%d', self.uid);
            
            % initial distance calculation to get the intial length of the
            % tether to help with calculating dl during the simulation
            v_init_dist = mbdyn.pre.variable ( 'real', tether_length_var_name, ...
                                               'Value', dist_calc_str, ...
                                               'LabelRepObjects', { self.node1, ...
                                                                    self.node2 } ...
                                              );
            
            d_instant_dist = mbdyn.pre.stringDrive ( dist_calc_str, ...
                                           'LabelRepObjects', { self.node1, self.node2 } ...
                                                    );
                                     
            dlength_calc_str = sprintf ( 'model::drive(UID:%d,Time) - %s + %.17g', ...
                                         d_instant_dist.uid, ...
                                         tether_length_var_name, ...
                                         self.preStrain );

            d_dlength_calc = mbdyn.pre.stringDrive ( dlength_calc_str, ...
                                                     'LabelRepObjects', { d_instant_dist } );
            
            tension_calc_str = sprintf ('%s = model::drive(UID:%d,Time); (%s > 0) * %.17g * %s', ...
                                        tether_dl_var_name, ...
                                        d_dlength_calc.uid, ...
                                        tether_dl_var_name, ...
                                        self.springConstant, ...
                                        tether_dl_var_name );
                                    
            d_tension_calc = mbdyn.pre.stringDrive ( tension_calc_str, ...
                                                     'LabelRepObjects', { d_dlength_calc } );
                                                 
            F_tension = tensionForces (self, d_tension_calc);
            
            sysinputs.Elements = { F_tension };
            sysinputs.Drives = { d_instant_dist, d_dlength_calc, d_tension_calc };
            sysinputs.Nodes = {};
            sysinputs.Variables = { v_dl, v_init_dist };
            
        end
        
    end
    
    methods (Access=private)
        
        function F_tension = tensionForces (self, d_tension_calc)

            Fx_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::xunitvec(UID:%d,UID:%d)', ...
                                    d_tension_calc.uid, ...
                                    self.node1.uid, ...
                                    self.node2.uid );

            Fy_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::yunitvec(UID:%d,UID:%d)', ...
                                    d_tension_calc.uid, ...
                                    self.node1.uid, ...
                                    self.node2.uid );

            Fz_calc_str = sprintf ( 'model::drive(UID:%d,Time) * model::zunitvec(UID:%d,UID:%d)', ...
                                    d_tension_calc.uid, ...
                                    self.node1.uid, ...
                                    self.node2.uid );

            d_fx_string_drive = mbdyn.pre.stringDrive ( Fx_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, self.node1, self.node2} );

            d_fy_string_drive = mbdyn.pre.stringDrive ( Fy_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, self.node1, self.node2} );

            d_fz_string_drive = mbdyn.pre.stringDrive ( Fz_calc_str, ...
                                                        'LabelRepObjects', {d_tension_calc, self.node1, self.node2} );

            d_force_drive = mbdyn.pre.componentTplDriveCaller ( { d_fx_string_drive, ...
                                                                  d_fy_string_drive, ...
                                                                  d_fz_string_drive } );

            F_tension = mbdyn.pre.structuralInternalForce ( self.node1, self.node2, 'absolute', d_force_drive, ...
                                      'Position1', 'null', ...
                                      'Position1Ref', 'node', ...
                                      'Position2', 'null', ...
                                      'Position2Ref', 'other node' );

        end 
        
    end
    
    
end