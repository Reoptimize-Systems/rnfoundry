classdef angularVelocity < mbdyn.pre.singleNodeJoint
    
    properties (GetAccess=public, SetAccess=protected)
        omegaDrive;
        rotAxis;
    end
    
    methods
        
        function self = angularVelocity (node, rot_axis, omega_drive)
            % angular velocity joint constructor
            %
            % Syntax
            %
            % av = mbdyn.pre.angularVelocity (node, rot_axis, omega_drive)
            %
            % Description
            %
            % Imposes the absolute angular velocity of a node about a given
            % axis.
            %
            % Input
            %
            %  node - node for which the angular velocity is to be
            %   prescribed
            %
            %  rot_axis - (3 x 1) vector representing the axis or rotation.
            %   Will be normalised to unity by MBDyn.
            %
            %  omega_drive - mbdyn.pre.drive object which sets the speed of
            %   rotation.
            %
            % Output
            %
            %  av - mbdyn.pre.angularVelocity object
            %
            %
            %
            % See Also: 
            %

            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.checkCartesianVector (rot_axis);
            
            assert ( isa (omega_drive, 'mbdyn.pre.drive'), ...
                'omega_drive must be an obect derived from mbdyn.pre.drive (i.e. it must be a drive of some type)' );
            
            self.type = 'angular velocity';
            self.omegaDrive = omega_drive;
            self.rotAxis = rot_axis;
            
        end
        
        function str = generateMBDynInputString (self)
            % generate an MBDyn input file string for the element
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList (self.rotAxis), ...
                                       2, ...
                                       true, ...
                                       'rotation axis' );
                                   
            str = self.addOutputLine ( str, ...
                                       self.omegaDrive.generateMBDynInputString (), ...
                                       2, ...
                                       false );
                                   
            str = self.addOutputLine ( str, ...
                                       ';', ...
                                       1, ...
                                       false, ...
                                       ['end ', self.type] );
            
        end
        
    end
    
end