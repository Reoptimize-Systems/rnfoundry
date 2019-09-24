classdef linearVelocity < mbdyn.pre.singleNodeJoint
    
    properties (GetAccess=public, SetAccess=protected)
        velocityDrive;
        relativeDirection;
    end
    
    methods
        
        function self = linearVelocity (node, relative_direction, velocity_drive)
            % linear velocity joint constructor
            %
            % Syntax
            %
            % av = mbdyn.pre.linearVelocity (node, rot_axis, velocity_drive)
            %
            % Description
            %
            % Imposes the absolute linear velocity of a node about a given
            % axis.
            %
            % Input
            %
            %  node - node for which the linear velocity is to be
            %   prescribed
            %
            %  relative_direction - (3 x 1) vector representing the
            %   relative direction of the velocity. Will be normalised to
            %   unity by MBDyn.
            %
            %  velocity_drive - mbdyn.pre.drive object which sets the value
            %   of the velocity.
            %
            % Output
            %
            %  av - mbdyn.pre.linearVelocity object
            %
            %
            %
            % See Also: 
            %

            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.checkCartesianVector (relative_direction);
            
            assert ( isa (velocity_drive, 'mbdyn.pre.drive'), ...
                'velocity_drive must be an obect derived from mbdyn.pre.drive (i.e. it must be a drive of some type)' );
            
            self.type = 'linear velocity';
            self.velocityDrive = velocity_drive;
            self.relativeDirection = relative_direction;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for linearVelocity joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (av)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  av - mbdyn.pre.linearVelocity object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            
            str = self.addOutputLine ( str, ...
                                       self.commaSepList (self.relativeDirection), ...
                                       2, ...
                                       true, ...
                                       'relative direction' );
                                   
            str = self.addOutputLine ( str, ...
                                       self.velocityDrive.generateMBDynInputString (), ...
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