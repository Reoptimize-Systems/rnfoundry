classdef loggingSettings
    
    properties (GetAccess=public, SetAccess = public)
        
        positions = true;
        velocities = true;
        accelerations = true;
        angularPositions = true;
        angularVelocities = true;
        angularAccelerations = true;
        
        nodeForces = true;
        nodeForcesUncorrected = true;
        forceHydro = true;
        forceExcitation = true;
        forceExcitationRamp = true;
        forceExcitationLin = true;
        forceExcitationNonLin = true;
        forceRadiationDamping = true;
        forceRestoring = true;
        forceMorrison = true;
        forceViscousDamping = true;
        forceAddedMassUncorrected = true;
        forceAddedMass = true;
        
        nodeMoments = true;
        nodeMomentsUncorrected = true;
        momentHydro = true;
        momentExcitation = true;
        momentExcitationRamp = true;
        momentExcitationLin = true;
        momentExcitationNonLin = true;
        momentRadiationDamping = true;
        momentRestoring = true;
        momentMorrison = true;
        momentViscousDamping = true;
        momentAddedMassUncorrected = true;
        momentAddedMass = true;
        
        powerTakeOffInternal = true;
        
    end
    
    methods
        
        function self = set.positions (self, val)
            check.isLogicalScalar (val, true, 'positions');
            self.positions = val;
        end
        
        function self = set.velocities (self, val)
            check.isLogicalScalar (val, true, 'velocities');
            self.velocities = val;
        end
        
        function self = set.accelerations (self, val)
            check.isLogicalScalar (val, true, 'accelerations');
            self.accelerations = val;
        end
        
        function self = set.angularPositions (self, val)
            check.isLogicalScalar (val, true, 'angularPositions');
            self.angularPositions = val;
        end
        
        function self = set.angularVelocities (self, val)
            check.isLogicalScalar (val, true, 'angularVelocities');
            self.angularVelocities = val;
        end
        
        function self = set.angularAccelerations (self, val)
            check.isLogicalScalar (val, true, 'angularAccelerations');
            self.angularAccelerations = val;
        end
        
        % forces logging settings
        
        function self = set.nodeForces (self, val)
            check.isLogicalScalar (val, true, 'nodeForces');
            self.nodeForces = val;
            
            self = self.setNodeForces (val);
            
        end
        
        function self = set.nodeForcesUncorrected (self, val)
            check.isLogicalScalar (val, true, 'nodeForcesUncorrected');
            self.nodeForcesUncorrected = val;
        end
        
        function self = set.forceHydro (self, val)
            check.isLogicalScalar (val, true, 'forceHydro');
            self.forceHydro = val;
        end
        
        function self = set.forceExcitation (self, val)
            check.isLogicalScalar (val, true, 'forceExcitation');
            self.forceExcitation = val;
        end
        
        function self = set.forceExcitationRamp (self, val)
            check.isLogicalScalar (val, true, 'forceExcitationRamp');
            self.forceExcitationRamp = val;
        end
        
        function self = set.forceExcitationLin (self, val)
            check.isLogicalScalar (val, true, 'forceExcitationLin');
            self.forceExcitationLin = val;
        end
        
        function self = set.forceExcitationNonLin (self, val)
            check.isLogicalScalar (val, true, 'forceExcitationNonLin');
            self.forceExcitationNonLin = val;
        end
        
        function self = set.forceRadiationDamping (self, val)
            check.isLogicalScalar (val, true, 'forceRadiationDamping');
            self.forceRadiationDamping = val;
        end
        
        function self = set.forceRestoring (self, val)
            check.isLogicalScalar (val, true, 'forceRestoring');
            self.forceRestoring = val;
        end
        
        function self = set.forceMorrison (self, val)
            check.isLogicalScalar (val, true, 'forceMorrison');
            self.forceMorrison = val;
        end
        
        function self = set.forceViscousDamping (self, val)
            check.isLogicalScalar (val, true, 'forceViscousDamping');
            self.forceViscousDamping = val;
        end
        
        function self = set.forceAddedMass (self, val)
            
            check.isLogicalScalar (val, true, 'forceAddedMass');
            
            self.forceAddedMass = val;
            
            self = self.setForceAddedMass (val);
            
        end
        
        % moments logging settings
        
        function self = set.nodeMoments (self, val)
            check.isLogicalScalar (val, true, 'nodeForces');
            self.nodeMoments = val;
            
            self = self.setNodeMoments (val);
            
        end
        
        function self = set.nodeMomentsUncorrected (self, val)
            check.isLogicalScalar (val, true, 'nodeMomentsUncorrected');
            self.nodeMomentsUncorrected = val;
        end
        
        function self = set.momentHydro (self, val)
            check.isLogicalScalar (val, true, 'momentHydro');
            self.momentHydro = val;
        end
        
        function self = set.momentExcitation (self, val)
            check.isLogicalScalar (val, true, 'momentExcitation');
            self.momentExcitation = val;
        end
        
        function self = set.momentExcitationRamp (self, val)
            check.isLogicalScalar (val, true, 'momentExcitationRamp');
            self.momentExcitationRamp = val;
        end
        
        function self = set.momentExcitationLin (self, val)
            check.isLogicalScalar (val, true, 'momentExcitationLin');
            self.momentExcitationLin = val;
        end
        
        function self = set.momentExcitationNonLin (self, val)
            check.isLogicalScalar (val, true, 'momentExcitationNonLin');
            self.momentExcitationNonLin = val;
        end
        
        function self = set.momentRadiationDamping (self, val)
            check.isLogicalScalar (val, true, 'momentRadiationDamping');
            self.momentRadiationDamping = val;
        end
        
        function self = set.momentRestoring (self, val)
            check.isLogicalScalar (val, true, 'momentRestoring');
            self.momentRestoring = val;
        end
        
        function self = set.momentMorrison (self, val)
            check.isLogicalScalar (val, true, 'momentMorrison');
            self.momentMorrison = val;
        end
        
        function self = set.momentViscousDamping (self, val)
            check.isLogicalScalar (val, true, 'momentViscousDamping');
            self.momentViscousDamping = val;
        end
        
        function self = set.momentAddedMass (self, val)
            
            check.isLogicalScalar (val, true, 'momentAddedMass');
            
            self.momentAddedMass = val;
            
            self = self.setMomentAddedMass (val);
            
        end
        
        
        function self = set.powerTakeOffInternal (self, val)
            check.isLogicalScalar (val, true, 'powerTakeOffInternal');
            self.powerTakeOffInternal = val;
        end
        
    end
    
    methods (Access=private)
        
        function self = setForceAddedMass (self, val)
            % function needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.forceAddedMassUncorrected = true;
                self.momentAddedMassUncorrected = true;
            end
        end
        
        function self = setNodeForces (self, val)
            % function needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.nodeForcesUncorrected = true;
                self.nodeMomentsUncorrected = true;
                self.forceAddedMassUncorrected = true;
                self.momentAddedMassUncorrected = true;
            end
        end
        
        function self = setMomentAddedMass (self, val)
            % function needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.forceAddedMassUncorrected = true;
                self.momentAddedMassUncorrected = true;
            end
        end
        
        function self = setNodeMoments (self, val)
            % function needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.nodeForcesUncorrected = true;
                self.nodeMomentsUncorrected = true;
                self.forceAddedMassUncorrected = true;
                self.momentAddedMassUncorrected = true;
            end
        end
        
    end
    
end