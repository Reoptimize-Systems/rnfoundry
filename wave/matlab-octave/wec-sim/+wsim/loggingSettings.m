classdef loggingSettings
    
    properties (GetAccess=public, SetAccess = public)
        
        positions = true;
        velocities = true;
        accelerations = true;
        nodeForcesAndMoments = true;
        nodeForcesAndMomentsUncorrected = true;
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
        
        function self = set.nodeForcesAndMoments (self, val)
            check.isLogicalScalar (val, true, 'nodeForcesAndMoments');
            self.nodeForcesAndMoments = val;
            
            self = self.setNodeForcesAndMoments (val);
            
        end
        
        function self = set.nodeForcesAndMomentsUncorrected (self, val)
            check.isLogicalScalar (val, true, 'nodeForcesAndMomentsUncorrected');
            self.nodeForcesAndMomentsUncorrected = val;
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
            
            self = self.setforceAddedMass (val);
            
        end
        
        function self = set.powerTakeOffInternal (self, val)
            check.isLogicalScalar (val, true, 'powerTakeOffInternal');
            self.powerTakeOffInternal = val;
        end
        
    end
    
    methods (Access=private)
        
        function self = setforceAddedMass (self, val)
            % funciton needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.forceAddedMassUncorrected = true;
            end
        end
        
        function self = setNodeForcesAndMoments (self, val)
            % funciton needed so we can set value of other property which
            % is otherwise prohibited
            
            if val == true
                % the uncorrected added mass force must be logged to get
                % corrected added mass force
                self.nodeForcesAndMomentsUncorrected = true;
                self.forceAddedMassUncorrected = true;
            end
        end
        
    end
    
end