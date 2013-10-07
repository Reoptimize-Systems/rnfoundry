function evaloptions = designandevaloptions_ACPMSM(evaloptions)
% parses the evaloptions structure for the function
% designandevaloptions_ACPMSM filling in missing values with defaults etc.
%
% Input:
%
%   evaloptions - optional structure containing various optional parameters for
%             the simultaion and evaluation. The following fields can be
%             specified:
%
%             E - (1 x 2) vector of youngs modulus values for the I-beams
%                 and the translator central section respectively (default:
%                 [200e9 151e9])
%
%             magCost - the cost per kg of magnet material (default: 30)
%
%             copperCost - the cost per kg of copper wire (default: 10)
%
%             backIronCost - the cost per kg of the back iron (default: 0)
%
%             armatureIronCost - the cost per kg of laminated iron core
%             material (default: 3)
%
%             magnetDensity - density of magnet material (default: 7500
%             kg/m^3)
%
%             copperDensity - density of copper wire (default: 8960 kg/m^3)
%
%             backIronDensity - density of back iron (default: 7800 kg/m^3)
%
%             armatureIronDensity - density of laminated core material 
%             (default: 7800 kg/m^3)
%
%             structMaterialDensity - density of structuaral material 
%             (default: 7800 kg/m^3)
%
%             sections - the number of sections into which the beam is
%             split for analysis of the deflections due to maxwell stresses
%             (default: 10)
%
%             minRMSemf - the minimum operating RMS voltage
%
%             maxJrms - the minimum current density
%
%             alphab - the beam overlap factor for the bearings, i.e. if ls
%             is 1m and the bearings are 1.2 m long
%
%             mode - specifies whether we are simulating a double sided or
%             single sided machine. If mode is 1, it is double sided, if 0,
%             it is single sided (default: 1)
%
%             v - operating speed for linear power rating (default: 2.2
%             m/s)
%
%             pointsPerPole - number of points at which the voltage and
%             hence power will be evaluated for each pole (default: 40)
%
%             control - Determines if some kind of stator switching control
%             is used, if 0, none is used and the total phase resistance
%             will be higher due to inactive coils. If 1, only aoverlapping
%             coils are active and the resistance will be lower.
%
%             targetPower - If not zero, this is a target power rating for
%             the machine. The number of field Poles will be multiplied to
%             the required number to achieve this output power. If this is
%             present the system will look for the mlength field and act 
%             accordingly depending on what is found. (default: 500e3 W)
%
%             mlength - If targetPower is present and is zero, and mlength
%             is supplied, mlength should contain a vector of two lengths,
%             the field and armature respectively). If targetPower is zero
%             and mlength is not present it is set to one pole pitch for
%             each part. If targetPower is present and greater than zero,
%             mlength should be a scalar value of the overlap between field
%             and armature in metres. If not supplied, it is set to zero.
%             If there is no targetPower supplied, mlength should contain
%             the length of the field and armature. 


    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = designandevaloptions_FM(evaloptions);
    
    if nargin == 0
        
        evaloptions.E = [200e9 151e9];
        evaloptions.mode = 1; % mode: 1 is double sided machine, 0 is single sided
        evaloptions.v = 2.2;
        evaloptions.pointsPerPole = 40;
        evaloptions.control = 1;
        evaloptions.presimfinfun = 'prescribedmotfinfun_ACPMSM';
        
    elseif nargin == 1
        % If the evaloptions structure is supplied, we check each possible
        % optional values and set any not supplied to the default value

        if ~isfield(evaloptions, 'E')
            evaloptions.E = [200e9 151e9];
        end
        
        if ~isfield(evaloptions, 'mode')
            evaloptions.mode = 1;
        end
        
        if ~isfield(evaloptions, 'v')
            evaloptions.v = 2.2;
        end
        
        if ~isfield(evaloptions, 'pointsPerPole')
            evaloptions.pointsPerPole = 40;
        end
                
        if ~isfield(evaloptions, 'control')
            evaloptions.control = 1;
        end
        
        evaloptions.presimfinfun = 'prescribedmotfinfun_ACPMSM';
        
    end  
    
end
