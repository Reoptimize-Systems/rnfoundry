function design = materialmasses_RADIAL_SLOTTED(design, simoptions)
% calculates the material masses and volumes of materails in a slotted
% radial flux permanent magnet machine design
%
% Syntax
%
% [design, simoptions] = materialmasses_RADIAL_SLOTTED(design, simoptions)
% 
% Input
%
%   design - a structure containing the design of the radial flux machine,
%     it must contain at least the following fields:
%
%     Rmi - the inner radius of the magnets
%     Rmo - the outer radius of the magnets
%     ls - the stack length of the machine
%     thetam - the magnet pitch in radians
%     thetap - the pole pitch in radians
%     Poles - the number of Poles in the machine
%     Dc - diameter of wire used in the coils
%     MTL - mean turn length in a coil
%     CoilTurns - number of turns in a coil
%     NPhaseCoils - number of coils per phase winding
%     ArmatureIronAreaPerPole - cross-sectional area of the armature yoke,
%       and teeth
%
%   simoptions - another structure, this structure is expected to contiain
%     the field 'evaloptions'. evaloptions is another structure which must
%     contain the following fields:
%
%     MagnetDensity - desity of magnet material used in the rotor
%     FieldIronDensity - density of the field iron used in the rotor
%     CopperDensity - density of the wire material used
% 
% Output
%
%   design - the input design structure with the following fields appended:
%     
%     MagnetVolume - Volume of magnet material in the design
%     MagnetMass - Mass of magnet material in the design
%     FieldIronVolume - Volume of Field back iron in the design
%     FieldIronMass - Mass of Field back iron in the design
%     CopperVol - Volume of copper used in the design
%     CopperMass - Mass of copper used in the design
%     ArmatureIronVol - Volume of armature iron used in the design
%     ArmatureIronMass - Mass of armature iron used in the design
%     ArmatureSupportMass - currently always zero
%     RotorSupportMass - currently always zero
%     RotorMass - Mass of the rotor
%     StatorMass - Mass of the stator
%     TotalMass - Total mass of the machine
%
%   simoptions is returned unchanged
%

    % get the magnet mass (the rotor structural mass should have been
    % calculated by integrating over the FE structure mesh)
    design = materialmasses_RADIAL_ROTOR(design, simoptions);
    
    % calculate the copper mass
    design.CopperVol = design.NCoilsPerPhase ... % number of coils per phase
                       * design.Phases ... % number of Phases
                       * design.MTL ... % mean turn length in coil
                       * design.CoilTurns ... % number of turns per coil
                       * (pi * (design.Dc/2)^2); % cross-sectional area of wire
    
    design.CopperMass = design.CopperVol * simoptions.evaloptions.CopperDensity;

    % calculate the armature iron amss
    if isfield(design, 'ArmatureIronArea')
        design.ArmatureIronVol = design.ArmatureIronArea * design.ls;
    else
        error('design.ArmatureIronArea must be available.');
%     design.ArmatureIronVol = design.ls ...
%                              * ( design.Poles*annularsecarea(design.Ryi, design.Ryo, design.thetap) ...
%                                 + design.Qs * annularsecarea(design.Ryo, design.Rtsb, (design.thetas - design.thetac)) ...
%                                 + design.Qs * annularsecarea(design.Rtsg, design.Rao, (design.thetas - design.thetasg)) ...
%                                 + ( design.Qs * (annularsecarea(design.Rtsg, design.Rtsb, (design.thetas - design.thetasg)) ...
%                                                   - 2*sectorarea() ...
%                                                   - 2*(0.5*(design.Rtsg - design.Rtsb)*(sqrt(h^2 - ))))));
    end
    
    design.ArmatureIronMass = design.ArmatureIronVol * simoptions.evaloptions.ArmatureIronDensity;
    
    % set the epoxy mass to zero for the moment
    design.EpoxyVol = 0;
    design.EpoxyMass = design.EpoxyVol * simoptions.evaloptions.EpoxyDensity;
    
    % TODO calculate the support mass
    design.ArmatureSupportMass = 0;
    design.RotorSupportMass = 0;
    
    % calculate the component and module masses
    design.RotorMass = design.FieldIronMass ...
                       + design.MagnetMass ...
                       + design.RotorSupportMass;
    
    design.StatorMass = design.CopperMass ...
                        + design.EpoxyMass ...
                        + design.ArmatureSupportMass ...
                        + design.ArmatureIronMass;
                    
	design.TotalMass = design.RotorMass + design.StatorMass;
    
end