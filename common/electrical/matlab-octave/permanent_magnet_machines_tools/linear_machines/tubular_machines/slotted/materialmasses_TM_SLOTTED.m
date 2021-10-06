function design = materialmasses_TM_SLOTTED (design, simoptions)
% calculates the material masses and volumes of materails in a slotted
% radial flux permanent magnet machine design
%
% Syntax
%
% [design, simoptions] = materialmasses_TM_SLOTTED(design, simoptions)
% 
% Input
%
%   design - a structure containing the design of the radial flux machine,
%     it must contain at least the following fields:
%
%     Rfo - the inner radius of the magnets
%     Rag - the outer radius of the magnets
%     zm - the magnet pitch in radians
%     zp - the pole pitch in radians
%     Poles - the number of Poles in the machine
%     Dc - diameter of wire used in the coils
%     MTL - mean turn length in a coil
%     CoilTurns - number of turns in a coil
%     NPhaseCoils - number of coils per phase winding
%     ArmatureIronVol -volume of the armature yoke and teeth
%
%   simoptions - another structure, this structure is expected to contiain
%     the field 'evaloptions'. evaloptions is another structure which must
%     contain the following fields:
%
%     MagnetDensity - desity of magnet material used in the Field
%     FieldIronDensity - density of the field iron used in the Field
%     CopperDensity - density of the wire material used
% 
% Output
%
%   design - the input design structure with the following fields appended:
%     
%     MagnetVolume - Volume of magnet material in the design
%     MagnetMass - Mass of magnet material in the design
%     FieldSpacerVol - Volume of Field back iron in the design
%     FieldIronMass - Mass of Field back iron in the design
%     CopperVol - Volume of copper used in the design
%     CopperMass - Mass of copper used in the design
%     ArmatureIronVol - Volume of armature iron used in the design
%     ArmatureIronMass - Mass of armature iron used in the design
%     ArmatureSupportMass - currently always zero
%     FieldSupportMass - currently always zero
%     FieldMass - Mass of the Field
%     ArmatureMass - Mass of the Armature
%     TotalMass - Total mass of the machine
%
%   simoptions is returned unchanged
%

    % calculate the armature iron mass
    if ~isfield(design, 'ArmatureIronVol')
        error('ArmatureIronVol must be available in the design structure.')
    end
    
    % get the magnet mass (the Field structural mass should have been
    % calculated by integrating over the FE structure mesh)
    design = materialmasses_TM (design, simoptions);
    
end