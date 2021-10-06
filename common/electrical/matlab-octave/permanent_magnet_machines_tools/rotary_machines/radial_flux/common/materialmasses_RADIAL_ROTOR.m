function [design, simoptions] = materialmasses_RADIAL_ROTOR(design, simoptions)
% calcualtes the mass of a radial flux permanent magnet machine rotor
%
% Syntax
%
% [design, simoptions] = materialmasses_RADIAL_ROTOR(design, simoptions)
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
%
%   simoptions - another structure, this structure is expected to contiain
%     the field 'evaloptions'. evaloptions is another structure which must
%     contain the following fields:
%
%     MagnetDensity - desity of magnet material used in the rotor
%     FieldIronDensity - density of the field iron used in the rotor
% 
% Output
%
%   design - the input design structure with the following fields appended:
%     
%     MagnetVolume - Volume of magnet material in the design
%     MagnetMass - Mass of magnet material in the design
%     FieldIronVolume - Volume of Field back iron in the design
%     FieldIronMass - Mass of Field back iron in the design
%
%   simoptions is returned unchanged
%




    design.MagnetVolume = annularsecarea(design.Rmi, design.Rmo, design.thetam) ...
                          * design.ls ...
                          * design.Poles;
                      
    design.MagnetMass = design.MagnetVolume * simoptions.Evaluation.MagnetDensity; 
    
    design.FieldIronVolume = annularsecarea(design.Rbi, design.Rbo, design.thetap) ...
                             * design.ls ...
                             * design.Poles;
    
    design.FieldIronMass = design.FieldIronVolume * simoptions.Evaluation.FieldIronDensity;
    
end