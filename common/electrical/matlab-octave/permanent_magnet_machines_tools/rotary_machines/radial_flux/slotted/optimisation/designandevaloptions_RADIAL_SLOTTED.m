function evaloptions = designandevaloptions_RADIAL_SLOTTED(evaloptions)
% parses the evaloptions structure for the function for slotted radial flux
% machines
%
% Syntax
%
% evaloptions = designandevaloptions_RADIAL_SLOTTED(evaloptions)
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
%             MagnetCost - the cost per kg of magnet material (default: 80)
%
%             CopperCost - the cost per kg of copper wire (default: 10)
%
%             FieldIronCost - the cost per kg of the back iron (default: 0)
%
%             ArmatureIronCost - the cost per kg of laminated iron core
%             material (default: 3)
%
%             MagnetDensity - density of magnet material (default: 7500
%             kg/m^3)
%
%             CopperDensity - density of copper wire (default: 8960 kg/m^3)
%
%             FieldIronDensity - density of back iron (default: 7800 kg/m^3)
%
%             ArmatureIronDensity - density of laminated core material 
%             (default: 7800 kg/m^3)
%
%             StructMaterialDensity - density of structuaral material 
%             (default: 7800 kg/m^3)
%

    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = setfieldifabsent(evaloptions, 'ArmatureIronCost', 8);
    
    % set the default structural evaluation function for slotted radial
    % machines
    evaloptions = setfieldifabsent(evaloptions, 'structevalfcn', @evaluatestructure_RADIAL_SLOTTED);
    
    evaloptions = designandevaloptions_RADIAL(evaloptions);
    
end