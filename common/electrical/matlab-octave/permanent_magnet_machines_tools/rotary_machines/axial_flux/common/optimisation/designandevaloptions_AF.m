function evaloptions = designandevaloptions_AF(evaloptions)
% parses the evaloptions structure for the function for axial flux machines
%
% Syntax
%
% evaloptions = designandevaloptions_AF(evaloptions)
%
% Input
%
%  evaloptions - optional structure containing various optional parameters
%   for the simultaion and evaluation. The following fields can be
%   specified:
%
%   E : (1 x 2) youngs modulus value for the structural material
%
%   MagnetCost : the cost per kg of magnet material (default: 80)
%
%   CopperCost : the cost per kg of copper wire (default: 10)
%
%   FieldIronCost : the cost per kg of the back iron (default: 0)
%
%   ArmatureIronCost : the cost per kg of laminated iron core material
%    (default: 3)
%
%   MagnetDensity : density of magnet material (default: 7500 kg/m^3)
%
%   CopperDensity : density of copper wire (default: 8960 kg/m^3)
%
%   FieldIronDensity : density of back iron (default: 7800 kg/m^3)
%
%   ArmatureIronDensity : density of laminated core material (default:
%    7800 kg/m^3)
%
%   StructMaterialDensity : density of structuaral material (default: 
%    7800 kg/m^3)
%
% Output
%
%  evaloptions - input structure with default options added where missing
%
%
% See also: designandevaloptions_ROTARY
%
%
    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = designandevaloptions_ROTARY(evaloptions);

    evaloptions = setfieldifabsent(evaloptions, 'E', 160e9);
    
    evaloptions = setfieldifabsent(evaloptions, 'structmeshoptions', struct());
    
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'ShaftAxialLayersPerM', 50);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'DiscAxialLayersPerM', 100);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'SupportAxialLayersPerM', 50);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'CircumPointsPerM', 20);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'BackIronRadialPointsPerM', 30);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'MagnetRadialPointsPerM', 30);
    evaloptions.structmeshoptions = setfieldifabsent(evaloptions.structmeshoptions, 'MeshType', 'H8');

end