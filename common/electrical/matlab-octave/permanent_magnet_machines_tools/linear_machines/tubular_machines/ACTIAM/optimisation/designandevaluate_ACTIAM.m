function [score, design, simoptions, T, Y, results] = designandevaluate_ACTIAM(design, simoptions)
% designandevaluate_ACTIAM: a function to evaluate a design of the linear
% tubular permanent magnet machine with ferromagnetic sheath
%
% Syntax:
%
%   [score, design, simoptions, T, Y, results] = designandevaluate_ACTIAM(design, simoptions)
%
% Input:
%
%   design, a structure containing the following members:
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio, the sheath outer radius to inner
%           radius ratio
%
%   RsoVRm - scalar value of the outer shaft radius to the translator
%            radius, Rso / Rm
%
%   RiVRm - scalar value of the Ri / Rm ratio, the inner coil radius to
%           translator radius
%
%   WcVWp - scalar ratio of the coil width to pole width in the machine 
%           (Wc / Wp)
%
%   AconVAcu - The conductor cross-sectional area to copper area ratio
%
%   Rm - Translator radius in m
%
%   RlVRp - The apparent grid resistance to coil resistance ratio
%
%   kfill - the copper fill factor in the coil slots
%
%   lm - magnet depth
%
%   simoptions.evaloptions - optional structure containing various optional parameters for
%             the simultaion and evaluation. The following fields can be
%             specified:
%
%             E - (1 x 3) Vector containing Young's modulus of the beam
%             materials. The first value is Young Modulus of the shaft, the
%             second the the young's modulus of the coil and the third the
%             young's modulus of the outer sheath [shaft coil sheath]
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
%
% Output:
%
%   CostPerW - cost per watt of the machine
%
%   design - a structure containing detailed information on the
%            machine design
%
    if nargin < 2
        simoptions.evaloptions = [];
    end
    
    simoptions.evaloptions = designandevaloptions_ACTIAM(simoptions.evaloptions);

    simoptions.FieldIronDensity = simoptions.evaloptions.FieldIronDensity;
    simoptions.MagnetDensity = simoptions.evaloptions.MagnetDensity;
    simoptions.ShaftDensity = simoptions.evaloptions.StructMaterialDensity;
    simoptions.CopperDensity = simoptions.evaloptions.CopperDensity;
    simoptions.ArmatureIronDensity = simoptions.evaloptions.ArmatureIronDensity;
    
    [T, Y, results, design, simoptions] = designandevaluate_TM(design, simoptions);
    
    % Choose appropriate sections to split up the coil
    %[Sr, Sz] = ChooseCoilSections_TM(design.RiVRm, design.RoVRm, design.PoleWidthVRm, design.WcVWp, design.CoilTurns);
    Sr = 5;
    Sz = 5;
    
    % determine the structure required for the machine
    design = designstructure_ACTIAM(design, simoptions.evaloptions, Sr, Sz);                                                 
        
    % score the machine design
    [score, design] = machinescore_ACTIAM(design, simoptions);
    
    if isfield(design,'i')
        fprintf(1, 'Ind %d, Score:  %f\n', design.i, score);
    end
    
end 
