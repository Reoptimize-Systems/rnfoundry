function [score, design, simoptions, T, Y, results] = designandevaluate_ACTM(design, simoptions)
% designandevaluate_ACTM: a function to evaluate a design of the linear
% tubular permanent magnet machine with ferromagnetic sheath
%
% Syntax:
%
%   score = designandevaluate_ACTIAM(design, simoptions, evaloptions)
%   [..., design] = designandevaluate_ACTIAM(design, simoptions, evaloptions)
%   [..., T] = designandevaluate_ACTIAM(design, simoptions, evaloptions)
%   [..., Y] = designandevaluate_ACTIAM(design, simoptions, evaloptions)
%   [..., results] = designandevaluate_ACTIAM(design, simoptions, evaloptions)
%
% Input
%
%  design - structure describing the air cored tubular machine design,
%   containing the following fields:
%
%   WmVWp : scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm : scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRi : scalar value of Ro/Ri Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo : scalar value of Ra/Ro Ratio for machine to be evaluated, the
%           supporting sheath outer to inner radius ratio 
%
%   RsoVRm : scalar value of the outer shaft radius to the translator
%            radius, Rso / Rm
%
%   RiVRm : scalar value of the Ri / Rm ratio, the inner coil radius to
%           translator radius
%
%   WcVWp : scalar ratio of the coil width to pole width in the machine 
%           (Wc / Wp)
%
%   Rm : Translator radius in m
%
%   RlVRp : The apparent grid resistance to coil resistance ratio
%
%   kfill : the copper fill factor in the coil slots
%
%   Dc : the conductor diameter
%
%   lm : magnet depth
%
%   Evaluation : optional structure containing various optional parameters
%    for the simultaion and evaluation. The following fields can be
%    specified:
%
%    E : (1 x 3) Vector containing Young's modulus of the beam materials.
%     The first value is Young Modulus of the shaft, the second the the
%     young's modulus of the coil and the third the young's modulus of the
%     outer sheath [shaft coil sheath]
%
%    MagnetCost : the cost per kg of magnet material (default: 30)
%
%    CopperCost : the cost per kg of copper wire (default: 10)
%
%    FieldIronCost : the cost per kg of the back iron (default: 0)
%
%    armatureIronCost : the cost per kg of laminated iron core material 
%     (default: 3)
%
%    buoyMassCost : the cost of the ballast in the buoy per kg, if
%     evaluating a machine and buoy system (default: 0.2)
%
%    MagnetDensity : density of magnet material (default: 7500 kg/m^3)
%
%    CopperDensity : density of copper wire (default: 8960 kg/m^3)
%
%    FieldIronDensity : density of back iron (default: 7800 kg/m^3)
%
%    ArmatureIronDensity : density of laminated core material 
%    (default: 7800 kg/m^3)
%
%    StructMaterialDensity : density of structuaral material (default: 
%     7800 kg/m^3)
%
%    sections : the number of sections into which the beam is split for
%     analysis of the deflections due to maxwell stresses (default: 10)
%
%    minRMSemf : the minimum operating RMS voltage
%
%    maxJrms : the minimum current density
%
%    v : operating speed for linear power rating (default: 2.2 m/s)
%
%    pointsPerPole : number of points at which the voltage and hence power
%     will be evaluated for each pole (default: 40)
%
%    control : Determines if some kind of stator switching control is used,
%     if 0, none is used and the total phase resistance will be higher due
%     to inactive coils. If 1, only aoverlapping coils are active and the
%     resistance will be lower.
%
%    targetPower : If not zero, this is a target power rating for the
%     machine. The number of field Poles will be multiplied to the required
%     number to achieve this output power. If this is present the system
%     will look for the mlength field and act accordingly depending on what
%     is found. (default: 500e3 W)
%
%    mlength : If targetPower is present and is zero, and mlength is
%     supplied, mlength should contain a vector of two lengths, the field
%     and armature respectively). If targetPower is zero and mlength is not
%     present it is set to one pole pitch for each part. If targetPower is
%     present and greater than zero, mlength should be a scalar value of
%     the overlap between field and armature in metres. If not supplied, it
%     is set to zero. If there is no targetPower supplied, mlength should
%     contain the length of the field and armature.
%
% Output
%
%  CostPerW - cost per watt of the machine
%
%   design - the input design structure with additional fields added
%
%
% See Also: designandevaluate_TM

% TODO: correct and complete help 

    simoptions = setfieldifabsent (simoptions, 'Evaluation', []);

    % parse the options structure filling in any missing items
    simoptions.Evaluation = designandevaloptions_ACTM( simoptions.Evaluation);
    
    simoptions.FieldIronDensity = simoptions.Evaluation.FieldIronDensity;
    simoptions.MagnetDensity = simoptions.Evaluation.MagnetDensity;
    simoptions.ShaftDensity = simoptions.Evaluation.StructMaterialDensity;
    simoptions.CopperDensity = simoptions.Evaluation.CopperDensity;
    simoptions.ArmatureIronDensity = simoptions.Evaluation.ArmatureIronDensity;
    
    [T, Y, results, design, simoptions] = designandevaluate_TM(design, simoptions);
    
    % design the ACTM structure
    design = designstructure_ACTM(design,  simoptions.Evaluation);

    % Finally we determine the machine score and costs etc.
    [score, design] = machinescore_ACTM(design, simoptions);

    if isfield(design,'i')
        fprintf(1, 'Ind %d, Score:  %f\n', design.i, score);
    end
    
end