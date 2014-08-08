% Validate_iron_losses_RADIAL_SLOTTED

% losses will be validated by comparison with results taken from the paper:
%
% MI et al., "MODELING OF IRON LOSSES OF PM SYNCHRONOUS MOTORS", IEEE
% TRANSACTIONS ON INDUSTRY APPLICATIONS, VOL. 39, NO. 3, MAY/JUNE 2003
%
%

%% set up design parameters
             
design.Poles = 4;
design.Qs = 36;
design.CoilLayers = 2;
design.Qc = design.Qs;
design.Phases = 3;

design.Ryo = 95e-3;
design.g = 2e-3;
design.ty = 17.4e-3;
design.tm = 6.3e-3;
design.Rai = 58.5e-3;
design.Rmo = design.Rai - design.g;
design.tbi = 29.9523e-3;
design.tsb = 2.604e-3;
design.tsg = 1.7364e-3;
design.tc = 19.1e-3 - design.tsb;
design.Rtsb = design.Rai + design.tsb;
design.Ryi = design.Rtsb + design.tc;
design.thetam = (tau / design.Poles) * 0.667;

chordanglefcn = @(r,l) 2 .* asin(l / 2.*r);

design.thetacg = feval(chordanglefcn, design.Rtsb, 5.6432e-3);
design.thetacy = feval(chordanglefcn, design.Ryi, 6.9455e-3);
design.thetasg = feval(chordanglefcn, design.Ryi, 3e-3);
design.ls = 88.9e-3;

design.ArmatureType = 'external';

design = completedesign_RADIAL_SLOTTED(design);

design.MagSimMaterials.Magnet = 'NdFeB 40 MGOe';
design.MagSimMaterials.FieldIron = '1117 Steel';
design.MagSimMaterials.ArmatureIron = design.MagSimMaterials.FieldIron;
design.MagSimMaterials.CoilWinding = '36 AWG';
design.MagSimMaterials.Gap = 'Air';
    
    
%% get the machine data

design.CoilFillFactor = 0.6;
design.CoilTurns = 1;

design.RlVRp = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_RADIAL_SLOTTED (design, simoptions);

[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);


