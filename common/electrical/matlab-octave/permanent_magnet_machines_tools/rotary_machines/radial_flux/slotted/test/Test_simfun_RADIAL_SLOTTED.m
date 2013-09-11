% Test_simfun_RADIAL_SLOTTED
%
%

clear design simoptions 

design.StatorType = 'so';
design.poles = 2;
design.phases = 3;
design.CoilLayers = 2;
design.Qc = design.phases * design.poles;
if design.CoilLayers == 1
    design.Qs = design.phases * 2 * design.Qc;
elseif design.CoilLayers == 2
    design.Qs = design.phases * 1 * design.Qc;
end
design.yd = 4;
design.thetap = 2*pi/design.poles;
design.thetam = design.thetap * 0.8;
design.thetac = (2*pi / design.Qs) * 0.85;
design.thetasg = design.thetac * 0.95;
design.tm = 0.01;
design.tbi = 0.01;
design.ty = 0.01;
design.tc = 0.03;
design.tsb = 0.01;
design.tsg = 0; %0.01;
design.g = 3/1000;
design.Rmo = 0.5;
design.Rmi = 0.5;
design.ls = 0.3;

if strcmp(design.StatorType, 'si')
    design.Rmo = 0.5;
    design.Rmi = design.Rmi - design.tm;
    design.Rmm = mean([design.Rmi, design.Rmo]);
    design.Rci = design.Rmo + design.g + design.tsb;
    design.Rco = design.Rci + design.tc;
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbo = design.Rmi;
    design.Rbi = design.Rbo - design.tbi;
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Ryi = design.Rco;
    design.Ryo = design.Rco + design.ty;
    design.Rym = mean([design.Ryi, design.Ryo]);
elseif strcmp(design.StatorType, 'so')
    design.Rmi = 0.5;
    design.Rmo = design.Rmi + design.tm;
    design.Rmm = mean([design.Rmi, design.Rmo]);
    design.Rco = design.Rmi - design.g - design.tsb;
    design.Rso = design.Rco + design.tsb;
    design.Rci = design.Rco - design.tc;
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbi = design.Rmo;
    design.Rbo = design.Rbi + design.tbi;
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Ryo = design.Rci;
    design.Ryi = design.Ryo - design.ty;
    design.Rym = mean([design.Ryi, design.Ryo]);
end

design.Dc = design.Rcm * design.thetac / 100;
design.fillfactor = 0.7;

design.Hc = design.tc / design.CoilLayers;
design.CoilTurns = 250;

design.NCoilsPerPhase = design.Qc / design.phases;

design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.YokeMaterial = design.BackIronMaterial;
design.CoilMaterial = '36 AWG';

design.RgVRc = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);
fprintf(1, 'done simfun\n');
[design, simoptions] = finfun_RADIAL_SLOTTED(design, simoptions);

