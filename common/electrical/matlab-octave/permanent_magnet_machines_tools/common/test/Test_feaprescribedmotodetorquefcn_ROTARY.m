% Test_feaprescribedmotodetorquefcn_ROTARY.m

clear design simoptions

% design = test_design_RADIAL_SLOTTED ('external');

design.Rbo = 0.161219144643046;

design.RmoVRbo = 0.972957708842189;

design.RmiVRmo = 0.965365222464475;

design.RaoVRmi = 0.990094211892778;

design.RtsbVRao = 0.977065387763626;

design.RyoVRtsb = 0.698015840417782;

design.RyiVRyo = 0.865558917014541;

design.tsgVtsb = 0.264918009501899;

design.thetamVthetap = 0.684728045138945;

design.thetacgVthetas = 0.513909963593437;

design.thetacyVthetas = 0.5;

design.thetasgVthetacg = 0.654923245696986;

design.lsVtm = 7.36269843608709;

design.ArmatureType = 'internal';

design.Phases = 3;

design.Poles = 52;

design.Qc = 39;

design.yd = 1;

design.CoilLayers = 1;

design.Branches = 1;

design.CoilsPerBranch = 13;

design.CoilTurns = 37;

design.Dc = 2.032411e-03;

design.MagnetSkew = 1.430928e-01;

design.NSkewMagnetsPerPole = 2;

design.CoilInsulationThickness = 2.000000e-04;

design.NStrands = 1;

design.NStages = 1;

design.MagFEASimMaterials.Magnet = 'NdFeB 40 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureYoke = design.MagFEASimMaterials.FieldBackIron;
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.HeatFEASimMaterials.Magnet = 'Iron, Pure';
design.HeatFEASimMaterials.FieldBackIron = 'Iron, Pure';
design.HeatFEASimMaterials.ArmatureYoke = design.HeatFEASimMaterials.FieldBackIron;
design.HeatFEASimMaterials.ArmatureCoil = 'Copper, Pure';
design.HeatFEASimMaterials.AirGap = 'Water';
design.HeatFEASimMaterials.CoilInsulation = 'Water';


%%
design.RlVRp = 10;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'InterpolateMissing', false);


simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_RADIAL_SLOTTED (design, simoptions);

% setup simulation options
%simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
%                                'Rpm', 30, ...
%                                'TSpan', [0,1]);

simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'feaprescribedmotfinfun_RADIAL_SLOTTED', ...
                                'torquefcn', 'torquefcn_RADIAL_SLOTTED', ...
                                'odeevfun', 'feaprescribedmotodetorquefcn_ROTARY', ...
                                'PoleCount', 30, ...
                                'RampPoles', 1, ...
                                'Rpm', 60 );

% simoptions.ODESim.OutputFcn = 'feaprescribedmotodeoutputfcn_ROTARY';

design.FEAFluxLinkageFCN = @feaode_RADIAL_SLOTTED;

simoptions.reltol = 1e-6;
simoptions.PhaseCurrentTols = repmat(0.001, 1, design.Phases);
% simoptions.maxstep = (simoptions.tspan(2) - simoptions.tspan(1)) / 10000;

simoptions.evaloptions = designandevaloptions_RADIAL_SLOTTED ();

simoptions.MagFEASim.UseParFor = true;

[design, simoptions] = feval(simoptions.simfun, design, simoptions);
simoptions.simfun = [];

[design, simoptions] = feval(simoptions.finfun, design, simoptions);
simoptions.finfun = [];

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions);

%%
%
%fsimoptions = simoptions;
%
%bp = 0.1;
%tmax = 1;
%
%fsimoptions.tspan = [ 0:0.01/100:bp, ...
%                      (bp+0.01/100):(bp+tmax)/1000:(bp+tmax)  ];
%
%fsimoptions.odesolver = 'odef1';
%
%[fT, fY, fresults, fdesign, fsimoptions] = simulatemachine_AM(design, ...
%                                                         fsimoptions);

