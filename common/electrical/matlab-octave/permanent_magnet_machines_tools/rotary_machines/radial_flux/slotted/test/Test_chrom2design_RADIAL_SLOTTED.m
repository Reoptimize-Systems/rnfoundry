% Test_chrom2design_RADIAL_SLOTTED_nova

clear design simoptions

% Set the simulation speed
simoptions.RPM = 17.0;
% choose the number of poles to cross during the fixed speed sim
simoptions.PoleCount = 300;
% choose the maximum machine outer radius
simoptions.MaxRbo = 0.9;

% the maximum allowed rms current density in a coil
simoptions.max_JCoilRms = 6e6;
simoptions.max_JCoilRms_penfactor = [100, 100];
% the maximum allowed peak current density in the coil
% simoptions.maxAllowedJpeak = 10e6;
% simoptions.maxAllowedJpeakPenFactor = 100;
% the minimum allowed rms voltage produced in a coil
simoptions.min_EMFPhaseRms = 600;
simoptions.min_EMFPhaseRms_penfactor = [100, 10];

simoptions.max_EMFPhaseRms = 700;
simoptions.max_EMFPhaseRms_penfactor = [100, 10];

% the target power rating
simoptions.target_PowerLoadMean = 100e3;
simoptions.target_PowerLoadMean_penfactor = [100, 10];
% the minimum allowed power
% simoptions.minPowerLoadMean = 90e3;
% simoptions.minPowerLoadMeanPenFactor = [1000, 100];
% simoptions.maxPowerLoadMean = 110e3;
% simoptions.maxPowerLoadMeanPenFactor = [10, 0];

% simoptions.maxAllowedDeflectionFactor = 0.1;
% set an absoltute value of the max deflection of 2mm, which supercedes the
% value from maxAllowedDeflectionFactor if it is larger
% simoptions.maxAllowedDeflection = 2/1000;

simoptions.addEfficiencyPenalty = true;
simoptions.EfficiencyPenFactor = [1, 0];

% limit torque ripple to 5% of torque value
simoptions.max_TorqueRippleFactor = 0.1; 
simoptions.max_TorqueRippleFactor_penfactor = [10, 1];

%% Simulation Parameters
% set the other simulation parameters 

% determines method used to calculate inductance
simoptions.Lmode = 1;
% the initial currents in the coils at t=0
simoptions.ODESim.InitialConditions = [0, 0, 0];
% the number of calculations to skip when producing output after the ode
% solver finishes
simoptions.ODESim.ResultsTSkip = 1;
%
simoptions.reltol = 1e-4;

% determines if a preliminary simulation at a fixed speed is performed
% prior to the real evaluation sim to determine if full sim should be
% skipped due to poorness of the design
simoptions.DoPreLinSim = false;

% Determines if multiple simulations are performed to determine the forces
% as the air gap closes, or a single sim, and forces based on this
simoptions.GetVariableGapForce = false;

% Core loss interpolation data, using material M-36 A-S Sheared 26 gage
% laminated magnetic steel
[fq, Bq, Pq] = m36assheared26gagecorelossdata(false);
simoptions.CoreLossData.fq = fq;
simoptions.CoreLossData.Bq = Bq;
simoptions.CoreLossData.Pq = Pq;

% common machine/cost options
% evaloptions.MagnetCost = 80;
% evaloptions.CopperCost = 10;
% evaloptions.FieldIronCost= 4;
% evaloptions.ArmatureIronCost = 6;
% evaloptions.StructMaterialCost = 3;
evaloptions.CapacityFactor = 0.3;
evaloptions.ProjectYears = 25;
evaloptions.DiscountRate = 0.08;

% Machine magnetics mesh sizes
simoptions.MagFEASim.MagnetRegionMeshSize = -1;
simoptions.MagFEASim.BackIronRegionMeshSize = -1;
simoptions.MagFEASim.OuterRegionsMeshSize = [-1, -1];
simoptions.MagFEASim.AirGapMeshSize = -1;

% Machine structural FEA mesh sizes
% simoptions.Evaluation.structmeshoptions.ShaftAxialLayersPerM = 10;
% simoptions.Evaluation.structmeshoptions.DiscAxialLayersPerM = 100;
% simoptions.Evaluation.structmeshoptions.SupportAxialLayersPerM = 10;
% simoptions.Evaluation.structmeshoptions.CircumPointsPerM = 20;
% simoptions.Evaluation.structmeshoptions.BackIronRadialPointsPerM = 10; 
% simoptions.Evaluation.structmeshoptions.MagnetRadialPointsPerM = 10;

% The simulation functions, these are suitible for a prescribed motion
% simulation of a coreless torus machine
simoptions.ODESim.PreProcFcn = 'simfun_RADIAL_SLOTTED';
simoptions.ODESim.PostPreProcFcn = 'prescribedmotfinfun_RADIAL_SLOTTED';
simoptions.ODESim.EvalFcn = 'prescribedmotodetorquefcn_ROTARY'; 
simoptions.ODESim.TorqueFcn = 'torquefcn_ROTARY'; 
simoptions.ODESim.PostSimFcn = 'prescribedmotresfun_ROTARY';
simoptions.MagFEASim.UseFemm = false;
simoptions.MagFEASim.QuietFemm = true;

% set some default spawning settings if not supplied
evaloptions.starttime = [];
evaloptions.endtime = [7,30,0];
evaloptions.maxslaves = 20;
evaloptions.matlicencebuffer = 10;
evaloptions.waitforotherfea = false;
evaloptions.waitforotherode = false;
evaloptions.spawnslaves = false;
evaloptions.MCoreFEADir = '';
evaloptions.MCoreODEDir = '';

simoptions.Evaluation = evaloptions;

%%

% design.Rbo = Chrom(1);
% tyVtm = Chrom(2);
% tcVMax_tc = Chrom(3);
% tsbVMax_tsb = Chrom(4);
% design.tsgVtsb = Chrom(5);
% g = Chrom(6);
% tmVMax_tm = Chrom(7);
% tbiVtm = Chrom(8);
% design.thetamVthetap = Chrom(9);
% design.thetacgVthetas = Chrom(10);
% design.thetacyVthetas = Chrom(11);
% design.thetasgVthetacg = Chrom(12);
% lsVMax_ls = Chrom(13);
% design.NBasicWindings = round(Chrom(14));
% design.DcAreaFac = Chrom(15);
% design.BranchFac = Chrom(16);
% 
% if numel(Chrom) > 16
%     design.MagnetSkew = Chrom(17);
% end

Chrom = [ 320e-3; % Rbo
          0.8; % tyVtm
          0.9; % tcVMax_tc
          0.5; % tsbVMax_tsb
          0.5; % tsgVtsb
          2e-3; % g
          0.1; % tmVMax_tm
          2.0; % tbiVtm
          0.8; % thetamVthetap
          0.8; % thetacgVthetas
          0.8; % thetacyVthetas
          0.4; % thetasgVthetacg
          0.5; % lsVMax_ls
          5; % NBasicWindings
          0.025; % DcAreaFac
          1; % BranchFac
          0.25; % MagnetSkew
          ]; 

[design, tsimoptions] = chrom2design_RADIAL_SLOTTED (simoptions, Chrom);

design.MagFEASimMaterials.Magnet = matstr2matstruct_mfemm ('NdFeB 48M@20C');                                         
ElectricalSteel = matstr2matstruct_mfemm ('M-19 Steel');
BackIronSteel = matstr2matstruct_mfemm ('Cold drawn carbon steel, annealed');
design.MagFEASimMaterials.FieldBackIron = BackIronSteel;
design.MagFEASimMaterials.ArmatureYoke = ElectricalSteel;
design.MagFEASimMaterials.ArmatureCoil = matstr2matstruct_mfemm('36 AWG');
design.MagFEASimMaterials.ArmatureCoil.WireD = design.Dc*1000;
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.CoilInsulation = 'Air';
design.CoilTurns = 1;

FemmProblem = slottedfemmprob_radial (design, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers );
                        
plotfemmproblem (FemmProblem)
                        
dispstruct(design, 25)

%%

Chrom = [0.1541;
    1.7273;
    0.0690 ;
    0.9160 ;
    0.8507  ;
    0.0027 ;
    0.9881;
    4.9468;
    0.7267 ;
    0.2819 ;
    0.2982 ;
    0.6774 ;
    0.1990 ;
    46.6448 ;
    0.9852 ;
    0.0001;
    0.0008 ;
    0.7163 ]

% set a sensible limit on the stack length
simoptions.Max_ls = 40e-3;

% choose the maximum possibe coil height and shoe base height
simoptions.Max_tc = 50e-3;
simoptions.Max_tsb = 10e-3;
simoptions.Max_tm = 10e-3;

simoptions.RlVRp = 1;

% choose the maximum machine outer radius
simoptions.Max_Rbo = 320e-3 / 2;
% choose the minimum air gap
simoptions.Min_g = 1.5e-3;

[design, tsimoptions] = ML0125_221115.chrom2design (simoptions, Chrom);

design.MagFEASimMaterials.Magnet = matstr2matstruct_mfemm ('NdFeB 48M@20C');                                         
ElectricalSteel = matstr2matstruct_mfemm ('M-19 Steel');
BackIronSteel = matstr2matstruct_mfemm ('Cold drawn carbon steel, annealed');
design.MagFEASimMaterials.FieldBackIron = BackIronSteel;
design.MagFEASimMaterials.ArmatureYoke = ElectricalSteel;
design.MagFEASimMaterials.ArmatureCoil = matstr2matstruct_mfemm('36 AWG');
design.MagFEASimMaterials.ArmatureCoil.WireD = design.Dc*1000;
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.CoilInsulation = 'Air';
design.CoilTurns = 1;

FemmProblem = slottedfemmprob_radial (design, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers );
                        
plotfemmproblem (FemmProblem)
                        
dispstruct(design, 25)


%%

Chrom = [128.000000000000e-003;
    4.00000000000000e+000;
    1.00000000000000e+000;
    0.00000000000000e+000;
    788.019787424054e-003;
    5.66342443257644e-003;
    1.00000000000000e+000 ;
    5.00000000000000e+000;
    950.000000000000e-003;
    900.000000000000e-003 ;
    200.000000000000e-003 ;
    103.370323711086e-003;
    11.2757874323161e-003  ;
    40.0188322962846e+000  ;
    0.00000000000000e+000 ;
    1.00000000000000e-003 ;
    121.077012075819e-006 ;
    103.144538484018e-003 ];


% set a sensible limit on the stack length
simoptions.Max_ls = 40e-3;

% choose the maximum possibe coil height and shoe base height
simoptions.Max_tc = 50e-3;
simoptions.Max_tsb = 10e-3;
simoptions.Max_tm = 10e-3;

simoptions.RlVRp = 1;

% choose the maximum machine outer radius
simoptions.Max_Rbo = 320e-3 / 2;
% choose the minimum air gap
simoptions.Min_g = 1.5e-3;

[design, tsimoptions] = ML0125_221115.chrom2design (simoptions, Chrom);

design.MagFEASimMaterials.Magnet = matstr2matstruct_mfemm ('NdFeB 48M@20C');                                         
ElectricalSteel = matstr2matstruct_mfemm ('M-19 Steel');
BackIronSteel = matstr2matstruct_mfemm ('Cold drawn carbon steel, annealed');
design.MagFEASimMaterials.FieldBackIron = BackIronSteel;
design.MagFEASimMaterials.ArmatureYoke = ElectricalSteel;
design.MagFEASimMaterials.ArmatureCoil = matstr2matstruct_mfemm('36 AWG');
design.MagFEASimMaterials.ArmatureCoil.WireD = design.Dc*1000;
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.CoilInsulation = 'Air';
design.CoilTurns = 1;

FemmProblem = slottedfemmprob_radial (design, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers );
                        
plotfemmproblem (FemmProblem)
                        
dispstruct(design, 25)   


%%

Chrom = [144.894567412053e-003 ;
    3.43935474440623e+000 ;
    894.473335060377e-003 ;
    629.715273760603e-003 ;
    1.00000000000000e+000  ;
    2.76375713124924e-003;
    881.657977413938e-003  ;
    3.06650072019637e+000 ;
    532.616862733833e-003 ;
    505.752285814302e-003 ;
    568.624333417970e-003 ;
    532.192351338307e-003;
    699.534337827985e-003 ;
    56.1154105686251e+000 ;
    320.860533716092e-003 ;
    252.622787288614e-006 ;
    309.830319257120e-006 ;
    480.246171651591e-003];


% set a sensible limit on the stack length
simoptions.Max_ls = 40e-3;

% choose the maximum possibe coil height and shoe base height
simoptions.Max_tc = 50e-3;
simoptions.Max_tsb = 10e-3;
simoptions.Max_tm = 10e-3;

simoptions.RlVRp = 1;

% choose the maximum machine outer radius
simoptions.Max_Rbo = 320e-3 / 2;
% choose the minimum air gap
simoptions.Min_g = 1.5e-3;

[design, tsimoptions] = ML0125_221115.chrom2design (simoptions, Chrom);

design.MagFEASimMaterials.Magnet = matstr2matstruct_mfemm ('NdFeB 48M@20C');                                         
ElectricalSteel = matstr2matstruct_mfemm ('M-19 Steel');
BackIronSteel = matstr2matstruct_mfemm ('Cold drawn carbon steel, annealed');
design.MagFEASimMaterials.FieldBackIron = BackIronSteel;
design.MagFEASimMaterials.ArmatureYoke = ElectricalSteel;
design.MagFEASimMaterials.ArmatureCoil = matstr2matstruct_mfemm('36 AWG');
design.MagFEASimMaterials.ArmatureCoil.WireD = design.Dc*1000;
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.CoilInsulation = 'Air';
design.CoilTurns = 1;

FemmProblem = slottedfemmprob_radial (design, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers );
                        
plotfemmproblem (FemmProblem)
                        
dispstruct(design, 25)   