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
simoptions.IC = [0, 0, 0];
% the number of calculations to skip when producing output after the ode
% solver finishes
simoptions.skip = 1;
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
simoptions.femmmeshoptions.MagnetRegionMeshSize = -1;
simoptions.femmmeshoptions.BackIronRegionMeshSize = -1;
simoptions.femmmeshoptions.OuterRegionsMeshSize = [-1, -1];
simoptions.femmmeshoptions.AirGapMeshSize = -1;

% Machine structural FEA mesh sizes
% simoptions.evaloptions.structmeshoptions.ShaftAxialLayersPerM = 10;
% simoptions.evaloptions.structmeshoptions.DiscAxialLayersPerM = 100;
% simoptions.evaloptions.structmeshoptions.SupportAxialLayersPerM = 10;
% simoptions.evaloptions.structmeshoptions.CircumPointsPerM = 20;
% simoptions.evaloptions.structmeshoptions.BackIronRadialPointsPerM = 10; 
% simoptions.evaloptions.structmeshoptions.MagnetRadialPointsPerM = 10;

% The simulation functions, these are suitible for a prescribed motion
% simulation of a coreless torus machine
simoptions.simfun = 'simfun_RADIAL_SLOTTED';
simoptions.finfun = 'prescribedmotfinfun_RADIAL_SLOTTED';
simoptions.odeevfun = 'prescribedmotodetorquefcn_ROTARY'; 
simoptions.torquefcn = 'torquefcn_ROTARY'; 
simoptions.resfun = 'prescribedmotresfun_ROTARY';
simoptions.usefemm = false;
simoptions.quietfemm = true;

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

simoptions.evaloptions = evaloptions;

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
