% Test_evaluatestructure_AF
% setup design
clear design simoptions

% number of Phases in the machine
design.Phases = 3; %
% number of parallel branches per phase
design.Branches = 1;
% number of series coils per branch
design.CoilsPerBranch = 44;
% number of stages in the machine
design.NStages = 4;

% coil epoxy thickness
design.tepoxy = 0.001; %
% total air gap (from magnetic perspective)
design.g = 5/1000 + design.tepoxy;  %
% coil height in air gap
design.tc = 0.022; %

% Mean radius of magnets
design.Rmm = 2.86;
% Outer radius of magnets 
design.Rmo = 3.016; %
% Inner radius of magnets
design.Rmi = 2.704;  %
% Inner shaft radius
design.Rs = 2.440;
% Outer shaft radius, make shaft radius plus half a bolt thickness
design.Rbi = design.Rs + 24/1000/2;
% Outer back iron radius
design.Rbo = 1.05*design.Rmo;
% Inner radial height of coil
design.hcoil = 326/1000;
% Inner radius of coil
design.Rci = design.Rmm - design.hcoil/2;
% Outer Radius of coil
design.Rco = design.Rmm + design.hcoil/2;
% total number of magnetic Poles
design.Poles = 176; %

% mean pole pitch
design.taupm = (pi * (design.Rmo + design.Rmi)) / design.Poles;
% mean magnet pitch
design.taumm = 0.076; %
% mean outer coil pitch
design.tauco = 0.133; %
% mean coil former (inner coil) pitch
design.tauci = 0.042; %
% actual coil spacing (total coil spacing is slightly larger that tauco)
design.taupcg = design.Phases * 0.136;


% complete standard design structure values used for calculations of some
% parameters
design.Hc = design.tc;
design.Wc = (design.tauco - design.tauci) / 2;

% area of one of the wires used in prototype made up of two strands of
% rectangular wire
strandarea = 2 * 0.02178 * 1.96E-03; %
design.CoilTurns = 20; %
% calculate coil wire diameter to give cross-sectional area same as actual
% non-round wire cross-sectional area
design.Dc = 2 * area2radius(strandarea); %
% determine the copper fill factor for completeness
design.CoilFillFactor = (design.CoilTurns * strandarea) / (design.Hc * design.Wc); %
% thickness of the magnet
design.tm = 0.015; %
% outer and inner plate thickness
design.tbi = [0.023, 0.012];% from spreadsheet [0.020, 0.012]; %
% Number of coils per phase
design.NPhaseCoils = 44; %
% number of modules
design.NModules = round(176 / 8);
% number of module ribs
design.NModuleSupports = 4;
% size of module rib base in z direction
design.tsuppb = 55/1000;
% width, or span, or a module ribs
design.tausupp = 40/1000;

% tdiscsep = (2*design.g) + design.tc;
% Ratio of load resistance to phase resistance
design.RlVRp = 10;
% Ratio of load inductance to phase inductance
design.LgVLc = 0;
% Resistivity of wire at 20 deg Celcius
design.Rho20Wire = 1.7e-8;
% Temperature coefficient of wire resistivity
design.AlphaRho20Wire = 3.93e-3;

% initially insert 32 MGOe magnet material from library
design.MagFEASimMaterials.Magnet = 'NdFeB 32 MGOe';
% but force a change to the coercive value to be the same as protype spec
% from data sheet (nominal Hcb for Sintered neodymium-iron-boron NdFeB
% magnet N38M)
design.HcMag = 915e3;
% design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
% don't know what material is back iron yet, use an electrical (Silicon
% Iron) steel M-27 Steel from FEMM library until known
design.MagFEASimMaterials.FieldBackIron = 'M-27 Steel';
% use any wire type, is not used in resistance calc, use wire diameter
% specified in design.Dc for that
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';

% setup simulation options
simoptions = simsetup_ROTARY(design, 'simfun_TORUS_CORELESS', 'finfun_TORUS_CORELESS', ...
                                'Velocity', 1, 'TSpan', [0,10]);
                            
simoptions.reltol = 1e-4;
simoptions.abstol = repmat(0.001, 1, design.Phases);
simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 10000;

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions);

simoptions = setfieldifabsent(simoptions, 'evaloptions', []);

simoptions.evaloptions = designandevaloptions_TORUS_CORELESS(simoptions.evaloptions);
    
[maxzdef, design] = evaluatestructure_AF(design, simoptions);


