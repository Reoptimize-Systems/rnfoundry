% Test_simulatemachine_TORUS_CORELESS


clear design simoptions

% setup design
design.g = 5/1000; 

design.tc = 0.1;

design.Rmo = 10;
design.Rmi = 9.5;  
Npoles = 400;
design.NPhaseCoils = Npoles;
design.RlVRp = 0.1;
design.LgVLc = 0;
design.Phases = 3;

design.taupm = (pi * (design.Rmo + design.Rmi)) / Npoles;
design.taumm = 0.85 * design.taupm;
design.tausm = design.taupm / 3;
design.tauco = 0.95 * design.tausm;
design.tauci = 0.05 * design.tauco;
design.taupcg = design.Phases * 2 * design.tausm;
design.Dc = design.taumm / 1000;
design.CoilFillFactor = 0.8;
design.Rco = design.Rmo + 0.01*design.Rmo;
design.Rci = 0.99*design.Rmi;
design.Hc = design.tc;
design.Wc = (design.tauco - design.tauci) / 2;
[design.CoilTurns, design.Dc] = CoilTurns(design.Hc * design.Wc, design.CoilFillFactor, design.Dc);

design.tm = 0.15 * design.taumm;
design.tbi = 0.05;

simoptions.GetVariableGapForce = false;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagSimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagSimMaterials.FieldIron = '1117 Steel';
design.MagSimMaterials.CoilWinding = '36 AWG';

% setup simulation options
simoptions = simsetup_ROTARY(design, 'simfun_TORUS_CORELESS', 'finfun_TORUS_CORELESS', ...
                                'Velocity', 1, ...
                                'TSpan', [0,10], ...
                                'odeevfun', 'prescribedmotodeforcefcn_linear', ...
                                'forcefcn', 'lossforces_TORUS_CORELESS');

                            
simoptions.reltol = 1e-4;
simoptions.PhaseCurrentTols = repmat(0.001, 1, design.Phases);
simoptions.maxstep = (simoptions.tspan(2) - simoptions.tspan(1)) / 10000;

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions, ...
                                                         simoptions.simfun, ...
                                                         simoptions.finfun, ...
                                                         simoptions.odeevfun, ...
                                                         simoptions.resfun);

                                                     