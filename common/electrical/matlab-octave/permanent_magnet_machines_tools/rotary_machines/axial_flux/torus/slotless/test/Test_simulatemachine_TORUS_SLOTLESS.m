% Test_simulatemachine_TORUS_SLOTLESS

clear design simoptions
% setup design
% setup design
design.g = 5/1000; 

design.tc = 0.02744;

design.hm = 0.3182;
design.Rmm = 1.005;
design.Rmo = design.Rmm + design.hm/2;
design.Rmi = design.Rmm - design.hm/2;  
Npoles = 28;

design.taupm = 0.226;
design.taumm = 0.15;
design.tauco = design.taupm * 0.3;
design.taupcg = design.taupm;
design.Wc = design.tauco;
design.Dc = design.taumm / 1000;
design.CoilFillFactor = 0.8;
design.Hc = design.tc;
design.CoilTurns =  75;
design.Dc = 6.04 / 1000;
design.tm = 0.15 * design.taumm;
design.tbi = 0.045;
design.ty = 2 * design.tbi;

[design.CoilTurns, design.Dc] = CoilTurns(design.Hc * design.Wc, design.CoilFillFactor, design.Dc);

design.NPhaseCoils = Npoles;
design.RgVRc = 0.1;
design.LgVLc = 0;
design.Phases = 3;
design.Branches = 7;
design.CoilsPerBranch = 4;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.CoilMaterial = '36 AWG';

% setup simulation options
simoptions = simsetup_ROTARY(design, 'simfun_TORUS_SLOTLESS', 'finfun_TORUS_SLOTLESS', ...
                                'Velocity', 1, ...
                                'TSpan', [0,10], ...
                                'odeevfun', 'prescribedmotodeforcefcn_linear', ...
                                'forcefcn', 'lossforces_TORUS_SLOTLESS', ...
                                'forcefcnargs', {});
                                    
simoptions.reltol = 1e-4;
simoptions.abstol = repmat(0.001, 1, design.Phases);
simoptions.maxstep = (simoptions.tspan(2) - simoptions.tspan(1)) / 10000;

% add core loss interpolation data
[fq, Bq, Pq] = m19corelossdata();
simoptions.CoreLossData.fq = fq;
simoptions.CoreLossData.Bq = Bq;
simoptions.CoreLossData.Pq = Pq;

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions, ...
                                                         simoptions.simfun, ...
                                                         simoptions.finfun, ...
                                                         simoptions.odeevfun, ...
                                                         simoptions.resfun);

                                                     