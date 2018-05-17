% Test_feaprescribedmotodeforcefcn_linear.m

clear design simoptions T Y results

% design = test_design_TM_SLOTLESS ('dims');

load ('/home/rcrozier/Sync/work/enercro/Projects/WES-Wavedrive/Outputs/design_004_converted_to_rnfoundry.mat');

simoptions = struct ();

simoptions.set_CoilResistance = 0.497091;

%%
design.RlVRp = 4;

% design.RlVRp = 1e9;

% [design.CoreLoss.kh, ...
%  design.CoreLoss.kc, ...
%  design.CoreLoss.ke, ...
%  design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'InterpolateMissing', false);

simoptions.GetVariableGapForce = false;
simoptions.MagFEASimType = 'multiple';
simoptions.UseParFor = true;
simoptions.MagFEASim.UseParFor = false;
simoptions.AddPhaseCurrentsComponents = true;

design = completedesign_TM_SLOTLESS (design, simoptions);

design.zs = design.zp / 3;

% setup simulation options
%simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
%                                'Rpm', 30, ...
%                                'TSpan', [0,1]);

simoptions = simsetup_linear (design, 'simfun_TM_SLOTLESS', 'feaprescribedmotfinfun_TM_SLOTLESS', ...
                                'ForceFcn', 'forcefcn_linear_pscbmot', ...
                                'EvalFcn', 'feaprescribedmotodeforcefcn_linear', ...
                                'simoptions', simoptions, ...
                                'PoleCount', 15, ...
                                'RampPoles', 5, ...
                                'Velocity', 1.0 );
                            


% simoptions.ODESim.OutputFcn = 'feaprescribedmotodeoutputfcn_ROTARY';

design.FEAFluxLinkageFCN = @feaode_TM_SLOTLESS;

simoptions.ODESim.RelTol = 1e-5;
% simoptions.PhaseCurrentTols = repmat (0.001, 1, design.Phases);
% simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 10000;

simoptions.Evaluation = designandevaloptions_TM_SLOTLESS ();

[design, simoptions] = feval (simoptions.ODESim.PreProcFcn, design, simoptions);
simoptions.ODESim.PreProcFcn = [];

[design, simoptions] = feval (simoptions.ODESim.PostPreProcFcn, design, simoptions);
% simoptions.ODESim.PostPreProcFcn = [];

steptol = 1e-5;

Tind = find (simoptions.xT >= 3*steptol);

initsteps = linspace (simoptions.drivetimes(1), simoptions.drivetimes(Tind(1)), 100);

[xE, vE] = prescribedmotvelpos (initsteps, simoptions);

Tind = find (xE >= 3*steptol);

simoptions.ODESim.InitialStep = min ([simoptions.ODESim.MaxStep, max([simoptions.ODESim.MaxStep/5, initsteps(Tind(1))])]);

%%
close all
simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.ResetFcn ()


% simoptions.ODESim.MaxStep = simoptions.ODESim.MaxStep / 10;

[T, Y, results, design, simoptions] = simulatemachine_AM ( design, ...
                                                           simoptions );

%%
%
%fsimoptions = simoptions;
%
%bp = 0.1;
%tmax = 1;
%
%fsimoptions.ODESim.TimeSpan = [ 0:0.01/100:bp, ...
%                      (bp+0.01/100):(bp+tmax)/1000:(bp+tmax)  ];
%
%fsimoptions.ODESim.Solver = 'odef1';
%
%[fT, fY, fresults, fdesign, fsimoptions] = simulatemachine_AM(design, ...
%                                                         fsimoptions);

