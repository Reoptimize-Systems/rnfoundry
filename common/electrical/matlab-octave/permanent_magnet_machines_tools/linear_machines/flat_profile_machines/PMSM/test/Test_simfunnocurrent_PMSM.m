% Test_simfunnocurrent_PMSM

clear design simoptions

design.Phases = 3;
design.Wp = 0.12;
design.Wm = 0.8*design.Wp;
design.hm = 0.015;
design.kw = 0.84;
%Ns = 470;
%Ns = 6; % ???
design.Hc = 979000;
%ht = 0.1;
design.CoilFillFactor = 0.585;
design.g = 0.003; 
design.ls = 1; 
design.Dc = 0.005; % gives Ns of 468
design.E = [200e9 151e9];
design.Ws=design.Wp/design.Phases;
design.Wt=design.Ws/2;
design.ht=5*design.Wt;
design.hbf = design.hm;
design.hba = design.hbf;

% Number of Poles
design.Poles = [10 40];
design.CoilTurns = 400;

design.RgVRc = 10;
design.LgVLc = 0;

% Get the dimensionless ratios from the parameters
design = dimensions2ratios_PMSM(design);

% set design mode, double-sided, double-layered (half of slot filled by
% coil)
design.mode = [1, 1, 0, 1];

[design, simoptions] = simfunnocurrent_PMSM(design, struct());

figure; plot([design.FEAFx(:), design.FEAFy(:)]);

figure; plot(design.gvar, design.gforce);

figure; plot(design.psi);


