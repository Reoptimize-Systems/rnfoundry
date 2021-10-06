%% Test_designAndEvaluate_ACTIAM

design.Phases = 3;         % Number of Phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.CoilFillFactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RlVRp = 10; % Ratio of machine resistance to grid resistance

design.Poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

E = [200e9 151e9 200e9];  % Young's modulus of elasticity for Structural steel and laminated steel respectively

options.targetPower = 1e5; % 100kW machine
options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

% % set up the functions
% simoptions.ODESim.PreProcFcn = @RunStructFEMMSimNew_ACTIAM;
% simoptions.xycoords = randMat([1.00000001*design.Rm; 0], [0.99999*design.Ro; 1.0], 0, 2000); 
% simoptions.ODESim.PostPreProcFcn = @finfun_ACTIAM;
% simoptions.odefun = @simplelinearmachineode_proscribedmotion; 
% simoptions.dpsidxfun = @polypsidot_ACTIAM; %@dpsidx_tubular; 
% simoptions.ODESim.PostSimFcn = @resfun_ACTM;
% 
% % Test with linear motion
% speed = 1;
% simoptions.ODESim.InitialConditions = 0;
% simoptions.ODESim.ResultsTSkip = 1;
% simoptions.ODESim.TimeSpan = [0, 10];
% simoptions.drivetimes = 0:simoptions.ODESim.TimeSpan(2);
% simoptions.vT = repmat(speed, size(simoptions.drivetimes));
% simoptions.xT = simoptions.vT .* simoptions.drivetimes;
% simoptions.Lmode = 0;

% set up the functions
simoptions.ODESim.PreProcFcn = 'simfun_ACTIAM';
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTIAM';
simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostSimFcn = 'systemresfun_linear';

% use buoy number 37, 4m diameter, 2m draft
design.buoynum = 37;

simoptions.BuoySim.SeaParameters.sigma = 2 * pi * 0.35;
simoptions.BuoySim.SeaParameters.phase = pi / 2;
simoptions.BuoySim.SeaParameters = defaultseaparamaters(simoptions.BuoySim.SeaParameters);
simoptions.BuoySim.tether_length = 5;
simoptions.water_depth = 40;

simoptions.ODESim.InitialConditions = [0,0,0];
simoptions.ODESim.ResultsTSkip = 1;
simoptions.ODESim.TimeSpan = [0, 45];
simoptions.Lmode = 0;

%%

[score, design, simoptions, T, Y, results] = designandevaluate_ACTIAM(design, simoptions, options);

%%
% if ~isfemmopen
%     openfemm
% end
% % First convert the ratios into actual dimensions
% [Wp, Wm, Ws, Ri, Ro, Ra, g, Rsi, Rso, Wc, Hc] = ratios2dimensions_ACTIAM(WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, 0, RsoVRm, WcVWp, Rm);
% % Determine the area of the coil
% coilArea = Hc * Wc;
% [Ntot, dc] = CoilTurns(coilArea, kfill, dc);
% Dc = RunFEMMSimWithCoils_ACTIAM(WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, RsoVRm, WcVWp, Rm, Ntot, kfill, [0 0 0], [0 1]);

%%

plotresultsbuoysys_linear(T, Y, results, 1)
