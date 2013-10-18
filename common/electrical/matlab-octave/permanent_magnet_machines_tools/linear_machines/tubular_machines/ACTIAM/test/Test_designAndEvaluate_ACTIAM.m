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
% simoptions.simfun = @RunStructFEMMSimNew_ACTIAM;
% simoptions.xycoords = randMat([1.00000001*design.Rm; 0], [0.99999*design.Ro; 1.0], 0, 2000); 
% simoptions.finfun = @finfun_ACTIAM;
% simoptions.odefun = @simplelinearmachineode_proscribedmotion; 
% simoptions.dpsidxfun = @polypsidot_ACTIAM; %@dpsidx_tubular; 
% simoptions.resfun = @resfun_ACTM;
% 
% % Test with linear motion
% speed = 1;
% simoptions.IC = 0;
% simoptions.skip = 1;
% simoptions.tspan = [0, 10];
% simoptions.drivetimes = 0:simoptions.tspan(2);
% simoptions.vT = repmat(speed, size(simoptions.drivetimes));
% simoptions.xT = simoptions.vT .* simoptions.drivetimes;
% simoptions.Lmode = 0;

% set up the functions
simoptions.simfun = 'simfun_ACTIAM';
simoptions.finfun = 'systemfinfun_ACTIAM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.resfun = 'systemresfun_linear';

% use buoy number 37, 4m diameter, 2m draft
design.buoynum = 37;

simoptions.SeaParameters.sigma = 2 * pi * 0.35;
simoptions.SeaParameters.phase = pi / 2;
simoptions.SeaParameters = defaultseaparamaters(simoptions.SeaParameters);
simoptions.tether_length = 5;
simoptions.water_depth = 40;

simoptions.IC = [0,0,0];
simoptions.skip = 1;
simoptions.tspan = [0, 45];
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
