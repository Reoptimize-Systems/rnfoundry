% Test_ddlinearsystemsim


% First set the machine physical variables
design.Phases = 3;         
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.75;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.05;
design.RsoVRm = 0.5;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.CoilFillFactor = 0.65;

design = ratios2dimensions_ACTIAM(design);

% set the number of turns or the wire diameter or both
%design.Dc = 1/1000;
design.CoilTurns = 500;
% set the mode of the design
design.mode = 2;
% Ratio of grid resistance to machine resistance
design.RlVRp = 10; 
% ratio of grid inductance to machine inductance
design.LgVLc = 0;
% set the number of Poles in each part to 1 as we will be multiplying up
% the Poles to get the required power specified in optins.targetpower
design.Poles = [1 1];

simoptions.Lmode = true;

simoptions.tether_length = 3;

simoptions.ODESim.InitialConditions = [0, 0, 0];

simoptions.tspan = [0, 30];

% set up the functions
simoptions.simfun = @RunStructFEMMSimNew_ACTIAM;
simoptions.finfun = @finfun_ACTIAM;
simoptions.resfun = @buoysysresfun_ACTIAM;

Bparams = load('buoyparams_d3103v4.mat');

simfunction = @simulatebuoysys_TM;
% simulate the system
[T, Y, results, design, params] = ddlinearsystemsim(design, simoptions, simfunction, 'BuoyParameters', Bparams);

%% Useful code snippets

polypsidot_TM(Mdesign, -Mdesign.Wp / 2)

RunFEMMSimWithCoils_ACTIAM(Mdesign.WmVWp, Mdesign.WpVRm, Mdesign.RiVRm, Mdesign.RoVRm, Mdesign.RaVRo,  ...
                                               Mdesign.RsoVRm, Mdesign.WcVWp, Mdesign.Rm, Mdesign.CoilTurns, Mdesign.CoilFillFactor, ...
                                               [0 1e6 0], [2, 1, 1]);
                                           
-yforce_TM(design, 1e6, 0.5*design.Wp);  

4.77517 * 2 * pi / (2 * Mdesign.Wp)   


RunFEMMSimWithCoils_ACTIAM(Mdesign.WmVWp, Mdesign.WpVRm, Mdesign.RiVRm, Mdesign.RoVRm, Mdesign.RaVRo,  ...
                                               Mdesign.RsoVRm, Mdesign.WcVWp, Mdesign.Rm, Mdesign.CoilTurns, Mdesign.CoilFillFactor, ...
                                               fliplr(Jz'), [2, 1, 1]);
                                           
                                           
scatter3(xycoords(:,1), xycoords(:,2), p(:,2));       


RunFEMMSimWithCoils_ACTIAM(design.WmVWp, design.WpVRm, design.RiVRm, design.RoVRm, design.RaVRo,  ...
                                               design.RsoVRm, design.WcVWp, design.Rm, design.CoilTurns, design.CoilFillFactor, ...
                                               [0 1e6 0], [2, 1, 0]);
                                           
 [tintBx, tintBy] = coilfluxdensityintegral_TM(design, 0.5);                                           
                                           


