% Test_simulatemachine_0

clear

design.WmVWp = 0.5; 
design.WpVRm = 0.3;
design.RsoVRm = 0.05;
design.RsiVRso = 0.8;
design.RsiVRm = design.RsiVRso * design.RsoVRm;
design.RiVRm = 1.02;
design.RoVRi = 1.2;
design.RoVRm = design.RiVRm * design.RoVRi;
design.WcVWp = 1/3;
design.RlVRp = 1;
design.LgVLc = 0.1;
design.CoilTurns = 1000;
design.CoilFillFactor = 0.55;
design.Rm = 0.15; 
design.Wp = design.WpVRm * design.Rm;
design.Ro = design.RoVRm .* design.Rm;
design.Poles = 50;
design.mode = 2; 

% Test with linear motion
speed = 1;
simoptions.IC = 0;
simoptions.skip = 1;
simoptions.tspan = [0, 10];
simoptions.drivetimes = 0:simoptions.tspan(2);
simoptions.vT = repmat(speed, size(simoptions.drivetimes));
simoptions.xT = simoptions.vT .* simoptions.drivetimes;
simoptions.Lmode = 0;
%simoptions.maxstep = (design.Wp / 4.1) / speed;

%% ACTM

% set up the functions
simfun = @RunStructFEMMSimNew_ACTM;
xycoords = randMat([1.00000001*design.Rm; 0], [2.9999999*design.Rm; 1.0], 0, 2000); 
finfun = @finfun_ACTM;
odefun = @simplelinearmachineode_proscribedmotion; 
dpsidxfun = @polypsidot_ACTM;%@dpsidx_tubular; 
resfun = @resfun_ACTM;

[T, Y, results, design] = simulatemachine_tubular(design, simoptions, simfun, xycoords,...
                                          finfun, odefun, dpsidxfun, resfun);

%% ACTIAM

design.RaVRo = 1.05;

% set up the functions
simfun = @RunStructFEMMSimNew_ACTIAM;
xycoords = randMat([1.00000001*design.Rm; 0], [0.99999*design.Ro; 1.0], 0, 2000); 
finfun = @finfun_ACTIAM;
odefun = @simplelinearmachineode_proscribedmotion; 
dpsidxfun = @dpsidx_tubular; 
resfun = @resfun_ACTM;

[T, Y, results, design] = simulatemachine_tubular(design, simoptions, simfun, xycoords,...
                                                  finfun, odefun, dpsidxfun, resfun);
                                      
%% ACPMSM

clear

design.dgVlm = 2; 
design.bpVlm = 5; 
design.taupVbp = 1.1; 
design.lsVbp = 5;
design.dbiVlm = 1;
design.WcVtaup = 1/3;
design.hcV2dg = 0.95;
design.lm = 0.02;
design.CoilTurns = 1000;
design.kcufill = 0.55;

[design.Taup, design.bp, design.dg, design.ls, design.dbi, design.Wc, design.Hc] = ratios2dimensions_ACPMSM(design.dgVlm, design.bpVlm, design.taupVbp, design.lsVbp, design.dbiVlm, design.lm, design.WcVtaup, design.hcV2dg);

design.Wp = design.Taup;

% set up the functions
simfun = @RunStructFEMMSim_ACPMSM;
xycoords = randMat([1.000001*(design.lm + design.dbi); 0], [0.99999*((2*design.dg) - design.lm); 1.0], 0, 2000); 
finfun = @finfun_ACPMSM;
odefun = @simplelinearmachineode_proscribedmotion; 
dpsidxfun = @dpsidx_tubular; 
resfun = @resfun_ACTM;

[T, Y, results, design] = simulatemachine_tubular(design, simoptions, simfun, xycoords,...
                                                  finfun, odefun, dpsidxfun, resfun);
     
                                      
                                      


