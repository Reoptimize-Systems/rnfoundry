% Test_objactm

path = fileparts(which('GetEMF_ACTM'));

cd(fullfile(path, 'Optimisations'));

% E = [200e9 151e9];  % Young's modulus of elasticity for Structural steel
% and laminated steel respectively

options.E = [207e9 100e5 207e9];
options.targetPower = 1e5; % 100kW machine
options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

% Test with linear motion
% speed = 1;
% simoptions.ODESim.InitialConditions = 0;
% simoptions.ODESim.ResultsTSkip = 1;
% simoptions.ODESim.TimeSpan = [0, 5];
% simoptions.drivetimes = 0:simoptions.ODESim.TimeSpan(2);
% simoptions.vT = repmat(speed, size(simoptions.drivetimes));
% simoptions.xT = simoptions.vT .* simoptions.drivetimes;
% simoptions.Lmode = 1;

% set up the functions
% simoptions.ODESim.PreProcFcn = @RunStructFEMMSimNew_ACTM;
% simoptions.xycoords = randMat([1.00000001*design.Rm; 0], [2.9999999*design.Rm; 1.0], 0, 2000); 
% simoptions.ODESim.PostPreProcFcn = @finfun_ACTM;
% simoptions.odefun = @simplelinearmachineode_proscribedmotion; 
% simoptions.dpsidxfun = @dpsidx_tubular; 
% simoptions.ODESim.PostSimFcn = @resfun_ACTM;

simoptions.desiredMeanPower = 100*3;
simoptions.maxAllowedJrms = 5e6;
simoptions.desiredEMF = 400;

mpgaoptions.OBJ_F = 'objactm';
mpgaoptions.TERMEXACT = 1e-8;
mpgaoptions.SEL_F = 'sus';
mpgaoptions.XOV_F = 'recint';
mpgaoptions.MUT_F = 'mutbga';
mpgaoptions.SUBPOP = 3;
mpgaoptions.NIND = 20;
mpgaoptions.MAXGEN = 600;
mpgaoptions.MUTR = 0.2;
mpgaoptions.DISPLAYMODE = 2;
mpgaoptions.SAVEMODE = 1;
mpgaoptions.filename = 'objobjactm_output_20100526.mat';
%mpgaoptions.resumefile = 'objobjactm_output.mat';

%%

[Best, IndAll, Chrom, ObjV] = mpgafun2(mpgaoptions, 'ObjectiveArgs', {simoptions, options});


%%

objactm(BestInd(1:end-1), [], simoptions, options)
