%% Test_buoys

clear 
% First we need some plausible machine variables for testing, we will use
% the AWS variables for this

design.phases = 3;         % Number of phases in machine
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.2;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.025;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.fillfactor = 0.65;
%design.Dc = 1/1000;  % 1 mm diameter wire 
design.Ntot = 500;
design.mode = 2; 
design.LgVLc = 0;
design.poles = [15 30];
% Acon = pi * (dc/2)^2;
% Acu = (Rm * WpVRm * WcVWp) * ((RoVRm * Rm) - Ri) * kfill;
% AconVAcu = Acon / Acu;

% E = [200e9 151e9];  % Young's modulus of elasticity for Structural steel
% and laminated steel respectively
design.RgVRc = 10; % Ratio of grid resistance to machine resistance

options.E = [207e9 100e5 207e9];
options.targetPower = 10e3; % 10kW machine
options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

% Chose the method used to evaluate inductance
simoptions.Lmode = 1;

simoptions.buoylibdir = buoylibdir ();

simoptions.BuoySim.SeaParameters.peak_freq = 1/9;
simoptions.BuoySim.SeaParameters.sigma_range = [2*pi*0.055, 0.1*(0.6)^0.5 + 2.24];
simoptions.BuoySim.SeaParameters.freqcount = 50;
simoptions.BuoySim.SeaParameters.water_depth = 50;

simoptions.BuoySim.SeaParameters = defaultseaparamaters(simoptions.BuoySim.SeaParameters);

simoptions.BuoySim.tether_length = 10;

simoptions.ODESim.TimeSpan = [0, 500];

% set up the functions
simoptions.ODESim.PreProcFcn = 'dummysimfun_ACTM';
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTM';
simoptions.ODESim.EvalFcn = 'systemode_ACTM';
simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
simoptions.ODESim.PostSimFcn = 'splitsystemresfun_ACTM';
simoptions.ODESim.Split = 5;
simoptions.ODESim.SplitPointFcn = 'splitodesystemres_ACTM';
simoptions.events = 'vevents_ACTM';

design = Ratios2Dimensions_ACTM(design);

design = systemsimfun_ACTM(design);

[design, simoptions] = finfun_ACTM(design, simoptions);

%%

numbuoys = buoynum2buoydata(simoptions.buoylibdir);
buoystrs = {};
buoyres = [];

%%
seatests = 5;
simoptions = buoysimsetup(1, simoptions);
simoptionsArr = repmat(simoptions, 1, seatests);

diary('buoy_test_output.log')

for j = 1:seatests
    
    simoptionsArr(j).SeaParameters = defaultseaparamaters(simoptions.BuoySim.SeaParameters);
    
    for i = 1:numbuoys

        % Set up buoy simulation
        simoptionsArr(j) = buoysimsetup(i, simoptionsArr(j));

        try
            
            fprintf(1, '\nBeginning test of buoy %d', i);
            [score, design] = designAndEvaluate_ACTM(design, simoptionsArr(j), options);
            fprintf(1, '\nBuoy %d test complete', i);
            buoystrs{i,j} = simoptionsArr(j).HeaveFile;
            buoyres(i,j) = 0;

        catch ME

            if ~isempty(strfind(ME.identifier, 'nomem'))
                %out of memory
                fprintf(1, '\nBuoy %d was dodgy, ran out of memory', i);
                buoystrs{i,j} = simoptionsArr(j).HeaveFile;
                buoyres(i,j) = 1;
            else
                %An unexpected error happened
                buoystrs{i,j} = simoptionsArr(j).HeaveFile;
                buoyres(i,j) = 2;
            end

        end

    end

end

diary off

%% Use multicore

actmrootpath = fileparts(which('GetEMF_ACTM'));
       
% Determines whether the master performs work or only coordinates
settings.masterIsWorker    = true;
% This is the number of function evaluations given to each worker
% in a batch
settings.nrOfEvalsAtOnce   = 1;
% The maximum time a single evaluation should take, determines
% the timeout for a worker
settings.maxEvalTimeSingle = 60*60;
% Determines whether a wait bar is displayed, 0 means no wait bar
settings.useWaitbar = 0;
% Post processing function info
settings.postProcessHandle   = '';
settings.postProcessUserData = {};

settings.debugMode = 0;
settings.showWarnings = 1;

settings.multicoreDir = fullfile(actmrootpath, 'Temp', 'ODE');

for i = 1:numbuoys
    
    parameterCell{i,1} = {i, design, simoptions, options};

end

out = startmulticoremaster2(@testbuoy, parameterCell, settings);
  
%%

temp = dir(fullfile(simoptions.buoylibdir, 'cyl_*'));

%%

outmat = cell2mat(out);

badbuoys = find(outmat==1);
goodbuoys = find(outmat==0);

badtemp = temp(badbuoys);

for i = 1:numel(badtemp)
    fprintf(1, '\n%s', badtemp(i).name);
end

%%

simoptions = rmfield(simoptions, 'splitode');
simoptions.ODESim.PreProcFcn = @systemsimfun_ACTM;
simoptions.ODESim.PostPreProcFcn = @systemfinfun_ACTM;
simoptions.odefun = @systemode_ACTM; 
simoptions.dpsidxfun = @polypsidot_ACTIAM; %@dpsidx_tubular; 
simoptions.ODESim.PostSimFcn = @systemresfun_ACTM;
simoptions.events = @vevents_ACTM;
buoynum = 1;
simoptions = buoysimsetup(buoynum, simoptions);

fprintf(1, '\nBeginning test of buoy %d', buoynum);
[score, design, T, Y, results] = designAndEvaluate_ACTM(design, simoptions, options);
fprintf(1, '\nBuoy %d test complete', buoynum);


%%

plotresultsbuoysys_linear(T, Y, results, 1)

%%

% fps = 30;
% tscale = 1;
% frames = tscale*round((max(T)-min(T))*fps);
% animatesytem_tubular(design, simoptions, T, Y, results, frames, fps, 'supergen_ACTM_mono_waves.avi')

%%
% if ~isfemmopen
%     openfemm;
%     %main_minimize;
% end
% 
% RunFEMMSimWithCoils_ACTM_Dc_Value = RunFEMMSimWithCoils_ACTM(design.WmVWp, design1.WpVRm, design1.RiVRm, design1.RoVRm, design1.RsoVRm, design1.WcVWp, design1.Rm, design1.Ntot, design1.fillfactor, [0 0 0], [0 1]);
% 

%% Using saved values

load splitodesystemres_ACTM.mat

simoptions.ODESim.TimeSpan = [0, 500];

% set up the functions
simoptions.ODESim.PreProcFcn = 'dummysimfun_ACTM';
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTM';
simoptions.ODESim.EvalFcn = 'systemode_ACTM';
simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
simoptions.ODESim.PostSimFcn = 'splitsystemresfun_ACTM';
simoptions.ODESim.Split = 5;
simoptions.ODESim.SplitPointFcn = 'splitodesystemres_ACTM';
simoptions.events = 'vevents_ACTM';

actmrootpath = fileparts(which('GetEMF_ACTM'));
       
% Determines whether the master performs work or only coordinates
settings.masterIsWorker    = true;
% This is the number of function evaluations given to each worker
% in a batch
settings.nrOfEvalsAtOnce   = 1;
% The maximum time a single evaluation should take, determines
% the timeout for a worker
settings.maxEvalTimeSingle = 600*60;
% Determines whether a wait bar is displayed, 0 means no wait bar
settings.useWaitbar = 0;
% Post processing function info
settings.postProcessHandle   = '';
settings.postProcessUserData = {};

settings.debugMode = 0;
settings.showWarnings = 1;

settings.multicoreDir = fullfile(actmrootpath, 'Temp', 'ODE');

numbuoys = buoynum2buoydata(simoptions.buoylibdir);

options.E = [207e9 100e5 207e9];
options.targetPower = 10e3; % 10kW machine
%options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

for i = 1:numbuoys
    
    simoptions = buoysimsetup(i, simoptions);
    
    parameterCell{i,1} = {i, design, simoptions, options};

end

out = startmulticoremaster2(@testbuoy, parameterCell, settings);

outmat = cell2mat(out);

badbuoys = find(outmat==3);
goodbuoys = find(outmat==0);

badtemp = temp(badbuoys);

for i = 1:numel(badtemp)
    fprintf(1, '\n%s', badtemp(i).name);
end

%% Now repeat with different sea

simoption = rmfield(simoptions, 'SeaParameters');
simoptions.BuoySim.SeaParameters.peak_freq = 1/9;
simoptions.BuoySim.SeaParameters.sigma_range = [2*pi*0.055, 0.1*(0.6)^0.5 + 2.24];
simoptions.BuoySim.SeaParameters.freqcount = 55;
simoptions.BuoySim.SeaParameters.water_depth = 50;

simoptions.BuoySim.SeaParameters = defaultseaparamaters(simoptions.BuoySim.SeaParameters);

for i = 1:numbuoys
    
    simoptions = buoysimsetup(i, simoptions);
    
    parameterCell{i,1} = {i, design, simoptions, options};

end

out2 = startmulticoremaster2(@testbuoy, parameterCell, settings);

outmat2 = cell2mat(out2);

badbuoys2 = find(outmat2~=0);
goodbuoys2 = find(outmat2==0);

badtemp2 = temp(badbuoys2);

for i = 1:numel(badtemp2)
    fprintf(1, '\n%s', badtemp2(i).name);
end

%%

outmat3 = cell2mat(out3);

badbuoys3 = find(outmat3~=0);
goodbuoys3 = find(outmat3==0);

badtemp3 = temp(badbuoys3);

for i = 1:numel(badtemp3)
    fprintf(1, '\n%s', badtemp3(i).name);
end
