%% ACTM

clear design

design.Phases = 3;         % Number of phases in machine
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.5;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.025;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.CoilFillFactor = 0.65;
%design.Dc = 1/1000;  % 1 mm diameter wire 
design.CoilTurns = 50;
design.mode = 2; 
design.LgVLc = 0;
design.Poles = [10 30];
design.CoilsPerBranch = 10;
design.Branches = 1;
% design.FieldDirection = 1;
% design.PowerPoles = poles(1);

design = ratios2dimensions_ACTM(design);



%% Set up Common Parameters

design.RlVRp = 10;

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;
simoptions.maxAllowedxT = 0.5;

%% Test with buoy in sinusoidal sea

% % Set up the buoy and sea data files, these are for the 2m buoy
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
% simoptions.HeaveFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'heave_coefficients_cyl_2di_1dr_d020610.mat');
% simoptions.SurgeFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'surge_coefficients_cyl_2di_1dr.mat');
% simoptions.HydroCoeffsFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.1');
% simoptions.ExcitationFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.2');
% simoptions.BuoySim.BuoyParameters = load(fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'buoyparams_d3103v4.mat'));
% simoptions.buoynum = -1;
% design.buoynum = simoptions.buoynum;
% 
% simoptions.ODESim.TimeSpan = [0, 60];
% % params.amp = 1;
% params.sigma = 2 * pi * 0.35;
% params.phase = pi/2;
% 
% simoptions.BuoySim.SeaParameters = defaultseaparamaters(params);
% 
% simoptions.BuoySim.tether_length = 4;
% 
% simoptions.ODESim.EvalFcn = 'systemode_linear'; 
% simoptions.ODESim.PostPreProcFcn = ['systemfinfun_', mname];
% simoptions.ODESim.PostSimFcn = 'systemresfun_linear';    
% 
% [T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.ODESim.PreProcFcn, ...
%                                                   simoptions.ODESim.PostPreProcFcn, simoptions.ODESim.EvalFcn, simoptions.ODESim.PostSimFcn); 
% 
% plotresultsbuoysys_linear(T, Y, results, design, 1)


%% Test with buoy in random sea

% Set up the buoy and sea data files, these are for the 2m buoy
%snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
% simoptions.HeaveFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'heave_coefficients_cyl_2di_1dr_d020610.mat');
% simoptions.SurgeFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'surge_coefficients_cyl_2di_1dr.mat');
% simoptions.HydroCoeffsFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.1');
% simoptions.ExcitationFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.2');
% simoptions.BuoySim.BuoyParameters = load(fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'buoyparams_d3103v4.mat'));
% simoptions.buoynum = -1;

simoptions.buoy = 37;
% design.buoynum = simoptions.buoynum;

simoptions.ODESim.TimeSpan = [0, 60];
% params.amp = 1;

simoptions.BuoySim.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.ODESim.PreProcFcn = 'systemsimfun_ACTM';
mname = 'ACTM';
simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = ['systemfinfun_', mname];
simoptions.ODESim.PostSimFcn = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

% params.peak_freq = 1/9; % centred at resonant frequency
% params.sigma_range = [0.345575191894877,2.31745966692415;];
% params.water_depth = 50;

simoptions.BuoySim.SeaParameters = seasetup ('PMPeakFreq', 1/9);

Nbuoys = 7;
circlerad = 30;
origphase = simoptions.BuoySim.SeaParameters.phase;

phi = cricwavephasediffs(simoptions.BuoySim.SeaParameters.L, circlerad, Nbuoys);

buoyangle = linspace(0, 2*pi, Nbuoys+1);

plot(circlerad .* cos(buoyangle(1:end)), circlerad .* sin(buoyangle(1:end)), 'x');

T = cell(1,5);

results = cell(1,5);

outdesign = cell(1,5);

interpdur = inf;

[design, simoptions] = feval (simoptions.ODESim.PreProcFcn, design, simoptions);
simoptions.ODESim.PreProcFcn = [];
for i = 1:Nbuoys

    disp(i)
    
    simoptions.BuoySim.SeaParameters.phase = origphase + phi(i, :);
    
    [T{i}, Y, results{i}, outdesign{i}, outsimoptions{i}] = simulatemachine_linear (design, simoptions); 
    
    interpdur = min(mean(T{i}(2:end) - T{i}(1:end-1))/2, interpdur);
    
end

%%
interpT = 0:interpdur:simoptions.ODESim.TimeSpan(2);

Forces = zeros(numel(interpT), Nbuoys);

for i = 1:Nbuoys
    
    Forces(:,i) = interp1(T{i}, results{i}.Ffea_heave, interpT);
    
end

plotresultsbuoysys_linear(T{end}, Y, results{end}, outdesign{end}, outsimoptions{end}, 5)

% figure; plot(interpT, Forces);

figure; plot(interpT, Forces, ':', interpT, sum(Forces, 2), '-');

 
 
 