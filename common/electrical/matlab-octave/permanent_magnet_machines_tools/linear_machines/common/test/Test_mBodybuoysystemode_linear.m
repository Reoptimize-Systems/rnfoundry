%% ACTM

clear design simoptions

design = test_design_ACTM ();

design = ratios2dimensions_ACTM(design);

% Set up Common Parameters

design.RlVRp = 10;

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;
simoptions.maxAllowedxT = 0.5;

% Test with buoy in sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy
simoptions.BuoySim.buoy = 'cyl_2dia_1dr';

simoptions.ODESim.TimeSpan = [0, 15];
% params.amp = 1;
% params.sigma = 2 * pi * 0.35;
% params.phase = pi/2;
simoptions.ODESim.MaxStep = 0.2;
simoptions.ODESim.InitialStep = 0.1;

% use a default sea
% simoptions.BuoySim.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);
simoptions.BuoySim.SeaParameters = seasetup ('Sigmas', 2 * pi * (1/5));
simoptions.BuoySim.NRadiationCoefs = 10;

% set the end stop position (the maximum allowed translator displacement)
simoptions.maxAllowedxT = inf;

% simoptions.BuoySim.tether_length = 4;

simoptions.ODESim.PreProcFcn = 'simfun_ACTM';
simoptions.ODESim.EvalFcn = 'mBodybuoysystemode_linear_nested'; 
simoptions.ODESim.PostPreProcFcn = 'mBodybuoysystemfinfun_ACTM';
simoptions.ODESim.PostAssemblyFcn = 'mBodybuoysyspostassembly_linear';
simoptions.ODESim.PostSimFcn = 'nestedsysresfun_AM'; 
% simoptions.events = 'systemevents_linear'; 

simoptions.ODESim.Solver = 'ode45';

%%
tic;
[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear(design, simoptions); 
toc;

%%
plotresultsbuoysys_linear(T, Y, results, design, outsimoptions, 1)

%%
mb = outdesign.MultiBodySystem;

axlims = [ min(results.xBs) - outsimoptions.BuoySim.BuoyParameters.a, ...
           max(results.xBs) + outsimoptions.BuoySim.BuoyParameters.a, ...
           -2*outsimoptions.BuoySim.BuoyParameters.a, ...
           2*outsimoptions.BuoySim.BuoyParameters.a, ...
           -10, ...
           max(results.xBh) + outsimoptions.ODESim.SolutionComponents.MultiBodySystem.Buoy.sz ];
           

mb.animate ('t', T, 'x', Y(:,21:end)', 'alpha', 0.1, 'PlotFcn', ...
                @(tind, t, x, mbod) plotvectors ( tind, t, x, mbod, ...
                                                  results.vRvec, ...
                                                  results.vRgenVec, ...
                                                  results.FptoVec, ...
                                                  results.FBDs, ...
                                                  results.FBDh, ...
                                                  axlims ));

%%

axlims = [ min(results.xBs) - outsimoptions.BuoySim.BuoyParameters.a, ...
           max(results.xBs) + outsimoptions.BuoySim.BuoyParameters.a, ...
           -2*outsimoptions.BuoySim.BuoyParameters.a, ...
           2*outsimoptions.BuoySim.BuoyParameters.a, ...
           -10, ...
           max(results.xBh) + outsimoptions.ODESim.SolutionComponents.MultiBodySystem.Buoy.sz ];
           

mb.animate ('t', T, ...
            'x', Y(:,21:end)', ...
            'alpha', 0.1, ...
            'filename', 'test_buoy_vid', ...
            'speed', 1, ...
            'PlotFcn', @(tind, t, x, mbod) plotvectors ( tind, t, x, mbod, ...
                                                         results.vRvec, ...
                                                         results.vRgenVec, ...
                                                         results.FptoVec, ...
                                                         results.FBDs, ...
                                                         results.FBDh, ...
                                                         axlims ) ...
           );
                                              
                                              
%%

figure;
maxT = 8; 
[ax, h1, h2] = plotyy (T(T<=maxT), [results.FBs(T<=maxT), results.FBh(T<=maxT)], T(T<=maxT), [results.vBs(T<=maxT), results.vBh(T<=maxT)]);
hFBs = h1(1); 
hFBh = h1(2); 
hvBs = h2(1); 
hvBh = h2(2);
legend ('FBs', 'FBh', 'vBs', 'vBh');

hvBh.Color = hFBh.Color;
hvBs.Color = hFBh.Color;
hvBs.LineStyle = '--';

hFBh.Color = hFBs.Color;
hFBs.LineStyle = '--';

%%

figure;
maxT = 8; 
[ax, h1, h2] = plotyy ( T(T<=maxT), ...
                        [ ...
                          results.buoyancy_force(T<=maxT), ...
                          results.excitation_force_heave(T<=maxT), ...
                          results.radiation_force_heave(T<=maxT), ...
                          results.FBDh(T<=maxT) ...
                          results.FBh(T<=maxT) ...
                          results.FptoVec(T<=maxT,3) ...
                          results.excitation_force_surge(T<=maxT), ...
                          results.radiation_force_surge(T<=maxT) ...
                          results.FBDs(T<=maxT) ...
                          results.FBs(T<=maxT) ...
                          results.FptoVec(T<=maxT,1) ...
                        ], ...
                        T(T<=maxT), ...
                        [ ...
                          results.xBs(T<=maxT), ...
                          results.xBh(T<=maxT) ...
                          results.vBs(T<=maxT), ...
                          results.vBh(T<=maxT) ...
                        ] ...
                      );
                  
h_buoyancy_force = h1(1); 
h_excitation_force_heave = h1(2); 
h_radiation_force_heave = h1(3); 
h_FBDh = h1(4);
h_FBh = h1(5);
h_Fptos = h1(6);
h_excitation_force_surge = h1(7);
h_radiation_force_surge = h1(8);
h_FBDs = h1(9);
h_FBs = h1(10);
h_Fptoh = h1(11);


h_buoyancy_force.Color = 'k';

h_excitation_force_surge.Color = h_excitation_force_heave.Color;
h_excitation_force_surge.LineStyle = '--';

h_radiation_force_surge.Color = h_radiation_force_heave.Color;
h_radiation_force_surge.LineStyle = '--';

h_FBDs.Color = h_FBDh.Color;
h_FBDs.LineStyle = '--';

h_FBs.Color = h_FBh.Color;
h_FBs.LineStyle = '--';

h_Fptos.Color = h_Fptoh.Color;
h_Fptos.LineStyle = '--';

legend ('buoyancy force', ...
        'excitation force heave', ...
        'radiation force heave', ...
        'FBDh', ...
        'FBh', ...
        'Fptoh', ...
        'excitation force surge', ...
        'radiation force surge', ...
        'FBDs', ...
        'FBh', ...
        'Fptos', ...
        'xBs', ...
        'xBh', ...
        'vBs', ...
        'vBh' );
    
xlabel ('Time, s')
ax1 = ax(1);
ax1.YLabel('Forces, N');
ax2 = ax(2);
ax2.YLabel('Velocity/Displacement');

%%

figure;
maxT = 8; 
[ax, h1, h2] = plotyy ( T(T<=maxT), ...
                        [ ...
                          results.buoyancy_force(T<=maxT), ...
                          ... results.excitation_force_heave(T<=maxT), ...
                          ... results.radiation_force_heave(T<=maxT), ...
                          ... results.excitation_force_surge(T<=maxT), ...
                          ... results.radiation_force_surge(T<=maxT) ...
                        ], ...
                        T(T<=maxT), ...
                        [ ...
                          ... results.vBs(T<=maxT), ...
                          results.xBh(T<=maxT) ...
                        ] ...
                      );

legend ('buoyancy force', 'xBh');

% h_buoyancy_force = h1(1); 
% h_excitation_force_heave = h1(2); 
% h_radiation_force_heave = h1(3); 
% h_excitation_force_surge = h1(4);
% h_radiation_force_surge = h1(5);
% 
% legend ('buoyancy force', 'excitation force heave', 'radiation force heave', 'excitation force surge', 'radiation force surge', 'vBs', 'vBh');


%% SImulate with buoy in random sea

% Set up the buoy and sea data files, these are for the 2m buoy
simoptions.buoy = 37;
% design.buoynum = simoptions.buoynum;

simoptions.ODESim.TimeSpan = [0, 60];
% params.amp = 1;

simoptions.BuoySim.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.ODESim.PreProcFcn = 'systemsimfun_ACTM';
simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTM';
simoptions.ODESim.PostSimFcn = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

% params.peak_freq = 1/9; % centred at resonant frequency
% params.sigma_range = [0.345575191894877,2.31745966692415;];
% params.water_depth = 50;

simoptions.BuoySim.SeaParameters = seasetup ('PMPeakFreq', 1/9);

[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear (design, simoptions); 
    
% produce some interesting plots of the output, but plotting only every 5th
% value in the output arrays for efficiency
skip = 5;
plotresultsbuoysys_linear(T, Y, results, outdesign, outsimoptions, skip)

%% plot
 
fnames = dir ([simoptions.ODESim.Nested.GeneratorSolver.solutionFileTemplate, '*.mat']);

S = load (fnames(1).name);
innerT = S.T';
innerY = S.Y;

for ind = 2:numel(fnames)
   
    S = load (fnames(ind).name);
    
    innerT = [innerT; S.T(2:end)'];
    
    innerY = [innerY; S.Y(2:end,:)];
    
end
              
figure; plotyy (T, Y(:,[simoptions.ODESim.SolutionComponents.MultiBodySystem.SolutionIndices(simoptions.ODESim.SolutionComponents.MultiBodySystem.BuoyPositionSolutionInd), ...
                        simoptions.ODESim.SolutionComponents.MultiBodySystem.SolutionIndices(simoptions.ODESim.SolutionComponents.MultiBodySystem.BuoyVelocitySolutionInd) ...
                       ] ...
                     ), ...
                innerT, innerY);
 