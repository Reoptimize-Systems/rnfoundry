% Test_simulatemachine_TORUS_CORELESS


clear design simoptions 

design = test_design_RADIAL_SLOTTED ('external');

%%
design.RlVRp = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_RADIAL_SLOTTED (design, simoptions);

% setup simulation options
simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
                                'Rpm', 30, ...
                                'TSpan', [0,1]);

                            
simoptions.reltol = 1e-4;
simoptions.PhaseCurrentTols = repmat(0.001, 1, design.Phases);
simoptions.maxstep = (simoptions.tspan(2) - simoptions.tspan(1)) / 10000;

simoptions.evaloptions = designandevaloptions_RADIAL_SLOTTED ();

[design, simoptions] = feval(simoptions.simfun, design, simoptions);
simoptions.simfun = [];

[design, simoptions] = feval(simoptions.finfun, design, simoptions);
simoptions.finfun = [];

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions, ...
                                                         simoptions.simfun, ...
                                                         simoptions.finfun, ...
                                                         simoptions.odeevfun, ...
                                                         simoptions.resfun);
                                                         
%%

fsimoptions = simoptions;

bp = 0.1;
tmax = 1;

fsimoptions.tspan = [ 0:0.01/100:bp, ...
                      (bp+0.01/100):(bp+tmax)/1000:(bp+tmax)  ];
                      
fsimoptions.odesolver = 'odef1';

[fT, fY, fresults, fdesign, fsimoptions] = simulatemachine_AM(design, ...
                                                         fsimoptions, ...
                                                         fsimoptions.simfun, ...
                                                         fsimoptions.finfun, ...
                                                         fsimoptions.odeevfun, ...
                                                         fsimoptions.resfun);

                                                     