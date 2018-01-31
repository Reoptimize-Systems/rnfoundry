% Test_simulatemachine_TORUS_CORELESS


clear design simoptions 

design = test_design_RADIAL_SLOTTED ('external');

%%
design.RlVRp = 10;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'InterpolateMissing', false);
 
 
simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_RADIAL_SLOTTED (design, simoptions);

% setup simulation options
%simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
%                                'Rpm', 30, ...
%                                'TSpan', [0,1]);

simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
                                'torquefcn', 'torquefcn_RADIAL_SLOTTED', ...
                                'PoleCount', 1000, ...
                                'RampPoles', 20 );
                            
simoptions.reltol = 1e-6;
%simoptions.PhaseCurrentTols = repmat(0.001, 1, design.Phases);
%simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 10000;

simoptions.Evaluation = designandevaloptions_RADIAL_SLOTTED ();

[design, simoptions] = feval(simoptions.ODESim.PreProcFcn, design, simoptions);
simoptions.ODESim.PreProcFcn = [];

[design, simoptions] = feval(simoptions.ODESim.PostPreProcFcn, design, simoptions);
simoptions.ODESim.PostPreProcFcn = [];

[T, Y, results, design, simoptions] = simulatemachine_AM(design, ...
                                                         simoptions);
                                                         
%%
%
%fsimoptions = simoptions;
%
%bp = 0.1;
%tmax = 1;
%
%fsimoptions.ODESim.TimeSpan = [ 0:0.01/100:bp, ...
%                      (bp+0.01/100):(bp+tmax)/1000:(bp+tmax)  ];
%                      
%fsimoptions.ODESim.Solver = 'odef1';
%
%[fT, fY, fresults, fdesign, fsimoptions] = simulatemachine_AM(design, ...
%                                                         fsimoptions);

                                                     