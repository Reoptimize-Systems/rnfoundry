% Test_simulatemachine_ACPMSM

clear design simoptions

design.bpVTaup = 0.85; 
design.lmVbp = 0.2; 
design.dgVlm = 2.0; 
design.lsVTaup = 3;
design.dbiVlm = 1;
design.WcVTaup = 1/3;
design.HcVgap = 0.95;
design.Taup = 0.2;
design.Ntot = 1000;
design.CoilFillFactor = 0.55;
design.J = 0;

design = ratios2dimensions_ACPMSM(design);

design.Phases = 3;
design.Poles = [10 30];
design.CoilTurns = 1000;
design.RlVRp = 10;
design.LgVLc = 0;


%% Test with linear motion

speed = 1;
simoptions.IC = zeros(1, design.Phases);
simoptions.skip = 1;
simoptions.tspan = [0, 5];
simoptions.drivetimes = 0:simoptions.tspan(2)/2:simoptions.tspan(2);
simoptions.vT = repmat(speed, size(simoptions.drivetimes));
simoptions.xT = simoptions.vT .* simoptions.drivetimes;
simoptions.tether_length = 0;
simoptions.NoOfMachines = 1;

simoptions.odeevfun = 'prescribedmotodeforcefcn_linear'; 
simoptions.forcefcn = 'forcefcn_linear_pscbmot'; 
simoptions.simfun = 'simfun_ACPMSM';
simoptions.resfun = 'prescribedmotresfun_linear';
simoptions.finfun = 'prescribedmotfinfun_ACPMSM';

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun); 

plotresultsproscribedmot_linear(T, Y, results, design, 1);   

