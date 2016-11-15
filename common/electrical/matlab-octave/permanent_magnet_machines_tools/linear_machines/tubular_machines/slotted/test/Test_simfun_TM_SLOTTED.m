% Test_simfun_TM_SLOTTED
%
%

clear design simoptions 

design = test_design_TM_SLOTTED ('dims');

%%
design.RlVRp = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_TM_SLOTTED (design, simoptions);

[design, simoptions] = simfun_TM_SLOTTED(design, simoptions);
fprintf(1, 'done simfun\n');

[design, simoptions] = finfun_TM_SLOTTED(design, simoptions);
fprintf(1, 'done finfun\n');