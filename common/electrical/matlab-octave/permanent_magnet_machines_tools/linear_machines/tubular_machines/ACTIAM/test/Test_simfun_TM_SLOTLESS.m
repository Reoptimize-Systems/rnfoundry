% Test_simfun_TM_SLOTTED
%

%%

design = test_design_TM_SLOTLESS ('dims');

design.RlVRp = 10;
simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_TM_SLOTLESS (design, simoptions);

[design, simoptions] = simfun_TM_SLOTLESS(design, simoptions);
fprintf(1, 'done simfun\n');

% [design, simoptions] = finfun_TM_SLOTLESS(design, simoptions);
% fprintf(1, 'done finfun\n');
[design, simoptions] = finfun_TM_SLOTLESS (design, simoptions);
fprintf(1, 'done finfun\n');

%%

design = test_design_TM_SLOTLESS ('dims');
design.RlVRp = 10;
simoptions = struct();
simoptions.GetVariableGapForce = false;
simoptions.MagFEASimType = 'multiple';

design = completedesign_TM_SLOTLESS (design, simoptions);

[design, simoptions] = simfun_TM_SLOTLESS (design, simoptions);
fprintf(1, 'done simfun\n');

[design, simoptions] = finfun_TM_SLOTLESS (design, simoptions);
fprintf(1, 'done finfun\n');


%%

load ('/home/rcrozier/Sync/work/enercro/Projects/WES-Wavedrive/Outputs/design_004_converted_to_rnfoundry.mat')

design.RlVRp = 8;
simoptions = struct();
simoptions.GetVariableGapForce = false;
simoptions.MagFEASimType = 'multiple';

design = completedesign_TM_SLOTLESS (design, simoptions);

design.zs = design.zp / 3;
% design.Phases = 6;
% design.WindingLayout.Phases = repmat ((1:design.Phases)', design.Poles(1), 1);

[design, simoptions] = simfun_TM_SLOTLESS (design, simoptions);
fprintf(1, 'done simfun\n');

[design, simoptions] = finfun_TM_SLOTLESS (design, simoptions);
fprintf(1, 'done finfun\n');

