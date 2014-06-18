% Test_simfun_RADIAL_SLOTTED
%
%

clear design simoptions 

design = test_design_RADIAL_SLOTTED ('external');

%%
design.RlVRp = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

design = completedesign_RADIAL_SLOTTED (design, simoptions);

[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);
fprintf(1, 'done simfun\n');

%%

simoptions.evaloptions = designandevaloptions_RADIAL_SLOTTED ();

[design, simoptions] = finfun_RADIAL_SLOTTED(design, simoptions);

