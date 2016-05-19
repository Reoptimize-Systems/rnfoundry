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



%% Internal armature
%
%

clear design simoptions 

design = test_design_RADIAL_SLOTTED ('internal');

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

%% Test parfor for fea

simoptions = struct();
simoptions.MagFEASim.UseParFor = true;
simoptions.GetVariableGapForce = false;

design.RlVRp = 10;
design = completedesign_RADIAL_SLOTTED (design, simoptions);

[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);
fprintf(1, 'done simfun\n');

