% Test_evaluatefestructure_PMSM


load PMSM_Test_Design.mat

% IVarsCell = {[0.1, 0.1, 0.098, 0.098]};

options = designandevaloptions_PMSM;

options.nu = 0.31;

options.IMethod = '1.3';
options.gfactor = 0.5;
options.sections = 10;

design.AngleFromHorizontal = pi/2;
design.OuterWebs = 1;
design.sides = 2;
design.GuideRailIVars = [0.1, 0.1, 0.095, 0.095];
design.GuideRailIMethod = '1.3';

design.tols = [0;0;0];

%% 

design1 = evaluateRectSecStructure_PMSM(design, options);

design.tols = [0;0;0];
options.IMethod = '1.6';
design = evaluateIBeamStructure_PMSM(design, options)

%% 

design.InnerStructureBeamVars = [];
listsearchevalfcn = @evaluatefestructure_PMSM;
structureevalfcn = @feairgapclosure_PMSM;

% Now we will design the machine structure
design2 = evaluatestructure_FM(design, options, listsearchevalfcn, structureevalfcn);

