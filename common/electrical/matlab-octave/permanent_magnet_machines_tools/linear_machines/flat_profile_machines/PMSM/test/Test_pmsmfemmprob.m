
design.Phases = 3;         % Number of Phases in machine
design.Wp = 0.1;         % Pole Pitch
design.Wm = 0.8*design.Wp;      % Magnet Pitch
design.hm = 15/1000;       % Magnet depth
design.CoilFillFactor = 0.585;      % Copper fill factor
design.g = 3/1000;         % Initial air-gap
design.Dc = 1/1000;        % 1 mm diameter wire for both machines
design.hbf = 30/1000;      % Translator/Field yoke (back iron) thickness
design.hba = design.hbf;          % Stator/Armature yoke (half central section) thickness
design.ht = 2 * design.hba;
design.ls = 0.4;
design.Ws = design.Wp / design.Phases;
design.Wt = 0.1 * design.Wp;
design.HcMag = 979000;
design.CoilTurns = 100;


design = dimensions2ratios_PMSM(design);

design.mode = 1;

pos = 0;

[FemmProblem, design] = pmsmfemmprob(design, pos, 'NWIndingLayers', 2);

openprobleminfemm_mfemm(FemmProblem, 'temp.fem');

RunStructFEMMSim_PMSM(design, pos)
