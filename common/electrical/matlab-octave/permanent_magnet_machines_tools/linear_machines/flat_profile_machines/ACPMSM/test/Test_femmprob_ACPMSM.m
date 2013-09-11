% Test_femmprob_ACPMSM

design.bpVTaup = 0.8;
design.lmVbp = 0.3;
design.dgVlm = 2;
design.lsVTaup = 2;
design.dbiVlm = 0.5;
design.Taup = 1;
design.WcVTaup = 1/3;
design.hcVgap = 0.95;
design.CoilTurns = 1;
design.kcufill = 0.55; 
design.HcMag = 979000;
design.Dc = 1.25e-3;


design.CoilLayers = 2;

design = ratios2dimensions_ACPMSM(design);

% pos = 0.5*design.Taup;

pos = -(design.bp + (design.Taup - design.bp)/2);

% pos = 3.5 * design.Taup;

FemmProblem = femmprob_ACPMSM(design, pos, 'NWindingLayers', design.CoilLayers);

% plotfemmproblem(FemmProblem);

openprobleminfemm_mfemm(FemmProblem);