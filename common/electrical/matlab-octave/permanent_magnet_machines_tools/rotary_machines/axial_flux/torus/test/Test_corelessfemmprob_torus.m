% Test_corelessfemmprob_torus

design.taupm = 1;
design.taumm = 0.8;

design.g = 0.1; 
design.tc = 0.3;
design.tauco = 0.333 * design.taupm;
design.tauci = 0.7 * design.tauco;

design.tm = 0.15;
design.tbi = 0.05;

design.Rmo = 10;
design.Rmi = 9;

design.CoilTurns = 100;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagFEASimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';


Inputs.MagnetRegionMeshSize = -1;
Inputs.BackIronRegionMeshSize = -1;
Inputs.OuterRegionsMeshSize = [-1, -1];
Inputs.AirGapMeshSize = -1;
Inputs.DrawCoils = true;
Inputs.NStages = 2;
FemmProblem = corelessfemmprob_torus(design, Inputs);

    
filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))