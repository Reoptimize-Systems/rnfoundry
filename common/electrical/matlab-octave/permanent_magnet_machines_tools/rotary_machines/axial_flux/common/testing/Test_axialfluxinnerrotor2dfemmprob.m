% Test_axialfluxinnerrotor2dfemmprob

FemmProblem = newproblem_mfemm('planar');

%   Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);

ypole = 3;
ymag = ypole * 0.8;
xmag = ymag / 3;
xbackiron = xmag / 2;
magsep = 2 * xmag;

FemmProblem = axialfluxinnerrotor2dfemmprob(ypole, ymag, xmag, xbackiron, magsep, ...
    'FemmProblem', FemmProblem, 'MagArrangement', 'NS', 'MagnetMaterial', 3, ...;
    'BackIronMaterial', 2, 'NInnerParts', 3);


filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))

%%

FemmProblem = newproblem_mfemm('planar');

%   Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);

ypole = 3;
ymag = ypole * 0.8;
xmag = ymag / 3;
xbackiron = 0;
magsep = 2 * xmag;

FemmProblem = axialfluxinnerrotor2dfemmprob(ypole, ymag, xmag, xbackiron, magsep, ...
    'FemmProblem', FemmProblem, 'MagArrangement', 'NS', 'MagnetMaterial', 3, ...;
    'BackIronMaterial', 2, 'NInnerParts', 3);


filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))


