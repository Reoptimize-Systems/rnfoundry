% Test_axialfluxouterrotor2dfemmprob

FemmProblem = newproblem_mfemm('planar');

% Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);


thetapole = 2*pi/20; 
thetamag = thetapole * 0.8;
rmag = 0.02;
rbackiron = rmag / 2;
drawnrotors = [1,1];
magsep = 0.2;
rrotor = [0.5, 0.5 - magsep];

FemmProblem = radialfluxrotor2dfemmprob(thetapole, thetamag, rmag, rbackiron, drawnrotors, rrotor, ...
                                        'FemmProblem', FemmProblem, ...
                                        'MagArrangement', 'NS', ...
                                        'MagnetMaterial', 3, ...;
                                        'BackIronMaterial', 2, ...
                                        'NPolePairs', 3);


plotfemmproblem (FemmProblem)
% openprobleminfemm_mfemm(FemmProblem);


