% Test_axialouterrotorfinsmesh3d

Rs = 1;
Rbi = Rs + 0.03;
Rmi = Rbi + 0.1;
Rmo = Rmi + 0.3;
Rbo = Rmo + 0.05;

tbi = 0.1;
tsuppb = 0.2;
tausupp = tsuppb/2;
nmodules = 6;
nmodulesupports = 2;
separation = 3 * tbi;
layers = 2;
circumpoints = 25;
maglabels = [1, 10];

[fens, gcells] = ...
    axialouterrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tausupp, ...
        nmodules, nmodulesupports, separation, 'layers', layers, 'circ', circumpoints, 'magl', maglabels)

gv = drawmesh({fens, gcells}, 'gcells', 'facecolor', 'red', 'label', maglabels, 'labelcolor', {'blue', 'green'});

% gv = headlight(gv);
% camlight right



