% Test_axialinnerrotorfinsmesh3d

Rs = 1;
Rbi = Rs + 0.03;
Rmi = Rbi + 0.1;
Rmo = Rmi + 0.3;
Rbo = Rmo + 0.05;

tbi = 0.1;
nfins = 3;
outersep = 10 * nfins * tbi;
maglabels = [];
layers = 2;
circumpoints = 25;
nmodules = 6;
nmodulesupports = 2;



[fens, gcells] = axialinnerrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, outersep, nfins, maglabels, layers, circumpoints);

% drawmesh({fens, gcells}, 'gcells', 'label', [1,2], 'labelcolor', {'red', 'blue'});

drawmesh({fens, gcells}, 'gcells');

