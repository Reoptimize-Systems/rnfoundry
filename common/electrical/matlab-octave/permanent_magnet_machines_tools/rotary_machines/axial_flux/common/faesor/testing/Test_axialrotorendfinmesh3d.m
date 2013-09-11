% Test_axialrotorendfinmesh3d

Rs = 1;
Rbi = Rs + 0.03;
Rmi = Rbi + 0.1;
Rmo = Rmi + 0.3;
Rbo = Rmo + 0.05;

tbi = 0.1;
separation = 3 * tbi;

tsuppb = tbi;
tausupp = 0.1;
nmodules = 6;
nmodulesupports = 1;
disclayers = 2;
supportlayers = 2;
circumpoints = 25;
maglabel = 1;

[fens,gcells,modulexs] = ...
    axialrotorendfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, ...
                           tausupp, nmodules, nmodulesupports, 'b', ...
                           disclayers, supportlayers, circumpoints, maglabel);

drawmesh({fens, gcells}, 'gcells', 'facecolor', 'red')


%%

E = 207e9;

nu = 0.31;

% create the finite element block
[feb,geom,u] = axialrotorfeb(fens, gcells, Rbi, separation, tbi, 2*tbi, E, nu);

% Assemble the system matrix
ems = stiffness(feb, geom, u)