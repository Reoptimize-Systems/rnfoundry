% Test_axialrotorfinmesh3d

Rs = 1;
Rbi = Rs + 0.03;
Rmi = Rbi + 0.1;
Rmo = Rmi + 0.3;
Rbo = Rmo + 0.05;

tbi = 0.1;

layers = 1;
circumpoints = 25;
nmodules = 2;

tdiscsep = 2 * tbi;
nstages = 1;
modulegap = 1e-2;
suppxs = [];
modulexs = sort(unique([suppxs, linspace(modulegap/2, 2*pi/nmodules - modulegap/2, ceil(circumpoints/nmodules)) ]));
outersep = ((nstages-1) * tbi) + (nstages * tdiscsep);

Inputs.bipointsperm = 100;
Inputs.magpointsperm = 100;
    
[fens,gcells] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, nmodules, layers, modulexs, 1, Inputs);

drawmesh({fens, gcells});

%%

E = 207e9;

nu = 0.31;

% create the finite element block
[feb,geom,u] = axialrotorfeb(fens, gcells, Rbi, outersep, tbi, 2*tbi, E, nu);

% Assemble the system matrix
ems = stiffness(feb, geom, u)
