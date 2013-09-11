% Test_axialrotormesh

% torus_rotor_faesor_script

Rs = 1;
Rbi = Rs + 0.03;
Rmi = Rbi + 0.1;
Rmo = Rmi + 0.3;
Rbo = Rmo + 0.05;
tbi = [0.1, 0.2];
tsuppb = 0.2;
tausupp = tsuppb/2;
nstages = 2;
nmodules = 6;
nmodulesupports = 2;

tdiscsep = 1.5*tbi(2);
% layers = 2;
circumpoints = 25;
maglabels = 1:nstages+1;

% tic
[fens,gcells,maglabels,shaftlabels] = axialrotormesh(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tdiscsep, ...
                    tausupp, nmodules, nmodulesupports, ...
                    'ShaftAxialLayersPerM', 15, ...
                    'DiscAxialLayersPerM', 20, ...
                    'SupportAxialLayersPerM', 30, ...
                    'CircumPoints', circumpoints, ...
                    'MagnetLabels', maglabels, ...
                    'NStages', nstages);
% toc
% tic
% gv = drawmesh({fens, gcells}, 'gcells', 'facecolor', 'red', ...
%     'label', [maglabels,shaftlabels], ...
%     'labelcolor', [repmat({'blue'}, size(maglabels)), repmat({'green'}, size(maglabels))]);
% toc
% view(0,0)

%%

celllist = gcell_select(fens, gcells, struct ('flood', true, 'startfen', 1) );
gv = drawmesh({fens, gcells}, 'gcells', 'facecolor', 'red', 'cell_list', celllist)

