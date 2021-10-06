% Test_innerstructuremesh_RF
%
% 

Rsii = 1;
Rsio = 1.1;
Rsoi = 3;
Rsoo = 3.1;
thetasp = 2*pi / 30;
ls = 1;
nspokes = 5;
% layers = 10;
% ncircumpoints = 50;
% ncircumsppoints = 3;
% nradsipoints = 3;
% nradsppoints = 3;
% nradsopoints = 3;
% labels.outersurface = 1;
faeprob = faesorprob;
[faeprob, labels] = innerstructuremesh_RF(faeprob, Rsii, Rsio, Rsoi, Rsoo, thetasp, ls, nspokes);

% faeprob.drawmesh('gcells', 'label', 1, 'labelcolor', 'red')

E = 200e9;
nu = 0.3;

[faeprob,mater] = innerstructurefeb_RF(faeprob, Rsii, ls, E, nu);

shearforce = 1000;
radialforce = 10000;
structdensity = 7500;
omega = 0;

[Kmat, Fmat] = innerstructurestresses_RF(faeprob, Rsoo, ls, shearforce, ...
                                         radialforce, structdensity, omega, labels);
                
x = Kmat \ Fmat;
% toc
% tic
% % use an iterative method to save memory
% x = bicgstab(Kmat, Fmat, [], 5000);
% toc
faeprob.u = scatter_sysvec(faeprob.u, x);

% draw the displacements as a colorfield
comp = 6;
% choose scale that makes max deflection appear as half gap between fins
% uscale = (tdiscsep/2) / maxdisp(3);
gv = graphic_viewer;
gv = reset (gv,[]);
nvals = getvals(faeprob.u);
nvals = nvals(:,3);
dcm = data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
colorfield = field(struct ('name', 'colorfield', 'data', map_data(dcm, nvals)));
% draw(feb, gv, struct ('x', geom, 'u', uscale*u, 'colorfield', colorfield, 'shrink',1));
draw(faeprob.feb, gv, struct ('x', faeprob.geom, 'u', 10000*faeprob.u, 'colorfield', colorfield, 'shrink',1));
