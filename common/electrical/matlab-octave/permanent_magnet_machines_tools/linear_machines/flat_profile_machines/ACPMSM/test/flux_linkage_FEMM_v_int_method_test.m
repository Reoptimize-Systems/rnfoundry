
design.J = 0;
design.mode = [1 0, 1];


pos = linspace(0.5, 1.5, 15);
circprops = zeros(numel(pos), 3);

for i = 1:numel(pos)
    
    % Draw the main sim to extract the vector potential
    RunStructFEMMSim_ACPMSM(design, pos(i)*design.Taup);

    mi_analyse(1);
    mi_loadsolution;
    
    circprops(i,:) = mo_getcircuitproperties('Coil A');
    
    mo_close;
    mi_close;
    
end

plot(design.psilookup(1,:), design.psilookup(2,:), pos, circprops(:,3))
legend('Agrid', 'FEMM')




%%

if isfield(design, 'Apoly')
    
    plotpolymodel(design.APoly, 'PlotType2D', 'scatter', 'GridSize', 50, 'PlotOpt', 'rx'); 
    
    hold on

    scatter3(xycoords(:,1), xycoords(:,2), p(:,1), 'b+')

    hold off

elseif isfield(design, 'Agrid')
    
    scatter3(design.X(:), design.Y(:), design.Agrid, 'rx')
    
end



%%

design.J = 0;
design.mode(1:2) = [1 0];
RunStructFEMMSim_ACPMSM(design, 0);

mi_analyse(1);
mi_loadsolution;

p = mo_getpointvalues(design.X(:), design.Y(:) * design.Taup);
    
    
scatter3(design.X(:), design.Y(:), design.Agrid(:), 'rx')

hold on

scatter3(design.X(:), design.Y(:), p(:,1), 'b+')

hold off


%%

tic
npoints = 45;
xrange = linspace(-design.dg+design.FEMMTol, design.dg-design.FEMMTol, npoints);
yrange = linspace(0, 1, npoints);
[design.X, design.Y] = meshgrid(xrange,yrange);

p = mo_getpointvalues(design.X(:), design.Y(:) * design.Taup);

toc

tic

npoints = 100;
xrange = linspace(-design.dg+design.FEMMTol, design.dg-design.FEMMTol, npoints);
yrange = linspace(0, 1, npoints);
[design.X, design.Y] = meshgrid(xrange,yrange);

p = mo_getpointvalues(design.X(:), design.Y(:) * design.Taup);

toc


%%
normpos = 0.5;
design.J = 0;
design.mode(1:2) = [1 0];
RunStructFEMMSim_ACPMSM(design, normpos*design.Taup);
int = integratehalfperiod2ddata(design.X, design.Y, design.Agrid, ...
    -0.558546, (normpos + 0.0861677) / design.Taup, 0, (normpos + 0.0967185) / design.Taup, 1e-10, 'max') / (design.Hc/2 * design.WcVTaup)
