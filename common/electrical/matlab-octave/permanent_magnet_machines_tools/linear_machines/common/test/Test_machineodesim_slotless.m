% Test_machineodesim_slotless

load actm_design_and_simoptions.mat

%%

design.slm_psidot = slmengine(design.psilookup(1,:), design.psilookup(2,:),...
        'plot','off',...
        'verbosity',0,...
        'knots',6,...
        'leftslope', 0,...
        'leftvalue', design.psilookup(2,1),...
        'rightslope', 0,...
        'rightvalue', design.psilookup(2,end),...
        'InteriorKnots', 'free');
    
plotslm(design.slm_psidot)   

hold on
plot(design.psilookup(1,:), slmeval(design.psilookup(1,:), design.slm_psidot, 1))
hold off

design.FieldDirection = -1;
design.PowerPoles = design.Poles(1);

xBh = 0;
xBs = 0;
vBh = 1;
vBs = 0;

Icoils = vBh * -[1, 0, -1]';
    
[dpsidxF, EMF, Force, ForceVec, xT, vT] = ...
    machineodesim_slotless(design, simoptions, Icoils, xBh, xBs, vBh, vBs)
