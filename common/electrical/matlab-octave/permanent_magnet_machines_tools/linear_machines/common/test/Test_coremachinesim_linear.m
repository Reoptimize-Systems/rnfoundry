% Test_coremachinesim_linear

x = 0:0.01:1;

design.slm_psidot = slmengine(x, cos(pi * x),...
        'plot','off',...
        'verbosity',0,...
        'knots',6,...
        'leftslope', 0,...
        'leftvalue', cos(pi * x(1)),...
        'rightslope', 0,...
        'rightvalue', cos(pi * x(end)),...
        'positiveinflection', 0.5, ...
        'InteriorKnots', 'fixed');
    
design.FieldDirection = 1;  
design.PoleWidth = 1;
design.PowerPoles  = 1;

design.Phases = 3;
design.tauco = 1/design.Phases;
design.taupcg = design.Phases * design.tauco;
% calculate the physical separation between adjacent coils in each phase
design.CoilPositions = ((0:design.Phases-1) .* (1/design.Phases));
% adjust positions to give correct polarity
design.CoilPositions(2:2:end) = (design.CoilPositions(2:2:end) + 1);
design.CoilPositions = design.CoilPositions * design.taupcg;
simoptions.NoOfMachines = 1;

design.CoilPositions = [0, 2/3, 4/3]

%%

xRE = 0;
vEF = 1;
vRE = 0;

xEF = linspace(0, 6, 100);

Icoils = zeros(numel(xEF), design.Phases);
EMF = Icoils;
dpsidxCRTF = Icoils;

for i = 1:numel(xEF)

    [junk, junk, tEMF, tdpsidxCRTF] = coremachinesim_linear(design, simoptions, xEF(i), xRE, vEF, vRE, 0);
    
    Icoils(i,:) = tEMF ./ 1;
    
    [FEF(i), FRE(i), EMF(i,:), dpsidxCRTF(i,:)] = coremachinesim_linear(design, simoptions, xEF(i), xRE, vEF, vRE, Icoils(i,:));
    
end

figure;
subplot(2, 1, 1); plot(xEF, EMF);

subplot(2, 1, 2); plot(xEF, FEF)
