% Test_machineodesim_linear_mvgarm

simoptions.BuoySim.tether_length = 1;
simoptions.NoOfMachines = 1;
design.Phases = 1;
Icoils=1;
xA=0;
vA=0;
xBh=design.PoleWidth/4;
xBs=0;
vBh=1;
vBs=0;

[dpsidxR, EMF, Force, ForceVec, xT, vT] = ...
    machineodesim_linear_mvgarm(design2, simoptions2, Icoils, xA, vA, xBh, xBs, vBh, vBs)


plotslm(design.slm_psidot, {'dy'})

tempx = [0:0.01:1];
figure;plot(tempx,slmpsidot_linear(design2, tempx, design.PoleWidth).* 641)