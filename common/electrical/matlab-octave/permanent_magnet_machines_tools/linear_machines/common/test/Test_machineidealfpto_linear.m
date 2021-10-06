% Test_machineidealfpto_linear

load generator_design_and_simoptions.mat

Fpto = 1000;
xBh = 0;
xBs = 0;
vBh = 1;
vBs = 0;
Icoils = [0,0,0];

[EMF, idealIcoils, limitIcoils, Force, Ploss] = ...
    machineidealfpto_linear(design, simoptions, Fpto, xBh, xBs, vBh, vBs, Icoils);