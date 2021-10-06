% Test_drawsquaremags

bpVTaup = 0.8;
lmVbp = 0.3;
dgVlm = 2;
lsVTaup = 2;
dbiVlm = 0.5;
Taup = 1;
WcVTaup = 1/3;
hcV2dg = 0.95;
Ntot = 1;
kcufill = 0.55; 
J = 0;

[bp, lm, dg, dbi, ls, Wc, hc] = ratios2dimensions_ACPMSM(bpVTaup, lmVbp, dgVlm, lsVTaup, dbiVlm, Taup, WcVTaup, hcV2dg);

% pos = 0.5*Taup;

pos = -(bp + (Taup - bp)/2);

% pos = 3.5 * Taup;

openfemm;

Dc = RunFEMMSim_ACPMSM(bpVTaup, lmVbp, dgVlm, lsVTaup, dbiVlm, Taup, WcVTaup, hcV2dg, Ntot, kcufill, J, pos);

