

dgVlm = 4;
bpVlm= 6; 
taupVbp = 1.5;
%lsVbp = numbers(i,4);
dbiVlm = 0.2;
lm = 0.02;

dg = dgVlm * lm;

bp = bpVlm * lm;

taup = taupVbp * bp;

ls = 1.0;

lsVbp = ls / bp;

dbi = dbiVlm * lm;

RunACPMSMFEMMForceSim(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm);

% Calculate force using integral of B^2
[closeForce, closeForceVA] = GetFEMMClosingForce_ACPMSM(dbi, bp, taup, ls, lm);

[magFixForce, magFixForceVA] = GetFEMMMagFixForce_ACPMSM(dbi, bp, taup, ls);

% Calculate force using weighted stress tensor integral 
MaxIntcloseForce = GetFEMMIntStressTensorForce_ACPMSM(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm);

% Get the energy and coenergy stored in the magnetic field in the air
% gap between the two magnets 
fEnergy = GetFEMMFieldEnergyIntegral_ACPMSM(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm);

% close the simulations
mo_close;
mi_close;

g = 2*(dg - lm);

% reduce the airgap by x %
displ = g*0.01;

dg = dg - (displ/2);

% Run the sim with reduced airgap
RunACPMSMFEMMForceSim(dg/lm, bpVlm, taupVbp, lsVbp, dbiVlm, lm);

% recalculate the energy integrals
fEnergy(:,3:4) = GetFEMMFieldEnergyIntegral_ACPMSM(dg/lm, bpVlm, taupVbp, lsVbp, dbiVlm, lm);

EnergyForce = (fEnergy(:,3)-fEnergy(:,1)) ./ displ

CoEnergyForce = (fEnergy(:,4)-fEnergy(:,2)) ./ displ

% close the simulations
mo_close;
mi_close;

