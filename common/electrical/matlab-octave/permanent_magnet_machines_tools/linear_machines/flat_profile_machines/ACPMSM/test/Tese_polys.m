

E = 200e9;
Taup = 100/1000; % pole pitch
hs = 44/1000; % armature thickness (the entire thickness of the coils)
hy = 25/1000;
g = 5/1000; % Air gap
lm = 23/1000; % Magnet thickness
dbi = hy; % back-iron thickness, we will use this as the minimum flange thickness?
ls = 1.0; % Stack length
bs = 33/1000; % Slot width
fieldSurfArea = 19; % Surface area of the field (stator in this case)
armSurfArea = fieldSurfArea + 12; % Surface area of the armature (translator in this case)
sections = 10;

bp = 0.9*Taup; % magnet pitch, will assume 90% of pole pitch (is probably less)
dg = lm + (hs/2) + g;
fieldPoles = ceil(fieldSurfArea / (Taup * ls)); % total number of poles in the field/stator of machine
fieldLength = fieldSurfArea / ls;
armPoles = ceil(armSurfArea / (Taup * ls));

dgVlm = dg / lm;
bpVlm = bp / lm;
taupVbp = Taup / bp;
lsVbp = ls / bp;
dbiVlm = dbi / lm;

dgThresh = 3/1000;

Force = gapclosingforce_ACPMSM(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm)

RunACPMSMFEMMForceSim(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm)

[Taup, bp, dg, ls, dbi] = ratios2dimensions_ACPMSM(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm)

[ppForce, ForceVArea] = GetFEMMClosingForce_ACPMSM(dbi, bp, Taup, ls, lm)