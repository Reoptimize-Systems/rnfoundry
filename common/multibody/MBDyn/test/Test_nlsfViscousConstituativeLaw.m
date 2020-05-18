
sn1 = mbdyn.pre.structuralNode6dof ('dynamic');

sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', [0;0;-0.1]);


cl = mbdyn.pre.clamp (sn1, 'node', 'node');

F_pto_max = 5;
c_pto = 10;

v_pto_max = F_pto_max / c_pto;

x = [-2*v_pto_max, -v_pto_max, v_pto_max, 2*v_pto_max];
y = [-F_pto_max,   -F_pto_max, F_pto_max, F_pto_max];

% plot (x, y);

nlsf = mbdyn.pre.multiLinearScalarFunction ( 'nlsf damping', x, y);

mass = 1;

bd1 = mbdyn.pre.sphericalBody (mass, 0.1, sn1);
bd2 = mbdyn.pre.sphericalBody (mass, 0.1, sn2);

% nlsf = mbdyn.pre.constScalarFunction ('nlsf damping', 0.5 * mass * 9.81);

nlsf_law = mbdyn.pre.nlsfViscousConstituativeLaw ('null', {'null', 'null', nlsf.name});

j_pto_damper = mbdyn.pre.deformableDisplacementJoint (sn1, sn2, nlsf_law, 'null', 'null', 'Offset2Reference', 'other node', 'Name', 'PTO_Damper');
% j_pto_damper=[];
f_g = mbdyn.pre.gravity ();

prob = mbdyn.pre.initialValueProblem (0, 5, 0.01, 'LinearSolver', mbdyn.pre.linearSolver ('naive') ); 

% assemble the system
mbsys = mbdyn.pre.system ( prob, ...
                           'Nodes', {sn1, sn2}, ...
                           'Elements', {bd1, bd2, cl, j_pto_damper, f_g}, ...
                           'DefaultOrientation', 'orientation matrix' );

mbsys.addScalarFunctions (nlsf);

thisfilepathwoext = mfilename('fullpath');

delete (fullfile (thisfilepathwoext, '*'));

inputfile = fullfile (thisfilepathwoext, 'input.mbd');

mbsys.generateMBDynInputFile (inputfile);

[mbstatus, cmdout, pid, inputfile] = mbdyn.mint.start_mbdyn (inputfile, ...
    'OutputPrefix', fullfile (thisfilepathwoext, 'input'), ...
    'MBDynOutputFile', fullfile (thisfilepathwoext, 'mbdyn_output.txt') );

mbpp = mbdyn.postproc ([inputfile(1:end-4)], mbsys);

mbpp.plotNodePositions ();

mbpp.plotNodeVelocities ();

mbpp.plotNetCDFVar (j_pto_damper, 'f');
