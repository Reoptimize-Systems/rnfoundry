
sn1 = mbdyn.pre.structuralNode6dof ('dynamic');

mass = 1;

bd1 = mbdyn.pre.sphericalBody (mass, 0.1, sn1);

damping_matrix = diag (ones (6,1));

law = mbdyn.pre.linearViscousGenericConstituativeLaw (damping_matrix);

j_viscous_body = mbdyn.pre.viscousBody (sn1, law, 'null');

f_g = mbdyn.pre.gravity ();

prob = mbdyn.pre.initialValueProblem (0, 5, 0.01, 'LinearSolver', mbdyn.pre.linearSolver ('naive') ); 

% assemble the system
mbsys = mbdyn.pre.system ( prob, ...
                           'Nodes', {sn1}, ...
                           'Elements', {bd1, j_viscous_body, f_g}, ...
                           'DefaultOrientation', 'orientation matrix' );

thisfilepathwoext = mfilename('fullpath');

mkdir (thisfilepathwoext);

delete (fullfile (thisfilepathwoext, '*'));

inputfile = fullfile (thisfilepathwoext, 'input.mbd');

mbsys.generateMBDynInputFile (inputfile);

[mbstatus, cmdout, pid, inputfile] = mbdyn.mint.start_mbdyn (inputfile, ...
    'OutputPrefix', fullfile (thisfilepathwoext, 'input'), ...
    'MBDynOutputFile', fullfile (thisfilepathwoext, 'mbdyn_output.txt') );

mbpp = mbdyn.postproc ([inputfile(1:end-4)], mbsys);

mbpp.plotNodePositions ();

mbpp.plotNodeVelocities ();

mbpp.plotNetCDFVar (j_viscous_body, 'f');
