%% example_basic_cosimulation.m
%
% 
% Make a very simple problem of a single body (a sphere) in a uniform
% gravitational field. Apply a force after some time, then turn it off
% after some more time

%% Problem Definition

% define the sphere
mass = 1;
radius = 0.05;
inertiamat = eye (3) * (2/5) * mass * radius^2;

% sim start time
itime = 0;
% sim end time
ftime = 10;
% sim time step
nsteps = 100;
tstep = (ftime - itime) / nsteps;

% set aceleration due to gravity
g = 9.81;

% size of force to apply to the sphere
force_magnitude = 2 * mass * g;
% angle of applied force in the X-Z plane
force_angle = deg2rad (45);
% time when to start applying the force
force_start_time = itime + 0.1 * (ftime - itime);
% time when to stop applying the force
force_stop_time = itime + 0.75 * (ftime - itime);

% now make the MBDyn system

% create a dynamic node at the origin
node = mbdyn.pre.structuralNode6dof ('dynamic');

% create the body attached to the node
cog = [0;0;0];
body = mbdyn.pre.body (mass, cog, inertiamat, node);

% add gravity
grav = mbdyn.pre.gravity ('GravityAcceleration', g);

% set up the communication methods between MBDyn and Matlab
if ispc
    % windows does not have 'local' sockets like Linux, so use inet socket
    communicator = mbdyn.pre.socketCommunicator ('Port', 5500, ...
                                                 'Coupling', 'tight', ...
                                                 'Create', 'yes', ...
                                                 'SleepTime', 0, ...
                                                 'SendAfterPredict', 'no');
else
    % add the socket forces
    communicator = mbdyn.pre.socketCommunicator ('Path', '/tmp/mbdyn.sock', ...
                                                 'Coupling', 'tight', ...
                                                 'Create', 'yes', ...
                                                 'SleepTime', 0, ...
                                                 'SendAfterPredict', 'no');
end

% create the external force element
matlab_force = mbdyn.pre.externalStructuralForce ({node}, [], communicator);

% set up the initial-value problem
pbm = mbdyn.pre.initialValueProblem ( itime, ftime, tstep, ...
                                      'MaxIterations', 10, ...
                                      'ResidualTolerance', 1e-9 );

% put it all together in an MBDyn system
mbsys = mbdyn.pre.system ( pbm, ...
                           'Nodes', {node}, ...
                           'Elements', {body, grav, matlab_force}, ...
                           'DefaultOrientation', 'orientation matrix' );

%% Start MBDyn

% put MBDyn input file in temporary directory, we don't need it once MBDyn
% has read it
mbdpath = fullfile (tempdir(), 'example_cosim.mbd' );

% put output in new folder in current directory
outputdir = fullfile ('.', 'Test_mbdyn_MBCNodal_output');

mkdir (outputdir);

% make output file be named Test_mbdyn_MBCNodal.mov etc. in the output
% directory
outputfile_prefix = fullfile (outputdir, 'Test_mbdyn_MBCNodal');

% create the MBCNodal object to handle communication with MBDyn
mb = mbdyn.mint.MBCNodal ('MBDynPreProc', mbsys, ...
                          'UseMoments', true, ...
                          'MBDynInputFile', mbdpath, ...
                          'OverwriteInputFile', true, ...
                          'OutputPrefix', outputfile_prefix ...
                          );
                      
mb.start ('Verbosity', 1);

nnodes = mb.GetNodes ();

%% Perform Simulation Loop

% to match example program 
maxiter = 3; 

% initialise some variables
status = 0;
time = 0;
ind = 1;
pos = zeros (3,1);
rot = zeros (3,3,1);

while status == 0
    
    % Apply the force depending on the specified start and end times
    if time(ind) < force_start_time
        fmag = 0;
    elseif time(ind) >= force_stop_time
        fmag = 0;
    else
        fmag = force_magnitude;
    end
    
    % Make 3D force vector to apply to node
    forces = [ fmag * cos(force_angle); 
               0; 
               fmag * sin(force_angle); ];
    
    for iter = (maxiter-1):-1:1

        status = mb.GetMotion ();
        
        if status ~= 0
            break;
        end
    
        fprintf (1, 'iter: %d, iterstatus: %d\n', iter, status);
        
        % set the forces
        mb.F (forces);
        mb.M (zeros (3,1));

        % apply the force but tell MBDyn that the simulation at the Matlab end
        % has not converged, so it cannot proceed to the next time step, but
        % should recalculate the motion based on the new forces and resend
        % the result back to the Matlab simulation.
        mb.applyForcesAndMoments (false)

    end
    
    if status ~= 0
        break;
    end
    
    status = mb.GetMotion ();
    
    if status ~= 0
        break;
    end

    pos(:,ind) = mb.X(1);
    
    rot(:,:,ind) = mb.GetRot();

    % set the forces
    mb.F (forces);
    mb.M (zeros (3,1));

    % apply the force and tell MBDyn that the simulation at the Matlab end
    % has converged, so it can proceed to the next time step
    mb.applyForcesAndMoments (true)
    
    % get next time step
    ind = ind + 1;
    time(ind) = time(ind-1) + tstep;

end

clear mb;

% strip the extra time step which will have been added
time(end) = [];

%% Perform Post-Processing

% we can plot the local Matlab data as normal
plot ( time, pos');
legend ('x position', 'y position', 'z position');
xlabel ('time [s]');
ylabel ('displacement [m]');

% load the results files
mbout = mbdyn.postproc ( outputfile_prefix, mbsys );

% plot the node trajectory, we replace the default title with our own
mbout.plotNodeTrajectories ('Title', 'Sphere Trajectory');

% plot all the node velocities
mbout.plotNodeVelocities ()







