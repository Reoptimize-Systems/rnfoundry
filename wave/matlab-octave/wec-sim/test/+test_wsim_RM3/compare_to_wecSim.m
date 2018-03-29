% +test_wsim_RM3/run.m.m

clear body ...
      constraint ...
      output ...
      pto ...
      simu ...
      thisdir ...
      waves ....
      wecSim_rm3_dir ...
      waves_wsim ...
      simu_wsim ...
      hsys ...
      mbsys ...
      float_hbody ...
      spar_hbody ...
      ref_float ...
      ref_spar ...
      mb ...
      problem_options ...
      initptodpos ...
      mbdpath ...
      outputfile_prefix ...
      F_ExcitLin ...
      F_ViscousDamping ...
      F_AddedMass ...
      F_Restoring ...
      F_RadiationDamping ...
      F_ExcitNonLin ...
      F_MorrisonElement ...
      F_Excit ...
      F_ExcitRamp ...
      FptoVec ...
      pos vel accel eul R ...
      newhydroforces out forces ...
      gf_* test_*
      

if exist ('wecSim', 'file')
    wecSim_rm3_dir = fullfile (fileparts (which ('wecSim')), '..', 'tutorials', 'RM3');
else
    error ('wecSim.m must be on the path');
end

% Run wecSim original
thisdir = pwd ();

cd (wecSim_rm3_dir);

wecSim ();

cd (thisdir);


%% Hydro Simulation Data

simu_wsim = wsim.simSettings (getmfilepath ('test_wsim_RM3.run'));  % Create the Simulation Variable
% simu.mode = 'normal';                 % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
% simu.explorer='on';                   % Turn SimMechanics Explorer (on/off)
simu_wsim.startTime = simu.startTime;   % Simulation Start Time [s]
simu_wsim.endTime = simu.endTime;       % Simulation End bdcloseTime [s]
simu_wsim.solver = 'ode4';              % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu_wsim.dt = simu.dt; 				% Simulation time-step [s]
simu_wsim.rampT = simu.rampT;           % Wave Ramp Time Length [s]
simu_wsim.multibodySolver = 'MBDyn';
simu_wsim.b2b = simu.b2b;
simu_wsim.ssCalc = simu.ssCalc;

%% Wave Information 

waves_wsim = wsim.waveSettings (waves.type);       % Create the Wave Variable and Specify Type                               
waves_wsim.H = waves.H;                            % Wave Height [m]
waves_wsim.T = waves.T;                            % Wave Period [s]
waves_wsim.spectrumType = waves.spectrumType;

%% Hydrodynamic body system

% Float
float_hbody = wsim.hydroBody('hydroData/rm3.h5', 'CaseDirectory', simu_wsim.caseDir);      
    %Create the wsim.hydroBody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
float_hbody.mass = 'equilibrium';                  
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
float_hbody.momOfInertia = body(1).momOfInertia;  %Moment of Inertia [kg*m^2]     
float_hbody.geometryFile = fullfile ('geometry', 'float.stl');    %Location of Geomtry File

% Spar/Plate
spar_hbody = wsim.hydroBody('hydroData/rm3.h5', 'CaseDirectory', simu_wsim.caseDir); 
spar_hbody.mass = 'equilibrium';                   
spar_hbody.momOfInertia = body(2).momOfInertia;
spar_hbody.geometryFile = fullfile ('geometry', 'plate.stl'); 

% make a hydrosys object for simulation
hsys = wsim.hydroSystem (waves_wsim, simu_wsim, [float_hbody, spar_hbody]);

% set up transient simulation
hsys.initialiseHydrobodies ();
hsys.odeSimSetup ();
[hydro_mbnodes, hydro_mbbodies, hydro_mbelements] = hsys.makeMBDynComponents ();

%% Multibody dynamics system specification (mbdyn)

problem_options.ResidualTol = 1e-6;
problem_options.MaxIterations = 200;
problem_options.Output = {}; % 'iterations', 'solution', 'jacobian matrix', 'matrix condition number', 'solver condition number'
% problem_options.NonLinearSolver = mbdyn.pre.newtonRaphsonSolver ();
problem_options.NonLinearSolver = [];
% problem_options.LinearSolver = mbdyn.pre.linearSolver ('naive');
problem_options.LinearSolver = [];
problem_options.DefaultElementOutput = {'none'};

[mbsys, initptodpos] = test_wsim_RM3.make_multibody_system (waves_wsim, simu_wsim, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options);

mbdpath = fullfile (simu_wsim.caseDir, 'RM3.mbd');

%% Set up PTO

k = pto(1).k;
c = pto(1).c;

forcefcn = @(time, xRpto, vRpto) -k*xRpto -c*vRpto;

pto = mbdyn.mint.twoNodeTranslationalForce ( hydro_mbnodes{2}, hydro_mbnodes{1}, 3, 'forcefcn', forcefcn);

%% Run the simulation

% start mbdyn
outputfile_prefix = fullfile (simu_wsim.caseDir, 'output', 'RM3');

delete ([outputfile_prefix, '.*']);

% create the communicator object
mb = mbdyn.mint.MBCNodal ('MBDynPreProc', mbsys, ...
                          'UseMoments', true, ...
                          'MBDynInputFile', mbdpath, ...
                          'OverwriteInputFile', true, ...
                          'OutputPrefix', outputfile_prefix ...
                          );
                      
mb.start ('Verbosity', 0);

%%
nnodes = mb.GetNodes ();

time = mbsys.problems{1}.initialTime;
ind = 1;

status = mb.GetMotion ();

if status ~= 0
    error ('mbdyn returned %d, aborting sim, check output file:\n%s\nfor clues at to why.', status, mb.MBDynOutputFile)
end

eul = zeros (3,nnodes);

R = mb.GetRot();
for Rind = 1:size (R,3)
    om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
    eul(1:3,Rind) = om.euler123();
end

pos = [ mb.NodePositions(); 
        eul ];
    
vel = [ mb.NodeVelocities(); 
        mb.NodeOmegas() ];
    
accel = [ mb.NodeAccelerations(); 
          mb.NodeAngularAccels() ];

[forces, out] = hsys.hydroForces (time, pos, vel, accel);

% set the forces
mb.F (forces(1:3,:));
mb.M (forces(4:6,:));

mbconv = mb.applyForcesAndMoments (false);

% preallocate data vectors
nsteps = (simu_wsim.endTime - simu_wsim.startTime) ./ simu_wsim.dt + 1;

time = repmat (time, [1, nsteps]);
pos = repmat (pos, [1, 1, nsteps]);
vel = repmat (vel, [1, 1, nsteps]);
accel = repmat (accel, [1, 1, nsteps]);
ptoforce = ones (1, nsteps) * nan;
xRpto = ptoforce;
vRpto = ptoforce;
xRptoVec = repmat ([0;0;0], [1, nsteps]);
vRptoVec = repmat ([0;0;0], [1, nsteps]);
FptoVec = repmat ([0;0;0], [1, 1, nsteps]);
% twoNodeTranslationalForce test variables
test_ptoforce = ptoforce;
test_xRpto = xRpto;
test_vRpto = vRpto;
test_FptoVec = FptoVec;
test_xRptoVec = xRptoVec;
test_vRptoVec = vRptoVec;

F_ExcitLin = repmat (out.F_ExcitLin, [1,1,nsteps]);
F_ViscousDamping = repmat (out.F_ViscousDamping, [1,1,nsteps]);
F_AddedMass = repmat (out.F_AddedMass, [1,1,nsteps]);
F_Restoring = repmat (out.F_Restoring, [1,1,nsteps]);
F_RadiationDamping = repmat (out.F_RadiationDamping, [1,1,nsteps]);
F_ExcitNonLin = repmat (out.F_ExcitNonLin, [1,1,nsteps]);
F_MorrisonElement = repmat (out.F_MorrisonElement, [1,1,nsteps]);
F_Excit = repmat (out.F_Excit, [1,1,nsteps]);
F_ExcitRamp = repmat (out.F_ExcitRamp, [1,1,nsteps]);

% accept the last data into the time history of solutions
hsys.advanceStep (time(1), vel(:,:,ind), accel(:,:,ind));
    
ind = 2;

plotvectors = false;
checkoutputs = false;
miniters = 0;
maxiters = mbsys.problems{1}.maxIterations;
absforcetol = 100;
relforcetol = 1e-5;

if plotvectors
    figure;
    hvectplotax = axes;
end

test_ptoforce(1) = 0;
xRpto(1) = 0;
vRpto(1) = 0;

tic
while status == 0
    
    status = mb.GetMotion ();
    
    if status ~= 0
        continue;
    end
    
    time(ind) = time(ind-1) + mbsys.problems{1}.timeStep;                              
    
    R = mb.GetRot();
    
    for Rind = 1:size (R,3)
        om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
        eul(1:3,Rind) = om.euler123();
    end

    pos(:,:,ind) = [ mb.NodePositions(); eul];
    vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
    accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

    [hydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));

    forces (:,:,ind) = hydroforces;
    
    % PTO force
    [test_Fpto, ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();
    
    FptoVec(1:3,1,ind) = test_Fpto(:,1);
    forces (1:3,1,ind) = forces (1:3,1,ind) + test_Fpto(:,1);
    forces (1:3,2,ind) = forces (1:3,2,ind) + test_Fpto(:,2);
    
	mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (false);
    
    status = mb.GetMotion ();
    
    if status ~= 0
                
        F_ExcitLin(:,:,ind) = out.F_ExcitLin;
        F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
        F_AddedMass(:,:,ind) = out.F_AddedMass;
        F_Restoring(:,:,ind) = out.F_Restoring;
        F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
        F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
        F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
        F_Excit(:,:,ind) = out.F_Excit;
        F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
        
        ind = ind + 1;
        
        continue;
    end
    
    R = mb.GetRot();
    
    for Rind = 1:size (R,3)
        om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
        eul(1:3,Rind) = om.euler123();
    end

    pos(:,:,ind) = [ mb.NodePositions(); eul];
    vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
    accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

    % Hydrodynamic forces
    [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));

    forces (:,:,ind) = newhydroforces;

    % PTO force
    [test_Fpto, test_ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();
    
    FptoVec(1:3,1,ind) = test_Fpto(:,1);
    forces (1:3,1,ind) = forces (1:3,1,ind) + test_Fpto(:,1);
    forces (1:3,2,ind) = forces (1:3,2,ind) + test_Fpto(:,2);
    
	mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (false);
    
    
    forcediff = abs (hydroforces - newhydroforces);

    maxforces = max(hydroforces, newhydroforces);
    relforcediff = abs(forcediff) ./ abs(maxforces);
    relforcediff(maxforces == 0) = 0;
    itercount = 1;
    while mbconv ~= 0 ...
        || itercount < miniters ...
        || (max (forcediff(:)) > absforcetol) ...
        || (ind > 3 && (max (relforcediff(:)) > relforcetol))
        
        % store the previously calculated hydrodynamic forces
        hydroforces = newhydroforces;
        
        status = mb.GetMotion ();
        
        if status ~= 0
            break;
        end
        
        R = mb.GetRot();
    
        for Rind = 1:size (R,3)
            om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
            eul(1:3,Rind) = om.euler123();
        end

        pos(:,:,ind) = [ mb.NodePositions(); eul];
        vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
        accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

        % Hydrodynamic forces
        [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));

        forces (:,:,ind) = newhydroforces;
        
        % PTO force
        [test_Fpto, ptoforce(ind), test_xRpto(ind) vRpto(ind)] = pto.force ();

        FptoVec(1:3,1,ind) = test_Fpto(:,1);
        forces (1:3,1,ind) = forces (1:3,1,ind) + test_Fpto(:,1);
        forces (1:3,2,ind) = forces (1:3,2,ind) + test_Fpto(:,2);

        mb.F (forces(1:3,:,ind));
        mb.M (forces(4:6,:,ind));

        mbconv = mb.applyForcesAndMoments (false);
        
        itercount = itercount + 1;
        
        if itercount > maxiters
            error ('mbdyn iterations exceeded max allowed');
        end
        
        forcediff = abs (hydroforces - newhydroforces);
        maxforces = max(hydroforces, newhydroforces);
        relforcediff = abs(forcediff) ./ abs(maxforces);
        relforcediff(maxforces == 0) = 0;
    
    end
    
    status = mb.GetMotion ();
        
    if status ~= 0
        
        F_ExcitLin(:,:,ind) = out.F_ExcitLin;
        F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
        F_AddedMass(:,:,ind) = out.F_AddedMass;
        F_Restoring(:,:,ind) = out.F_Restoring;
        F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
        F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
        F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
        F_Excit(:,:,ind) = out.F_Excit;
        F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
        
        ind = ind + 1;
    
        break;
    end

    mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (true);
    
    F_ExcitLin(:,:,ind) = out.F_ExcitLin;
    F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
    F_AddedMass(:,:,ind) = out.F_AddedMass;
    F_Restoring(:,:,ind) = out.F_Restoring;
    F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
    F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
    F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
    F_Excit(:,:,ind) = out.F_Excit;
    F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
    
    % accept the last data into the time history of solutions
    hsys.advanceStep (time(ind), vel(:,:,ind), accel(:,:,ind));
    
    ind = ind + 1;
    
end

[F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys, forces, F_AddedMass, accel);
toc;
clear mb;

fprintf (1, 'Reached time %f, in %d steps\n', time(end), ind-1);

%%
figure;
tmin = 0;
tmax = 400;
plotinds =  time>=tmin & time<=tmax;
plotyy (time(plotinds), [ squeeze(forces(1:3,1,plotinds))',  test_ptoforce(plotinds)'], time(plotinds), vRpto(plotinds));
legend ('fx', 'fy', 'fz', 'ptoforce', 'vRpto');

%%

if ~exist ('output', 'var')
    warning ('Not comparing output to original WEC-Sim as ''output'' variable is not in the workspace.')
else

    %% plot force comparison
    figure;
    tmin = 0;
    tmax = 400;

    bodyind = 2;

    tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end) ]);
    plotinds =  time>=tmin & time<=tmax;

    plotyy (time(plotinds), ...
          [ squeeze(F_ExcitRamp(:,bodyind,plotinds))', ...
            ... squeeze(F_ViscousDamping(3,bodyind,plotinds)), ...
            squeeze(F_AddedMassCorrected(:,bodyind,plotinds))', ...
            squeeze(F_Restoring(:,bodyind,plotinds))', ...
            squeeze(F_RadiationDamping(:,bodyind,plotinds))', ...
            ...squeeze(F_ExcitNonLin(3,bodyind,plotinds)), ...
            ... squeeze(F_MorrisonElement(3,bodyind,plotinds)), ...
            ...squeeze(F_Excit(3,bodyind,plotinds)), ...
            ...squeeze(F_ExcitRamp(3,bodyind,plotinds)), ...
            ...ptoforce(plotinds)', ...
            output.bodies(bodyind).forceExcitation(plotinds,:), ...
            output.bodies(bodyind).forceAddedMass(plotinds,:), ...
            output.bodies(bodyind).forceRestoring(plotinds,:), ...
            output.bodies(bodyind).forceRadiationDamping(plotinds,:), ...
            ], ...
            time(plotinds), vRpto(plotinds) );

    legend ( 'wsim F ExcitRamp x', ...
             'wsim F ExcitRamp y', ...
             'wsim F ExcitRamp z', ...
             'wsim M ExcitRamp x', ...
             'wsim M ExcitRamp y', ...
             'wsim M ExcitRamp z', ...
             ... 'F ViscousDamping',  ...
             'wsim F addedmass x',  ...
             'wsim F addedmass y',  ...
             'wsim F addedmass z',  ...
             'wsim M addedmass x',  ...
             'wsim M addedmass y',  ...
             'wsim M addedmass z',  ...
             ...
             'wsim F Restoring x',  ...
             'wsim F Restoring y',  ...
             'wsim F Restoring z',  ...
             'wsim M Restoring x',  ...
             'wsim M Restoring y',  ...
             'wsim M Restoring z',  ...
             ...
             'wsim F RadiationDamping x', ...
             'wsim F RadiationDamping y', ...
             'wsim F RadiationDamping z', ...
             'wsim M RadiationDamping x', ...
             'wsim M RadiationDamping y', ...
             'wsim M RadiationDamping z', ...
             ...
             'wecSim F ExcitRamp x', ...
             'wecSim F ExcitRamp y', ...
             'wecSim F ExcitRamp z', ...
             'wecSim M ExcitRamp x', ...
             'wecSim M ExcitRamp y', ...
             'wecSim M ExcitRamp z', ...
             ... 'F ViscousDamping',  ...
             'wecSim F addedmass x',  ...
             'wecSim F addedmass y',  ...
             'wecSim F addedmass z',  ...
             'wecSim M addedmass x',  ...
             'wecSim M addedmass y',  ...
             'wecSim M addedmass z',  ...
             ...
             'wecSim F Restoring x',  ...
             'wecSim F Restoring y',  ...
             'wecSim F Restoring z',  ...
             'wecSim M Restoring x',  ...
             'wecSim M Restoring y',  ...
             'wecSim M Restoring z',  ...
             ...
             'wecSim F RadiationDamping x', ...
             'wecSim F RadiationDamping y', ...
             'wecSim F RadiationDamping z', ...
             'wecSim M RadiationDamping x', ...
             'wecSim M RadiationDamping y', ...
             'wecSim M RadiationDamping z', ...
             ... 'ptoforce',  ...
             'vRpto' );


    %% PTO force comparison
    figure;
    tmin = 0;
    tmax = 400;

    tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end) ]);
    plotinds =  time>=tmin & time<=tmax;
    
    plotyy (time(plotinds), ...
          [ squeeze(FptoVec(:,1,plotinds))', ...
            output.ptos(1).forceTotalWorld(plotinds,:) - output.ptos(1).forceConstraintWorld(plotinds,:), ...
            ], ...
            time(plotinds), vRpto(plotinds) );
        
    legend ( 'wsim Fpto x', ...
             'wsim Fpto y', ...
             'wsim Fpto z', ...
             'wecSim Fpto x', ...
             'wecSim Fpto y', ...
             'wecSim Fpto z', ...
             'vRpto' ...
            );
        
    %% plot force comparison
    figure;
    tmin = 0;
    tmax = 400;

    bodyind = 1;

    tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end) ]);
    plotinds =  time>=tmin & time<=tmax;

    plotyy (time(plotinds), ...
          [ squeeze(pos(1,bodyind,plotinds)), ...
            squeeze(pos(2,bodyind,plotinds)), ...
            squeeze(pos(3,bodyind,plotinds)), ...
            squeeze(pos(4,bodyind,plotinds)), ...
            squeeze(pos(5,bodyind,plotinds)), ...
            squeeze(pos(6,bodyind,plotinds)), ...
            squeeze(vel(1,bodyind,plotinds)), ...
            squeeze(vel(2,bodyind,plotinds)), ...
            squeeze(vel(3,bodyind,plotinds)), ...
            squeeze(vel(4,bodyind,plotinds)), ...
            squeeze(vel(5,bodyind,plotinds)), ...
            squeeze(vel(6,bodyind,plotinds)), ...
            output.bodies(bodyind).position(plotinds,1),   ...
            output.bodies(bodyind).position(plotinds,2),   ...
            output.bodies(bodyind).position(plotinds,3),   ...
            output.bodies(bodyind).position(plotinds,4),   ...
            output.bodies(bodyind).position(plotinds,5),   ...
            output.bodies(bodyind).position(plotinds,6),   ...
            output.bodies(bodyind).velocity(plotinds,1),   ...
            output.bodies(bodyind).velocity(plotinds,2),   ...
            output.bodies(bodyind).velocity(plotinds,3),   ...
            output.bodies(bodyind).velocity(plotinds,4),   ...
            output.bodies(bodyind).velocity(plotinds,5),   ...
            output.bodies(bodyind).velocity(plotinds,6),   ...
          ], ...
            time(plotinds), vRpto(plotinds) );

    legend ( 'wsim x', ...
             'wsim y', ...
             'wsim z', ...
             'wsim a', ...
             'wsim b', ...
             'wsim c', ...
             'wsim vx', ...
             'wsim vy', ...
             'wsim vz', ...
             'wsim omega a', ...
             'wsim omega b', ...
             'wsim omega c', ...
             'wecSim x', ...
             'wecSim y', ...
             'wecSim z', ...
             'wecSim a', ...
             'wecSim b', ...
             'wecSim c', ...
             'wecSim vx', ...
             'wecSim vy', ...
             'wecSim vz', ...
             'wecSim omega a', ...
             'wecSim omega b', ...
             'wecSim omega c', ...
             'vRpto' ...
            );
    
    
    %%
% 
%     mbout = mbdyn.postproc ( outputfile_prefix, mbsys ); 
% 
%     mbout.plotNodeTrajectories ('AxLims', [-1.5, 1.5; -1.5, 1.5; -25, 5]);


    %%
    % mbout.animate ( 'DrawMode', 'solid', ...
    %                 'Light', true, ...
    %                 'skip', 5, ...
    %                 'AxLims', [-30, 30; -30, 30; -35, 35])


    %% Compare to origninal WEC-Sim

    %
    % Requires to have run the RM3 example first using original WEC-Sim

    doplot = false;

    if doplot
        for bodyind = 1:numel (output.bodies)

            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceTotal, ...
                  time, squeeze(forces(:,bodyind,:)));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceTotal,  ...
                  time, squeeze(F_Total(:,bodyind,:)));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            title (sprintf ('forceTotal vs F\\_Total for body %d', bodyind));
    %         figure;
    %         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceExcitation,  output.bodies(bodyind).time, squeeze(F_ExcitRamp(:,bodyind,:))); 
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceAddedMass,  ...
                  output.bodies(bodyind).time, squeeze(F_AddedMass(:,bodyind,:)));
            title (sprintf ('forceAddedMass vs F\\_addedmass for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceAddedMass,  ...
                  output.bodies(bodyind).time, squeeze(F_AddedMassCorrected(:,bodyind,:)));
            title (sprintf ('forceAddedMass vs F\\_AddedMassCorrected for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceRadiationDamping,  ...
                  time, squeeze(F_RadiationDamping(:,bodyind,:)));
            title (sprintf ('forceRadiationDamping vs F\\_RadiationDamping for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
    %         figure;
    %         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceRestoring,  output.bodies(bodyind).time, squeeze(F_Restoring(:,bodyind,:)));
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).position,  ...
                  time, squeeze(pos(:,bodyind,:)));
            title (sprintf ('output.bodies(%d).position vs pos for body %d', bodyind, bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
        end

    end

    tmin = 0;
    tmax = 400;

    bodyind = 1;

    tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end)]);
    plotinds =  time>=tmin & time<=tmax;

%     figure;
%     plotyy ( output.bodies(bodyind).time(plotinds), ...
%              [output.ptos(1).forceTotal(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%              time(plotinds), ...
%              [vRptoVec(plotinds)', [ output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3)] ] );
%     title (sprintf ('output.ptos(%d).forceTotal vs FptoVec for body %d', bodyind, bodyind));
%     legend ('1', '2', '3', '4', '5', '6', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
%     figure;
%     plotyy ( output.bodies(bodyind).time(plotinds), ...
%             [output.ptos(1).forceTotalWorld(plotinds,:) - output.ptos(1).forceConstraintWorld(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%             time(plotinds), ...
%             [vRptoVec(:,plotinds)', output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3) ] );
%     title (sprintf ('output.ptos(1).forceTotalWorld - output.ptos(1).forceConstraintWorld vs FptoVec for body %d', bodyind, bodyind));
%     legend ('1', '2', '3', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
%     figure;
%     plot ( output.bodies(bodyind).time(plotinds), ...
%            output.ptos(1).forceInternalMechanics(plotinds,3),  ...
%            time(plotinds), ...
%            squeeze(ptoforce(plotinds)));
%     title (sprintf ('output.ptos(%d).forceInternalMechanics(:,3) vs ptoforce for body %d', bodyind, bodyind));
%     legend ('forceInternalMechanics(:,3)', 'ptoforce');

    %% Statistical comparison

    data = [];
    rowheadings = {};

    stats = {'3', '5', '8'};
    gf_forceTotal = nan * ones (numel (output.bodies), numel(stats));
    gf_momentTotal = nan * ones (numel (output.bodies), numel(stats));
    gf_F_Total = gf_forceTotal;
    gf_M_Total = gf_F_Total;
    gf_forceExcitation = gf_forceTotal;
    gf_momentExcitation = gf_forceTotal;
    gf_forceAddedMass = gf_forceTotal;
    gf_F_AddedMassCorrected = gf_forceTotal;
    gf_M_AddedMassCorrected = gf_forceTotal;
    gf_forceRadiationDamping = gf_forceTotal;
    gf_momentRadiationDamping= gf_forceTotal;
    gf_forceRestoring = gf_forceTotal;
    gf_momentRestoring = gf_forceTotal;

    gf_pos = gf_forceTotal;
    gf_theta = gf_forceTotal;
    gf_vel = gf_forceTotal;
    gf_omega = gf_forceTotal;
    gf_accel = gf_forceTotal;
    gf_omegaaccel = gf_forceTotal;

    for bodyind = 1:numel (output.bodies)

        % calculate some proper stats
        gf_xforceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,1), squeeze(forces(1,bodyind,:))', stats);
        gf_yforceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,2), squeeze(forces(2,bodyind,:))', stats);
        gf_zforceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,3), squeeze(forces(3,bodyind,:))', stats);
        
        gf_amomentTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,4), squeeze(forces(4,bodyind,:))', stats);
        gf_bmomentTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,5), squeeze(forces(5,bodyind,:))', stats);
        gf_cmomentTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,6), squeeze(forces(6,bodyind,:))', stats);

    %     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, squeeze(F_Total(:,bodyind,:))', stats);
        gf_F_X_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,1), squeeze(F_Total(1,bodyind,:))', stats);
        gf_F_Y_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,2), squeeze(F_Total(2,bodyind,:))', stats);
        gf_F_Z_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,3), squeeze(F_Total(3,bodyind,:))', stats);
        
        gf_M_A_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,4), squeeze(F_Total(4,bodyind,:))', stats);
        gf_M_B_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,5), squeeze(F_Total(5,bodyind,:))', stats);
        gf_M_C_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,6), squeeze(F_Total(6,bodyind,:))', stats);

        gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
        gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);
        
        gf_F_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,1:3), squeeze(F_AddedMassCorrected(1:3,bodyind,:))', stats);
        gf_M_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,4:6), squeeze(F_AddedMassCorrected(4:6,bodyind,:))', stats);

        gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,1:3), squeeze(F_RadiationDamping(1:3,bodyind,:))', stats);
        gf_momentRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,4:6), squeeze(F_RadiationDamping(4:6,bodyind,:))', stats);

        gf_forceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,1:3),  squeeze(F_Restoring(1:3,bodyind,:))', stats); 
        gf_a_momentRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,4),  squeeze(F_Restoring(4,bodyind,:))', stats); 
        gf_b_momentRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,5),  squeeze(F_Restoring(5,bodyind,:))', stats); 
        gf_c_momentRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,6),  squeeze(F_Restoring(6,bodyind,:))', stats); 

        gf_x_pos(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,1),  squeeze(pos(1,bodyind,:))', stats);
        gf_y_pos(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,2),  squeeze(pos(2,bodyind,:))', stats);
        gf_z_pos(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,3),  squeeze(pos(3,bodyind,:))', stats);
        
        gf_a_theta(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,4),  squeeze(pos(4,bodyind,:))', stats); 
        gf_b_theta(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,5),  squeeze(pos(5,bodyind,:))', stats); 
        gf_c_theta(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,6),  squeeze(pos(6,bodyind,:))', stats); 

        gf_vel(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,1:3),  squeeze(vel(1:3,bodyind,:))', stats); 
        gf_omega(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,4:6),  squeeze(vel(4:6,bodyind,:))', stats); 

        gf_accel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,1:3),  squeeze(accel(1:3,bodyind,:))', stats); 
        gf_omegaaccel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,4:6),  squeeze(accel(4:6,bodyind,:))', stats); 

        rowheadings = [rowheadings, {
                   ...['forceTotal_body_', int2str(bodyind)], ...
                   sprintf('Body %d Total X Force Uncorrected', bodyind), ...
                   sprintf('Body %d Total Y Force Uncorrected', bodyind), ...
                   sprintf('Body %d Total Z Force Uncorrected', bodyind), ...
                   sprintf('Body %d Total A Moment Uncorrected', bodyind), ...
                   sprintf('Body %d Total B Moment Uncorrected', bodyind), ...
                   sprintf('Body %d Total C Moment Uncorrected', bodyind), ...
                   sprintf('Body %d Total X Force', bodyind), ...
                   sprintf('Body %d Total Y Force', bodyind), ...
                   sprintf('Body %d Total Z Force', bodyind), ...
                   sprintf('Body %d Total A Moment', bodyind), ...
                   sprintf('Body %d Total B Moment', bodyind), ...
                   sprintf('Body %d Total C Moment', bodyind), ...
                   sprintf('Body %d Excitation Force', bodyind), ...
                   sprintf('Body %d Excitation Moment', bodyind), ...
                   ...['forceAddedMass_body_', int2str(bodyind)], ...
                   sprintf('Body %d Added Mass Force', bodyind), ...
                   sprintf('Body %d Added Mass Moment', bodyind), ...
                   sprintf('Body %d Radiation and Damping Force', bodyind), ...
                   sprintf('Body %d Radiation and Damping Moment', bodyind), ...
                   sprintf('Body %d Hydrostatic Restoring Force', bodyind), ...
                   sprintf('Body %d A Hydrostatic Restoring Moment', bodyind), ...
                   sprintf('Body %d B Hydrostatic Restoring Moment', bodyind), ...
                   sprintf('Body %d C Hydrostatic Restoring Moment', bodyind), ...
                   sprintf('Body %d X Position', bodyind), ...
                   sprintf('Body %d Y Position', bodyind), ...
                   sprintf('Body %d Z Position', bodyind), ...
                   sprintf('Body %d A Angular Position', bodyind), ...
                   sprintf('Body %d B Angular Position', bodyind), ...
                   sprintf('Body %d C Angular Position', bodyind), ...
                   sprintf('Body %d Velocity', bodyind), ...
                   sprintf('Body %d Angular Velocity', bodyind), ...
                   sprintf('Body %d Acceleration', bodyind), ...
                   sprintf('Body %d Angular Acceleration', bodyind) }];


        data = [ data;
                 gf_xforceTotal(bodyind,:);
                 gf_yforceTotal(bodyind,:);
                 gf_zforceTotal(bodyind,:);
                 gf_amomentTotal(bodyind,:);
                 gf_bmomentTotal(bodyind,:);
                 gf_cmomentTotal(bodyind,:);
                 gf_F_X_Total(bodyind,:);
                 gf_F_Y_Total(bodyind,:);
                 gf_F_Z_Total(bodyind,:);
                 gf_M_A_Total(bodyind,:);
                 gf_M_B_Total(bodyind,:);
                 gf_M_C_Total(bodyind,:);
                 gf_forceExcitation(bodyind,:);
                 gf_momentExcitation(bodyind,:);
                 ...gf_forceAddedMass(bodyind,:);
                 gf_F_AddedMassCorrected(bodyind,:);
                 gf_M_AddedMassCorrected(bodyind,:);
                 gf_forceRadiationDamping(bodyind,:);
                 gf_momentRadiationDamping(bodyind,:); 
                 gf_forceRestoring(bodyind,:);
                 gf_a_momentRestoring(bodyind,:);
                 gf_b_momentRestoring(bodyind,:);
                 gf_c_momentRestoring(bodyind,:);
                 gf_x_pos(bodyind,:);
                 gf_y_pos(bodyind,:);
                 gf_z_pos(bodyind,:);
                 gf_a_theta(bodyind,:);
                 gf_b_theta(bodyind,:);
                 gf_c_theta(bodyind,:);
                 gf_vel(bodyind,:);
                 gf_omega(bodyind,:);
                 gf_accel(bodyind,:);
                 gf_omegaaccel(bodyind,:); ];
    end

    % display table
    colheadings = { 'RMSE', 'MAE', 'R2' };

    wid = 16;
    fms = {};
    fileID = 1;

    displaytable (data,colheadings,wid,fms,rowheadings,fileID);

end
