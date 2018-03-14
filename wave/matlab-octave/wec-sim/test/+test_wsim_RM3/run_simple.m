% +test_wsim_RM3/run.m.m

clear waves simu hsys mbsys float_hbody spar_hbody ref_float ref_spar mb

%% Hydro Simulation Data
simu = wsim.simSettings (getmfilepath ('test_wsim_RM3.run'));  % Create the Simulation Variable
% simu.mode = 'normal';                 %Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
% simu.explorer='on';                   %Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                   %Simulation Start Time [s]
simu.endTime=400;                       %Simulation End bdcloseTime [s]
simu.solver = 'ode4';                   %simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.1; 							%Simulation time-step [s]
simu.rampT = 100;                       %Wave Ramp Time Length [s]
simu.multibodySolver = 'MBDyn';
simu.b2b = true;

%% Wave Information 

% Regular Waves
waves = wsim.waveSettings ('regularCIC');        %Create the Wave Variable and Specify Type                               
waves.H = 2.5;                          %Wave Height [m]
waves.T = 8;                            %Wave Period [s]
%%%%%%%%%%%%%%%%%%

%% Hydrodynamic body system

% Float
float_hbody = wsim.hydroBody('hydroData/rm3.h5', 'CaseDirectory', simu.caseDir);      
    %Create the wsim.hydroBody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
float_hbody.mass = 'equilibrium';                   
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
float_hbody.momOfInertia = [20907301, 21306090.66, 37085481.11];  %Moment of Inertia [kg*m^2]     
float_hbody.geometryFile = fullfile ('geometry', 'float.stl');    %Location of Geomtry File

% Spar/Plate
spar_hbody = wsim.hydroBody('hydroData/rm3.h5', 'CaseDirectory', simu.caseDir); 
spar_hbody.mass = 'equilibrium';                   
spar_hbody.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody.geometryFile = fullfile ('geometry', 'plate.stl'); 

% make a hydrosys object for simulation
hsys = wsim.hydroSystem (waves, simu, [float_hbody, spar_hbody]);

% set up transient simulation
hsys.initialiseHydrobodies ();
hsys.odeSimSetup ();
[hydro_mbnodes, hydro_mbbodies, hydro_mbelements] = hsys.makeMBDynComponents ();

%% Multibody dynamics system specification (mbdyn)

problem_options.ResidualTol = 1e-5;
problem_options.MaxIterations = 200;
problem_options.Output = {}; 
problem_options.NonLinearSolver = [];
problem_options.LinearSolver = [];

[mbsys, initptodpos] = test_wsim_RM3.make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options);
                     
% draw it
% mbsys.draw ('Mode', 'wireghost', 'Light', false);

mbsys.draw ( 'Mode', 'solid', ...
             'Light', true, ...
             'AxLims', [-30, 30; -30, 30; -35, 35], ...
             'Joints', false, ...
             'StructuralNodes', false)

mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

%% Set up PTO

k = 0;
c = 1200000;

forcefcn = @(xRpto, vRpto) -k*xRpto -c*vRpto;

pto = mbdyn.mint.twoNodeTranslationalForce ( hydro_mbnodes{2}, hydro_mbnodes{1}, 3, 'forcefcn', forcefcn);

%% Run the simulation
%

% set some tolerance used to determine if forces have converged
absforcetol = 100;
relforcetol = 1e-3;

% all of the following will be wrapped up into wecsim class in future to
% hide details of implementation (instead you will just make a wsim.wecSim
% object and call wsim.wecSim.run () to run the simulation).
%

% start mbdyn
outputfile_prefix = fullfile (simu.caseDir, 'output', 'RM3');

delete ([outputfile_prefix, '.*']);

% create the communicator object
mb = mbdyn.mint.MBCNodal ('MBDynPreProc', mbsys, ...
                          'UseMoments', true, ...
                          'MBDynInputFile', mbdpath, ...
                          'OverwriteInputFile', true, ...
                          'OutputPrefix', outputfile_prefix ...
                          );
                      
mb.start ('Verbosity', 0);

nnodes = mb.GetNodes ();

time = mbsys.problems{1}.initialTime;

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

F_ExcitLin = out.F_ExcitLin;
F_ViscousDamping = out.F_ViscousDamping;
F_AddedMass = out.F_AddedMass;
F_Restoring = out.F_Restoring;
F_RadiationDamping = out.F_RadiationDamping;
F_ExcitNonLin = out.F_ExcitNonLin;
F_MorrisonElement = out.F_MorrisonElement;
F_Excit = out.F_Excit;
F_ExcitRamp = out.F_ExcitRamp;
FptoVec = [0;0;0];

% accept the last data into the time history of solutions
hsys.advanceStep (time(end), vel, accel);
    
ind = 2;

miniters = 0;
maxiters = mbsys.problems{1}.maxIterations;

ptoforce = 0;
xRpto = 0;
vRpto = 0;
xRptoVec = [0;0;0];
vRptoVec = [0;0;0];

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
    
    [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] =  pto.force ();
    
    forces (1:3,1,ind) = forces (1:3,1,ind) + Fpto(:,1);
    forces (1:3,2,ind) = forces (1:3,2,ind) + Fpto(:,2);
    
    FptoVec(1:3,1,ind) = Fpto(:,1); % store for comparison later
    
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

    [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));

    forces (:,:,ind) = newhydroforces;

    [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] =  pto.force ();
    
    forces (1:3,1,ind) = forces (1:3,1,ind) + Fpto(:,1);
    forces (1:3,2,ind) = forces (1:3,2,ind) + Fpto(:,2);
    
    FptoVec(1:3,1,ind) = Fpto(:,1); % store for comparison later
    
	mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (false);
    
    forcediff = abs (hydroforces - newhydroforces);

    maxforces = max(hydroforces, newhydroforces);
    relforcediff = abs(forcediff) ./ abs(maxforces);
    relforcediff(maxforces == 0) = 0;
%     disp(relforcediff)
    itercount = 1;
    while mbconv ~= 0 ...
        || itercount < miniters ...
        || (max (forcediff(:)) > absforcetol) ...
        || (ind > 3 && (max (relforcediff(:)) > relforcetol))
            
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

        [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));

        forces (:,:,ind) = newhydroforces;


        [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] = pto.force ();

        forces (1:3,1,ind) = forces (1:3,1,ind) + Fpto(:,1);
        forces (1:3,2,ind) = forces (1:3,2,ind) + Fpto(:,2);

        FptoVec(1:3,1,ind) = Fpto(:,1); % store for comparison later

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
    hsys.advanceStep (time(end), vel(:,:,ind), accel(:,:,ind));
    
    ind = ind + 1;
    
end

[F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys, forces, F_AddedMass, accel);
toc;
clear mb;

fprintf (1, 'Reached time %f, in %d steps\n', time(end), ind-1);

% return

%%
figure;
tmin = 0;
tmax = 400;
plotinds =  time>=tmin & time<=tmax;
plotyy (time(plotinds), [ squeeze(forces(1:3,1,plotinds))',  ptoforce(plotinds)'], time(plotinds), vRpto(plotinds));
legend ('fx', 'fy', 'fz', 'ptoforce', 'vRpto');

