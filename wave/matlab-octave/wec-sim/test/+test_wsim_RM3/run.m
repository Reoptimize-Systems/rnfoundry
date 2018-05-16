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
simu.b2b = true;

%% Wave Information 

% noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       %Create the Wave Variable and Specify Type      
%%%%%%%%%%%%%%%%%%%

% Regular Waves
waves = wsim.waveSettings ('regularCIC');        %Create the Wave Variable and Specify Type                               
waves.H = 2.5;                          %Wave Height [m]
waves.T = 8;                            %Wave Period [s]
%%%%%%%%%%%%%%%%%%

% Regular Waves
% waves = wsim.waveSettings ('regular');        %Create the Wave Variable and Specify Type                               
% waves.H = 2.5;                          %Wave Height [m]
% waves.T = 8;                            %Wave Period [s]
% simu.ssCalc = 1;
%%%%%%%%%%%%%%%%%%%

% Irregular Waves using PM Spectrum with Convolution Integral Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'PM';
%%%%%%%%%%%%%%%%%%%

% % Irregular Waves using BS Spectrum with Convolution Integral Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'BS';
% %%%%%%%%%%%%%%%%%%%

% Irregular Waves using BS Spectrum with State Space Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'BS';
% simu.ssCalc = 1;	
%Control option to use state space model 
%%%%%%%%%%%%%%%%%%%

% Irregular Waves using User-Defined Spectrum
% waves = wsim.waveSettings ('irregularImport');  %Create the Wave Variable and Specify Type
% waves.spectrumDataFile = 'ndbcBuoyData.txt';  %Name of User-Defined Spectrum File [2,:] = [omega, Sf]
%%%%%%%%%%%%%%%%%%%

% User-Defined Time-Series
% waves = wsim.waveSettings ('userDefined');     %Create the Wave Variable and Specify Type
% waves.etaDataFile = 'umpqua46229_6_2008.mat'; % Name of User-Defined Time-Series File [:,2] = [time, wave_elev]
%%%%%%%%%%%%%%%%%%%

%% Hydrodynamic body system

% Float
float_hbody = wsim.hydroBody('rm3.h5');      
    %Create the wsim.hydroBody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
float_hbody.mass = 'equilibrium';                   
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
float_hbody.momOfInertia = [20907301, 21306090.66, 37085481.11];  %Moment of Inertia [kg*m^2]     
float_hbody.geometryFile = 'float.stl';    %Location of Geomtry File

% Spar/Plate
spar_hbody = wsim.hydroBody('rm3.h5'); 
spar_hbody.mass = 'equilibrium';                   
spar_hbody.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody.geometryFile = 'plate.stl'; 

% make a hydrosys object for simulation
hsys = wsim.hydroSystem (waves, simu, [float_hbody, spar_hbody]);

% set up transient simulation
hsys.initialiseHydrobodies ();
hsys.timeDomainSimSetup ();
[hydro_mbnodes, hydro_mbbodies, hydro_mbelements] = hsys.makeMBDynComponents ();

%% Multibody dynamics system specification (mbdyn)

problem_options.ResidualTol = 1e-5;
problem_options.MaxIterations = 200;
problem_options.Output = {}; % 'iterations', 'solution', 'jacobian matrix', 'matrix condition number', 'solver condition number'
% problem_options.NonLinearSolver = mbdyn.pre.newtonRaphsonSolver ();
problem_options.NonLinearSolver = [];
% problem_options.LinearSolver = mbdyn.pre.linearSolver ('naive');
problem_options.LinearSolver = [];

[mbsys, initptodpos] = test_wsim_RM3.make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options);
                     
% draw it
% mbsys.draw ('Mode', 'wireghost', 'Light', false);

mbsys.draw ( 'Mode', 'ghost', ...
             'Light', true, ...
             'AxLims', [-30, 30; -30, 30; -35, 35], ...
             'Joints', false, ...
             'StructuralNodes', true)

mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

%% Set up PTO

k = 0;
c = 1200000;

forcefcn = @(time, xRpto, vRpto) -k*xRpto -c*vRpto;

pto = mbdyn.mint.twoNodeTranslationalForce ( hydro_mbnodes{2}, hydro_mbnodes{1}, 3, 'forcefcn', forcefcn);

%% Run the simulation

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

%%
nnodes = mb.GetNodes ();

time = mbsys.problems{1}.initialTime;

% % set the forces
% mb.F (zeros(3,nnodes));
% mb.M (zeros(3,nnodes));
% 
% mbconv = mb.applyForcesAndMoments (true);

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
                               
% [forces, out] = hsys.hydroForces ( time, ...
%                                    [output.bodies(1).position(ind,:)', output.bodies(2).position(ind,:)'], ...
%                                    [output.bodies(1).velocity(ind,:)', output.bodies(2).velocity(ind,:)'], ...
%                                    [output.bodies(1).acceleration(ind,:)', output.bodies(2).acceleration(ind,:)'], ...
%                                    );

% set the forces
mb.F (forces(1:3,:));
mb.M (forces(4:6,:));

mbconv = mb.applyForcesAndMoments (false);

% fprintf (1, '    time: %f, ind %d, mbconv: %d\n', time, 1, mbconv);
% 
% status = mb.GetMotion ();
% 
% pos = [ mb.NodePositions(); 
%         mb.GetRot() ];
%     
% vel = [ mb.NodeVelocities(); 
%         mb.NodeOmegas() ];
%     
% accel = [ mb.NodeAccelerations(); 
%           mb.NodeAngularAccels() ];
%       
% % % set the forces
% mb.F (forces(1:3,:));
% mb.M (forces(4:6,:));
% 
% mbconv = mb.applyForcesAndMoments (true);
% 
% fprintf (1, '    time: %f, ind %d, mbconv: %d\n', time, 1, mbconv);

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

plotvectors = false;
checkoutputs = false;
miniters = 0;
maxiters = mbsys.problems{1}.maxIterations;
absforcetol = 100;
relforcetol = 1e-3;

if plotvectors
    figure;
    hvectplotax = axes;
end

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
    
    if checkoutputs
        
        % get the last line of output file to see if time diverges
        fid = fopen ( [ outputfile_prefix, '.out'], 'rt');

        lastline = '';
        L = fgetl (fid);
        while L ~= -1
            lastline = L;
            L = fgetl (fid);
        end
        fclose (fid)
        C = textscan (lastline, '%s %f %f %f %f %f %f %f %f');
        mbtime = C{3};

        if round2 (mbtime, 0.001) ~= round2 (time(end), 0.001)
            fprintf (1, 'times diverged at t = %f\n', mbtime);
            ind = ind - 1;
            convflag = true;
        else
            convflag = false;
            time(ind) = time(ind-1) + mbsys.problems{1}.timeStep;
        end
    else
        time(ind) = time(ind-1) + mbsys.problems{1}.timeStep;
    end
    
%     for ind = 1:nnodes
%     
%         fprintf (1, 'Node %d has label %d\n', ind, mb.KinematicsLabel (ind));
%     
%     end

%     [forces, out] = hsys.hydroForces ( time, ...
%                                        [output.bodies(1).position(ind,:)', output.bodies(2).position(ind,:)'], ...
%                                        [output.bodies(1).velocity(ind,:)', output.bodies(2).velocity(ind,:)'], ...
%                                        [output.bodies(1).acceleration(ind,:)', output.bodies(2).acceleration(ind,:)'], ...
%                                        output.wave.elevation(ind) );
%                                    
    
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
    
    [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] =  pto.forceFromFcn ();
    
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

    [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] =  pto.forceFromFcn ();
    
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


        [Fpto, ptoforce(ind), xRpto(ind) vRpto(ind)] = pto.forceFromFcn ();

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

    if plotvectors && time(ind) > 150
        
        if exist ('hvectplotax', 'var') && isvalid (hvectplotax)
            cla (hvectplotax);
        else
            figure;
            hvectplotax = axes;
        end
        
        olen = 2;
        
        vect.plotvec3 (unit (xRptoVec(:,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (unit (vRptoVec(:,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (olen * unit (FptoVec(1:3,1,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (unit (vRvec), [], 'PlotAxes', hvectplotax);
        
        v1 = vel(1:3,1,ind);
        v2 = vel(1:3,2,ind);
        v1magn = magn(v1);
        v2magn = magn(v2);
        
        v1vec = unit (v1) * v1magn / max(v1magn, v2magn);
        v2vec = unit (v2) * v2magn / max(v1magn, v2magn);
        vect.plotvec3 (v1vec, [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (v2vec, [], 'PlotAxes', hvectplotax);
        
        
        vox = [0.8*olen;0;0]' * om.orientationMatrix;
        voy = [0;0.8*olen;0]' * om.orientationMatrix;
        voz = [0;0;0.8*olen]' * om.orientationMatrix;
        
        vect.plotvec3 (vox', [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'k'});
        vect.plotvec3 (voy', [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'k'});
        vect.plotvec3 (voz', [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'k'});
        
        vgx = [0.3;0;0];
        vgy = [0;0.3;0];
        vgz = [0;0;0.3];
        
        vect.plotvec3 (vgx, [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'b'});
        vect.plotvec3 (vgy, [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'b'});
        vect.plotvec3 (vgz, [], 'PlotAxes', hvectplotax, 'Properties', {'Color', 'b'});
        
        axis square;
        view (3);
        xlabel ('x'); ylabel ('y'); zlabel ('z');
        legend ('xRptoVec', 'vRptoVec', 'FptoVec', 'vRvec', 'v1vec', 'v2vec', 'vox', 'voy', 'voz');
        
        set (hvectplotax, 'Xlim', [-2,2], 'Ylim', [-2,2], 'Zlim', [-2,2]);
        view (0,-1);
        drawnow;
    end
    
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

%%

if ~exist ('output', 'var')
    warning ('Not comparing output to original WEC-Sim as ''output'' variable is not in the workspace.')
else

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

    legend ( 'my F ExcitRamp x', ...
             'my F ExcitRamp y', ...
             'my F ExcitRamp z', ...
             'my M ExcitRamp x', ...
             'my M ExcitRamp y', ...
             'my M ExcitRamp z', ...
             ... 'F ViscousDamping',  ...
             'my F addedmass x',  ...
             'my F addedmass y',  ...
             'my F addedmass z',  ...
             'my M addedmass x',  ...
             'my M addedmass y',  ...
             'my M addedmass z',  ...
             ...
             'my F Restoring x',  ...
             'my F Restoring y',  ...
             'my F Restoring z',  ...
             'my M Restoring x',  ...
             'my M Restoring y',  ...
             'my M Restoring z',  ...
             ...
             'my F RadiationDamping x', ...
             'my F RadiationDamping y', ...
             'my F RadiationDamping z', ...
             'my M RadiationDamping x', ...
             'my M RadiationDamping y', ...
             'my M RadiationDamping z', ...
             ...
             'their F ExcitRamp x', ...
             'their F ExcitRamp y', ...
             'their F ExcitRamp z', ...
             'their M ExcitRamp x', ...
             'their M ExcitRamp y', ...
             'their M ExcitRamp z', ...
             ... 'F ViscousDamping',  ...
             'their F addedmass x',  ...
             'their F addedmass y',  ...
             'their F addedmass z',  ...
             'their M addedmass x',  ...
             'their M addedmass y',  ...
             'their M addedmass z',  ...
             ...
             'their F Restoring x',  ...
             'their F Restoring y',  ...
             'their F Restoring z',  ...
             'their M Restoring x',  ...
             'their M Restoring y',  ...
             'their M Restoring z',  ...
             ...
             'their F RadiationDamping x', ...
             'their F RadiationDamping y', ...
             'their F RadiationDamping z', ...
             'their M RadiationDamping x', ...
             'their M RadiationDamping y', ...
             'their M RadiationDamping z', ...
             ... 'ptoforce',  ...
             'vRpto' );


    %%

    mbout = mbdyn.postproc ( outputfile_prefix, mbsys ); 

    mbout.plotNodeTrajectories ('AxLims', [-1.5, 1.5; -1.5, 1.5; -25, 5]);


    %%
    mbout.animate ( 'DrawMode', 'solid', ...
                    'Light', true, ...
                    'skip', 5, ...
                    'AxLims', [-30, 30; -30, 30; -35, 35])


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

    % % figure;
    % % fcomp = 3;
    % % plot (output.bodies(1).time, [body1_F_AddedMass_Simulink.signals.values(:,fcomp), body2_F_AddedMass_Simulink.signals.values(:,fcomp), squeeze(F_AddedMass(:,fcomp,1:2))]);
    % % title (sprintf ('comp %d of body2\\_F\\_AddedMass\\_Simulink vs hydrobody F\\_addedmass for body 2', fcomp));
    % % legend (sprintf ('body 1 simulink comp %d', fcomp), ...
    % %     sprintf ('body 2 simulink comp %d', fcomp),...
    % %     sprintf ('body 1 hydrobody F_AddedMass comp %d', fcomp), ...
    % %     sprintf ('body 2 hydrobody F_AddedMass comp %d', fcomp))
    % % 
    % % figure;
    % % fcomp = 3;
    % % plot (output.bodies(1).time, [output.bodies(1).forceAddedMass(:,fcomp), output.bodies(2).forceAddedMass(:,fcomp), squeeze(F_AddedMassCorrected(:,fcomp,1:2))]);
    % % title (sprintf ('comp %d forceAddedMass vs F\\_AddedMassCorrected for both bodies', fcomp));
    % % legend (sprintf ('body 1 forceAddedMass post-processed comp %d', fcomp), ...
    % %     sprintf ('body 2 forceAddedMass post-processed comp %d', fcomp),...
    % %     sprintf ('body 1 hydrobody F_AddedMassCorrected comp %d', fcomp), ...
    % %     sprintf ('body 2 hydrobody F_AddedMassCorrected comp %d', fcomp))
    % 
    % data = [];
    % rowheadings = {};
    % 
    % stats = {'3', '5', '8'}
    % gf_forceTotal = nan * ones (numel (output.bodies), numel(stats));
    % gf_F_Total = gf_forceTotal;
    % gf_M_Total = gf_F_Total;
    % gf_forceExcitation = gf_forceTotal;
    % gf_forceAddedMass = gf_forceTotal;
    % gf_F_AddedMassCorrected = gf_forceTotal;
    % gf_forceRadiationDamping = gf_forceTotal;
    % gf_forceRestoring = gf_forceTotal;
    % 
    % for bodyind = 1:numel (output.bodies)
    %     
    %     % calculate some proper stats
    %     ...gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, forceTotal(:,:,bodyind));
    %     
    %     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,1:3), squeeze(F_Total(1:3,bodyind,:))', stats);
    %     
    %     gf_M_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,4:6), squeeze(F_Total(4:6,bodyind,:))', stats);
    %     
    %     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation, squeeze(F_ExcitRamp(:,bodyind,:))', stats);
    %     
    % %     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMass(:,:,bodyind));
    %     
    %     gf_F_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, squeeze(F_AddedMassCorrected(:,bodyind,:))', stats);
    %     
    %     gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping, squeeze(F_RadiationDamping(:,bodyind,:))', stats);
    %     
    %     gf_forceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring,  squeeze(F_Restoring(:,bodyind,:))', stats); 
    % 
    %     rowheadings = [rowheadings, {
    %                ...['forceTotal_body_', int2str(bodyind)], ...
    %                sprintf('Body %d Total Force', bodyind), ...
    %                sprintf('Body %d Total Moments', bodyind), ...
    %                sprintf('Body %d Excitation Force', bodyind), ...
    %                ...['forceAddedMass_body_', int2str(bodyind)], ...
    %                sprintf('Body %d Added Mass Force', bodyind), ...
    %                sprintf('Body %d Radiation and Damping Force', bodyind), ...
    %                sprintf('Body %d Hydrostatic Restoring Force', bodyind) }];
    % 
    %     data = [ data;
    %              ...gf_forceTotal(bodyind,:); 
    %              gf_F_Total(bodyind,:);
    %              gf_M_Total(bodyind,:);
    %              gf_forceExcitation(bodyind,:); 
    %              ...gf_forceAddedMass(bodyind,:);
    %              gf_F_AddedMassCorrected(bodyind,:);
    %              gf_forceRadiationDamping(bodyind,:); 
    %              gf_forceRestoring(bodyind,:)];
    % end
    % 
    % % display table
    % colheadings = { 'RMSE', 'MAE', 'R2' };
    % 
    % wid = 16;
    % % fms = {'d','.4f','.5E'};
    % fms = {};
    % fileID = 1;
    % 
    % displaytable (data,colheadings,wid,fms,rowheadings,fileID);
    % 
    % colheadings = [{'Force Description'}, colheadings];
    % fms = {'.3g','.3g','.2e'};
    % 
    % latextable (data, ...
    %             'ColumnHeadings', colheadings, ...
    %             'NumberWidth', wid, ...
    %             'FormatSpec', fms, ...
    %             'RowHeadings', rowheadings, ...
    %             'BookTabs', true)

    data = [];
    rowheadings = {};

    stats = {'3', '5', '8'}
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
        gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,1:3), squeeze(forces(1:3,bodyind,:))', stats);
        gf_momentTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,4:6), squeeze(forces(4:6,bodyind,:))', stats);

    %     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, squeeze(F_Total(:,bodyind,:))', stats);
        gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,1:3), squeeze(F_Total(1:3,bodyind,:))', stats);
        gf_M_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,4:6), squeeze(F_Total(4:6,bodyind,:))', stats);

    %     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
    %     gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);

        gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
        gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);

    %     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMass(:,:,bodyind));

        gf_F_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,1:3), squeeze(F_AddedMassCorrected(1:3,bodyind,:))', stats);
        gf_M_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,4:6), squeeze(F_AddedMassCorrected(4:6,bodyind,:))', stats);

        gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,1:3), squeeze(F_RadiationDamping(1:3,bodyind,:))', stats);
        gf_momentRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,4:6), squeeze(F_RadiationDamping(4:6,bodyind,:))', stats);

        gf_forceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,1:3),  squeeze(F_Restoring(1:3,bodyind,:))', stats); 
        gf_momentRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,4:6),  squeeze(F_Restoring(4:6,bodyind,:))', stats); 

        gf_pos(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,1:3),  squeeze(pos(1:3,bodyind,:))', stats); 
        gf_theta(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,4:6),  squeeze(pos(4:6,bodyind,:))', stats); 

        gf_vel(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,1:3),  squeeze(vel(1:3,bodyind,:))', stats); 
        gf_omega(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,4:6),  squeeze(vel(4:6,bodyind,:))', stats); 

        gf_accel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,1:3),  squeeze(accel(1:3,bodyind,:))', stats); 
        gf_omegaaccel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,4:6),  squeeze(accel(4:6,bodyind,:))', stats); 

        rowheadings = [rowheadings, {
                   ...['forceTotal_body_', int2str(bodyind)], ...
                   sprintf('Body %d Total Force Uncorrected', bodyind), ...
                   sprintf('Body %d Total Uncorrected', bodyind), ...
                   sprintf('Body %d Total Force', bodyind), ...
                   sprintf('Body %d Total Moments', bodyind), ...
                   sprintf('Body %d Excitation Force', bodyind), ...
                   sprintf('Body %d Excitation Moment', bodyind), ...
                   ...['forceAddedMass_body_', int2str(bodyind)], ...
                   sprintf('Body %d Added Mass Force', bodyind), ...
                   sprintf('Body %d Added Mass Moment', bodyind), ...
                   sprintf('Body %d Radiation and Damping Force', bodyind), ...
                   sprintf('Body %d Radiation and Damping Moment', bodyind), ...
                   sprintf('Body %d Hydrostatic Restoring Force', bodyind), ...
                   sprintf('Body %d Hydrostatic Restoring Moment', bodyind), ...
                   sprintf('Body %d Position', bodyind), ...
                   sprintf('Body %d Angular Position', bodyind), ...
                   sprintf('Body %d Velocity', bodyind), ...
                   sprintf('Body %d Angular Velocity', bodyind), ...
                   sprintf('Body %d Acceleration', bodyind), ...
                   sprintf('Body %d Angular Acceleration', bodyind) }];


        data = [ data;
                 ...gf_forceTotal(bodyind,:); 
                 gf_forceTotal(bodyind,:)
                 gf_momentTotal(bodyind,:)
                 gf_F_Total(bodyind,:);
                 gf_M_Total(bodyind,:);
                 gf_forceExcitation(bodyind,:);
                 gf_momentExcitation(bodyind,:);
                 ...gf_forceAddedMass(bodyind,:);
                 gf_F_AddedMassCorrected(bodyind,:);
                 gf_M_AddedMassCorrected(bodyind,:);
                 gf_forceRadiationDamping(bodyind,:);
                 gf_momentRadiationDamping(bodyind,:); 
                 gf_forceRestoring(bodyind,:);
                 gf_momentRestoring(bodyind,:);
                 gf_pos(bodyind,:); 
                 gf_theta(bodyind,:);
                 gf_vel(bodyind,:);
                 gf_omega(bodyind,:);
                 gf_accel(bodyind,:);
                 gf_omegaaccel(bodyind,:); ];
    end

    % display table
    colheadings = { 'RMSE', 'MAE', 'R2' };

    wid = 16;
    % fms = {'d','.4f','.5E'};
    fms = {};
    fileID = 1;

    displaytable (data,colheadings,wid,fms,rowheadings,fileID);

    colheadings = [{'Force Description'}, colheadings];
    fms = {'.2g','.2g','.2f'};

    latextable (data, ...
                'ColumnHeadings', colheadings, ...
                'NumberWidth', wid, ...
                'FormatSpec', fms, ...
                'RowHeadings', rowheadings, ...
                'BookTabs', true)
end
