% +test_wsim_RM3/run.m.m

%% first run the reference sim

test_wsim_RM3.run 

clear hsys mbsys float_hbody spar_hbody initptodpos mb

%% Hydrodynamic body system

% Float
float_hbody = wsim.hydrobody('RM3_NEMOH_output/RM3_NEMOH_output.h5', 'CaseDirectory', simu.caseDir);      
    %Create the wsim.hydrobody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
float_hbody.mass = 'equilibrium';                   
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
float_hbody.momOfInertia = [20907301, 21306090.66, 37085481.11];  %Moment of Inertia [kg*m^2]     
float_hbody.geometryFile = fullfile ('geometry', 'float.stl');    %Location of Geomtry File

% Spar/Plate
spar_hbody = wsim.hydrobody('RM3_NEMOH_output/RM3_NEMOH_output.h5', 'CaseDirectory', simu.caseDir); 
spar_hbody.mass = 'equilibrium';                   
spar_hbody.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody.geometryFile = fullfile ('geometry', 'plate.stl'); 

% make a hydrosys object for simulation
hsys = wsim.hydrosys (waves, simu, [float_hbody, spar_hbody]);

% set up transient simulation
hsys.initialiseHydrobodies ();
hsys.odeSimSetup ();
[hydro_mbnodes, hydro_mbbodies, hydro_mbelements] = hsys.makeMBDynComponents ();

%% Multibody dynamics system specification (mbdyn)

% note problem_options comes from test_wsim_RM3.run 
[mbsys, initptodpos] = test_wsim_RM3.make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options);
                     
% draw it
% mbsys.draw ('Mode', 'wireghost', 'Light', false);

mbsys.draw ( 'Mode', 'solid', ...
             'Light', true, ...
             'AxLims', [-30, 30; -30, 30; -35, 35], ...
             'Joints', false, ...
             'StructuralNodes', false)

mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

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

nemoh_eul = zeros (3,nnodes);

R = mb.GetRot();
for Rind = 1:size (R,3)
    om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
    nemoh_eul(1:3,Rind) = om.euler123();
end

nemoh_pos = [ mb.NodePositions(); 
        nemoh_eul ];
    
nemoh_vel = [ mb.NodeVelocities(); 
        mb.NodeOmegas() ];
    
nemoh_accel = [ mb.NodeAccelerations(); 
          mb.NodeAngularAccels() ];

[nemoh_forces, out] = hsys.hydroForces (time, nemoh_pos, nemoh_vel, nemoh_accel);
                               
% [forces, out] = hsys.hydroForces ( time, ...
%                                    [output.bodies(1).position(ind,:)', output.bodies(2).position(ind,:)'], ...
%                                    [output.bodies(1).velocity(ind,:)', output.bodies(2).velocity(ind,:)'], ...
%                                    [output.bodies(1).acceleration(ind,:)', output.bodies(2).acceleration(ind,:)'], ...
%                                    );

% set the forces
mb.F (nemoh_forces(1:3,:));
mb.M (nemoh_forces(4:6,:));

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

NEMOH_F_ExcitLin = out.F_ExcitLin;
NEMOH_F_ViscousDamping = out.F_ViscousDamping;
NEMOH_F_AddedMass = out.F_AddedMass;
NEMOH_F_Restoring = out.F_Restoring;
NEMOH_F_RadiationDamping = out.F_RadiationDamping;
NEMOH_F_ExcitNonLin = out.F_ExcitNonLin;
NEMOH_F_MorrisonElement = out.F_MorrisonElement;
NEMOH_F_Excit = out.F_Excit;
NEMOH_F_ExcitRamp = out.F_ExcitRamp;
FptoVec = [0;0;0];

% accept the last data into the time history of solutions
hsys.advanceStep (time(end), nemoh_vel, nemoh_accel);
    
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
nemoh_xRpto = 0;
nemoh_vRpto = 0;
nemoh_xRptoVec = [0;0;0];
nemoh_vRptoVec = [0;0;0];

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
        nemoh_eul(1:3,Rind) = om.euler123();
    end

    nemoh_pos(:,:,ind) = [ mb.NodePositions(); nemoh_eul];
    nemoh_vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
    nemoh_accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

    [hydroforces, out] = hsys.hydroForces (time(ind), nemoh_pos(:,:,ind), nemoh_vel(:,:,ind), nemoh_accel(:,:,ind));

    nemoh_forces (:,:,ind) = hydroforces;
    
	% calculate spring damping PTO force here
    xRvec = nemoh_pos(1:3,1,ind) - nemoh_pos(1:3,2,ind) - initptodpos;
    vRvec = nemoh_vel(1:3,1,ind) - nemoh_vel(1:3,2,ind);

    nemoh_xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
    nemoh_vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;
    
    % velocities and displacements are the z components of the vectors in
    % the pto coordinate system
    nemoh_xRpto(ind) = nemoh_xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
    nemoh_vRpto(ind) = nemoh_vRptoVec(3,ind); % magn (vRgenvec) ;
    
    ptoforce(ind) = -k*nemoh_xRpto(ind) -c*nemoh_vRpto(ind);
    
%     FptoVec = om.orientationMatrix * [0; 0; -ptoforce(ind)];
    FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * om.orientationMatrix)' ;
    
    nemoh_forces (1:3,1,ind) = nemoh_forces (1:3,1,ind) + FptoVec(1:3,1,ind);
    nemoh_forces (1:3,2,ind) = nemoh_forces (1:3,2,ind) - FptoVec(1:3,1,ind);
    
	mb.F (nemoh_forces(1:3,:,ind));
    mb.M (nemoh_forces(4:6,:,ind));
    
%     f = dir([ outputfile_prefix, '.out']);
%     outsize = f.bytes;
    
    mbconv = mb.applyForcesAndMoments (false);
    
%     pause (1);
%     
%     f = dir([ outputfile_prefix, '.out']);
%     
%     if outsize == f.bytes
%         fprintf (1, 'outfile size did not change size at time t = %f\n', time(ind));
%     end
%     
%     fprintf (1, '    time: %f, ind %d, mbconv: %d\n', time(ind), ind, mbconv);

    status = mb.GetMotion ();
    
    if status ~= 0
                
        NEMOH_F_ExcitLin(:,:,ind) = out.F_ExcitLin;
        NEMOH_F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
        NEMOH_F_AddedMass(:,:,ind) = out.F_AddedMass;
        NEMOH_F_Restoring(:,:,ind) = out.F_Restoring;
        NEMOH_F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
        NEMOH_F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
        NEMOH_F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
        NEMOH_F_Excit(:,:,ind) = out.F_Excit;
        NEMOH_F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
        
        ind = ind + 1;
        
        continue;
    end
    
    R = mb.GetRot();
    
    for Rind = 1:size (R,3)
        om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
        nemoh_eul(1:3,Rind) = om.euler123();
    end

    nemoh_pos(:,:,ind) = [ mb.NodePositions(); nemoh_eul];
    nemoh_vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
    nemoh_accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

    [newhydroforces, out] = hsys.hydroForces (time(ind), nemoh_pos(:,:,ind), nemoh_vel(:,:,ind), nemoh_accel(:,:,ind));

    nemoh_forces (:,:,ind) = newhydroforces;
    
	% calculate spring damping PTO force here
    xRvec = nemoh_pos(1:3,1,ind) - nemoh_pos(1:3,2,ind) - initptodpos;
    vRvec = nemoh_vel(1:3,1,ind) - nemoh_vel(1:3,2,ind);

    nemoh_xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
    nemoh_vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;
    
    % velocities and displacements are the z components of the vectors in
    % the pto coordinate system
    nemoh_xRpto(ind) = nemoh_xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
    nemoh_vRpto(ind) = nemoh_vRptoVec(3,ind); % magn (vRgenvec) ;
    
    ptoforce(ind) = -k*nemoh_xRpto(ind) - c*nemoh_vRpto(ind);
    
    FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * R(:,:,1))';
    
    nemoh_forces (1:3,1,ind) = nemoh_forces (1:3,1,ind) + FptoVec(1:3,1,ind);
    nemoh_forces (1:3,2,ind) = nemoh_forces (1:3,2,ind) - FptoVec(1:3,1,ind);
    
	mb.F (nemoh_forces(1:3,:,ind));
    mb.M (nemoh_forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (false);
    
%     fprintf (1, '    time: %f, mbconv: %d\n', time(ind), mbconv);
    
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
        
%         fprintf (1, '    time: %f, ind: %d, iterating, mbconv: %d, iteration: %d\n', time(ind), ind, mbconv, itercount);
            
        hydroforces = newhydroforces;
        
        status = mb.GetMotion ();
        
        if status ~= 0
            break;
        end
        
        R = mb.GetRot();
    
        for Rind = 1:size (R,3)
            om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
            nemoh_eul(1:3,Rind) = om.euler123();
        end

        nemoh_pos(:,:,ind) = [ mb.NodePositions(); nemoh_eul];
        nemoh_vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
        nemoh_accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];

        [newhydroforces, out] = hsys.hydroForces (time(ind), nemoh_pos(:,:,ind), nemoh_vel(:,:,ind), nemoh_accel(:,:,ind));

        nemoh_forces (:,:,ind) = newhydroforces;

        % calculate spring damping PTO force here
        xRvec = nemoh_pos(1:3,1,ind) - nemoh_pos(1:3,2,ind) - initptodpos;
        vRvec = nemoh_vel(1:3,1,ind) - nemoh_vel(1:3,2,ind);
        
        nemoh_xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
        nemoh_vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;

        % velocities and displacements are the z components of the vectors in
        % the pto coordinate system
        nemoh_xRpto(ind) = nemoh_xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
        nemoh_vRpto(ind) = nemoh_vRptoVec(3,ind); % magn (vRgenvec) ;

        ptoforce(ind) = -k*nemoh_xRpto(ind) -c*nemoh_vRpto(ind);

        FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * R(:,:,1))' ;

        nemoh_forces (1:3,1,ind) = nemoh_forces (1:3,1,ind) + FptoVec(1:3,1,ind);
        nemoh_forces (1:3,2,ind) = nemoh_forces (1:3,2,ind) - FptoVec(1:3,1,ind);

        mb.F (nemoh_forces(1:3,:,ind));
        mb.M (nemoh_forces(4:6,:,ind));

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
        
        NEMOH_F_ExcitLin(:,:,ind) = out.F_ExcitLin;
        NEMOH_F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
        NEMOH_F_AddedMass(:,:,ind) = out.F_AddedMass;
        NEMOH_F_Restoring(:,:,ind) = out.F_Restoring;
        NEMOH_F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
        NEMOH_F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
        NEMOH_F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
        NEMOH_F_Excit(:,:,ind) = out.F_Excit;
        NEMOH_F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
        
        ind = ind + 1;
    
        break;
    end

    mb.F (nemoh_forces(1:3,:,ind));
    mb.M (nemoh_forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (true);
    
%     fprintf (1, 'time: %f, tind: %d, final status: %d\n', time(ind), ind, status);

    if plotvectors && time(ind) > 150
        
        if exist ('hvectplotax', 'var') && isvalid (hvectplotax)
            cla (hvectplotax);
        else
            figure;
            hvectplotax = axes;
        end
        
        olen = 2;
        
        vect.plotvec3 (unit (nemoh_xRptoVec(:,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (unit (nemoh_vRptoVec(:,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (olen * unit (FptoVec(1:3,1,ind)), [], 'PlotAxes', hvectplotax);
        vect.plotvec3 (unit (vRvec), [], 'PlotAxes', hvectplotax);
        
        v1 = nemoh_vel(1:3,1,ind);
        v2 = nemoh_vel(1:3,2,ind);
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
    
    NEMOH_F_ExcitLin(:,:,ind) = out.F_ExcitLin;
    NEMOH_F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
    NEMOH_F_AddedMass(:,:,ind) = out.F_AddedMass;
    NEMOH_F_Restoring(:,:,ind) = out.F_Restoring;
    NEMOH_F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
    NEMOH_F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
    NEMOH_F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
    NEMOH_F_Excit(:,:,ind) = out.F_Excit;
    NEMOH_F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
    
    % accept the last data into the time history of solutions
    hsys.advanceStep (time(end), nemoh_vel(:,:,ind), nemoh_accel(:,:,ind));
    
    ind = ind + 1;
    
end

[NEMOH_F_Total, NEMOH_F_AddedMassCorrected] = correctAddedMassForce (hsys, nemoh_forces, NEMOH_F_AddedMass, nemoh_accel);
toc;
clear mb;

fprintf (1, 'Reached time %f, in %d steps\n', time(end), ind-1);

%%
figure;
tmin = 0;
tmax = 400;
plotinds =  time>=tmin & time<=tmax;
plotyy (time(plotinds), [ squeeze(forces(1:3,1,plotinds))',  ptoforce(plotinds)'], time(plotinds), vRpto(plotinds));
legend ('fx', 'fy', 'fz', 'ptoforce', 'vRpto');

%%

figure;
tmin = 0;
tmax = 400;

bodyind = 2;

tmax = min ( [tmax, time(end), time(end) ]);
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
        ...output.bodies(bodyind).forceExcitation(plotinds,:), ...
        ...output.bodies(bodyind).forceAddedMass(plotinds,:), ...
        ...output.bodies(bodyind).forceRestoring(plotinds,:), ...
        ...output.bodies(bodyind).forceRadiationDamping(plotinds,:), ...
        squeeze(NEMOH_F_ExcitRamp(:,bodyind,plotinds))', ...
        ... squeeze(F_ViscousDamping(3,bodyind,plotinds)), ...
        squeeze(NEMOH_F_AddedMassCorrected(:,bodyind,plotinds))', ...
        squeeze(NEMOH_F_Restoring(:,bodyind,plotinds))', ...
        squeeze(NEMOH_F_RadiationDamping(:,bodyind,plotinds))', ...
        ], ...
        time(plotinds), vRpto(plotinds) );

legend ( 'ref F ExcitRamp x', ...
         'ref F ExcitRamp y', ...
         'ref F ExcitRamp z', ...
         'ref M ExcitRamp x', ...
         'ref M ExcitRamp y', ...
         'ref M ExcitRamp z', ...
         ... 'F ViscousDamping',  ...
         'ref F addedmass x',  ...
         'ref F addedmass y',  ...
         'ref F addedmass z',  ...
         'ref M addedmass x',  ...
         'ref M addedmass y',  ...
         'ref M addedmass z',  ...
         ...
         'ref F Restoring x',  ...
         'ref F Restoring y',  ...
         'ref F Restoring z',  ...
         'ref M Restoring x',  ...
         'ref M Restoring y',  ...
         'ref M Restoring z',  ...
         ...
         'ref F RadiationDamping x', ...
         'ref F RadiationDamping y', ...
         'ref F RadiationDamping z', ...
         'ref M RadiationDamping x', ...
         'ref M RadiationDamping y', ...
         'ref M RadiationDamping z', ...
         ...
         'nemoh F ExcitRamp x', ...
         'nemoh F ExcitRamp y', ...
         'nemoh F ExcitRamp z', ...
         'nemoh M ExcitRamp x', ...
         'nemoh M ExcitRamp y', ...
         'nemoh M ExcitRamp z', ...
         ... 'F ViscousDamping',  ...
         'nemoh F addedmass x',  ...
         'nemoh F addedmass y',  ...
         'nemoh F addedmass z',  ...
         'nemoh M addedmass x',  ...
         'nemoh M addedmass y',  ...
         'nemoh M addedmass z',  ...
         ...
         'nemoh F Restoring x',  ...
         'nemoh F Restoring y',  ...
         'nemoh F Restoring z',  ...
         'nemoh M Restoring x',  ...
         'nemoh M Restoring y',  ...
         'nemoh M Restoring z',  ...
         ...
         'nemoh F RadiationDamping x', ...
         'nemoh F RadiationDamping y', ...
         'nemoh F RadiationDamping z', ...
         'nemoh M RadiationDamping x', ...
         'nemoh M RadiationDamping y', ...
         'nemoh M RadiationDamping z', ...
         ... 'ptoforce',  ...
         'vRpto' );


%%

mbout = mbdyn.postproc ( outputfile_prefix, mbsys ); 

mbout.plotNodeTrajectories ('AxLims', [-1.5, 1.5; -1.5, 1.5; -25, 5]);


%%
% mbout.animate ( 'DrawMode', 'solid', ...
%                 'Light', true, ...
%                 'skip', 5, ...
%                 'AxLims', [-30, 30; -30, 30; -35, 35])


%% Compare to reference WEC-Sim

%
% Requires to have run the RM3 example first using original WEC-Sim

doplot = false;

if doplot
    for bodyind = 1:numel (output.bodies)

        figure;
        plot ( time, squeeze(forces(:,bodyind,:)), ...
               time, squeeze(nemoh_forces(:,bodyind,:)) );
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
        figure;
        plot ( time, squeeze(F_Total(:,bodyind,:)),  ...
               time, squeeze(NEMOH_F_Total(:,bodyind,:)) );
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
        title (sprintf ('ref F\\_Total vs NEMOH\\_F\\_Total for body %d', bodyind));
%         figure;
%         plot (time, output.bodies(bodyind).forceExcitation,  time, squeeze(F_ExcitRamp(:,bodyind,:))); 
        figure;
        plot (time, squeeze(F_AddedMass(:,bodyind,:)),  ...
              time, squeeze(NEMOH_F_AddedMass(:,bodyind,:)) );
        title (sprintf ('F\\_AddedMass vs NEMOH\\_F\\_AddedMass for body %d', bodyind));
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
        figure;
        plot (time, output.bodies(bodyind).forceAddedMass,  ...
              time, squeeze(NEMOH_F_AddedMassCorrected(:,bodyind,:)));
        title (sprintf ('forceAddedMass vs F\\_AddedMassCorrected for body %d', bodyind));
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
        figure;
        plot (time, squeeze(F_RadiationDamping(:,bodyind,:)),  ...
              time, squeeze(NEMOH_F_RadiationDamping(:,bodyind,:)) );
        title (sprintf ('F\\_RadiationDamping vs NEMOH\\_F\\_RadiationDamping for body %d', bodyind));
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%         figure;
%         plot (time, output.bodies(bodyind).forceRestoring,  time, squeeze(F_Restoring(:,bodyind,:)));
        figure;
        plot ( time, squeeze(pos(:,bodyind,:)),  ...
               time, squeeze(nemoh_pos(:,bodyind,:)) );
        title (sprintf ('pos vs pnemoh_os for body %d', bodyind, bodyind));
        legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
    end

end

tmin = 0;
tmax = 400;

bodyind = 1;

tmax = min ( [tmax, time(end), time(end)]);
plotinds =  time>=tmin & time<=tmax;

% figure;
% plotyy ( time(plotinds), ...
%          [output.ptos(1).forceTotal(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%          time(plotinds), ...
%          [vRptoVec(plotinds)', [ output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3)] ] );
% title (sprintf ('output.ptos(%d).forceTotal vs FptoVec for body %d', bodyind, bodyind));
% legend ('1', '2', '3', '4', '5', '6', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
% figure;
% plotyy ( time(plotinds), ...
%         [output.ptos(1).forceTotalWorld(plotinds,:) - output.ptos(1).forceConstraintWorld(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%         time(plotinds), ...
%         [vRptoVec(:,plotinds)', output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3) ] );
% title (sprintf ('output.ptos(1).forceTotalWorld - output.ptos(1).forceConstraintWorld vs FptoVec for body %d', bodyind, bodyind));
% legend ('1', '2', '3', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
% figure;
% plot ( time(plotinds), ...
%        output.ptos(1).forceInternalMechanics(plotinds,3),  ...
%        time(plotinds), ...
%        squeeze(ptoforce(plotinds)));
% title (sprintf ('output.ptos(%d).forceInternalMechanics(:,3) vs ptoforce for body %d', bodyind, bodyind));
% legend ('forceInternalMechanics(:,3)', 'ptoforce');

%% Statistical comparison

data = [];
rowheadings = {};

nbodies = 2;

stats = {'3', '5', '8'}
gf_forceTotal = nan * ones (nbodies, numel(stats));
gf_momentTotal = nan * ones (nbodies, numel(stats));
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

for bodyind = 1:nbodies

    % calculate some proper stats
    gf_forceTotal(bodyind,:) = gfit2 (squeeze(forces(1:3,bodyind,:))', squeeze(nemoh_forces(1:3,bodyind,:))', stats);
    gf_momentTotal(bodyind,:) = gfit2 (squeeze(forces(4:6,bodyind,:))', squeeze(nemoh_forces(4:6,bodyind,:))', stats);

%     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, squeeze(F_Total(:,bodyind,:))', stats);
    gf_F_Total(bodyind,:) = gfit2 (squeeze(F_Total(1:3,bodyind,:))', squeeze(NEMOH_F_Total(1:3,bodyind,:))', stats);
    gf_M_Total(bodyind,:) = gfit2 (squeeze(F_Total(4:6,bodyind,:))', squeeze(NEMOH_F_Total(4:6,bodyind,:))', stats);

%     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
%     gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);

    gf_forceExcitation(bodyind,:) = gfit2 (squeeze(F_ExcitRamp(1:3,bodyind,:))', squeeze(NEMOH_F_ExcitRamp(1:3,bodyind,:))', stats);
    gf_momentExcitation(bodyind,:) = gfit2 (squeeze(F_ExcitRamp(4:6,bodyind,:))', squeeze(NEMOH_F_ExcitRamp(4:6,bodyind,:))', stats);

%     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMass(:,:,bodyind));

    gf_F_AddedMassCorrected(bodyind,:) = gfit2 ( squeeze(F_AddedMassCorrected(1:3,bodyind,:))', squeeze(NEMOH_F_AddedMassCorrected(1:3,bodyind,:))', stats);
    gf_M_AddedMassCorrected(bodyind,:) = gfit2 (squeeze(F_AddedMassCorrected(4:6,bodyind,:))', squeeze(NEMOH_F_AddedMassCorrected(4:6,bodyind,:))', stats);

    gf_forceRadiationDamping(bodyind,:) = gfit2 (squeeze(F_RadiationDamping(1:3,bodyind,:))', squeeze(NEMOH_F_RadiationDamping(1:3,bodyind,:))', stats);
    gf_momentRadiationDamping(bodyind,:) = gfit2 (squeeze(F_RadiationDamping(4:6,bodyind,:))', squeeze(NEMOH_F_RadiationDamping(4:6,bodyind,:))', stats);

    gf_forceRestoring(bodyind,:) = gfit2 (squeeze(F_Restoring(1:3,bodyind,:))',  squeeze(NEMOH_F_Restoring(1:3,bodyind,:))', stats); 
    gf_momentRestoring(bodyind,:) = gfit2 (squeeze(F_Restoring(4:6,bodyind,:))',  squeeze(NEMOH_F_Restoring(4:6,bodyind,:))', stats); 

    gf_pos(bodyind,:) = gfit2 (squeeze(pos(1:3,bodyind,:))',  squeeze(nemoh_pos(1:3,bodyind,:))', stats); 
    gf_theta(bodyind,:) = gfit2 (squeeze(pos(4:6,bodyind,:))',  squeeze(nemoh_pos(4:6,bodyind,:))', stats); 

    gf_vel(bodyind,:) = gfit2 (squeeze(vel(1:3,bodyind,:))',  squeeze(nemoh_vel(1:3,bodyind,:))', stats); 
    gf_omega(bodyind,:) = gfit2 (squeeze(vel(4:6,bodyind,:))',  squeeze(nemoh_vel(4:6,bodyind,:))', stats); 

    gf_accel(bodyind,:) = gfit2 (squeeze(accel(1:3,bodyind,:))',  squeeze(nemoh_accel(1:3,bodyind,:))', stats); 
    gf_omegaaccel(bodyind,:) = gfit2 (squeeze(accel(4:6,bodyind,:))',  squeeze(nemoh_accel(4:6,bodyind,:))', stats); 

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
            'BookTabs', true);
