% +test_wsim_RM3/run.m.m

clear waves simu hsys mbsys float_hbody spar_hbody ref_float ref_spar mb

%% Hydro Simulation Data
simu = wsim.simSettings (getmfilepath ('test_wsim_RM3.run'));    %Create the Simulation Variable
% simu.mode = 'normal';                 %Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
% simu.explorer='on';                   %Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                   %Simulation Start Time [s]
simu.endTime=400;                       %Simulation End bdcloseTime [s]
simu.solver = 'ode4';                   %simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.25; % 0.1; 							%Simulation time-step [s]
simu.rampT = 100;                       %Wave Ramp Time Length [s]
simu.multibodySolver = 'MBDyn';
simu.b2b = 1;

%% Wave Information 

% noWaveCIC, no waves with radiation CIC  
% waves = wsim.waveSettings ('noWaveCIC');       %Create the Wave Variable and Specify Type      
%%%%%%%%%%%%%%%%%%%

% Regular Waves
waves = wsim.waveSettings ('regularCIC');        %Create the Wave Variable and Specify Type                               
waves.H = 2.5;                          %Wave Height [m]
waves.T = 8;                            %Wave Period [s]
%%%%%%%%%%%%%%%%%%%

% Irregular Waves using PM Spectrum with Convolution Integral Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'PM';
%%%%%%%%%%%%%%%%%%%

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

% set up ode simulation
hsys.initialiseHydrobodies ();
hsys.timeDomainSimSetup ();

%% Multibody dynamics system specification (mbdyn)

gref = mbdyn.pre.globalref;

ref_seabed = mbdyn.pre.reference ( [0;0;-waves.waterDepth], [], [], [], 'Parent', gref); 

% Float
ref_float = mbdyn.pre.reference ( float_hbody.cg, [], [], [], 'Parent', gref);

float_node = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', ref_float.pos);

float_mb_body = mbdyn.pre.body ( float_hbody.mass,  ...
                                 [0;0;0], ...
                                 diag (float_hbody.momOfInertia), ...
                                 float_node, ...
                                 'STLFile', fullfile (simu.caseDir, float_hbody.geometryFile) );

float_mb_body.setSize (20, 20, 5);

% make it a yellow colour
float_mb_body.setColour ([247, 244, 17]./255);

% Spar
ref_spar = mbdyn.pre.reference ( spar_hbody.cg, [], [], [], 'Parent', gref);

spar_node = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', ref_spar.pos);

spar_mb_body = mbdyn.pre.body ( spar_hbody.mass, ...
                                [0;0;0], ...
                                diag (spar_hbody.momOfInertia), ...
                                spar_node, ...
                                'STLFile', fullfile (simu.caseDir, spar_hbody.geometryFile) );

spar_mb_body.setSize (6, 6, 38);

% make it a grey colour
spar_mb_body.setColour ([181, 181, 181]./255);


% create a node to clamp to the global frame which we can use to attach a
% planar joint (actually made up of two joints) to the spar.
%
% Since  we are going to clamp the node, we can make it static (getting rid
% of the degrees of freedom associated with it)
clamped_node = mbdyn.pre.structuralNode6dof ('static', 'AbsolutePosition', ref_seabed.pos);
% clamp it
jclamp = mbdyn.pre.clamp (clamped_node, 'node', 'node');

% add the joints
j1 = mbdyn.pre.prismatic (float_node, spar_node);

om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1., 0., 0.], 'ib', 3, 'vecB', [0., 0., 1.]));
j2 = mbdyn.pre.inline (float_node, spar_node, 'RelativeLinePosition', [0;0;0], ...
                                              'RelativeOrientation', om, ...
                                              'OrientationReference', 'global');

j3 = mbdyn.pre.inPlane (  clamped_node, spar_node, [0;1;0], ...
                        'RelativeOffset', 'null', ...
                        'RelativeOffsetReference', 'other node');
                    
j4 = mbdyn.pre.revoluteRotation ( spar_node, clamped_node, ...
                                 'RelativeOffset1', 'null', ...
                                 'RelativeOffset1Reference', 'node', ...
                                 'RelativeOffset2', 'null', ...
                                 'RelativeOffset2Reference', 'other node', ...
                                 'RelativeOrientation1', mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1,0,0], 'ib', 2, 'vecB', [0,0,1])), ...
                                 'RelativeOrientation2', mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1,0,0], 'ib', 2, 'vecB', [0,0,1])));

% add the socket forces
if ispc
    % add the socket forces
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

socket_force = mbdyn.pre.externalStructuralForce ({float_node, spar_node}, [], communicator, ...
                                                  'Labels', 'no', ...
                                                  'Sorted', 'yes', ...
                                                  'Orientation', 'orientation matrix', ...
                                                  'Echo', fullfile (simu.caseDir, 'output', 'socket_echo.txt'));

prob = mbdyn.pre.initialValueProblem (simu.startTime, simu.endTime, simu.dt, ...
                                'ResidualTol', 1e-6, 'MaxIterations', 200);

% assemble the system
mbsys = mbdyn.pre.system ( prob, ...
                           'Nodes', {float_node, spar_node, clamped_node}, ...
                           'Elements', {float_mb_body, spar_mb_body, jclamp, j1, j2, j3, j4 socket_force} );
                       
% draw it
setSize (float_node, 10, 10, 10);
setSize (spar_node, 10, 10, 10);
mbsys.draw ('Mode', 'solid', 'Light', true, 'StructuralNodes', false, 'AxLims', [-35, 35; -35, 35; -35, 35]);

mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

%% Set up PTO
% 
% k = 0;
% c = 1200000;
initptodpos = ref_float.pos - ref_spar.pos; 

% set up the nested generator simulation

load ('/home/rcrozier/Sync/work/enercro/Projects/WES-Wavedrive/Outputs/design_004_converted_to_rnfoundry.mat');

simoptions = struct ();

simoptions.set_CoilResistance = 0.497091;

simoptions.GetVariableGapForce = false;
simoptions.MagFEASimType = 'multiple';
simoptions.UseParFor = false;
simoptions.MagFEASim.UseParFor = false;
% simoptions.AddPhaseCurrentsComponents = true;

design.RlVRp = 4;
design = completedesign_TM_SLOTLESS (design, simoptions);
design.zs = design.zp / 3;

% Set up Common Parameters
simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;

[design, simoptions] = simfun_TM_SLOTLESS(design, simoptions);
[design, simoptions] = finfun_TM_SLOTLESS(design, simoptions);

simoptions.ODESim.ForceFcn = 'forcefcn_linear_pscbmot';
simoptions.ODESim.ForceFcnArgs = {};

innert0 = 0;
innerx0 = [0,0,0];

innerodeevfcn = @(t, y, mc, flag) nestedodeforcefcn_linear (t, y, mc, flag, design, simoptions);

% initial displacement and velocity of the generator
interpdat = [0; 0];

machinesolver = ode.odesolver ( innerodeevfcn, innert0, innerx0, interpdat, odeset (), ...
                    'Solver', 'ode15s', ...
                    'SaveSolutions', false ...
                    ... 'SplitFcn', @(flag, results, sol, mc, evalfcn) nestedsysresults_linear (flag, results, sol, mc, evalfcn, design, simoptions) ...
                              );
                              
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

nnodes = mb.GetNodes ();

time = prob.initialTime;

status = mb.GetMotion ();

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

verbose = false;
plotvectors = false;
checkoutputs = false;
miniters = 0;
maxiters = prob.maxIterations;
forcetol = 1;

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
            if verbose
                fprintf (1, 'times diverged at t = %f\n', mbtime);
            end
            ind = ind - 1;
            convflag = true;
        else
            convflag = false;
            time(ind) = time(ind-1) + prob.timeStep;
        end
    else
        time(ind) = time(ind-1) + prob.timeStep;
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
    
	% calculate spring damping PTO force here
    xRvec = pos(1:3,1,ind) - pos(1:3,2,ind) - initptodpos;
    vRvec = vel(1:3,1,ind) - vel(1:3,2,ind);

    xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
    vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;
    
    % velocities and displacements are the z components of the vectors in
    % the pto coordinate system
    xRpto(ind) = xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
    vRpto(ind) = vRptoVec(3,ind); % magn (vRgenvec) ;
    
%     ptoforce(ind) = -k*xRpto(ind) -c*vRpto(ind);
    
    machinesolver.solve (time(end), [xRpto(ind); vRpto(ind)]);
    
    results = machinesolver.getOutputs ();
    
    ptoforce(ind) = results.Fpto;
    
%     FptoVec = om.orientationMatrix * [0; 0; -ptoforce(ind)];
    FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * om.orientationMatrix)' ;
    
    forces (1:3,1,ind) = forces (1:3,1,ind) + FptoVec(1:3,1,ind);
    forces (1:3,2,ind) = forces (1:3,2,ind) - FptoVec(1:3,1,ind);
    
	mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));
    
    mbconv = mb.applyForcesAndMoments (false);
    
    status = mb.GetMotion ();
    
    if status ~= 0
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
    
	% calculate spring damping PTO force here
    xRvec = pos(1:3,1,ind) - pos(1:3,2,ind) - initptodpos;
    vRvec = vel(1:3,1,ind) - vel(1:3,2,ind);

    xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
    vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;
    
    % velocities and displacements are the z components of the vectors in
    % the pto coordinate system
    xRpto(ind) = xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
    vRpto(ind) = vRptoVec(3,ind); % magn (vRgenvec) ;
    
    machinesolver.solve (time(end), [xRpto(ind); vRpto(ind)]);
    
    results = machinesolver.getOutputs ();
    
    ptoforce(ind) = results.Fpto;
%     ptoforce(ind) = -k*xRpto(ind) - c*vRpto(ind);
    
    FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * R(:,:,1))';
    
    forces (1:3,1,ind) = forces (1:3,1,ind) + FptoVec(1:3,1,ind);
    forces (1:3,2,ind) = forces (1:3,2,ind) - FptoVec(1:3,1,ind);
    
	mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (false);
    
%     fprintf (1, '    time: %f, mbconv: %d\n', time(ind), mbconv);
    forcediff = abs (hydroforces - newhydroforces);
    
    itercount = 1;
    while mbconv ~= 0 || itercount < miniters || (max (forcediff(:)) > forcetol)
        
        if verbose
            fprintf (1, '    time: %f, ind: %d, iterating, mbconv: %d, iteration: %d, max forcediff: %g\n', time(ind), ind, mbconv, itercount, max (forcediff(:)));
        end
        
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

        % calculate spring damping PTO force here
        xRvec = pos(1:3,1,ind) - pos(1:3,2,ind) - initptodpos;
        vRvec = vel(1:3,1,ind) - vel(1:3,2,ind);
        
        xRptoVec(1:3,ind) = R(:,:,1).' * xRvec;
        vRptoVec(1:3,ind) = R(:,:,1).' * vRvec;

        % velocities and displacements are the z components of the vectors in
        % the pto coordinate system
        xRpto(ind) = xRptoVec(3,ind); % magn (pos(:,2) - pos(:,1));
        vRpto(ind) = vRptoVec(3,ind); % magn (vRgenvec) ;

        machinesolver.solve (time(end), [xRpto(ind); vRpto(ind)]);
    
        results = machinesolver.getOutputs ();
        
        ptoforce(ind) = results.Fpto;
        
%         ptoforce(ind) = -k*xRpto(ind) -c*vRpto(ind);

        FptoVec(1:3,1,ind) = ([0; 0; ptoforce(ind)]' * R(:,:,1))' ;

        forces (1:3,1,ind) = forces (1:3,1,ind) + FptoVec(1:3,1,ind);
        forces (1:3,2,ind) = forces (1:3,2,ind) - FptoVec(1:3,1,ind);

        mb.F (forces(1:3,:,ind));
        mb.M (forces(4:6,:,ind));

        mbconv = mb.applyForcesAndMoments (false);
        
        itercount = itercount + 1;
        
        if itercount > maxiters
            error ('mbdyn iterations exceeded max allowed');
        end
        
        forcediff = abs (hydroforces - newhydroforces);
    
    end
    
    status = mb.GetMotion ();
        
    if status ~= 0
        break;
    end

    mb.F (forces(1:3,:,ind));
    mb.M (forces(4:6,:,ind));

    mbconv = mb.applyForcesAndMoments (true);
    
    if verbose
        fprintf (1, 'time: %f, tind: %d, final status: %d\n', time(ind), ind, status);
    end

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
    hsys.advanceStep (time(end),  vel(:,:,ind), accel(:,:,ind));
    
    % get the full output from the last generator sim
    
    results = machinesolver.getOutputs ({}, true);
    
    fnames = fieldnames (results);
    
    for fnind = 1:numel (fnames)
        if ind == 2
            genresults = results;
        else
            genresults.(fnames{fnind}) = [ genresults.(fnames{fnind}); results.(fnames{fnind})(2:end,:) ];        
        end
    end
    
    if ind == 2
        genresults.time = machinesolver.sol.x';
        genresults.PhaseCurents = machinesolver.sol.y';
    else
        genresults.time = [ genresults.time;
                            machinesolver.sol.x(2:end)'];
                        
        genresults.PhaseCurents = [ genresults.PhaseCurents; 
                                    machinesolver.sol.y(:,2:end)' ];
    end
    
    machinesolver.acceptState ([xRpto(ind); vRpto(ind)]);
    
    ind = ind + 1;
    
end

[F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys, forces, F_AddedMass, accel);

toc

% close the socket connection
clear mb;

%% generator outputs

figure;
tmin =0;
tmax = 150;
plotinds =  time>=tmin & time<=tmax;
genplotinds = genresults.time>=tmin & genresults.time<=tmax;
[ax, hI, hvRpto] = plotyy ( genresults.time(genplotinds), genresults.PhaseCurents(genplotinds,:), ...
                            time(plotinds), vRpto(plotinds));
                        
ylabel(ax(1), 'Phase Current (A)');
ylabel(ax(2), 'Relative Velocity (ms^{-1})');
xlabel (ax(1), 'Time (s)');

legend ('I_A', 'I_b', 'I_c', 'vRpto');


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


tmax = min ( tmax, time(end));

if exist ('output', 'var')
    tmax = min (tmax, output.bodies(bodyind).time(end)); 
end

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
         ... 'ptoforce',  ...
         'vRpto' );

     

%%

mbout = mbdyn.postproc (outputfile_prefix, mbsys);


mbout.loadResultsFromFiles ( outputfile_prefix );

[hfig, hax] = mbout.plotNodeTrajectories ('OnlyNodes', [1,2]); 

%%
mbout.animate ( 'DrawMode', 'solid', ...
                'Light', true, ...
                'skip', 2, ...
                ...'AxLims', [-30, 30; -30, 30; -35, 35], ...
                'VideoFile', 'rm3_with_linear_gen_ds.avi', ...
                'VideoSpeed', 3, ...
                'OnlyNodes', [1,2])


%% Compare to origninal WEC-Sim

%
% Requires to have run the RM3 example first using original WEC-Sim

% doplot = false;
% 
% if doplot
%     for bodyind = 1:numel (output.bodies)
% 
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceTotal, ...
%               time, squeeze(forces(:,bodyind,:)));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceTotal,  ...
%               time, squeeze(F_Total(:,bodyind,:)));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%         title (sprintf ('forceTotal vs F\\_Total for body %d', bodyind));
% %         figure;
% %         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceExcitation,  output.bodies(bodyind).time, squeeze(F_ExcitRamp(:,bodyind,:))); 
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceAddedMass,  ...
%               output.bodies(bodyind).time, squeeze(F_AddedMass(:,bodyind,:)));
%         title (sprintf ('forceAddedMass vs F\\_addedmass for body %d', bodyind));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceAddedMass,  ...
%               output.bodies(bodyind).time, squeeze(F_AddedMassCorrected(:,bodyind,:)));
%         title (sprintf ('forceAddedMass vs F\\_AddedMassCorrected for body %d', bodyind));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceRadiationDamping,  ...
%               time, squeeze(F_RadiationDamping(:,bodyind,:)));
%         title (sprintf ('forceRadiationDamping vs F\\_RadiationDamping for body %d', bodyind));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
% %         figure;
% %         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceRestoring,  output.bodies(bodyind).time, squeeze(F_Restoring(:,bodyind,:)));
%         figure;
%         plot (output.bodies(bodyind).time, output.bodies(bodyind).position,  ...
%               time, squeeze(pos(:,bodyind,:)));
%         title (sprintf ('output.bodies(%d).position vs pos for body %d', bodyind, bodyind));
%         legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
%     end
% 
% end

% tmin = 0;
% tmax = 400;
% 
% bodyind = 1;
% 
% tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end)]);
% plotinds =  time>=tmin & time<=tmax;
% 
% figure;
% plotyy ( output.bodies(bodyind).time(plotinds), ...
%          [output.ptos(1).forceTotal(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%          time(plotinds), ...
%          [vRptoVec(plotinds)', [ output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3)] ] );
% title (sprintf ('output.ptos(%d).forceTotal vs FptoVec for body %d', bodyind, bodyind));
% legend ('1', '2', '3', '4', '5', '6', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
% figure;
% plotyy ( output.bodies(bodyind).time(plotinds), ...
%         [output.ptos(1).forceTotalWorld(plotinds,:) - output.ptos(1).forceConstraintWorld(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
%         time(plotinds), ...
%         [vRptoVec(:,plotinds)', output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3) ] );
% title (sprintf ('output.ptos(1).forceTotalWorld - output.ptos(1).forceConstraintWorld vs FptoVec for body %d', bodyind, bodyind));
% legend ('1', '2', '3', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');
% 
% figure;
% plot ( output.bodies(bodyind).time(plotinds), ...
%        output.ptos(1).forceInternalMechanics(plotinds,3),  ...
%        time(plotinds), ...
%        squeeze(ptoforce(plotinds)));
% title (sprintf ('output.ptos(%d).forceInternalMechanics(:,3) vs ptoforce for body %d', bodyind, bodyind));
% legend ('forceInternalMechanics(:,3)', 'ptoforce');
        
% 
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
% gf_forceTotal = nan * ones (numel (output.bodies), 11);
% gf_F_Total = gf_forceTotal;
% gf_forceExcitation = gf_forceTotal;
% gf_forceAddedMass = gf_forceTotal;
% gf_F_AddedMassCorrected = gf_forceTotal;
% gf_forceRadiationDamping = gf_forceTotal;
% gf_forceRestoring = gf_forceTotal;
% 
% for bodyind = 1:numel (output.bodies)
%     
%     % calculate some proper stats
%     gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, forceTotal(:,:,bodyind));
%     
%     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, F_Total(:,:,bodyind));
%     
%     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation, F_ExcitRamp(:,:,bodyind));
%     
%     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMass(:,:,bodyind));
%     
%     gf_F_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMassCorrected(:,:,bodyind));
%     
%     gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping, F_RadiationDamping(:,:,bodyind));
%     
%     gf_forceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,3:5),  F_Restoring(:,3:5,bodyind)); 
% 
%     rowheadings = [rowheadings, {
%                ['forceTotal_body_', int2str(bodyind)], ...
%                ['F_Total_body_', int2str(bodyind)], ...
%                ['forceExcitation_body_', int2str(bodyind)], ...
%                ['forceAddedMass_body_', int2str(bodyind)], ...
%                ['F_AddedMassCorrected_body_', int2str(bodyind)], ...
%                ['forceRadiationDamping_body_', int2str(bodyind)], ...
%                ['forceRestoring_body_', int2str(bodyind)] }];
% 
%     data = [ data;
%              gf_forceTotal(bodyind,:); 
%              gf_F_Total(bodyind,:);
%              gf_forceExcitation(bodyind,:); 
%              gf_forceAddedMass(bodyind,:);
%              gf_F_AddedMassCorrected(bodyind,:);
%              gf_forceRadiationDamping(bodyind,:); 
%              gf_forceRestoring(bodyind,:)];
% end
% 
% % display table
% colheadings = { 'mse', 'nmse', 'rmse', 'nrmse', 'mae', 'mare', 'r', 'r2', 'e', 'maxae', 'mxare' };
% 
% wid = 16;
% % fms = {'d','.4f','.5E'};
% fms = {};
% fileID = 1;
% 
% displaytable (data,colheadings,wid,fms,rowheadings,fileID);
% 
% 
