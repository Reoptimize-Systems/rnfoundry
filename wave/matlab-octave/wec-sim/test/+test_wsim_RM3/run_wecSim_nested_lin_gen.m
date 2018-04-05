% +test_wsim_RM3/run.m.m

clear waves simu hsys mbsys float_hbody spar_hbody hydro_mbnodes hydro_mbbodies hydro_mbelements lssett wsobj

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

% generate the nodes and elements for simulation of the hydrodynamic system
% in MBDyn. One node and one body element are created for each hydrodynamic
% body interacting with the waves.
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

% mbsys.draw ( 'Mode', 'solid', ...
%              'Light', true, ...
%              'AxLims', [-30, 30; -30, 30; -35, 35], ...
%              'Joints', false, ...
%              'StructuralNodes', false)

mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

%% Set up Power Take-Off (PTO)

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

% create a power take-off object attached to the two hydro nodes, with the
% force being based on the relative velocity and displacement along axis 3
% of the first (reference) node (in it's coordinate frame).
pto = wsim.linearPMMachinePTO ( hydro_mbnodes{2}, hydro_mbnodes{1}, 3, design, simoptions );

%% Run the simulation

lssett = wsim.loggingSettings ();

lssett.positions = true;
lssett.velocities = true;
lssett.accelerations = true;
lssett.nodeForcesAndMoments = true;
lssett.nodeForcesAndMomentsUncorrected = true;
lssett.forceHydro = true;
lssett.forceExcitation = true;
lssett.forceExcitationRamp = true;
lssett.forceExcitationLin = true;
lssett.forceExcitationNonLin = true;
lssett.forceRadiationDamping = true;
lssett.forceRestoring = true;
lssett.forceMorrison = true;
lssett.forceViscousDamping = true;
% lssett.ForceAddedMassUncorrected = false;
lssett.forceAddedMass = true;
        
% create the wesim object
wsobj = wsim.wecSim ( hsys, mbsys, ...
                      'PTO', pto, ... % PTO(s) could also be added later using the 
                      'LoggingSettings', lssett );

% initialise the simulation
wsobj.prepare ();

% run it and get the output data
datalog = wsobj.run ('TimeExecution', true);

%%
% figure;
% tmin = 0;
% tmax = 400;
% plotinds =  time>=tmin & time<=tmax;
% plotyy (time(plotinds), [ squeeze(forces(1:3,1,plotinds))',  ptoforce(plotinds)'], time(plotinds), vRpto(plotinds));
% legend ('fx', 'fy', 'fz', 'ptoforce', 'vRpto');

%%

if ~exist ('output', 'var')
    warning ('Not comparing output to original WEC-Sim as ''output'' variable is not in the workspace.')
else

    figure;
    tmin = 0;
    tmax = 400;

    bodyind = 2;
    
    time = datalog.data.Time;
    
    tmax = min ( [tmax, time(end), output.bodies(bodyind).time(end) ]);
    plotinds =  time>=tmin & time<=tmax;

    plotyy (time(plotinds), ...
          [ squeeze(datalog.data.ForceExcitationRamp(:,bodyind,plotinds))', ...
            ... squeeze(F_ViscousDamping(3,bodyind,plotinds)), ...
            squeeze(datalog.data.ForceAddedMass(:,bodyind,plotinds))', ...
            squeeze(datalog.data.ForceRestoring(:,bodyind,plotinds))', ...
            squeeze(datalog.data.ForceRadiationDamping(:,bodyind,plotinds))', ...
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
            time(plotinds), datalog.data.PTO_1_RelativeVelocity(plotinds) );

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
             'PTO\_1\_RelativeVelocity' );


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
                  time, squeeze(datalog.data.NodeForcesAndMomentsUncorrected(:,bodyind,:)));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceTotal,  ...
                  time, squeeze(datalog.data.NodeForcesAndMoments(:,bodyind,:)));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            title (sprintf ('forceTotal vs F\\_Total for body %d', bodyind));
    %         figure;
    %         plot (output.bodies(bodyind).time, output.bodies(bodyind).forceExcitation,  output.bodies(bodyind).time, squeeze(F_ExcitRamp(:,bodyind,:))); 
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).ForceAddedMass,  ...
                  output.bodies(bodyind).time, squeeze(datalog.data.ForceAddedMass(:,bodyind,:)));
            title (sprintf ('ForceAddedMass vs F\\_addedmass for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).ForceAddedMass,  ...
                  output.bodies(bodyind).time, squeeze(datalog.data.ForceAddedMass(:,bodyind,:)));
            title (sprintf ('ForceAddedMass vs F\\_AddedMassCorrected for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).forceRadiationDamping,  ...
                  time, squeeze(datalog.data.ForceRadiationDamping(:,bodyind,:)));
            title (sprintf ('forceRadiationDamping vs F\\_RadiationDamping for body %d', bodyind));
            legend ('1', '2', '3', '4', '5', '6', '1', '2', '3', '4', '5', '6');
    %         figure;
    %         plot (output.bodies(bodyind).time, output.bodies(bodyind).ForceRestoring,  output.bodies(bodyind).time, squeeze(data.data.ForceRestoring(:,bodyind,:)));
            figure;
            plot (output.bodies(bodyind).time, output.bodies(bodyind).position,  ...
                  time, squeeze(datalog.data.Positions(:,bodyind,:)));
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
    % % plot (output.bodies(1).time, [output.bodies(1).ForceAddedMass(:,fcomp), output.bodies(2).ForceAddedMass(:,fcomp), squeeze(data.data.ForceAddedMass(:,fcomp,1:2))]);
    % % title (sprintf ('comp %d ForceAddedMass vs F\\_AddedMassCorrected for both bodies', fcomp));
    % % legend (sprintf ('body 1 ForceAddedMass post-processed comp %d', fcomp), ...
    % %     sprintf ('body 2 ForceAddedMass post-processed comp %d', fcomp),...
    % %     sprintf ('body 1 hydrobody data.data.ForceAddedMass comp %d', fcomp), ...
    % %     sprintf ('body 2 hydrobody data.data.ForceAddedMass comp %d', fcomp))
    % 
    % data = [];
    % rowheadings = {};
    % 
    % stats = {'3', '5', '8'}
    % gf_forceTotal = nan * ones (numel (output.bodies), numel(stats));
    % gf_F_Total = gf_forceTotal;
    % gf_M_Total = gf_F_Total;
    % gf_forceExcitation = gf_forceTotal;
    % gf_ForceAddedMass = gf_forceTotal;
    % gf_data.data.ForceAddedMass = gf_forceTotal;
    % gf_forceRadiationDamping = gf_forceTotal;
    % gf_ForceRestoring = gf_forceTotal;
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
    % %     gf_ForceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).ForceAddedMass, F_AddedMass(:,:,bodyind));
    %     
    %     gf_data.data.ForceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).ForceAddedMass, squeeze(data.data.ForceAddedMass(:,bodyind,:))', stats);
    %     
    %     gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping, squeeze(data.data.ForceRadiationDamping(:,bodyind,:))', stats);
    %     
    %     gf_ForceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).ForceRestoring,  squeeze(data.data.ForceRestoring(:,bodyind,:))', stats); 
    % 
    %     rowheadings = [rowheadings, {
    %                ...['forceTotal_body_', int2str(bodyind)], ...
    %                sprintf('Body %d Total Force', bodyind), ...
    %                sprintf('Body %d Total Moments', bodyind), ...
    %                sprintf('Body %d Excitation Force', bodyind), ...
    %                ...['ForceAddedMass_body_', int2str(bodyind)], ...
    %                sprintf('Body %d Added Mass Force', bodyind), ...
    %                sprintf('Body %d Radiation and Damping Force', bodyind), ...
    %                sprintf('Body %d Hydrostatic Restoring Force', bodyind) }];
    % 
    %     data = [ data;
    %              ...gf_forceTotal(bodyind,:); 
    %              gf_F_Total(bodyind,:);
    %              gf_M_Total(bodyind,:);
    %              gf_forceExcitation(bodyind,:); 
    %              ...gf_ForceAddedMass(bodyind,:);
    %              gf_data.data.ForceAddedMass(bodyind,:);
    %              gf_forceRadiationDamping(bodyind,:); 
    %              gf_ForceRestoring(bodyind,:)];
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

    stattabledata = [];
    rowheadings = {};

    stats = {'3', '5', '8'}
    gf_forceTotal = nan * ones (numel (output.bodies), numel(stats));
    gf_momentTotal = nan * ones (numel (output.bodies), numel(stats));
    gf_F_Total = gf_forceTotal;
    gf_M_Total = gf_F_Total;
    gf_ForceExcitation = gf_forceTotal;
    gf_momentExcitation = gf_forceTotal;
    gf_ForceAddedMass = gf_forceTotal;
    gf_ForceAddedMassCorrected = gf_forceTotal;
    gf_M_AddedMassCorrected = gf_forceTotal;
    gf_ForceRadiationDamping = gf_forceTotal;
    gf_momentRadiationDamping= gf_forceTotal;
    gf_ForceRestoring = gf_forceTotal;
    gf_momentRestoring = gf_forceTotal;

    gf_pos = gf_forceTotal;
    gf_theta = gf_forceTotal;
    gf_vel = gf_forceTotal;
    gf_omega = gf_forceTotal;
    gf_accel = gf_forceTotal;
    gf_omegaaccel = gf_forceTotal;

    for bodyind = 1:numel (output.bodies)

        % calculate some proper stats
        gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,1:3), squeeze(datalog.data.NodeForcesAndMoments(1:3,bodyind,:))', stats);
        gf_momentTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotalOrig(:,4:6), squeeze(datalog.data.NodeForcesAndMoments(4:6,bodyind,:))', stats);

    %     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, squeeze(F_Total(:,bodyind,:))', stats);
        gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,1:3), squeeze(datalog.data.NodeForcesAndMomentsUncorrected(1:3,bodyind,:))', stats);
        gf_M_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal(:,4:6), squeeze(datalog.data.NodeForcesAndMomentsUncorrected(4:6,bodyind,:))', stats);

    %     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
    %     gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);

        gf_ForceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(datalog.data.ForceExcitation(1:3,bodyind,:))', stats);
        gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(datalog.data.ForceExcitation(4:6,bodyind,:))', stats);

    %     gf_ForceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).ForceAddedMass, F_AddedMass(:,:,bodyind));

        gf_ForceAddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,1:3), squeeze(datalog.data.ForceAddedMass (1:3,bodyind,:))', stats);
        gf_M_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass(:,4:6), squeeze(datalog.data.ForceAddedMass(4:6,bodyind,:))', stats);

        gf_ForceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,1:3), squeeze(datalog.data.ForceRadiationDamping(1:3,bodyind,:))', stats);
        gf_momentRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping(:,4:6), squeeze(datalog.data.ForceRadiationDamping(4:6,bodyind,:))', stats);

        gf_ForceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,1:3),  squeeze(datalog.data.ForceRestoring(1:3,bodyind,:))', stats); 
        gf_momentRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,4:6),  squeeze(datalog.data.ForceRestoring(4:6,bodyind,:))', stats); 

        gf_pos(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,1:3),  squeeze(datalog.data.Positions(1:3,bodyind,:))', stats); 
        gf_theta(bodyind,:) = gfit2 (output.bodies(bodyind).position(:,4:6),  squeeze(datalog.data.Positions(4:6,bodyind,:))', stats); 

        gf_vel(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,1:3),  squeeze(datalog.data.Velocities(1:3,bodyind,:))', stats); 
        gf_omega(bodyind,:) = gfit2 (output.bodies(bodyind).velocity(:,4:6),  squeeze(datalog.data.Velocities(4:6,bodyind,:))', stats); 

        gf_accel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,1:3),  squeeze(datalog.data.Accelerations(1:3,bodyind,:))', stats); 
        gf_omegaaccel(bodyind,:) = gfit2 (output.bodies(bodyind).acceleration(:,4:6),  squeeze(datalog.data.Accelerations(4:6,bodyind,:))', stats); 

        rowheadings = [rowheadings, {
                   ...['forceTotal_body_', int2str(bodyind)], ...
                   sprintf('Body %d Total Force Uncorrected', bodyind), ...
                   sprintf('Body %d Total Uncorrected', bodyind), ...
                   sprintf('Body %d Total Force', bodyind), ...
                   sprintf('Body %d Total Moments', bodyind), ...
                   sprintf('Body %d Excitation Force', bodyind), ...
                   sprintf('Body %d Excitation Moment', bodyind), ...
                   ...['ForceAddedMass_body_', int2str(bodyind)], ...
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


        stattabledata = [ stattabledata;
                 ...gf_forceTotal(bodyind,:); 
                 gf_forceTotal(bodyind,:)
                 gf_momentTotal(bodyind,:)
                 gf_F_Total(bodyind,:);
                 gf_M_Total(bodyind,:);
                 gf_ForceExcitation(bodyind,:);
                 gf_momentExcitation(bodyind,:);
                 ...gf_ForceAddedMass(bodyind,:);
                 gf_ForceAddedMassCorrected(bodyind,:);
                 gf_M_AddedMassCorrected(bodyind,:);
                 gf_ForceRadiationDamping(bodyind,:);
                 gf_momentRadiationDamping(bodyind,:); 
                 gf_ForceRestoring(bodyind,:);
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

    displaytable (stattabledata,colheadings,wid,fms,rowheadings,fileID);

    colheadings = [{'Force Description'}, colheadings];
    fms = {'.2g','.2g','.2f'};

    latextable (stattabledata, ...
                'ColumnHeadings', colheadings, ...
                'NumberWidth', wid, ...
                'FormatSpec', fms, ...
                'RowHeadings', rowheadings, ...
                'BookTabs', true)
end
