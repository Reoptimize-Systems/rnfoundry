% +test_wsim_OSWEC/run.m.m

if exist ('hmbsysfig', 'var')
    if ishghandle (hmbsysfig)
        close (hmbsysfig)
    end
end
clear waves simu hsys mbsys flap_hbody base_hbody mb hmbsysfig


%% Hydro Simulation Data
simu = wsim.simSettings (getmfilepath ('test_wsim_OSWEC.run'));  % Create the Simulation Variable
% simu.mode = 'normal';                       % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
% simu.explorer='on';                         % Turn SimMechanics Explorer (on/off)
% simu.startTime = 0;                         % Simulation Start Time [s]
simu.endTime=400;
simu.dt = 0.1; 							      %Simulation time-step [s] for a convolution function in the radiation force calculation 
simu.rampT = 100;
simu.CITime = 30;

%% Wave Information 

%% Regular Waves 
waves = wsim.waveSettings ('regular');
waves.H = 2.5;
waves.T = 8;

%% Hydrodynamic body system

%% Flap
flap_hbody = wsim.hydroBody ('oswec.h5');   % Initialize bodyClass for Flap
flap_hbody.mass = 127000;                         % User-Defined mass [kg]
flap_hbody.momOfInertia = [1.85e6 1.85e6 1.85e6]; % Moment of Inertia [kg-m^2]
flap_hbody.geometryFile = 'flap.stl';    % Geometry File

%% Base
base_hbody = wsim.hydroBody ('oswec.h5');   % Initialize bodyClass for Base
base_hbody.geometryFile = 'base.stl';    % Geometry File
base_hbody.mass = 'fixed';     

% make a hydrosys object for simulation
hsys = wsim.hydroSystem (waves, simu, [flap_hbody, base_hbody]);

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

[mbsys, pto_hinge] = test_wsim_OSWEC.make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options);
                     
% draw it
% mbsys.draw ('Mode', 'wireghost', 'Light', false);

[~, hmbsysfig] = mbsys.draw ( 'Mode', 'ghost', ...
             'Light', true, ...
             'AxLims', [-5, 5; -15, 15; -12, 5], ...
             'Joints', true, ...
             'StructuralNodes', true);

mbdpath = fullfile (simu.caseDir, 'OSWEC.mbd');

%% Set up PTO

k = 0;
% c = 6.2e6 / 0.2;

c = 1;

torquefcn = @(time, thetaRpto, omegaRpto) -k*thetaRpto -c*omegaRpto;

pto = wsim.rotaryPowerTakeOff ( pto_hinge, torquefcn );

%% Run the simulation


lssett = wsim.loggingSettings ();

lssett.positions = true;
lssett.velocities = true;
lssett.accelerations = true;
lssett.nodeForces = true;
lssett.nodeForcesUncorrected = true;
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
lssett.momentAddedMass = true;
lssett.nodeMoments = true;
lssett.nodeMomentsUncorrected = true;
lssett.momentHydro = true;
lssett.momentExcitation = true;
lssett.momentExcitationRamp = true;
lssett.momentExcitationLin = true;
lssett.momentExcitationNonLin = true;
lssett.momentRadiationDamping = true;
lssett.momentRestoring = true;
lssett.momentMorrison = true;
lssett.momentViscousDamping = true;
% lssett.momentAddedMassUncorrected = false;
lssett.momentAddedMass = true;

% create the wesim object
wsobj = wsim.wecSim ( hsys, mbsys, ...
                      'PTO', pto, ... % PTO(s) could also be added later using the 
                      'LoggingSettings', lssett );

% initialise the simulation
wsobj.prepare ();

% run it and get the output data
[ datalog, mbout ] = wsobj.run ('TimeExecution', true);

%% 

% figure;
datalog.plotVar ('PTO_1_InternalTorque');

datalog.plotVar ('PTO_1_RelativeAngularDisplacement');

datalog.plotVar ('PTO_1_RelativeAngularVelocity');

%% animate the sim
% wsobj.animate ( 'DrawMode', 'solid', ...
%                 'Light', true, ...
%                 'skip', 5, ...
%                 'AxLims', [-5, 5; -15, 15; -12, 5])

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

    figure;
    plotyy ( output.bodies(bodyind).time(plotinds), ...
             [output.ptos(1).forceTotal(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
             time(plotinds), ...
             [vRptoVec(plotinds)', [ output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3)] ] );
    title (sprintf ('output.ptos(%d).forceTotal vs FptoVec for body %d', bodyind, bodyind));
    legend ('1', '2', '3', '4', '5', '6', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');

    figure;
    plotyy ( output.bodies(bodyind).time(plotinds), ...
            [output.ptos(1).forceTotalWorld(plotinds,:) - output.ptos(1).forceConstraintWorld(plotinds,:), squeeze(FptoVec(:,bodyind,plotinds))'],  ...
            time(plotinds), ...
            [vRptoVec(:,plotinds)', output.bodies(1).velocity(plotinds,1:3) - output.bodies(2).velocity(plotinds,1:3) ] );
    title (sprintf ('output.ptos(1).forceTotalWorld - output.ptos(1).forceConstraintWorld vs FptoVec for body %d', bodyind, bodyind));
    legend ('1', '2', '3', 'FptoVec 1', 'FptoVec 2', 'FptoVec 3', 'myvR1', 'myvR2', 'myvR3', 'wsimvR1', 'wsimvR2', 'wsimvR3');

    figure;
    plot ( output.bodies(bodyind).time(plotinds), ...
           output.ptos(1).forceInternalMechanics(plotinds,3),  ...
           time(plotinds), ...
           squeeze(ptoforce(plotinds)));
    title (sprintf ('output.ptos(%d).forceInternalMechanics(:,3) vs ptoforce for body %d', bodyind, bodyind));
    legend ('forceInternalMechanics(:,3)', 'ptoforce');

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
