% test_hydrosys.m

clear hsys_mex float_hbody_mex spar_hbody_mex

test_wsim_RM3.run_wecsim_in_octave;

%% Hydrodynamic body system

% hydro data fiels assumed to be in assumed to be in folder 
% <case_directory>/hydroData

% Float, input is the name of the h5 file containing the hydrodynamic data,
% expected to be found in the <case_directory>/hydroData folder
float_hbody_mex = wsim.hydroBodyMexForces('float.mat');
% Body Mass. The 'equilibrium' Option Sets it to the Displaced Water Weight.
float_hbody_mex.mass = 'equilibrium';
% Moment of Inertia [kg*m^2]
float_hbody_mex.momOfInertia = [20907301, 21306090.66, 37085481.11];  
% Geomtry File Name (expected to be in <case_directory>/geometry)
float_hbody_mex.geometryFile = 'float.stl'; 

% Spar/Plate
spar_hbody_mex = wsim.hydroBodyMexForces('spar.mat'); 
spar_hbody_mex.mass = 'equilibrium';                   
spar_hbody_mex.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody_mex.geometryFile = 'plate.stl';

bodies = float_hbody_mex;
bodies(2) = spar_hbody_mex;

% make a hydrosys object for simulation
hsys_mex = wsim.hydroSystem (waves, simu, bodies);

% set up hydro system bodies (loading hydro data etc.)
hsys_mex.initialiseHydrobodies ();
% get it ready to do a transient simulation
hsys_mex.timeDomainSimSetup ();


%% Compare to origninal hydroBody

F_ExcitLin = [];
F_ViscousDamping = [];
F_addedmass = [];
F_RadiationDamping = [];
F_ExcitLinNonLin = [];
F_MorrisonElement = [];
F_Excit = [];
F_ExcitRamp = [];
F_RealAddedMass = [];
accel = [];

tic
for ind = 1:numel (datalog.data.Time)
    
    
    pos = [datalog.data.Positions(:,:,ind);  datalog.data.AngularPositions(:,:,ind)];
    vel = [datalog.data.Velocities(:,:,ind);  datalog.data.AngularVelocities(:,:,ind)];
    accel = [datalog.data.Accelerations(:,:,ind);  datalog.data.AngularAccelerations(:,:,ind)];

    [forces, out] = hsys_mex.hydroForces ( datalog.data.Time(ind), ...
                                       pos, ...
                                       vel, ...
                                       accel );


    forceTotal(1:6,1:2,ind) = forces;
    F_ExcitLin(1:6,1:2,ind) = out.F_ExcitLin;
    F_ViscousDamping(1:6,1:2,ind) = out.F_ViscousDamping;
    F_addedmass(1:6,1:2,ind) = out.F_AddedMass;
    F_Restoring(1:6,1:2,ind) = out.F_Restoring;
    F_RadiationDamping(1:6,1:2,ind) = out.F_RadiationDamping;
%     F_ExcitLinNonLin(1:6,1:2,ind) = out.F_ExcitLinNonLin;
    F_MorrisonElement(1:6,1:2,ind) = out.F_MorrisonElement;
    F_Excit(1:6,1:2,ind) = out.F_Excit;
    F_ExcitRamp(1:6,1:2,ind) = out.F_ExcitRamp;

    hsys_mex.advanceStep ( datalog.data.Time(ind), ...
                       vel, ...
                       accel );
   
end
toc

% for bodyind = 1:numel (output.bodies)
%     accel = cat (3, accel, output.bodies(bodyind).acceleration);
% end

% hsys.restoreMassMatrices ();
[F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys_mex, forceTotal, F_addedmass,  [datalog.data.Accelerations;  datalog.data.AngularAccelerations]);



%% Compare to WEC-Sim
%
% Requires to have run the RM3 example first using original WEC-Sim

% doplot = false;
% 
% if doplot
%     for ind = 1:numel (output.bodies)
% 
% %         figure;
% %         plot (output.bodies(ind).time, output.bodies(ind).forceTotal,  output.bodies(ind).time, forceTotal(:,:,ind));
% %         figure;
% %         plot (output.bodies(ind).time, output.bodies(ind).forceTotal,  output.bodies(ind).time, F_Total(:,:,ind)); 
% %         figure;
% %         plot (output.bodies(ind).time, output.bodies(ind).forceExcitation,  output.bodies(ind).time, F_ExcitRamp(:,:,ind)); 
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceAddedMass,  output.bodies(ind).time, F_addedmass(:,:,ind));
%         title (sprintf ('forceAddedMass vs F\\_addedmass for body %d', ind));
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceAddedMass,  output.bodies(ind).time, F_AddedMassCorrected(:,:,ind));
%         title (sprintf ('forceAddedMass vs F\\_AddedMassCorrected for body %d', ind));
% %         figure;
% %         plot (output.bodies(ind).time, output.bodies(ind).forceRadiationDamping,  output.bodies(ind).time, F_RadiationDamping(:,:,ind)); 
% %         figure;
% %         plot (output.bodies(ind).time, output.bodies(ind).forceRestoring,  output.bodies(ind).time, F_Restoring(:,:,ind));
% 
%     end
% 
% end
% 
% % figure;
% % fcomp = 3;
% % plot (output.bodies(1).time, [body1_F_AddedMass_Simulink.signals.values(:,fcomp), body2_F_AddedMass_Simulink.signals.values(:,fcomp), squeeze(F_addedmass(:,fcomp,1:2))]);
% % title (sprintf ('comp %d of body2\\_F\\_AddedMass\\_Simulink vs hydrobody F\\_addedmass for body 2', fcomp));
% % legend (sprintf ('body 1 simulink comp %d', fcomp), ...
% %     sprintf ('body 2 simulink comp %d', fcomp),...
% %     sprintf ('body 1 hydrobody F_addedmass comp %d', fcomp), ...
% %     sprintf ('body 2 hydrobody F_addedmass comp %d', fcomp))
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
%     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_addedmass(:,:,bodyind));
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


data = [];
rowheadings = {};

stats = {'3', '5', '8'}
gf_forceTotal = nan * ones (hsys_mex.nHydroBodies, numel(stats));
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

for bodyind = 1:numel (hsys_mex.nHydroBodies)
    
    % calculate some proper stats
    ...gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, forceTotal(:,:,bodyind));
    
%     gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, squeeze(F_Total(:,bodyind,:))', stats);
    gf_F_Total(bodyind,:) = gfit2 (forceTotal(1:3,bodyind,:), datalog.data.ForceHydro(:,bodyind,:), stats);
    gf_M_Total(bodyind,:) = gfit2 (forceTotal(4:6,bodyind,:), datalog.data.MomentHydro(:,bodyind,:), stats);
    
%     gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,1:3), squeeze(F_ExcitRamp(1:3,bodyind,:))', stats);
%     gf_momentExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation(:,4:6), squeeze(F_ExcitRamp(4:6,bodyind,:))', stats);
    
    gf_forceExcitation(bodyind,:) = gfit2 (F_ExcitRamp(1:3,bodyind,:), datalog.data.ForceExcitationRamp(:,bodyind,:), stats);
    gf_momentExcitation(bodyind,:) = gfit2 (F_ExcitRamp(4:6,bodyind,:), datalog.data.MomentExcitationRamp(:,bodyind,:), stats);
    
%     gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_addedmass(:,:,bodyind));
    
    gf_F_AddedMassCorrected(bodyind,:) = gfit2 (F_AddedMassCorrected(1:3,bodyind,:), datalog.data.ForceAddedMass(:,bodyind,:), stats);
    gf_M_AddedMassCorrected(bodyind,:) = gfit2 (F_AddedMassCorrected(4:6,bodyind,:), datalog.data.MomentAddedMass(:,bodyind,:), stats);
    
    gf_forceRadiationDamping(bodyind,:) = gfit2 (F_RadiationDamping(1:3,bodyind,:), datalog.data.ForceRadiationDamping(:,bodyind,:), stats);
    gf_momentRadiationDamping(bodyind,:) = gfit2 (F_RadiationDamping(4:6,bodyind,:), datalog.data.MomentRadiationDamping(:,bodyind,:), stats);
    
    gf_forceRestoring(bodyind,:) = gfit2 (F_Restoring(1:3,bodyind,:), datalog.data.ForceRestoring(:,bodyind,:), stats); 
    gf_momentRestoring(bodyind,:) = gfit2 (F_Restoring(4:6,bodyind,:), datalog.data.MomentRestoring(:,bodyind,:), stats); 

    rowheadings = [rowheadings, {
               ...['forceTotal_body_', int2str(bodyind)], ...
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
               sprintf('Body %d Hydrostatic Restoring Moment', bodyind) }];
           

    data = [ data;
             ...gf_forceTotal(bodyind,:); 
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
             gf_momentRestoring(bodyind,:) ];
end

% display table
colheadings = { 'RMSE', 'MAE', 'R2' };

wid = 16;
% fms = {'d','.4f','.5E'};
fms = {};
fileID = 1;

displaytable ( data, ...
               'ColHeadings', colheadings, ...
               'ColWidth', wid, ...
               'Format', fms, ...
               'RowHeadings', rowheadings);

colheadings = [{'Force Description'}, colheadings];
fms = {'.2g','.2g','.2f'};

latextable (data, ...
            'ColumnHeadings', colheadings, ...
            'NumberWidth', wid, ...
            'FormatSpec', fms, ...
            'RowHeadings', rowheadings, ...
            'BookTabs', true)

        
%% Compare speed of original and mex
% hydro data fiels assumed to be in assumed to be in folder 
% <case_directory>/hydroData

% Float, input is the name of the h5 file containing the hydrodynamic data,
% expected to be found in the <case_directory>/hydroData folder
float_hbody_orig = wsim.hydroBody('float.mat');
% Body Mass. The 'equilibrium' Option Sets it to the Displaced Water Weight.
float_hbody_orig.mass = 'equilibrium';
% Moment of Inertia [kg*m^2]
float_hbody_orig.momOfInertia = [20907301, 21306090.66, 37085481.11];  
% Geomtry File Name (expected to be in <case_directory>/geometry)
float_hbody_orig.geometryFile = 'float.stl'; 

% Spar/Plate
spar_hbody_orig = wsim.hydroBody('spar.mat'); 
spar_hbody_orig.mass = 'equilibrium';                   
spar_hbody_orig.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody_orig.geometryFile = 'plate.stl';

bodies = float_hbody_orig;
bodies(2) = spar_hbody_orig;

% make a hydrosys object for simulation
hsys_orig = wsim.hydroSystem (waves, simu, bodies);

% set up hydro system bodies (loading hydro data etc.)
hsys_orig.initialiseHydrobodies ();
% get it ready to do a transient simulation
hsys_orig.timeDomainSimSetup ();

% Float, input is the name of the h5 file containing the hydrodynamic data,
% expected to be found in the <case_directory>/hydroData folder
float_hbody_mex = wsim.hydroBodyMexForces('float.mat');
% Body Mass. The 'equilibrium' Option Sets it to the Displaced Water Weight.
float_hbody_mex.mass = 'equilibrium';
% Moment of Inertia [kg*m^2]
float_hbody_mex.momOfInertia = [20907301, 21306090.66, 37085481.11];  
% Geomtry File Name (expected to be in <case_directory>/geometry)
float_hbody_mex.geometryFile = 'float.stl'; 

% Spar/Plate
spar_hbody_mex = wsim.hydroBodyMexForces('spar.mat'); 
spar_hbody_mex.mass = 'equilibrium';                   
spar_hbody_mex.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody_mex.geometryFile = 'plate.stl';

bodies = float_hbody_mex;
bodies(2) = spar_hbody_mex;

% make a hydrosys object for simulation
hsys_mex = wsim.hydroSystem (waves, simu, bodies);

% set up hydro system bodies (loading hydro data etc.)
hsys_mex.initialiseHydrobodies ();
% get it ready to do a transient simulation
hsys_mex.timeDomainSimSetup ();

fprintf (1, 'Running mex\n');
tic
for ind = 1:numel (datalog.data.Time)
    
    
    pos = [datalog.data.Positions(:,:,ind);  datalog.data.AngularPositions(:,:,ind)];
    vel = [datalog.data.Velocities(:,:,ind);  datalog.data.AngularVelocities(:,:,ind)];
    accel = [datalog.data.Accelerations(:,:,ind);  datalog.data.AngularAccelerations(:,:,ind)];

    [forces, out] = hsys_mex.hydroForces ( datalog.data.Time(ind), ...
                                       pos, ...
                                       vel, ...
                                       accel );

    hsys_mex.advanceStep ( datalog.data.Time(ind), ...
                       vel, ...
                       accel );
   
end
toc

fprintf (1, 'Running orig\n');
tic
for ind = 1:numel (datalog.data.Time)
    
    
    pos = [datalog.data.Positions(:,:,ind);  datalog.data.AngularPositions(:,:,ind)];
    vel = [datalog.data.Velocities(:,:,ind);  datalog.data.AngularVelocities(:,:,ind)];
    accel = [datalog.data.Accelerations(:,:,ind);  datalog.data.AngularAccelerations(:,:,ind)];

    [forces, out] = hsys_orig.hydroForces ( datalog.data.Time(ind), ...
                                       pos, ...
                                       vel, ...
                                       accel );


    hsys_orig.advanceStep ( datalog.data.Time(ind), ...
                       vel, ...
                       accel );
   
end
toc
