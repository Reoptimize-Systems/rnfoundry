% test_hydrosys.m

clear waves simu hydrobodies hsys

%% Simulation Data
simu = simulationClass('/home/rcrozier/src/WEC-Sim-git/tutorials/RM3');    %Create the Simulation Variable
% simu.mode = 'normal';                 %Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
% simu.explorer='on';                   %Turn SimMechanics Explorer (on/off)
% simu.startTime = 0;                   %Simulation Start Time [s]
simu.endTime=400;                       %Simulation End bdcloseTime [s]
simu.solver = 'ode4';                   %simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.1; 							%Simulation time-step [s]
simu.rampT = 100;                       %Wave Ramp Time Length [s]

simu.b2b = 1;

%% Wave Information 
%% noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       %Create the Wave Variable and Specify Type      

%% Regular Waves  
waves = waveClass('regularCIC');        %Create the Wave Variable and Specify Type                               
waves.H = 2.5;                          %Wave Height [m]
waves.T = 8;                            %Wave Period [s]

%% Irregular Waves using PM Spectrum with Convolution Integral Calculation
% waves = waveClass('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'PM';

%% Irregular Waves using BS Spectrum with State Space Calculation
% waves = waveClass('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'BS';
% simu.ssCalc = 1;						%Control option to use state space model 

%% Irregular Waves using User-Defined Spectrum
% waves = waveClass('irregularImport');  %Create the Wave Variable and Specify Type
% waves.spectrumDataFile = 'ndbcBuoyData.txt';  %Name of User-Defined Spectrum File [2,:] = [omega, Sf]

%% User-Defined Time-Series
% waves = waveClass('userDefined');     %Create the Wave Variable and Specify Type
% waves.etaDataFile = 'umpqua46229_6_2008.mat'; % Name of User-Defined Time-Series File [:,2] = [time, wave_elev]

%% Body Data
%% Float
hydrobodies(1) = wsim.hydrobody('hydroData/rm3.h5', simu.caseDir);      
    %Create the wsim.hydrobody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
hydrobodies(1).mass = 'equilibrium';                   
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
hydrobodies(1).momOfInertia = [20907301 21306090.66 37085481.11];  %Moment of Inertia [kg*m^2]     
hydrobodies(1).geometryFile = 'geometry/float.stl';    %Location of Geomtry File

%% Spar/Plate
hydrobodies(2) = wsim.hydrobody('hydroData/rm3.h5', simu.caseDir); 
hydrobodies(2).mass = 'equilibrium';                   
hydrobodies(2).momOfInertia = [94419614.57 94407091.24 28542224.82];
hydrobodies(2).geometryFile = 'geometry/plate.stl'; 

% make a hydrosys object for simulation
hsys = wsim.hydrosys (waves, simu, hydrobodies);

% set up ode simulation (runs 
hsys.initialiseHydrobodies ();
hsys.odeSimSetup ();
% hsys.restoreMassMatrices (); % should be called *after* simulation

%% Compare to origninal WEC-Sim

F_ExcitLin = [];
F_ViscousDamping = [];
F_addedmass = [];
F_RadiationDamping = [];
F_ExcitLinNonLin = [];
F_MorrisonElement = [];
F_Excit = [];
F_ExcitRamp = [];
F_RealAddedMass = [];

for ind = 1:numel (output.bodies(1).time)

    [forces, out] = hsys.hydroForces ( output.bodies(1).time(ind), ...
                                       [output.bodies(1).position(ind,:)', output.bodies(2).position(ind,:)'], ...
                                       [output.bodies(1).velocity(ind,:)', output.bodies(2).velocity(ind,:)'], ...
                                       [output.bodies(1).acceleration(ind,:)', output.bodies(2).acceleration(ind,:)'], ...
                                       output.wave.elevation(ind) );

    for bodyind = 1:numel (output.bodies)
        
       forceTotal(ind,1:6,bodyind) = forces(:,bodyind)';
       F_ExcitLin(ind,1:6,bodyind) = out(bodyind).F_ExcitLin';
       F_ViscousDamping(ind,1:6,bodyind) = out(bodyind).F_ViscousDamping';
       F_addedmass(ind,1:6,bodyind) = out(bodyind).F_addedmass';
       F_Restoring(ind,1:6,bodyind) = out(bodyind).F_Restoring';
       F_RadiationDamping(ind,1:6,bodyind) = out(bodyind).F_RadiationDamping';
       F_ExcitLinNonLin(ind,1:6,bodyind) = out(bodyind).F_ExcitLinNonLin';
       F_MorrisonElement(ind,1:6,bodyind) = out(bodyind).F_MorrisonElement';
       F_Excit(ind,1:6,bodyind) = out(bodyind).F_Excit';
       F_ExcitRamp(ind,1:6,bodyind) = out(bodyind).F_ExcitRamp';

    end
   
end

accel = [];
for bodyind = 1:numel (output.bodies)
    accel = cat (3, accel, output.bodies(bodyind).acceleration);
end

% hsys.restoreMassMatrices ();
[F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys, forceTotal, F_addedmass, accel);

%% Compare to WEC-Sim
%
% Requires to have run the RM3 example first using original WEC-Sim

doplot = false;

if doplot
    for ind = 1:numel (output.bodies)

%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceTotal,  output.bodies(ind).time, forceTotal(:,:,ind));
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceTotal,  output.bodies(ind).time, F_Total(:,:,ind)); 
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceExcitation,  output.bodies(ind).time, F_ExcitRamp(:,:,ind)); 
        figure;
        plot (output.bodies(ind).time, output.bodies(ind).forceAddedMass,  output.bodies(ind).time, F_addedmass(:,:,ind));
        title (sprintf ('forceAddedMass vs F\\_addedmass for body %d', ind));
        figure;
        plot (output.bodies(ind).time, output.bodies(ind).forceAddedMass,  output.bodies(ind).time, F_AddedMassCorrected(:,:,ind));
        title (sprintf ('forceAddedMass vs F\\_AddedMassCorrected for body %d', ind));
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceRadiationDamping,  output.bodies(ind).time, F_RadiationDamping(:,:,ind)); 
%         figure;
%         plot (output.bodies(ind).time, output.bodies(ind).forceRestoring,  output.bodies(ind).time, F_Restoring(:,:,ind));

    end

end

% figure;
% fcomp = 3;
% plot (output.bodies(1).time, [body1_F_AddedMass_Simulink.signals.values(:,fcomp), body2_F_AddedMass_Simulink.signals.values(:,fcomp), squeeze(F_addedmass(:,fcomp,1:2))]);
% title (sprintf ('comp %d of body2\\_F\\_AddedMass\\_Simulink vs hydrobody F\\_addedmass for body 2', fcomp));
% legend (sprintf ('body 1 simulink comp %d', fcomp), ...
%     sprintf ('body 2 simulink comp %d', fcomp),...
%     sprintf ('body 1 hydrobody F_addedmass comp %d', fcomp), ...
%     sprintf ('body 2 hydrobody F_addedmass comp %d', fcomp))
% 
% figure;
% fcomp = 3;
% plot (output.bodies(1).time, [output.bodies(1).forceAddedMass(:,fcomp), output.bodies(2).forceAddedMass(:,fcomp), squeeze(F_AddedMassCorrected(:,fcomp,1:2))]);
% title (sprintf ('comp %d forceAddedMass vs F\\_AddedMassCorrected for both bodies', fcomp));
% legend (sprintf ('body 1 forceAddedMass post-processed comp %d', fcomp), ...
%     sprintf ('body 2 forceAddedMass post-processed comp %d', fcomp),...
%     sprintf ('body 1 hydrobody F_AddedMassCorrected comp %d', fcomp), ...
%     sprintf ('body 2 hydrobody F_AddedMassCorrected comp %d', fcomp))

data = [];
rowheadings = {};

gf_forceTotal = nan * ones (numel (output.bodies), 11);
gf_F_Total = gf_forceTotal;
gf_forceExcitation = gf_forceTotal;
gf_forceAddedMass = gf_forceTotal;
gf_F_AddedMassCorrected = gf_forceTotal;
gf_forceRadiationDamping = gf_forceTotal;
gf_forceRestoring = gf_forceTotal;

for bodyind = 1:numel (output.bodies)
    
    % calculate some proper stats
    gf_forceTotal(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, forceTotal(:,:,bodyind));
    
    gf_F_Total(bodyind,:) = gfit2 (output.bodies(bodyind).forceTotal, F_Total(:,:,bodyind));
    
    gf_forceExcitation(bodyind,:) = gfit2 (output.bodies(bodyind).forceExcitation, F_ExcitRamp(:,:,bodyind));
    
    gf_forceAddedMass(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_addedmass(:,:,bodyind));
    
    gf_F_AddedMassCorrected(bodyind,:) = gfit2 (output.bodies(bodyind).forceAddedMass, F_AddedMassCorrected(:,:,bodyind));
    
    gf_forceRadiationDamping(bodyind,:) = gfit2 (output.bodies(bodyind).forceRadiationDamping, F_RadiationDamping(:,:,bodyind));
    
    gf_forceRestoring(bodyind,:) = gfit2 (output.bodies(bodyind).forceRestoring(:,3:5),  F_Restoring(:,3:5,bodyind)); 

    rowheadings = [rowheadings, {
               ['forceTotal_body_', int2str(bodyind)], ...
               ['F_Total_body_', int2str(bodyind)], ...
               ['forceExcitation_body_', int2str(bodyind)], ...
               ['forceAddedMass_body_', int2str(bodyind)], ...
               ['F_AddedMassCorrected_body_', int2str(bodyind)], ...
               ['forceRadiationDamping_body_', int2str(bodyind)], ...
               ['forceRestoring_body_', int2str(bodyind)] }];

    data = [ data;
             gf_forceTotal(bodyind,:); 
             gf_F_Total(bodyind,:);
             gf_forceExcitation(bodyind,:); 
             gf_forceAddedMass(bodyind,:);
             gf_F_AddedMassCorrected(bodyind,:);
             gf_forceRadiationDamping(bodyind,:); 
             gf_forceRestoring(bodyind,:)];
end

% display table
colheadings = { 'mse', 'nmse', 'rmse', 'nrmse', 'mae', 'mare', 'r', 'r2', 'e', 'maxae', 'mxare' };

wid = 16;
% fms = {'d','.4f','.5E'};
fms = {};
fileID = 1;

displaytable (data,colheadings,wid,fms,rowheadings,fileID);


