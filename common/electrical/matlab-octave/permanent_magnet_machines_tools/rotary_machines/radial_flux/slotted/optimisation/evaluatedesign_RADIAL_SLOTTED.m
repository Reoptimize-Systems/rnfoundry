function [score, design, simoptions, T, Y, results] = evaluatedesign_RADIAL_SLOTTED(design, simoptions)
% evaluates a slotted radial flux permanent magnet machine design and
% assigns it a score
%
% Syntax
%
% [score, design, simoptions, T, Y, results] = evaluatedesign_RADIAL_SLOTTED(design, simoptions)
%
% 

    simoptions = setfieldifabsent(simoptions, 'evaloptions', []);

    simoptions.evaloptions = designandevaloptions_RADIAL_SLOTTED(simoptions.evaloptions);
    
    % pre-screen the design to see if a full simulation is worth it
    [sdesign, ssimoptions] = screendesign_RADIAL_SLOTTED(design, simoptions);
    
    simoptions = setfieldifabsent(simoptions, 'ForceFullSim', false);
    simoptions = setfieldifabsent(simoptions, 'DoStructEval', false);
    
    if ~simoptions.ForceFullSim && ...
         ( ...
            (sdesign.EMFPhaseRms > 2*ssimoptions.max_EMFPhaseRms) ...
            || (sdesign.EMFPhaseRms < 0.25*ssimoptions.min_EMFPhaseRms) ...
            || (sdesign.PowerLoadMean > 2*ssimoptions.max_PowerLoadMean) ...
            || (sdesign.PowerLoadMean < 0.5*ssimoptions.min_PowerLoadMean) ...
         )
        
        T = [];
        Y = [];
        results = [];
        
        % estimate the masses of the components
        sdesign = materialmasses_RADIAL_SLOTTED(sdesign, ssimoptions);
        
        % generate a score for the machine
        [prescore, sdesign, ssimoptions] = machinescore_RADIAL_SLOTTED(sdesign, ssimoptions);
        
        % make the score bigger to encourage getting into the full test zone
        score = prescore * 3; 
        design = sdesign;
        simoptions = ssimoptions;
        
    else
        % run the simulations and return the results using the generic axial
        % flux evaluation function
        [design, simoptions, T, Y, results] = evaluatedesign_RADIAL(design, simoptions);
        
        if simoptions.DoStructEval
            % evaluate the design structurally
            [design, simoptions] = feval(simoptions.evaloptions.structevalfcn, design, simoptions);
        end
        
        % estimate the masses of the components
        design = materialmasses_RADIAL_SLOTTED(design, simoptions);
        
        % generate a score for the machine
        [score, design, simoptions] = machinescore_RADIAL_SLOTTED(design, simoptions);
    end
    
end


function [sdesign, ssimoptions] = screendesign_RADIAL_SLOTTED(design, simoptions)
% pre-screens a design to assess if it is promising enough to 

    sdesign = design;
    ssimoptions = simoptions;
    
%     sdesign.Hc = sdesign.tc(1);
%         
%     % number of slots per pole and phase
%     if ~isfield(sdesign, 'qc')
%         sdesign.qc = fr(sdesign.Qs, sdesign.Phases * sdesign.Poles);
%     else
%         [sdesign.Qs,~] = rat(sdesign.qc * sdesign.Phases * sdesign.Poles);
%     end
%     
%     % slot pitch at the mean slot height
%     sdesign.tausm = sdesign.thetas * sdesign.Rcm;
%     
%     % get the numerator and denominator of qc
%     [sdesign.qcn,sdesign.qcd] = rat(sdesign.qc);
%     
%     % Average coil pitch as defined by (Qs/Poles)
%     sdesign.yp = fr(sdesign.Qs, sdesign.Poles);
%     
%     % get the numerator and denominator of the coil pitch in slots
%     [sdesign.ypn,sdesign.ypd] = rat(sdesign.yp);
%     
%     % calculate the actual coil pitch in slots if not supplied
%     if ~isfield(sdesign, 'yd')
%         if sdesign.ypd == 1
%             % the coil pitch in slots will be the same as the numerator of
%             % yp, being an integral slot winding
%             sdesign.yd = sdesign.ypn;
%         else
%             error('You must specify the coil pitch in fractional slot windings.')
%         end
%     end
%     
%     if sdesign.ypd ~= 1 && sdesign.ypd ~= 2
%     	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
%     end
% 
%     % \Tau_{cs} is the thickness of the winding, i.e. the pitch of a
%     % winding slot
%     sdesign.Wc = sdesign.thetac * sdesign.Rcm;
    
    if ~isfield(sdesign, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [sdesign.CoreLoss.fq, ...
         sdesign.CoreLoss.Bq, ...
         sdesign.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
     
        sdesign.CoreLoss(2) = sdesign.CoreLoss(1);
    end
    
    % do some common design processing things, but skipping calculations of
    % coil properties
    ssimoptions.SkipCheckCoilProps = true;
    [sdesign, ssimoptions] = simfun_RADIAL(sdesign, ssimoptions);
    
    % we are going to perform a single FEA simulation of the design to
    % extract some basic info about it to make a simple evaluation of the
    % performance
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'YokeRegionMeshSize', -1);
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'CoilRegionMeshSize', -1);

    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];
    
    armirongroup = 2;
    
    rmcoilturns = false;
    if ~isfield(sdesign, 'CoilTurns')
        sdesign.CoilTurns = 1;
        rmcoilturns = true;
    end
    
    % Draw the sim at position 0
    [sdesign.FemmProblem, sdesign.RotorDrawingInfo, sdesign.StatorDrawingInfo] = ...
                        slottedfemmprob_radial(sdesign, ...
                            'ArmatureType', sdesign.ArmatureType, ...
                            'NWindingLayers', sdesign.CoilLayers, ...
                            'Position', 0, ...
                            'MagnetRegionMeshSize', ssimoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                            'BackIronRegionMeshSize', ssimoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                            'RotorOuterRegionsMeshSize', ssimoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                            'StatorOuterRegionsMeshSize', ssimoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                            'AirGapMeshSize', ssimoptions.femmmeshoptions.AirGapMeshSize, ...
                            'ShoeGapRegionMeshSize', ssimoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                            'YokeRegionMeshSize', ssimoptions.femmmeshoptions.YokeRegionMeshSize, ...
                            'CoilRegionMeshSize', ssimoptions.femmmeshoptions.CoilRegionMeshSize);

    if rmcoilturns
        sdesign = rmfield(sdesign, 'CoilTurns');
    end
    
    % write the fem file to disk
    writefemmfile(femfilename, sdesign.FemmProblem);
    % analyse the problem
    try
        ansfilename = analyse_mfemm(femfilename, ...
                                    simoptions.usefemm, ...
                                    simoptions.quietfemm);
    catch err
        if strncmp (err.message, 'Material properties have not been defined', 41)
            warning (err.message);
            [sdesign, ssimoptions] = badscore (sdesign, ssimoptions);
            return;
        else
            rethrow (err);
        end
    end
	
    if (exist('fpproc_interface_mex', 'file')==3) && ~ssimoptions.usefemm

        solution = fpproc(ansfilename);
        solution.smoothon();

        % get the cross-sectional area of the armature iron for
        % calcuation of material masses later
        solution.clearblock();
        solution.groupselectblock(sdesign.FemmProblem.Groups.ArmatureBackIron);
        sdesign.ArmatureIronAreaPerPole = solution.blockintegral(5)/2;
        
        % get the peak flux linkage
        temp = solution.getcircuitprops('1');
        peakfl = temp(3);
        
        % get the cross-sectional area of the coil winding bundle
        sdesign.CoilArea = ...
            solution.blockintegral ( 5, ...
                                     sdesign.StatorDrawingInfo.CoilLabelLocations(1,1), ...
                                     sdesign.StatorDrawingInfo.CoilLabelLocations(1,2) );
                                 
        % estimate the coil resistance
%         [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
        
        % get the peak flux density in the armature back iron along
        % center line of a tooth
        NBpnts = 100;
        switch design.ArmatureType
            case 'external'
                [x, y] = pol2cart (repmat (sdesign.thetas, 1, NBpnts), ...
                                   linspace (sdesign.Rai, sdesign.Ryo, NBpnts));
            case 'internal'
                [x, y] = pol2cart (repmat (sdesign.thetas, 1, NBpnts), ...
                                   linspace (sdesign.Ryi, sdesign.Rao, NBpnts));
        end
        Bmag = magn (solution.getb (x, y));

        sdesign.ArmatureToothFluxDensityPeak = max (Bmag);
        
    else
        % open the solution in FEMM
        error('Not yet implemented')

    end
    % tidy up the fea files
    delete (femfilename); delete (ansfilename); 
    
    sdesign = checkcoilprops_AM(sdesign);
    
    if rmcoilturns
       peakfl = peakfl * sdesign.CoilTurns;
    end
    
    % now calculate coil resistance
    sdesign.MTL = rectcoilmtl( sdesign.ls, ...
                               sdesign.yd * sdesign.thetas * sdesign.Rcm, ...
                               mean([sdesign.thetacg, sdesign.thetacy] * sdesign.Rcm) );

    sdesign.CoilResistance = 1.7e-8 * sdesign.CoilTurns * sdesign.MTL ./ (pi * (sdesign.Dc/2)^2);

%     sdesign.PhaseResistance = sdesign.CoilResistance * sdesign.CoilsPerBranch / sdesign.Branches;

%     ssimoptions = simsetup_ROTARY(design, ssimoptions.simfun, ssimoptions.finfun, ...
%                                 'RPM', ssimoptions.RPM, ...
%                                 'PoleCount', ssimoptions.PoleCount, ...
%                                 'odeevfun', 'prescribedmotode_linear', ...
%                                 'simoptions', ssimoptions);
%                             
    omega = rpm2omega(ssimoptions.RPM);
    
    sdesign.CoilInductance = 0;
    
    sdesign.psilookup = [0,1;0,0];
    
    [sdesign, ssimoptions] = finfun_RADIAL(sdesign, ssimoptions);

    sdesign.EMFCoilPeak = peakemfest_ROTARY(abs(peakfl), omega, sdesign.Poles / 2);
    
    sdesign.EMFCoilRms = sdesign.EMFCoilPeak / sqrt(2);
    
    sdesign.EMFPhasePeak = sdesign.CoilsPerBranch * sdesign.EMFCoilPeak;
    
    sdesign.EMFPhaseRms = sdesign.EMFPhasePeak  / sqrt(2);
    
    sdesign.IPhasePeak = sdesign.EMFPhasePeak / (sdesign.PhaseResistance(1) + sdesign.LoadResistance(1));
    sdesign.ICoilPeak = sdesign.IPhasePeak / sdesign.Branches;
    
    sdesign.IPhaseRms = sdesign.IPhasePeak / sqrt(2);
    sdesign.ICoilRms = sdesign.IPhaseRms / sdesign.Branches;
    
    sdesign.PowerLoadMean = sdesign.IPhaseRms.^2 * sdesign.LoadResistance * sdesign.Phases;
    
    sdesign.JCoilRms = sdesign.ICoilRms / sdesign.ConductorArea;
    
    sdesign.JCoilPeak = sdesign.ICoilPeak / sdesign.ConductorArea;
    
    sdesign.Efficiency = 0.9 * (sdesign.LoadResistance / (sdesign.CoilResistance + sdesign.LoadResistance));
    sdesign.TorqueRippleFactor = 0.2;
    sdesign.VoltagePercentTHD = 50;
    if isfield (ssimoptions, 'max_TemperaturePeak')
        sdesign.TemperaturePeak = ssimoptions.max_TemperaturePeak - 1;
    end
    sdesign.CoggingTorquePeak = 0;
    
end

function [sdesign, ssimoptions] = badscore (sdesign, ssimoptions)

	sdesign.JCoilRms = 2 * ssimoptions.max_JCoilRms;
    sdesign.EMFPhaseRms = 5 * ssimoptions.max_EMFPhaseRms;
    sdesign.PowerLoadMean = 2 * ssimoptions.max_PowerLoadMean;
    if isfield (ssimoptions, 'max_TemperaturePeak')
        sdesign.TemperaturePeak = 2 * ssimoptions.max_TemperaturePeak;
    end
    sdesign.TorqueRippleFactor = 2 * ssimoptions.max_TorqueRippleFactor;
    sdesign.VoltagePercentTHD = 2 * ssimoptions.max_VoltagePercentTHD;
    sdesign.CoggingTorquePeak = 2 * ssimoptions.max_CoggingTorquePeak;
    sdesign.Efficiency = 0.5;
    sdesign.ArmatureToothFluxDensityPeak = 2 * ssimoptions.max_ArmatureToothFluxDensityPeak;
    
    % now calculate coil resistance
    sdesign.MTL = rectcoilmtl( sdesign.ls, ...
                               sdesign.yd * sdesign.thetas * sdesign.Rcm, ...
                               mean([sdesign.thetacg, sdesign.thetacy] * sdesign.Rcm) );
                           
    sdesign.CoilTurns = 1;
    sdesign.ArmatureIronAreaPerPole = 1;
    ssimoptions = setfieldifabsent(ssimoptions, 'basescorefcn', 'costscore_AM');
    
end
