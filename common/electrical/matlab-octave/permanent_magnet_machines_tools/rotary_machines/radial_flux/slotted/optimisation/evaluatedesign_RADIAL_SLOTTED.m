function [score, design, simoptions, T, Y, results] = evaluatedesign_RADIAL_SLOTTED(design, simoptions)
% evaluates a slotted radial flux permanent magnet machine design assigns
% it a score
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
    
    % estimate the masses of the components
    sdesign = materialmasses_RADIAL_SLOTTED(sdesign, ssimoptions);
    
    % generate a score for the machine
    [prescore, sdesign, ssimoptions] = machinescore_RADIAL_SLOTTED(sdesign, ssimoptions);
    
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
        
        % make the score bigger to encourage getting into the full test zone
        score = prescore * 2; 
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
    
    sdesign.Hc = sdesign.tc;
        
    % number of slots per pole and phase
    if ~isfield(sdesign, 'qc')
        sdesign.qc = fr(sdesign.Qs, sdesign.Phases * sdesign.Poles);
    else
        [sdesign.Qs,~] = rat(sdesign.qc * sdesign.Phases * sdesign.Poles);
    end
    
    % get the pitch of a whole slot in radians
    sdesign.thetas = (2*pi / sdesign.Qs);
    
    % slot pitch at the mean slot height
    sdesign.tausm = sdesign.thetas * sdesign.Rcm;
    
    % get the numerator and denominator of qc
    [sdesign.qcn,sdesign.qcd] = rat(sdesign.qc);
    
    % Average coil pitch as defined by (Qs/Poles)
    sdesign.yp = fr(sdesign.Qs, sdesign.Poles);
    
    % get the numerator and denominator of the coil pitch in slots
    [sdesign.ypn,sdesign.ypd] = rat(sdesign.yp);
    
    % calculate the actual coil pitch in slots if not supplied
    if ~isfield(sdesign, 'yd')
        if sdesign.ypd == 1
            % the coil pitch in slots will be the same as the numerator of
            % yp, being an integral slot winding
            sdesign.yd = sdesign.ypn;
        else
            error('You must specify the coil pitch in fractional slot windings.')
        end
    end
    
    if sdesign.ypd ~= 1 && sdesign.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    % \Tau_{cs} is the thickness of the winding, i.e. the pitch of a
    % winding slot
    sdesign.Wc = sdesign.thetac * sdesign.Rcm;
    
    if ~isfield(sdesign, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [sdesign.CoreLoss.fq, ...
         sdesign.CoreLoss.Bq, ...
         sdesign.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
     
        sdesign.CoreLoss(2) = sdesign.CoreLoss(1);
    end
    
    [sdesign, ssimoptions] = simfun_RADIAL(sdesign, ssimoptions);
    
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'YokeRegionMeshSize', -1);
    ssimoptions.femmmeshoptions = setfieldifabsent(ssimoptions.femmmeshoptions, 'CoilRegionMeshSize', -1);

    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];
    
    armirongroup = 2;
        
    % Draw the sim 
    [sdesign.FemmProblem, sdesign.outermagsep, sdesign.coillabellocs, sdesign.yokenodeids] = ...
                        slottedfemmprob_radial(sdesign, ...
                            'StatorType', sdesign.StatorType, ...
                            'NWindingLayers', sdesign.CoilLayers, ...
                            'Position', 0, ...
                            'ArmatureBackIronGroup', armirongroup, ...
                            'MagnetRegionMeshSize', ssimoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                            'BackIronRegionMeshSize', ssimoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                            'OuterRegionsMeshSize', ssimoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                            'AirGapMeshSize', ssimoptions.femmmeshoptions.AirGapMeshSize, ...
                            'ShoeGapRegionMeshSize', ssimoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                            'YokeRegionMeshSize', ssimoptions.femmmeshoptions.YokeRegionMeshSize, ...
                            'CoilRegionMeshSize', ssimoptions.femmmeshoptions.CoilRegionMeshSize);

    % write the fem file to disk
    writefemmfile(femfilename, sdesign.FemmProblem);
    % analyse the problem
    ansfilename = analyse_mfemm(femfilename, ...
                                simoptions.usefemm, ...
                                simoptions.quietfemm);

    if (exist('fpproc_interface_mex', 'file')==3) && ~ssimoptions.usefemm

        solution = fpproc(ansfilename);
        solution.smoothon();

        % get the peak flux density in the airgap close to the stator
        % surface
        [x, y] = pol2cart(linspace(0, sdesign.thetap, 100), sdesign.Rso+sdesign.FEMMTol);
        p = solution.getpointvalues(x, y);

%         maxB = max( sqrt( p(2,:).^2 + p(3,:).^2 ) );
        
        % get the cross-sectional area of the armature iron for
        % calcuation of material masses later
        solution.clearblock();
        solution.groupselectblock(armirongroup);
        sdesign.ArmatureIronAreaPerPole = solution.blockintegral(5)/2;
        
        temp = solution.getcircuitprops('1');
        
        peakfl = temp(3);
        
        % estimate the coil resistance
%         [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
        
    else
        % open the solution in FEMM
        error('Not yet implemented')

    end
    
    % now calculate coil resistance
    sdesign.MTL = rectcoilmtl( sdesign.ls, ...
                               sdesign.yd * sdesign.thetas * sdesign.Rcm, ...
                               sdesign.thetac * sdesign.Rcm );

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
    
end
