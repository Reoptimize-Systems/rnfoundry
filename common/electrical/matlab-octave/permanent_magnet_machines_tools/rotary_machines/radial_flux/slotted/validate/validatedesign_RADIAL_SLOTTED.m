function design = validatedesign_RADIAL_SLOTTED(design, simoptions)
% calculates various machine properties using an alternative method to
% validate the results of 
    

    fprintf (1, [textwrap2(['Beginning validation of radial flux slotted ' ...
                 'permanent magnet machine design. Currently only tooth-wound ' ...
                 'coils are valid.']), '\n\n']);

    % we are going to perform a single FEA simulation of the design to
    % extract some basic info about it to make a simple evaluation of the
    % performance
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'YokeRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'CoilRegionMeshSize', -1);

    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];

    if isfield (design, 'MagSimFEAPeakFluxLinkagePosition')
        pos = design.MagSimFEAPeakFluxLinkagePosition*design.thetap;
    elseif isfield (design, 'FemmDirectFluxLinkage')
        [~,maxflind] = max (abs (design.FemmDirectFluxLinkage));
        pos = (design.feapos(maxflind)+1) * design.thetap + design.FirstSlotCenter;
    else
        pos = 0;
    end

    fprintf (1, [textwrap2(['One set of validation calculations will be ' ...
                 'based on a single FEMM machine simulation, which ' ...
                 'will be combined with analytical calculations. Beginning ' ...
                 'FEMM simulation:']), '\n\n']);
    
    % Draw the sim at position 0
    design.Validation.FemmSingleFLAltInductance.FemmProblem = slottedfemmprob_radial (design, ...
                            'NPolePairs', 1, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers, ...
                            'Position', pos, ...
                            'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                            'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                            'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                            'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
                            'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                            'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
                            'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);
    
    % write the fem file to disk
    writefemmfile (femfilename, design.Validation.FemmSingleFLAltInductance.FemmProblem);
    % analyse the problem
    ansfilename = analyse_mfemm (femfilename, ...
                                 simoptions.usefemm, ...
                                 simoptions.quietfemm);
	
	fprintf (1, '\n');

    if (exist ('fpproc_interface_mex', 'file')==3) && ~simoptions.usefemm

        solution = fpproc (ansfilename);
        solution.smoothon ();

        % get the flux linkage in this position
        temp = solution.getcircuitprops('1');
        design.Validation.FemmSingleFLAltInductance.FluxLinkage = temp(3);

    else
        % open the solution in FEMM
        error('Not yet implemented')
    end
    
    % tidy up the fea files
    delete (femfilename); delete (ansfilename); 
	
    design.Validation.Omega = rpm2omega (1:500:5000);
    
    if isfield (simoptions, 'RPM')
        design.Validation.Omega = sort ([design.Validation.Omega, rpm2omega(simoptions.RPM)]);
    end
    
    design.Validation.RPM = omega2rpm (design.Validation.Omega);
    design.Validation.ElectricalOmega = design.Validation.Omega * (design.Poles / 2);
    
    design.Validation.CoilInductance = inductancefromreluctnet (design);
    design.Validation.PhaseInductance = design.Validation.CoilInductance * design.CoilsPerBranch / design.Branches;
    
    fprintf (1, [textwrap2(['Comparison of FEMM and reluctance network Inductnaces:']), '\n\n']);
    fprintf (1, 'FEMM Coil Inductance: %4.3f mH, Reluctance Coil Inductance %4.3f mH\n', ...
        design.CoilInductance(1)*1000, design.Validation.CoilInductance*1000);
    fprintf (1, 'FEMM Phase Inductance: %4.3f mH, Reluctance Phase Inductance %4.3f mH\n\n', ...
        design.PhaseInductance(1)*1000, design.Validation.PhaseInductance*1000);

    % calculate the phase impedance
    Zphase = design.PhaseResistance(1) + design.LoadResistance(1) ...
                + (1i.*design.Validation.ElectricalOmega.*design.Validation.PhaseInductance);
            
    design.Validation.FemmSingleFLAltInductance.PhaseImpedance = Zphase;

    design.Validation.FemmSingleFLAltInductance.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FemmSingleFLAltInductance.FluxLinkage), design.Validation.Omega, design.Poles / 2);
    
    design.Validation.FemmSingleFLAltInductance.EMFCoilRms = design.Validation.FemmSingleFLAltInductance.EMFCoilPeak / sqrt(2);
    
    design.Validation.FemmSingleFLAltInductance.EMFPhasePeak = design.CoilsPerBranch * design.Validation.FemmSingleFLAltInductance.EMFCoilPeak;
    
    design.Validation.FemmSingleFLAltInductance.EMFPhaseRms = design.Validation.FemmSingleFLAltInductance.EMFPhasePeak ./ sqrt(2);
    
    design.Validation.FemmSingleFLAltInductance.IPhasePeak = design.Validation.FemmSingleFLAltInductance.EMFPhasePeak ./ abs (Zphase);
    design.Validation.FemmSingleFLAltInductance.ICoilPeak = design.Validation.FemmSingleFLAltInductance.IPhasePeak ./ design.Branches;
    
    design.Validation.FemmSingleFLAltInductance.IPhaseRms = design.Validation.FemmSingleFLAltInductance.EMFPhaseRms ./ abs (Zphase);
    design.Validation.FemmSingleFLAltInductance.ICoilRms = design.Validation.FemmSingleFLAltInductance.IPhaseRms ./ design.Branches;
    
    design.Validation.FemmSingleFLAltInductance.PowerLoadMean = design.Validation.FemmSingleFLAltInductance.IPhaseRms.^2 * design.LoadResistance * design.Phases;
    
    design.Validation.FemmSingleFLAltInductance.JCoilRms = design.Validation.FemmSingleFLAltInductance.ICoilRms / design.ConductorArea;
    
    design.Validation.FemmSingleFLAltInductance.JCoilPeak = design.Validation.FemmSingleFLAltInductance.ICoilPeak / design.ConductorArea;
    
    design.Validation.FemmSingleFLAltInductance.Efficiency = design.LoadResistance / (design.CoilResistance + design.LoadResistance);
    
    
    % calculate the phase impedance
    Zphase = design.PhaseResistance(1) + design.LoadResistance(1) ...
                + (1i.*design.Validation.ElectricalOmega.*design.PhaseInductance(1));
            
	design.Validation.FemmSingleFLFEMMInductance.FluxLinkage = design.Validation.FemmSingleFLAltInductance.FluxLinkage;
    
    design.Validation.FemmSingleFLFEMMInductance.PhaseImpedance = Zphase;
    
    design.Validation.FemmSingleFLFEMMInductance.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FemmSingleFLFEMMInductance.FluxLinkage), design.Validation.Omega, design.Poles / 2);
    
    design.Validation.FemmSingleFLFEMMInductance.EMFCoilRms = design.Validation.FemmSingleFLFEMMInductance.EMFCoilPeak / sqrt(2);
    
    design.Validation.FemmSingleFLFEMMInductance.EMFPhasePeak = design.CoilsPerBranch * design.Validation.FemmSingleFLFEMMInductance.EMFCoilPeak;
    
    design.Validation.FemmSingleFLFEMMInductance.EMFPhaseRms = design.Validation.FemmSingleFLFEMMInductance.EMFPhasePeak ./ sqrt(2);
    
    design.Validation.FemmSingleFLFEMMInductance.IPhasePeak = design.Validation.FemmSingleFLFEMMInductance.EMFPhasePeak ./ abs (Zphase);
    design.Validation.FemmSingleFLFEMMInductance.ICoilPeak = design.Validation.FemmSingleFLFEMMInductance.IPhasePeak ./ design.Branches;
    
    design.Validation.FemmSingleFLFEMMInductance.IPhaseRms = design.Validation.FemmSingleFLFEMMInductance.EMFPhaseRms ./ abs (Zphase);
    design.Validation.FemmSingleFLFEMMInductance.ICoilRms = design.Validation.FemmSingleFLFEMMInductance.IPhaseRms ./ design.Branches;
    
    design.Validation.FemmSingleFLFEMMInductance.PowerLoadMean = design.Validation.FemmSingleFLFEMMInductance.IPhaseRms.^2 * design.LoadResistance * design.Phases;
    
    design.Validation.FemmSingleFLFEMMInductance.JCoilRms = design.Validation.FemmSingleFLFEMMInductance.ICoilRms / design.ConductorArea;
    
    design.Validation.FemmSingleFLFEMMInductance.JCoilPeak = design.Validation.FemmSingleFLFEMMInductance.ICoilPeak / design.ConductorArea;
    
    design.Validation.FemmSingleFLFEMMInductance.Efficiency = design.LoadResistance / (design.CoilResistance + design.LoadResistance);
    
    
    
    if isfield (design, 'FemmDirectFluxLinkage')
        
        design.Validation.FemmOrigDirectFL.FluxLinkage = max (abs (design.FemmDirectFluxLinkage));

        design.Validation.FemmOrigDirectFL.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FemmOrigDirectFL.FluxLinkage), design.Validation.Omega, design.Poles / 2);

        design.Validation.FemmOrigDirectFL.EMFCoilRms = design.Validation.FemmOrigDirectFL.EMFCoilPeak / sqrt(2);

        design.Validation.FemmOrigDirectFL.EMFPhasePeak = design.CoilsPerBranch * design.Validation.FemmOrigDirectFL.EMFCoilPeak;

        design.Validation.FemmOrigDirectFL.EMFPhaseRms = design.Validation.FemmOrigDirectFL.EMFPhasePeak / sqrt(2);

        design.Validation.FemmOrigDirectFL.IPhasePeak = design.Validation.FemmOrigDirectFL.EMFPhasePeak ./ abs (Zphase);
        design.Validation.FemmOrigDirectFL.ICoilPeak = design.Validation.FemmOrigDirectFL.IPhasePeak / design.Branches;

        design.Validation.FemmOrigDirectFL.IPhaseRms = design.Validation.FemmOrigDirectFL.EMFPhaseRms ./ abs (Zphase);
        design.Validation.FemmOrigDirectFL.ICoilRms = design.Validation.FemmOrigDirectFL.IPhaseRms / design.Branches;

        design.Validation.FemmOrigDirectFL.PowerLoadMean = design.Validation.FemmOrigDirectFL.IPhaseRms.^2 * design.LoadResistance * design.Phases;

        design.Validation.FemmOrigDirectFL.JCoilRms = design.Validation.FemmOrigDirectFL.ICoilRms / design.ConductorArea;

        design.Validation.FemmOrigDirectFL.JCoilPeak = design.Validation.FemmOrigDirectFL.ICoilPeak / design.ConductorArea;

        design.Validation.FemmOrigDirectFL.Efficiency = design.LoadResistance(1) / (design.CoilResistance + design.LoadResistance(1));
    
    end

    fprintf (1, '\n');
    reportvalidation (design, simoptions)

end

function L = inductancefromreluctnet (design)
% calculate the inductance using a reluctance networ approach, only works
% for single tooth wound coils at present
%
% Syntax
%
% L = inductancefromreluctnet (design)
%
% Description
%
% This function calculates the inductance of a single machine coil using a
% reluctance network approach. The circuit used is shown below, which is
% valid only for a tooth-wound single coil.

    fprintf (1, [textwrap2(['Note that the inductance calculated by this ' ...
                 'validation routine is based on a reluctance network ' ...
                 'only valid for a tooth-wound coil.']), '\n']);
    
    mu = mu_0 () * 3000;
    
    % using ari gap reluctance calculation method from [1]
    %
    % [1] A Balakrishnan A, "Air-Gap Reluctance and Inductance Calculations
    % for Mag Ckts Using a Schwarz–Christoffel Trans"
    %
    wc = design.Rai * (design.thetas - design.thetasg);
    la = design.g + design.tm;
    h = design.tsb + sum(design.tc);
    
    R_gm = 1 / (mu_0 * ((wc / la) + (4/pi)*(1 + log (pi * h / (4 * la))))) / design.ls;
    
%     R_gm = (design.g + design.tm) ...
%         / (mu_0 * design.ls * design.Rai * (design.thetas - design.thetasg + (design.thetasg/4)));
    
    R_TH = (design.Ryo - design.Rai) ...
        / (2 * mu_0 * design.MagFEASimMaterials.ArmatureYoke.Mu_x * (design.ls * design.Rcm * (design.thetas - mean([design.thetacy, design.thetacg]))));
    
    R_YH = (design.Rym * design.thetas) ...
        / (mu_0 * design.MagFEASimMaterials.ArmatureYoke.Mu_x * design.ls * design.ty);
    
    R_biH = (design.Rbm * design.thetas) / (mu * design.ls * design.tbi);
    
    L = design.CoilTurns^2 / (R_TH + 1.25*R_gm + 0.5*R_YH + 0.5*R_biH);

end


function reportvalidation (design, simoptions)

    

    colheadings = { 'RPM', ...
                    'Omega', ...
                    'EMFCoilPeak', ...
                    'EMFCoilRms', ...
                    'EMFPhasePeak', ...
                    'EMFPhaseRms', ...
                    'IPhasePeak', ...
                    'ICoilPeak', ...
                    'IPhaseRms', ...
                    'ICoilRms', ...
                    'JCoilRms',  ...
                    'JCoilPeak', ...
                    'PowerLoadMean' };
                
	veldata = [ design.Validation.RPM(:), design.Validation.Omega(:) ];
    
	if isfield (design.Validation, 'FemmSingleFLAltInductance')
        
        fprintf (1, textwrap2(['Validation results based on single new FEMM ', ...
            'simulation with the flux linkage extracted directly from this ', ...
            'single simulation, using the Inductance calculated ', ...
            'using a reluctance network approach:']));
        fprintf (1, '\n\n');
        
        data = [];
        for ind = 3:numel (colheadings)
            data = [ data, design.Validation.FemmSingleFLAltInductance.(colheadings{ind})(:) ];
        end

        wid = 10;
        fms = {'.2g'};

        displaytable([veldata, data], colheadings, wid, fms);

        if isfield (simoptions, 'RPM')

            fprintf (1, '\nAt simulation RPM of %6.2f:\n\n', simoptions.RPM);

            for ind = 3:numel (colheadings)
                if isfield (design, colheadings{ind})
                    simpleval = design.Validation.FemmSingleFLAltInductance.(colheadings{ind})(design.Validation.RPM <= simoptions.RPM*1.001 & design.Validation.RPM > simoptions.RPM*0.999);
                    fprintf (1, 'Full %s: %g, Simple %s: %g, Percent Difference: %g\n', ...
                        colheadings{ind}, ...
                        design.(colheadings{ind}), ...
                        colheadings{ind}, ...
                        simpleval(1), ...
                        100 * (simpleval(1) - design.(colheadings{ind}))/ design.(colheadings{ind}) );
                end
            end
        end
        
    else
        warning ('Field design.Validation.FemmSingleFLAltInductance not found');
    end
    
    
    if isfield (design.Validation, 'FemmSingleFLFEMMInductance')
        
        fprintf (1, '\n');
        
        fprintf (1, textwrap2(['Validation results based on single new FEMM ', ...
            'simulation with the flux linkage extracted directly from this ', ...
            'single simulation, using the Inductance previously calculated ', ...
            'using FEA:']));
        
        fprintf (1, '\n\n');
        
        data = [];
        for ind = 3:numel (colheadings)
            data = [ data, design.Validation.FemmSingleFLFEMMInductance.(colheadings{ind})(:) ];
        end
        
        wid = 10;
        fms = {'.2g'};
        
        displaytable([veldata, data], colheadings, wid, fms);
        
        if isfield (simoptions, 'RPM')
            
            fprintf (1, '\nAt simulation RPM of %6.2f:\n\n', simoptions.RPM);
            
            for ind = 3:numel (colheadings)
                if isfield (design, colheadings{ind})
                    simpleval = design.Validation.FemmSingleFLFEMMInductance.(colheadings{ind})(design.Validation.RPM <= simoptions.RPM*1.001 & design.Validation.RPM > simoptions.RPM*0.999);
                    fprintf (1, 'Full %s: %g, Simple %s: %g, Percent Difference: %g\n', ...
                        colheadings{ind}, ...
                        design.(colheadings{ind}), ...
                        colheadings{ind}, ...
                        simpleval(1), ...
                        100 * (simpleval(1) - design.(colheadings{ind}))/ design.(colheadings{ind}) );
                end
            end
        end
        
    else
        warning ('Field design.Validation.FemmSingleFLAltInductance not found');
    end
    
    
    if isfield (design.Validation, 'FemmOrigDirectFL')
        
        fprintf (1, '\n');
        
        fprintf (1, textwrap2('Validation using maximum value of flux linkage values taken direct from FEMM simulations done in main simfun:'));
        fprintf (1, '\n\n');
        
        data = [];
        for ind = 3:numel (colheadings)
            data = [ data, design.Validation.FemmOrigDirectFL.(colheadings{ind})(:) ];
        end

        wid = 10;
        fms = {'.2G'};

        displaytable([veldata, data], colheadings, wid, fms);

        if isfield (simoptions, 'RPM')

            fprintf (1, '\nAt simulation RPM of %6.2f:\n\n', simoptions.RPM);

            for ind = 3:numel (colheadings)
                if isfield (design, colheadings{ind})
                    simpleval = design.Validation.FemmOrigDirectFL.(colheadings{ind})(design.Validation.RPM <= simoptions.RPM*1.001 & design.Validation.RPM > simoptions.RPM*0.999);
                    fprintf (1, 'Full %s: %g, FEMM Direct FL %s: %g, Percent Difference: %g\n', ...
                        colheadings{ind}, ...
                        design.(colheadings{ind}), ...
                        colheadings{ind}, ...
                        simpleval(1), ...
                        100 * (simpleval(1) - design.(colheadings{ind}))/ design.(colheadings{ind}) );
                end
            end
        end
    else
        warning ('Field design.Validation.FemmOrigDirectFL not found');
    end

end


