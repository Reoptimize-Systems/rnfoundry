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

    if isfield (design, 'FemmDirectFluxLinkage')
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
    design.Validation.FemmSingleFL.FemmProblem = slottedfemmprob_radial (design, ...
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
    writefemmfile (femfilename, design.Validation.FemmSingleFL.FemmProblem);
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
        design.Validation.FemmSingleFL.FluxLinkage = temp(3);

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

    % calculate the phase impedance
    Zphase = design.PhaseResistance(1) + design.LoadResistance(1) ...
                + (1i.*design.Validation.ElectricalOmega.*design.Validation.PhaseInductance);
    design.Validation.PhaseImpedance = Zphase;

    design.Validation.FemmSingleFL.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FemmSingleFL.FluxLinkage), design.Validation.Omega, design.Poles / 2);
    
    design.Validation.FemmSingleFL.EMFCoilRms = design.Validation.FemmSingleFL.EMFCoilPeak / sqrt(2);
    
    design.Validation.FemmSingleFL.EMFPhasePeak = design.CoilsPerBranch * design.Validation.FemmSingleFL.EMFCoilPeak;
    
    design.Validation.FemmSingleFL.EMFPhaseRms = design.Validation.FemmSingleFL.EMFPhasePeak ./ sqrt(2);
    
    design.Validation.FemmSingleFL.IPhasePeak = design.Validation.FemmSingleFL.EMFPhasePeak ./ abs (Zphase);
    design.Validation.FemmSingleFL.ICoilPeak = design.Validation.FemmSingleFL.IPhasePeak ./ design.Branches;
    
    design.Validation.FemmSingleFL.IPhaseRms = design.Validation.FemmSingleFL.EMFPhaseRms ./ real (Zphase);
    design.Validation.FemmSingleFL.ICoilRms = design.Validation.FemmSingleFL.IPhaseRms ./ design.Branches;
    
    design.Validation.FemmSingleFL.PowerLoadMean = design.Validation.FemmSingleFL.IPhaseRms.^2 * design.LoadResistance * design.Phases;
    
    design.Validation.FemmSingleFL.JCoilRms = design.Validation.FemmSingleFL.ICoilRms / design.ConductorArea;
    
    design.Validation.FemmSingleFL.JCoilPeak = design.Validation.FemmSingleFL.ICoilPeak / design.ConductorArea;
    
    design.Validation.FemmSingleFL.Efficiency = design.LoadResistance / (design.CoilResistance + design.LoadResistance);
    
    
    if isfield (design, 'FemmDirectFluxLinkage')
        
        design.Validation.FemmOrigDirectFL.FluxLinkage = max (abs (design.FemmDirectFluxLinkage));

        design.Validation.FemmOrigDirectFL.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FemmOrigDirectFL.FluxLinkage), design.Validation.Omega, design.Poles / 2);

        design.Validation.FemmOrigDirectFL.EMFCoilRms = design.Validation.FemmOrigDirectFL.EMFCoilPeak / sqrt(2);

        design.Validation.FemmOrigDirectFL.EMFPhasePeak = design.CoilsPerBranch * design.Validation.FemmOrigDirectFL.EMFCoilPeak;

        design.Validation.FemmOrigDirectFL.EMFPhaseRms = design.Validation.FemmOrigDirectFL.EMFPhasePeak / sqrt(2);

        design.Validation.FemmOrigDirectFL.IPhasePeak = design.Validation.FemmOrigDirectFL.EMFPhasePeak ./ abs (Zphase);
        design.Validation.FemmOrigDirectFL.ICoilPeak = design.Validation.FemmOrigDirectFL.IPhasePeak / design.Branches;

        design.Validation.FemmOrigDirectFL.IPhaseRms = design.Validation.FemmOrigDirectFL.EMFPhaseRms ./ real (Zphase);
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
    
    mu = mu_0 () * 1000;
    
    R_gm = (design.g + design.tm) / (mu * design.ls * design.Rai * (design.thetas - design.thetasg));
    
    R_TH = (design.Ryo - design.Rai) / (2 * mu * (design.ls * design.Rcm * mean([design.thetacy, design.thetacg])));
    
    R_YH = (design.Rym * design.thetas) / (mu * design.ls * design.ty);
    
    R_biH = (mu * design.Rbm * design.thetas) / (mu * design.ls * design.tbi);
    
    L = design.CoilTurns / (4* R_TH + 3*R_gm + 2*R_YH + 2*R_biH);

end


function reportvalidation (design, simoptions)

    

    colheadings = { 'RPM', ...
                    'Omega', ...
                    'EMFCoilPeak', ...
                    'EMFCoilRms', ...
                    'EMFPhasePeak', ...
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
                
	if isfield (design.Validation, 'FemmSingleFL')
        
        fprintf (1, textwrap2('Validation results based on single new FEMM simulation with the flux linkage extracted directly from this single simulation:'));
        fprintf (1, '\n\n');
        
        data = [];
        for ind = 3:numel (colheadings)
            data = [ data, design.Validation.FemmSingleFL.(colheadings{ind})(:) ];
        end

        wid = 10;
        fms = {'.2g'};

        displaytable([veldata, data], colheadings, wid, fms);

        if isfield (simoptions, 'RPM')

            fprintf (1, '\nAt simulation RPM of %6.2f:\n\n', simoptions.RPM);

            for ind = 3:numel (colheadings)
                if isfield (design, colheadings{ind})
                    simpleval = design.Validation.FemmSingleFL.(colheadings{ind})(design.Validation.RPM <= simoptions.RPM*1.001 & design.Validation.RPM > simoptions.RPM*0.999);
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
        warning ('Field design.Validation.FemmSingleFL not found');
    end
    
    
    if isfield (design.Validation, 'FemmOrigDirectFL')
        
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


