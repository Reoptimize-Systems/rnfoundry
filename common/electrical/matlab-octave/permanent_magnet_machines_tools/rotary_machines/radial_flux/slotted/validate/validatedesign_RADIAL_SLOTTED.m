function design = validatedesign_RADIAL_SLOTTED(design, simoptions)
% calculates various machine properties using an alternative method to
% validate the results of 
    
    % we are going to perform a single FEA simulation of the design to
    % extract some basic info about it to make a simple evaluation of the
    % performance
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'YokeRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent (simoptions.femmmeshoptions, 'CoilRegionMeshSize', -1);

    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];
   
    
    % Draw the sim at position 0
    design.Validation.FemmProblem = slottedfemmprob_radial (design, ...
                            'ArmatureType', design.ArmatureType, ...
                            'NWindingLayers', design.CoilLayers, ...
                            'Position', 0, ...
                            'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                            'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                            'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                            'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
                            'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                            'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
                            'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);
    
    % write the fem file to disk
    writefemmfile (femfilename, design.Validation.FemmProblem);
    % analyse the problem
    ansfilename = analyse_mfemm (femfilename, ...
                                 simoptions.usefemm, ...
                                 simoptions.quietfemm);
	
    if (exist ('fpproc_interface_mex', 'file')==3) && ~simoptions.usefemm

        solution = fpproc (ansfilename);
        solution.smoothon ();
        
        % get the flux linkage in this position
        temp = solution.getcircuitprops('1');
        design.Validation.FluxLinkage = temp(3);        
        
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
    
    design.Validation.CoilInductance = NaN; % TODO: add manual calc of coil inductance for validation

    design.Validation.EMFCoilPeak = peakemfest_ROTARY(abs(design.Validation.FluxLinkage), design.Validation.Omega, design.Poles / 2);
    
    design.Validation.EMFCoilRms = design.Validation.EMFCoilPeak / sqrt(2);
    
    design.Validation.EMFPhasePeak = design.CoilsPerBranch * design.Validation.EMFCoilPeak;
    
    design.Validation.EMFPhaseRms = design.Validation.EMFPhasePeak  / sqrt(2);
    
    design.Validation.IPhasePeak = design.Validation.EMFPhasePeak / (design.PhaseResistance(1) + design.LoadResistance(1));
    design.Validation.ICoilPeak = design.Validation.IPhasePeak / design.Branches;
    
    design.Validation.IPhaseRms = design.Validation.IPhasePeak / sqrt(2);
    design.Validation.ICoilRms = design.Validation.IPhaseRms / design.Branches;
    
    design.Validation.PowerLoadMean = design.Validation.IPhaseRms.^2 * design.LoadResistance * design.Phases;
    
    design.Validation.JCoilRms = design.Validation.ICoilRms / design.ConductorArea;
    
    design.Validation.JCoilPeak = design.Validation.ICoilPeak / design.ConductorArea;
    
    design.Validation.Efficiency = design.LoadResistance / (design.CoilResistance + design.LoadResistance);
    
    reportvalidation (design, simoptions)
    
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
                  
    data = [];
    for ind = 1:numel (colheadings)
        data = [ data, design.Validation.(colheadings{ind})(:) ];
    end
    
    wid = 10;
    fms = {'.2E'};
  
  	displaytable(data,colheadings,wid,fms);
    
    if isfield (simoptions, 'RPM')
        
        fprintf (1, '\nAt simulation RPM of %f:\n\n', simoptions.RPM);
        
        for ind = 3:numel (colheadings)
            if isfield (design, colheadings{ind})
                simpleval = design.Validation.(colheadings{ind})(design.Validation.RPM <= simoptions.RPM*1.001 & design.Validation.RPM > simoptions.RPM*0.999);
                fprintf (1, 'Full %s: %g, Simple %s: %g, Percent Difference: %g\n', ...
                    colheadings{ind}, ...
                    design.(colheadings{ind}), ...
                    colheadings{ind}, ...
                    simpleval(1), ...
                    100 * (simpleval(1) - design.(colheadings{ind}))/ design.(colheadings{ind}) );
            end
        end
    end

end


