function design = DesignStructure_PMSM(design, options)
% DesignStructure_PMSM: finds the most appropriate beam for use in the
% support structure of the PMSM and reports any penalties if no adequate
% beam type is available

    if isfield(design, 'fens')
        design = rmfield(design, 'fens');
    end
    
    if isfield(design, 'gcells')
        design = rmfield(design, 'gcells');
    end

    beaminfofcn = @completebeaminfo_PMSM;
    AirGapClosurefcn = @feairgapclosure_PMSM;
    structvolfcn = @structvol_PMSM;
    
    % Now we will design the machine structure
    design = evaluatestructure_FM(design, options, structvolfcn, beaminfofcn, AirGapClosurefcn);
    
    if all(design.new_g > 0)
        def = design.g - design.new_g;
    else
        % get the deflection for smallest -ve g. Subtract the biggest -ve
        % air gap from the original air gap to get the total deflection
        def = design.g - design.new_g(design.new_g<=0);
    end
    
    maxdef = max(def(:));

    if ~isfield(options, 'gfactor')

        options.gfactor = 0.75;

    end
    
    design.defThresh = (1-options.gfactor) * design.g;

    % If the structure fails then there is no configuration which can
    % succeed so we add a design.StructPenalty factor. 
    if maxdef > design.defThresh
        design.StructPenalty = (maxdef - design.defThresh) ./ (design.g - design.defThresh);
    else
        design.StructPenalty = 0;
    end

end