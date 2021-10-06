function design = designstructure_ACPMSM(design, options)
% designstructure_ACPMSM: finds the most appropriate beam for use in the
% support structure of the snapper machine and reports any penalties if no
% adequate beam type is available

    listsearchevalfcn = @evaluatefestructure_ACPMSM;
    structureevalfcn = @feairgapclosure_ACPMSM;
    structvolfcn = @structvol_ACPMSM;
    
    % Now we will design the machine structure
    design = evaluatestructure_FM(design, options, structvolfcn, listsearchevalfcn, structureevalfcn);
    
    if all(design.new_g > 0)
        maxdef = max(design.g - design.new_g);
    else
        % get the deflection for smallest -ve g. Subtract the biggest -ve
        % air gap from the original air gap to get the total deflection
        maxdef = design.g - min(design.new_g(design.new_g<=0));
    end

    if ~isfield(design, 'defThresh')

        design.defThresh = 0.25 * design.g;

    end

    design.defthresh = (1-options.gfactor) * design.g;

    % If the structure fails then there is no configuration which can
    % succeed so we add a design.StructPenalty factor. 
    if maxdef > design.defthresh
        design.StructPenalty = (maxdef - design.defThresh) ./ (design.g - design.defThresh);
    else
        design.StructPenalty = 0;
    end

end