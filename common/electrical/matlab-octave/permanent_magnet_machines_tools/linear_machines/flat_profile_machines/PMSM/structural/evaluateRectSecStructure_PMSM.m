function design = evaluateRectSecStructure_PMSM(design, options)
% evaluateRectSecStructure_PMSM: evaluates whether any of a given selection
% of beams with rectangular sections can withstand the forces in a linear
% permanent magnet synchronous machine design, and which of these achieves
% this
%
% 
    % Load the database of available beams
    % beams columns: t, d, b
    beams = rectangularsections;
    
    % Change to metres from mm
    beams = beams ./ 1000;
    
    t = beams(:,1);
    d = beams(:,2);
    di = d - (2*t);
    b = beams(:,3);
    bi = b - (2*t);

    % First calculate the cross-sectional area of all the beams
    % Area is d*b - (d-2*t)*(b-2*t)
    area = CSArea([d, b, di, bi], options.IMethod);
    
    % Next determine the volume of a single support beam
    volume = area .* design.ls .* options.alphab;
    
    % Now calculate the total volume of beam material that would be
    % required if that beam were used. We do this by determining the number
    % of beams required, as they are not distributed with regard for pole
    % width.
    n_beams = zeros(size(beams,1), 1);
    for i = 1:size(beams,1)
        n_beams(i,1) = numbeams(span1([d(i), b(i), di(i), bi(i)], options.IMethod), design.Wp, design.poles(1) * design.Wp, design.BeamSpreadFactor);
    end
    
    totalVol = design.sides .* volume .* n_beams;

    % Now we sort the options.sections by the total volume required if that beam
    % were used for this machine (i.e. by ascending cost). So we can use
    % listsearch to find the smallest necessary
    beams(:,4) = n_beams;
    beams(:,5) = totalVol;
    beams(:,6) = n_beams .* MomentOfInertiaY1([d, b, di, bi], '1.3');
    % sort first by total volume of material required
    beams = sortrows(beams,5);
    % add a total volume index
    beams(:,7) = 1:size(beams,1);
    % resort by total second moment of area
    beams = sortrows(beams,6);
    
    % check if design tolerances have been specified
    if isfield(design, 'tols')
        if design.sides == 1 && size(design.tols,1) ~= 2 
            error('Tolerances have been specified for three components when only two are required for a single-sided machine')
        elseif design.sides == 2 && size(design.tols,1) ~= 3
            error('An incorrect number of tolerances have been specified for a double-sided machine')
        end
    else
        if design.sides == 1
            design.tols = [0;0];
        elseif design.sides == 2
            design.tols = [0;0;0];
        end
    end
    
    % Use listsearch to find the smallest beam total moment of area that
    % will do the job
%     [r, bestbeam] = listsearch(@evaluateStructure_PMSM, ...
%         [beams(:,2:3), beams(:,2)-(2.*beams(:,1)), beams(:,3)-(2.*beams(:,1))], ...
%         {options.IMethod, design, options.E, options.sections, design.mode(1), design.tols, options.alphab, options.gfactor});
    
    [r, bestbeam] = listsearch(@evaluatefestructure_PMSM, ...
                               [beams(:,2:3), beams(:,2)-(2.*beams(:,1)), beams(:,3)-(2.*beams(:,1))], ...
                               {design, options});
    
    if r == -1
        % no beam worked, so return the largest
        design.beam = beams(end,:);
    else
        % get all the working beams and sort them by their total volume, the
        % smallest of these will be the desired beam
        workingBeams = beams(r:end,:);
        workingBeams = sortrows(workingBeams, 5);
        design.r = workingBeams(1, 7);
        design.beam = workingBeams(1,:);
    end
    
    % Get the actual results using airgapclosure_PMSM directly as
    % listsearch only returns the winning list member     
    IVarsCell = {[design.beam(2), design.beam(3), design.beam(2)-(2.*design.beam(1)), design.beam(3)-(2.*design.beam(1))], ...
                 [design.hba, design.Wp, design.ht, design.Wt]};
         
%     forceRatio = forceratio_linear(design, span1(IVarsCell{1}, design.IMethod));
    
%     [design.new_g, design.actForce] = airgapclosure_PMSM(design, ...
%         options.E, IVarsCell, options.IMethod, options.sections, ...
%         design.mode(1), forceRatio, design.tols, options.alphab);

    [design.new_g, design.actForce] = feairgapclosure_PMSM(design, options, IVarsCell);
    
    design.supportStructIVars = [design.beam(2), design.beam(3), design.beam(2)-(2.*design.beam(1)), design.beam(3)-(2.*design.beam(1))];
    design.supportStructIMethod = options.IMethod;
    
    design = ratios2dimensions_PMSM(design);

end