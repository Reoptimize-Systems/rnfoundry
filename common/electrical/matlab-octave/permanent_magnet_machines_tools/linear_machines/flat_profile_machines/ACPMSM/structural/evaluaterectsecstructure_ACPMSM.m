function design = evaluaterectsecstructure_ACPMSM(design, options)
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
    b = beams(:,3);

    % First calculate the cross-sectional area of all the beams
    % Area is d*b - (d-2*t)*(b-2*t)
    area = d .* b;
    area = area - ((d - (2.*t)) .* (b - (2.*t)));
    
    % Next determine the volume of a single support beam
    volume = area .* design.ls .* options.alphab;
    
    % Now calculate the total volume of beam material that would be
    % required if that beam were used. We do this by determining the number
    % of beams required, as they are not distributed with regard for pole
    % width.
    n_beams = zeros(size(beams,1), 1);
    for i = 1:size(beams,1)
        n_beams(i,1) = numbeams(b(i,1), design.Taup, design.poles(1) * design.Taup, design.BeamSpreadFactor);
    end
    
    totalVol =  volume .* n_beams .* 2;

    % Now we sort the options.sections by the total volume required if that beam
    % were used for this machine (i.e. by ascending cost). So we can use
    % listsearch to find the smallest necessary
    beams(:,4) = n_beams;
    beams(:,5) = totalVol;
    beams(:,6) = n_beams .* MomentOfInertiaY1([beams(:,2:3), beams(:,2)-(2.*beams(:,1)), beams(:,3)-(2.*beams(:,1))], '1.3');
    % sort first by total volume of material required
    beams = sortrows(beams,5);
    % add a total volume index
    beams(:,7) = 1:size(beams,1);
    % resort by total second moment of area
    beams = sortrows(beams,6);

    % check if design tolerances have been specified
    if ~isfield(design, 'tols')
        design.tols = [0;0];
    end

    % Use listsearch to find the smallest beam total moment of area that
    % will do the job
    r = listsearch(@evaluatestructure_ACPMSM, ...
        [beams(:,2:3), beams(:,2)-(2.*beams(:,1)), beams(:,3)-(2.*beams(:,1))], ...
        {options.IMethod, design, options.E, options.sections, design.tols, options.alphab, options.gfactor});

    % get all the working beams and sort them by their total volume, the
    % smallest of these will be the desired beam
    workingBeams = beams(r:end,:);
    workingBeams = sortrows(workingBeams, 5);
    design.r = workingBeams(1, 7);
    design.beam = workingBeams(1,:);
    
    % Get the actual results using airgapclosure_ACPMSM directly as
    % listsearch only returns the winning list member     
    IVars = [design.beam(2),   design.beam(3), design.beam(2)-(2.*design.beam(1)), design.beam(3)-(2.*design.beam(1))];
         
    forceRatio = forceRatio_linear(design, design.beam(3));
    
    [design.new_g, design.actForce] = airgapclosure_ACPMSM(design, options.E, IVars, options.sections, options.IMethod, forceRatio, [0;0;0], options.alphab);
    
    design = ratios2dimensions_ACPMSM(design);

end