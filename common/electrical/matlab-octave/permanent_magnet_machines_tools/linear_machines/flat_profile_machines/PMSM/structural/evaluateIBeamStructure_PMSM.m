function design = evaluateIBeamStructure_PMSM(design, options)
% evaluateIBeamStructure_PMSM: evaluates whether any of a given selection
% of I-beams can withstand the forces in a linear permanent magnet
% synchronous machine design, and which of these achieves this
%
% 
    % Load the database of available I-Beams
    % beams columns: h, b, tw, t, r, d
    beams = Ibeamsections;
    
    % Change to metres from mm
    beams = beams ./ 1000;
    
    b = beams(:,2);
    t = beams(:,4);
    tw = beams(:,3);
    d = beams(:,1)-beams(:,4);

    % First calculate the cross-sectional area of all the beams
    % Area is 2*b*t + tw*d
    area = 2.*b.*t + tw.*d;
    
    % Next determine the volume of a single support beam
    volume = area .* design.ls .* options.alphab;
    
    % Now calculate the total volume of beam material that would be
    % required if that beam were used. We do this by determining the number
    % of beams required, as they are not distributed with regard for pole
    % width.
    n_beams = zeros(size(beams,1), 1);
    for i = 1:size(beams,1)
        n_beams(i,1) = numbeams(b(i,1), design.Wp, design.Poles(2) * design.Wp, design.BeamSpreadFactor);
    end
    
    totalVol =  design.sides .* volume .* n_beams;

    % Now we sort the options.sections by the total volume required if that beam
    % were used for this machine (i.e. by ascending cost). So we can use
    % listsearch to find the smallest necessary
    beams(:,7) = n_beams;
    beams(:,8) = totalVol;
    beams(:,9) = n_beams .* MomentOfInertiaY1([b, t, tw, d], '1.6');
    % sort first by total volume of material required
    beams = sortrows(beams,8);
    % add a total volume index
    beams(:,10) = 1:size(beams,1);
    % resort by total second moment of area
    beams = sortrows(beams,9);

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

    % Use listsearch to find the smallest beam moment of area that will do
    % the job
%     [r, bestbeam] = listsearch(@evaluateStructure_PMSM, ...
%                         [beams(:,2), beams(:,4), beams(:,3), beams(:,1)-beams(:,4)], ...
%                         {options.IMethod, design, options.E, options.sections, design.mode(1), ...
%                         design.tols, options.alphab, options.gfactor});
                    
    [r, bestbeam] = listsearch(@evaluatefestructure_PMSM, ...
                    [beams(:,2), beams(:,4), beams(:,3), beams(:,1)-beams(:,4)], ...
                    {design, options});

                    
    if r == -1
        % no beam worked, so return the largest
        design.beam = beams(end,:);
    else
        % get all the working beams and sort them by their total volume, the
        % smallest of these will be the desired beam
        workingBeams = beams(r:end,:);
        workingBeams = sortrows(workingBeams, 8);
        design.r = workingBeams(1, 10);
        design.beam = workingBeams(1,:);
    end
    
    % Get the actual results using airgapclosure_PMSM directly as
    % listsearch only returns the winning list member     
    IVarsCell = {[design.beam(2),   design.beam(4),     design.beam(3),    (design.beam(1)-design.beam(4))],;...
                 [design.hba,         design.Wp,         design.ht,           design.Wt              ]};
         
%     forceRatio = forceratio_linear(design, design.beam(1));
%     
%     [design.new_g, design.actForce] = airgapclosure_PMSM(design, options.E, IVarsCell, options.IMethod, ...
%                                                          options.sections, design.mode(1), forceRatio, ...
%                                                          design.tols, options.alphab);

    [design.new_g, design.actForce] = feairgapclosure_PMSM(design, options, IVarsCell);
    
    design.supportStructIVars = [design.beam(2),   design.beam(4),     design.beam(3),    (design.beam(1)-design.beam(4))];
    design.supportStructIMethod = options.IMethod;
    
    design = ratios2dimensions_PMSM(design);

end