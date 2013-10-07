function design = evaluateIbeamstructure_ACPMSM(design, options, initDef)
% DesignIBeamStructure_PMSM: evaluates whether a given I-beam can withstand
% the forces in a linear permanent magnet synchronous machine design
%
% 
    design = ratios2dimensions_ACPMSM(design);
    
    beams = Ibeamsections ./ 1000;
                                                  
    % We will calculatte the deflection of the I-Beams, stopping when we reach
    % an I-Beam that meets the specifications, as this will be the smallest
    % volume with which this can be achieved

    %           IVars(1,1): b, width of the I-beam flanges
    %           IVars(1,2): t, the thickness of the flanges
    %           IVars(1,3): tw, the thickness of the I-Beam vertical part
    %           IVars(1,4): d, height of the I-Beam vertical part
    %
    %           IVars(2,1): bc, thickness of central section
    %           IVars(2,2): Taup, the pole width
    %           IVars(2,3): dt, the tooth height
    %           IVars(2,4): bt, the tooth width
    
    b = beams(:,2);
    t = beams(:,4);
    tw = beams(:,3);
    d = beams(:,6);
    
    area = 2.*b.*t + tw.*d;
    
    volume = area .* design.ls .* options.alphab;
    
    % Now calculate the total volume of beam material that would be
    % required if that beam were used. We do this by determining the number
    % of beams required, as they are not distributed with regard for pole
    % width.
    n_beams = zeros(size(beams,1), 1);
    for i = 1:size(beams,1)
        n_beams(i,1) = numbeams(b(i,1), design.Taup, design.Poles(1) * design.PoleWidth, design.BeamSpreadFactor);
    end
    
    totalVol =  volume .* n_beams .* 2;

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
    if ~isfield(design, 'tols')
        design.tols = [0;0];
    end
    
    % Use listsearch to find the smallest beam total moment of area that
    % will do the job
    r = listsearch(@evaluatestructure_ACPMSM, [b, t, tw, d], ...
        {options.IMethod, design, options.E, options.sections, design.tols, options.alphab, options.gfactor});
    
    % get all the working beams and sort them by their total volume, the
    % smallest of these will be the desired beam
    workingBeams = beams(r:end,:);
    workingBeams = sortrows(workingBeams, 8);
    design.r = workingBeams(1, 10);
    design.beam = workingBeams(1,:);
    
    % Get the actual results using airgapclosure_ACPMSM directly as
    % listsearch only returns the winning list member              
    forceRatio = forceRatio_linear(design, design.beam(2));
    
    [design.new_g, design.actForce] = airgapclosure_ACPMSM(design, options.E, [b, t, tw, d], ...
                            options.sections, options.IMethod, forceRatio, [0;0;0], options.alphab);
    
    design = ratios2dimensions_ACPMSM(design);

%     
%     IVars = [IBeam(2),   IBeam(4),     IBeam(3),    (IBeam(1)-IBeam(4)) ];
%               
%     % Determine what the air-gap is with the structure loaded by the maxwell
%     % stresses
%     new_g = airgapclosure_ACPMSM(design, options.E, IVars, options.sections, initDef, options.alphab);
% 
%     % Check if beam is successful or a failure
%     if min(min(new_g)) >= options.gfactor * design.g
%         varargout{1} = 1;
%         varargout{2} = IBeam;
%     else
%         varargout{1} = 0;
%         varargout{2} = IBeam;
%     end

end