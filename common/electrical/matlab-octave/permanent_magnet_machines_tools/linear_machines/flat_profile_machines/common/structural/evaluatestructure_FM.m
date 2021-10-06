function design = evaluatestructure_FM(design, options, structvolfcn, airgapclosurefcn)
% evaluatestructure_FM: evaluates whether any of a given selection of beams
% can withstand the forces in a linear permanent magnet machine design of
% flat cross-section, and which of these achieves this
%
%

    fprintf(1, '\nBeginning Structural Evaluation: %s', datestr(now));
    
    design.OuterStructureBeamIMethod = options.IMethod;
    
    % Load the database of available beams base don the beam type stored in
    % the options structure
    IVars = beamvars(options.IMethod);
    
    IVarsLength = size(IVars, 2);

    beams = zeros(size(IVars, 1), IVarsLength + 4);
    
    beams(:, 1:IVarsLength) = IVars;
    
    beams(:, IVarsLength + 1) = MomentOfInertiaY1(IVars, design.OuterStructureBeamIMethod);
    
    % sort by moment of inertia
    beams = sortrows(beams, IVarsLength + 1);
    
    % create an index to beam with next biggest Moment of inertia
    beams(1:end-1, IVarsLength + 2) = 2:size(beams,1);
    
    % No stronger beam than last beam, so it should point to itself
    beams(end, IVarsLength + 2) = size(beams,1);
    IVarsLinkCol = IVarsLength + 2;
    
    for i = 1:size(IVars, 1)
        
        % get the number of pole support beams used in this design
        beams(i, IVarsLength + 4) = numbeams(span1(beams(i, 1:size(IVars, 2)), options.IMethod), ...
                                             design.PoleWidth, ...
                                             design.PowerPoles * design.PoleWidth, ...
                                             design.BeamSpreadFactor);
        
        % Set the outer structure beam variables, the first row of these
        % will be the outer pole support variables, while the second is the
        % variables for the members to which these outer pole supports are
        % fixed (typically referrred to with the variable name
        % 'SupportBeams' elsewhere. The support beams are the next
        % strongest available beam to the beam used by the pole supports
        design.OuterStructureBeamVars = [beams(i, 1:IVarsLength); 
                                         beams( beams(i,IVarsLinkCol), 1:IVarsLength )];                            
        
        % Determine the total volume of material that would be required in
        % each case
        beams(i, IVarsLength + 3) = structvolfcn(design, options);
        
    end
% keyboard
    % remove all structures which result in too many beams (usually done
    % for computational speed reasons)
    beams(:,IVarsLength + 5) = beams(:,IVarsLength + 4) <= max(options.maxPoleSupportBeams, design.PowerPoles);

    % sort by the total volume of the structural material
    beams = sortrows(beams, IVarsLength + 3);
    
    BeamInfo = [];
    
    % Now test each beam combination in turn, stopping at the combinations
    % with the lowest amount of required material which works
    for i = 1:size(beams, 1)
        
        if beams(i,IVarsLength + 5) 
            
            design.OuterStructureBeamVars = [beams(i, 1:IVarsLength); 
                                             beams( beams(i,IVarsLinkCol), 1:IVarsLength )];   

            design.OuterStructureNBeams = beams(i, IVarsLength + 4);

            % the first time the evaluation function is called it must create a
            % structure (BeamInfo) which is subsequently passes back on the
            % next calls so information can be maintained from one call to the
            % next
            if isempty(BeamInfo)
                [result, design, BeamInfo] = evaluatefestructure_FM(design, options, airgapclosurefcn);
            else
                [result, design, BeamInfo] = evaluatefestructure_FM(design, options, airgapclosurefcn, BeamInfo);
            end

            if result
                break;
            end
        
        end

    end
    
    % Ensure we have the details of the last beams that worked
    design.OuterStructureBeamVars = [beams(i, 1:IVarsLength); 
                                     beams( beams(i,IVarsLinkCol), 1:IVarsLength )];  
                                 
    design.OuterStructureBeamRow = i;
    
    design.StructVol = beams(i, IVarsLength + 3);
    
    % Get the actual results using the air gap closure function directly as
    % listsearch only returns the winning list member     
    [design.new_g, design.actForce] = airgapclosurefcn(design, options, BeamInfo);

    fprintf(1, '\nCompleted Structural Evaluation: %s\n', datestr(now));
    
end