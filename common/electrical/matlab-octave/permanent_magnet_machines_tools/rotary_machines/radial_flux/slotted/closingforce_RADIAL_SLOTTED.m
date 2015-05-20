function force = closingforce_RADIAL_SLOTTED (design, pos, method)
% gets the net force between rotor and stator of a slotted radial flux
% machine with a give deflection of the rotor
%
% Syntax
%
% force = closingforce_RADIAL_SLOTTED (design, pos, method)
%
% 
    if anyalldims (pos >= design.g)
        error ('Displacements cannot be greater than the initial air gap size.');
    end
    
    if nargin < 3
        method = 'full';
    end
    
    if strcmpi (method, 'full')
        
        yshift = pos;
        force = ones (size (yshift)) * nan;

        % parfor ind = 1:numel (yshift)

        for ind = 1:numel (yshift)

            FemmProblem = slottedfemmprob_radial(design, 'DrawingType', 'Full');

            FemmProblem = translategroups_mfemm(FemmProblem, ...
                                         [ FemmProblem.Groups.Magnet, ...
                                           FemmProblem.Groups.BackIron ], ...
                                         0, yshift(ind));

            [ansfilename, femfilename] = analyse_mfemm (FemmProblem);

            solution = fpproc (ansfilename);

            delete (ansfilename); delete (femfilename);

            % get the forces
            solution.groupselectblock ([FemmProblem.Groups.Magnet, FemmProblem.Groups.BackIron], true);
            force(ind) = solution.blockintegral (19);

        end
	
    else
        error ('Only ''full'' method specification is currently supported.');
    end

end