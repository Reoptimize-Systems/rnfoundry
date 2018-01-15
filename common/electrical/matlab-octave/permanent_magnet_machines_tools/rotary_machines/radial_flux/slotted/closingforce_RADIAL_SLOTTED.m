function [force, filenames] = closingforce_RADIAL_SLOTTED (design, pos, varargin)
% gets the net force between rotor and stator of a slotted radial flux
% machine with a give deflection of the rotor
%
% Syntax
%
% force = closingforce_RADIAL_SLOTTED (design, pos, method)
%
% Input
%
% 

    options.Method = 'full';
    options.UseFemm = false;
    options.Verbose = false;
    options.FilePrefix = tempname;
    options.DeleteFEMFiles = true;
    options.DeleteANSFiles = true;
    options.OnlyGenerateFEMFiles = false;
    
    options = parse_pv_pairs (options, varargin);
    
    % input checking
    assert (ischar (options.Method), 'Method must be a char array');
    check.isLogicalScalar (options.UseFemm, true, 'UseFemm');
    check.isLogicalScalar (options.Verbose, true, 'Verbose');
    assert (ischar (options.FilePrefix), 'FilePrefix must be a char array');
    check.isLogicalScalar (options.DeleteFEMFiles, true, 'DeleteFEMFiles');
    check.isLogicalScalar (options.DeleteANSFiles, true, 'DeleteANSFiles');
    
    assert (anyalldims (pos < design.g), ...
        'Rotor displacements must be smaller than the initial air gap size.');
    
    assert (anyalldims(pos >=0), ...
        'All air gap diplacements for force calculation must be greater than zero');
    
    nshifts = numel (pos);
    
    filenames = cell (nshifts, 2);
    
    if strcmpi (options.Method, 'full')
        
        xshift = pos;
        force = nan (size (xshift));

        % parfor ind = 1:numel (xshift)

        for ind = 1:nshifts
            
            fprintf ('Obtaining closing force for position %d of %d %g\n', ind, nshifts);

            FemmProblem = slottedfemmprob_radial(design, 'DrawingType', 'Full');

            FemmProblem = translategroups_mfemm(FemmProblem, ...
                                         [ FemmProblem.Groups.Magnet, ...
                                           FemmProblem.Groups.RotorBackIron ], ...
                                         -xshift(ind), 0);

            filename = [ options.FilePrefix, sprintf('_closing_force_pos_%03d',  ind), '.fem'];
            
            writefemmfile (filename, FemmProblem);
            
            if ~options.OnlyGenerateFEMFiles
                [ansfilename, femfilename] = analyse_mfemm (filename, options.UseFemm, ~options.Verbose);
            
                if nargout > 1
                    filenames(ind,:) = {femfilename, ansfilename};
                end

                solution = fpproc (ansfilename);

                if options.DeleteFEMFiles
                    delete (femfilename);
                end

                if options.DeleteANSFiles
                    delete (ansfilename); 
                end
                
                % get the forces
                solution.groupselectblock ([FemmProblem.Groups.Magnet, FemmProblem.Groups.RotorBackIron], true);
                force(ind) = abs (solution.blockintegral (18));
            
            else
                if nargout > 1
                    filenames{ind,1} = femfilename;
                end
            end

        end
	
    else
        error ('Only ''full'' method specification is currently supported.');
    end

end


% function [force, filenames] = closingforce_RADIAL (FemmProblems, xshift, yshift, shiftgroups, forcecalcgroups, varargin)
% 
%     options.UseFemm = false;
%     options.Verbose = false;
%     options.FilePrefix = tempname;
%     options.DeleteFEMFiles = true;
%     options.DeleteANSFiles = true;
%     
%     options = parse_pv_pairs (options, varargin);
%     
%     % input checking
%     assert (ischar (options.Method), 'Method must be a char array');
%     check.isLogicalScalar (options.UseFemm, true, 'UseFemm');
%     check.isLogicalScalar (options.Verbose, true, 'Verbose');
%     assert (ischar (options.FilePrefix), 'FilePrefix must be a char array');
%     check.isLogicalScalar (options.DeleteFEMFiles, true, 'DeleteFEMFiles');
%     check.isLogicalScalar (options.DeleteANSFiles, true, 'DeleteANSFiles');
%     
%     assert (samesize (FemmProblems, xshift, yshift), ...
%         'FemmProblems, xshift and yshift must all be the same size');
% 
%     nshifts = numel (FemmProblems);
%     
%     filenames = cell (numel(pos), 2);
%     force = nan (nshifts, 2);
%     
%     for ind = 1:nshifts
% 
%         fprintf ('Obtaining closing force for position %d of %d\n', ind, nshifts);
% 
%         FemmProblem = translategroups_mfemm ( FemmProblems(ind), ...
%                                               shiftgroups, ...
%                                               xshift(ind), ...
%                                               yshift(ind) );
% 
%         filename = [ options.FilePrefix, sprinft('_closing_force_pos_%d_of_%d', ind, nshifts)];
% 
%         writefemmfile (filename, FemmProblems(ind));
% 
%         [ansfilename, femfilename] = analyse_mfemm (filename, options.UseFemm, ~options.Verbose);
% 
%         if nargout > 1
%             filenames(ind,:) = {femfilename, ansfilename};
%         end
% 
%         solution = fpproc (ansfilename);
% 
%         if options.DeleteFEMFiles
%             delete (femfilename);
%         end
% 
%         if options.DeleteANSFiles
%             delete (ansfilename);
%         end
% 
%         % get the forces
%         solution.groupselectblock (forcecalcgroups, true);
%         force(ind,1) = abs (solution.blockintegral (18));
%         force(ind,2) = abs (solution.blockintegral (19));
% 
%     end
% 
% end