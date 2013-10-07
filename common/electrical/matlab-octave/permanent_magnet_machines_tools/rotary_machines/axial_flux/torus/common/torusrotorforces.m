function [force] = torusrotorforces(FemmProblem, nstages, nsimpoles, pos, halfdispgroupnos)

    currentpos = 0;

    force = zeros(size(pos));
    
    % generate a temporary file name location
    femfilename = [tempname(), '.fem'];
    
    for i = 1:numel(pos)

        % work out how far we need to move from the current position to get to
        % the next position
        disp = pos(i) - currentpos;

        % move the appropriate group
        FemmProblem = translategroups_mfemm(FemmProblem, nstages+1, disp, 0);

        if nargin > 4
            % move some other supplied groups half the displacement we
            % moved the rotor part (this can be useful for moving air gap
            % lables etc.)
            FemmProblem = translategroups_mfemm(FemmProblem, halfdispgroupnos, disp/2, 0);
        end
        
        currentpos = currentpos + disp;
        
        % write the fem file to disk
        writefemmfile(femfilename, FemmProblem);

    %     % solve the problem 
    %     % blah blah with mexfsolver??    
    %     % mesh the problem
    %     mexfmesher(femfilename);
    %     

        % for now, with femm
        openfemm;
        main_minimize;

        opendocument(femfilename);
        mi_analyse(1);
        mi_loadsolution;

        % get the forces
        mo_clearblock();
        mo_groupselectblock(nstages+1)
        
        % we divide by the number of Poles in the sim to get the per-pole
        % force
        force(i) = mo_blockintegral(18) / nsimpoles;
        
        mo_close;
        mi_close;

    end

end