function [force] = rotorforces_TORUS_SLOTLESS(design, ndrawnstages, nsimpoles, pos)
% calculates the air-gap forces in the slotless torus type machine with
% varying air gap
%
% Syntax
%
% [force] = rotorforces_TORUS_SLOTLESS(design, ndrawnstages, nsimpoles, pos)
%
% Input 
%
%  design - a slotless torus machine design structure
%
%  ndrawnstages - the number of stages to be drawn in the FEA simulation
%    used to create the forces. This is used to reduce the simulation size
%    by drawing less stages than are actually present in the machine
%
%  nsimpoles - the number of poles actually drawn in the simulation. This
%    is used to determine the per-pole forces from the force integral
%
%  pos - vector of rotor displacements positions at which the forces will
%    be determined. The end rotor will be displaced by each position in pos
%    and the forces calculated in each case.
%
% Output
%
%  force - vector of air-gap closing forces of the same size as pos
%    containing the forces at each position in pos
%
% 

    % get the id of the air gap label in the last stage
    aglabelid = findblocklabel_mfemm(design.FemmProblem, [(design.outermagsep/2)-(design.g/2),design.taupm]);
    
    currentpos = 0;

    force = zeros(size(pos));
    
    % generate a temporary file name location
    femfilename = [tempname(), '.fem'];
    
    previousg = design.g;
    originalsurfpos = design.FemmProblem.BlockLabels(aglabelid+1).Coords(1) - design.g/2;
    
    for i = 1:numel(pos)

        % work out how far we need to move from the current position to get to
        % the next position
        disp = pos(i) - currentpos;

        % move the appropriate group
        design.FemmProblem = translategroups_mfemm(design.FemmProblem, ndrawnstages+1, disp, 0);
        
        % destermine the new air gap
        previousg = (previousg + disp);
        
        % move the air gap block label to the point half-way across the new
        % air gap
        design.FemmProblem.BlockLabels(aglabelid+1).Coords(1) = ...
            originalsurfpos + (previousg / 2);
        
        currentpos = currentpos + disp;
        
        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);

        ansfilename = analyse_mfemm(femfilename);
        
        if exist('fpproc_interface_mex', 'file') == 3
            
            solution = fpproc();
            solution.opendocument(ansfilename);
            % get the forces
            solution.clearblock();
            solution.groupselectblock(ndrawnstages+1)
            
            % we divide by the number of poles in the sim to get the per-pole
            % force
            force(i) = solution.blockintegral(18) / nsimpoles;

        else
            %opendocument(ansfilename);
            % get the forces
            mo_clearblock();
            mo_groupselectblock(ndrawnstages+1)
            
            % we divide by the number of poles in the sim to get the per-pole
            % force
            force(i) = mo_blockintegral(18) / nsimpoles;
            
            mo_close;
            mi_close;
        end

    end

end