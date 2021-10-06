function [force] = rotorforces_TORUS_CORELESS(design, ndrawnstages, nsimpoles, pos)
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
%  nsimpoles - the number of Poles actually drawn in the simulation. This
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

    currentpos = 0;

    force = zeros(size(pos));
    
    % generate a temporary file name location
    femfilename = [tempname(), '.fem'];
    
    for i = 1:numel(pos)

        % work out how far we need to move from the current position to get to
        % the next position
        disp = pos(i) - currentpos;

        % move the appropriate group
        design.FemmProblem = translategroups_mfemm(design.FemmProblem, [ndrawnstages+1, 2*(ndrawnstages+1)], disp, 0);
        
        currentpos = currentpos + disp;
        
        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);

        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);

        ansfilename = analyse_mfemm(femfilename);
        
        if exist('fpproc_interface_mex', 'file') == 3
            
            solution = fpproc();
            solution.opendocument(ansfilename);
            % get the forces
            solution.clearblock();
            solution.groupselectblock(ndrawnstages+1)
            
            % we divide by the number of Poles in the sim to get the per-pole
            % force
            force(i) = solution.blockintegral(18) / nsimpoles;
        else
            %opendocument(ansfilename);
            % get the forces
            mo_clearblock();
            mo_groupselectblock(ndrawnstages+1)
            
            % we divide by the number of Poles in the sim to get the per-pole
            % force
            force(i) = mo_blockintegral(18) / nsimpoles;
            
            mo_close;
            mi_close;
        end

    end

end