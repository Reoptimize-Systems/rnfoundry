function [output, design] = getfeadata_PMSM(design, Jcoil, xRVWp)
% getfeadata_PMSM: performs FEA simulations for a PMSM design at a number
% of different positions and coil current densities and extracts relevent
% data from the sims.
%
% PMSM machine diagram with dimensions shown below:
%
%           hba        ht           g      hm  hbf
%         <-----><------------><---------><--><---->
%         |      |_____________            ___|    |  ^
%         |                    | ^     ^  |   |    |  :
%         |                    | : Wt  :  |   |    |  :
%         |       _____________| v     :  |   |    |  :
%         |      |                     :  |   |    |  :
%         |    ^ |_____________        :  |   |    |  :
%         |    : |             |       :  |   |    |  :
%         | Wc : |             |    Wm :  |   |    |  : Wp
%         |    : |_____________|       :  |   |    |  :
%         |    v |               ^Ws   :  |   |    |  :
%         |      |_____________  v     :  |   |    |  :
%         |                    |       :  |   |    |  :
%         |                    |       :  |   |    |  :
%         |       _____________|       v  |___|    |  :
%         |      |                            |    |  v   
%
%
% Syntax
%
% output = getfeadata_PMSM(design, Jcoil, xRVWp)
%
% Input
%
% design structure containingt the following members:
%
%   WmVWp - magnet width to pole width ratio
%
%   WtVWc - tooth width to combined tooth and slot width ratio
%
%   hmVWm - magnet height to magnet width ratio
% 
%   htVWt - tooth height to tooth width ratio
% 
%   hbaVht - armature back iron height to tooth height ratio
% 
%   hbfVhm - field back iron height to magnet height ratio
% 
%   lsVWp - stack length to pole width ratio
% 
%   gVhm - air gap to magnet height ratio
% 
%   Wp - pole width
%
%   Ntot - number of turns in the winding
%
%   fillfactor - winding fill factor
%
% Output:
%
%     output - A cell array containing the following members: 
%
%          {FEAFx, FEAFy, psi, intBx, intBy, indepvar}
%
%          Where:
%
%          FEAFx is a (k x 1) vector of force values in the x axis, i.e.
%          air-gap closing forces at each combination of position and
%          current density values specified in xRVWp and Jcoil
%
%          FEAFy is a (k x 1) vector of force values in the y axis, i.e.
%          shear forces, at each combination of position and current
%          density values specified in xRVWp and Jcoil
%
%          psi is a (k x 1) vector of flux linkage values in coil C,
%          at each combination of position and current density
%          values specified in xRVWp and Jcoil. If mode(3) == 1 this is
%          the actual flux linkage for the coil circuit as reported by
%          FEMM. Otherwise this is the integral of the vector potential
%          divided by the coil surface area to be used for finding the flux
%          linkage for any number of turns given the correctly calculated
%          current density value over the coil cross-section.
%    
%          intBx is a (k x 2) matrix of the integrals of the flux density 
%          in the coil cross-sections in the x direction. The first column
%          is the flux in the top part of the coil (with positive turn
%          direction). The second column is the flux in the bottom part
%          (with negative turn direction)
%
%          intBy is a (k x 2) matrix of the integrals of the flux density 
%          in the y direction in the coil cross-sections. The first column
%          is the flux in the top part of the coil (with positive turn
%          direction). The second column is the flux in the bottom part
%          (with negative turn direction)
%
%          indepvar, (k x 2) matrix, each row of which corresponds to the
%          pair of xRVWp anf Jcoil values which resulted in the values
%          stored in FEAFx(k,1), FEAFy(k,1), psi(k,1), intBx(k,:) and
%          intBy(k,:). This format is suitible for use with polyfitn as the
%          independent variables in a polynomial fit.



    % check the size of the current density vector
    if isvector(Jcoil)
        % if a vector is passed in, we will apply the same current density
        % to all coils, so reshape into a column vector and then
        % replicate
        Jcoil = reshape(Jcoil, [], 1);
        Jcoil = [Jcoil, Jcoil, Jcoil];
    else
        % otherwise we should have a matrix of pairs of current density
        % values, two in each row
        if size(Jcoil,2) ~= 3
            error('Jcoil, the current density input must be either a vector or an (n x 3) matrix')
        end
    end
    
%     design.CoilArea = design.ht * design.Ws;
%     if design.mode(4)
%         design.CoilArea = design.CoilArea / 2;
%     end

    % initialise some variables
    psi = zeros(size(Jcoil,1) * length(xRVWp),1);
    intBx = zeros(size(Jcoil,1) * length(xRVWp),2);
    intBy = zeros(size(Jcoil,1) * length(xRVWp),2);
    FEAFx = psi;
    FEAFy = psi;
    indepvar = zeros(size(Jcoil,1) * length(xRVWp),2);
    
    % CoreLoss(1) will contain the data for the yoke while
    % CoreLoss(2) will contain the data for the teeth
    design.CoreLoss(1).dy = design.hba / 10;
    design.CoreLoss(1).coreycoords = ((design.CoreLoss(1).dy/2):design.CoreLoss(1).dy:(design.hba-design.CoreLoss(1).dy/2))';

    design.CoreLoss(1).dx = min([0.01 / 5, design.Wp / 4]);
    design.CoreLoss(1).corexcoords = ((design.CoreLoss(1).dx/2):design.CoreLoss(1).dx:(design.Wp-design.CoreLoss(1).dx/2))';
    [yokemeshx, yokemeshy] = meshgrid(design.CoreLoss(1).coreycoords,design.CoreLoss(1).corexcoords);
    % The values along the first dimension (down the columns) of the arrays
    % will contain values sampled in the x-direction in the fea simulation.
    % In this case values along the dimension hbi, or ht in the case of the
    % teeth. The values along the second dimension are the values sampled
    % at each value xRVWp. The values along the third dimension of the
    % arrays will contain values which are sampled from the y direction of
    % the fea simulation, in this case along the dimension Wp, or Ws and Wt
    % in the case of the teeth.
    design.CoreLoss(1).Bx = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    design.CoreLoss(1).By = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    design.CoreLoss(1).Bz = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    design.CoreLoss(1).Hy = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    design.CoreLoss(1).Hx = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    design.CoreLoss(1).Hz = zeros([size(yokemeshx,2), length(xRVWp), size(yokemeshx,1)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(1).dz = design.ls;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(1).xstep = (xRVWp(2) - xRVWp(1)) * design.Wp;
            
    design.CoreLoss(2).dy = design.ht / 10;
    design.CoreLoss(2).coreycoords = (design.hba + (design.CoreLoss(1).dy/2):design.CoreLoss(1).dy:(design.hba+design.ht-design.CoreLoss(1).dy/2))';
    
    design.CoreLoss(2).dx = min([0.01 / 5, design.Wt / 4]);
    design.CoreLoss(2).corexcoords = [(design.CoreLoss(2).dx/2:design.CoreLoss(2).dx:design.Wt/2-design.CoreLoss(2).dx/2)'; ...
                                      ((design.Ws - design.Wt/2 + design.CoreLoss(2).dx/2):(design.CoreLoss(2).dx):(design.Ws + design.Wt/2-design.CoreLoss(2).dx/2))'; ...
                                      ((2*design.Ws - design.Wt/2 + design.CoreLoss(2).dx/2):(design.CoreLoss(2).dx):(2*design.Ws + design.Wt/2-design.CoreLoss(2).dx/2))'; ...
                                      ((3*design.Ws - design.Wt/2 + design.CoreLoss(2).dx/2):(design.CoreLoss(2).dx):(design.Wp-design.CoreLoss(2).dx/2))' ];
    [teethmeshx, teethmeshy] = meshgrid(design.CoreLoss(2).coreycoords,design.CoreLoss(2).corexcoords);
    
    design.CoreLoss(2).Bx = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    design.CoreLoss(2).By = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    design.CoreLoss(2).Bz = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    design.CoreLoss(2).Hy = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    design.CoreLoss(2).Hx = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    design.CoreLoss(2).Hz = zeros([size(teethmeshx,2), length(xRVWp), size(teethmeshx,1)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(2).dz = design.ls;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(2).xstep = (xRVWp(2) - xRVWp(1)) * design.Wp;
    
    design = setfieldifabsent(design, 'intBdata', struct());
    design.intBdata = setfieldifabsent(design.intBdata, 'pos', []);
    design.intBdata = setfieldifabsent(design.intBdata, 'intB1', []);
    design.intBdata = setfieldifabsent(design.intBdata, 'intB2', []);
    
    % initialise loop counter for filling indepvar
    k = 0;
    
%     % open FEMM if it's not already open
%     if ~isfemmopen
%         openfemm;
%     end
    
    % Now perform the simulations and gather the data
    for i = 1:size(Jcoil,1)
        for j = 1:length(xRVWp)

            k = k + 1;

            [design.FemmProblem] = pmsmfemmprob(design, xRVWp(j) .* design.Wp, ...
                                                'NWindingLayers', design.CoilLayers, ...
                                                'CoilCurrents', Jcoil(i,:) .* design.ConductorArea / design.fillfactor);

            % run the analysis in FEMM and load the simulation output
            femfilename = [tempname, '.fem'];
            writefemmfile(femfilename, design.FemmProblem);
            ansfilename = analyse_mfemm(femfilename);
            solution = fpproc(ansfilename);
            
            if i == 1 && j == 1
                % get the max flux density in the air gap for the
                % calculation of iron losses
                x = repmat(design.hba + design.ht + design.g - design.FEMMTol, 50, 1);
                y = linspace(0, design.Wp, 50)';
                
                p = solution.getpointvalues(x, y);
                
                BgPeak = max(abs(p(2,:)));
                
            end

            if design.mode(3) == 1

                % returned circuit properties are [I, V, flux linkage]
                temp = solution.getcircuitprops('Coil C');
                psi(k,1) = temp(3);

                % may want to extract some other properties but this will do
                % for now

            else
                % Determine the per-turn flux linkage from the magnetic
                % vector potential. This will give the correct flux linkage
                % when multiplied by the number of turns provided the
                % correct current density is used

                % coil C
                solution.clearblock();
                solution.groupselectblock(32);
                intA(1) = solution.blockintegral(1);
                solution.clearblock;
                solution.groupselectblock(31);
                intA(2) = solution.blockintegral(1);
                psi(k,1) = (intA(1) - intA(2)) / design.CoilArea;
                solution.clearblock;

            end

            % Next get the flux in the viscinity of the coil C for
            % calculating the Lorentz forces later
            solution.clearblock();
            % positive turns part
            solution.groupselectblock(32); 
            intBx(k,1) = solution.blockintegral(8);
            intBy(k,1) = solution.blockintegral(9);
            solution.clearblock();
            % negative turns part
            solution.groupselectblock(31);
            intBx(k,2) = solution.blockintegral(8);
            intBy(k,2) = solution.blockintegral(9);
            solution.clearblock();
            
            design = slotintBdata(design, solution, i, j, xRVWp);

            % get the forces, these forces are the force on one side for
            % two poles. Therefore as we want the total force between the
            % sides for one pole we leave them as they are
            solution.clearblock();
            solution.groupselectblock(3);

            if design.mode(1) == 1
                
                % mode(1) == 1 denotes a double-sided machine

                % Get the integral of the weighted maxwell stress tensor
                % over the translator. This is the per-pole force, as the
                % machine is double-sided, but the simulation consists of
                % two poles but with only one side of the machine
                FEAFx(k,1) = -solution.blockintegral(18)/2;

                % round off the force to 5 decimal places so that error in
                % the sim at the point where forces are theoretically zero
                % are in fact zero
                if roundoff(rem(xRVWp(j),1),5) == 0
                    FEAFy(k,1) = 0;
                else
                    FEAFy(k,1) = solution.blockintegral(19);
                end

            else
                % mode(3) ~= 1 denotes a single-sided machine

                % Get the integral of the weighted maxwell stress tensor
                % over the translator. This is the per-pole force, as the
                % machine is double-sided, but the simulation consists of
                % two poles but with only one side of the machine and we
                % only wish to know the air-gap closing force on one side
                FEAFx(k,1) = -solution.blockintegral(18)/2;

                % round off the force to 5 decimal places so that error in
                % the sim at the point where forces are theoretically zero
                % are in fact zero
                if roundoff(rem(xRVWp(j),1),5) == 0
                    FEAFy(k,1) = 0;
                else
                    FEAFy(k,1) = solution.blockintegral(19)/2;
                end

            end

            % Store the xRVWp and Jcoil pair which produced the data in
            % this pass of the loop in indepvar
            indepvar(k,:) = [xRVWp(j)-0.5, Jcoil(i)];
            
            % core material, note that the directions of the field are
            % changed when stored, as the FEA sim y-direction is not the
            % direction of motion of the translator.
            
            p = solution.getpointvalues(yokemeshx(:), yokemeshy(:));

            design.CoreLoss(1).By(:,j,:) = reshape(p(2,:)', size(yokemeshx))';
            design.CoreLoss(1).Bx(:,j,:) = reshape(p(3,:)', size(yokemeshx))';
            design.CoreLoss(1).Hy(:,j,:) = reshape(p(6,:)', size(yokemeshx))';
            design.CoreLoss(1).Hx(:,j,:) = reshape(p(7,:)', size(yokemeshx))';
            
            p = solution.getpointvalues(teethmeshx(:), teethmeshy(:));

            design.CoreLoss(2).By(:,j,:) = reshape(p(2,:)', size(teethmeshx))';
            design.CoreLoss(2).Bx(:,j,:) = reshape(p(3,:)', size(teethmeshx))';
            design.CoreLoss(2).Hy(:,j,:) = reshape(p(6,:)', size(teethmeshx))';
            design.CoreLoss(2).Hx(:,j,:) = reshape(p(7,:)', size(teethmeshx))';

            % delete the ans and .fem file
            delete(femfilename);
            delete(ansfilename);

        end
    end
    
    %closefemm;
    
    [design.intBdata.pos, ix] = sort(design.intBdata.pos);
    design.intBdata.intB1 = design.intBdata.intB1(ix,:);
    design.intBdata.intB2 = design.intBdata.intB2(ix,:);
    
    output = {FEAFx, FEAFy, psi, intBx, intBy, indepvar, BgPeak};

end


function design = slotintBdata(design, solution, Jind, xRVWpind, xRVWp)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 2
        
        slotypos = linspace(design.Wt/2 + design.Wc/2, 2*design.Wp - (design.Wt/2 + design.Wc/2), 6);

        sloty1pos = [slotypos(2:6), slotypos(1)];
        
        sloty2pos = [slotypos(5:6), slotypos(1:4)];

        slotxpos = [design.hba + design.ht/4, design.hba + 3*design.ht/4 ];

        dxR = 0;
        for i = 1:numel(slotypos)

            design.intBdata.pos(end+1) = xRVWp(xRVWpind) + dxR;

            solution.clearblock();
            solution.selectblock(slotxpos(1), sloty1pos(i));
            design.intBdata.intB1(end+1,1) = solution.blockintegral(8);
            design.intBdata.intB1(end,2) = solution.blockintegral(9);

            solution.clearblock();
            solution.selectblock(slotxpos(2), sloty2pos(i));
            design.intBdata.intB2(end+1,1) = solution.blockintegral(8);
            design.intBdata.intB2(end,2) = solution.blockintegral(9);
            
            dxR = dxR - design.Ws/design.Wp;
            
        end
    
    end

end
