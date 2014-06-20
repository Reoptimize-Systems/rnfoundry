function [design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions)
% generates simulation data for a radial flux slotted machine in
% preparation for a time series ode simulation
%
% 
%
%

% Stator slot/coil/winding terminology is based on that presented in:
%
% J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
% Single layer
%   q = 2qc
%   Qs = 2Qc
% Double layer
%   q = qc
%   Qs = Qc
%
% yp - Average coil pitch as defined by (Qs/Poles)
% yd - Actual coil pitch as defined by round(yp) +/- k
% Qs  -  number of stator slots
% Qc  -  number of winding coils
% q  -  number of slots per pole and phase, (1)
% qn  -  numerator of q
% qd  -  denominator of q
% qc - number of coils per pole and phase, (2)
% qcn  -  numerator of qc
% qcd  -  denominator of qc

%     design.Hc = design.tc;
%     
%     % \Tau_{cs} is the thickness of the winding, i.e. the pitch of a
%     % winding slot
%     design.Wc = design.thetac * design.Rcm;
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    if ~isfield(design, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
    end
    
    % placeholder for coil cross-sectional area, this is extracted from the
    % FEA below, this avoids checkcoilprops_AM changing it
    design.CoilArea = nan;
    
    % call the common radial simulation function
    [design, simoptions] = simfun_RADIAL(design, simoptions);
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20));
        else
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        end
    else
        if design.tsb > 1e-5
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        else
            simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
        end
    end
    
    % now get the flux linkage in the coils, we will do this by getting the
    % integral of the vector potential in a slot over a half pole of
    % displacement, giving a quarter of the wave period
    
    nfeapos = 10;
%     design.FirstSlotCenter = design.thetap + (design.thetap / slotsperpole / 2);
    design.FirstSlotCenter = 0;
    design.feapos = linspace(0, 1, nfeapos);

    design.intAdata.slotPos = [];
    design.intAdata.slotIntA = [];
    design.intBdata.slotPos = [];
    design.intBdata.slotIntB = [];
    
    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];
    
    % determine a unit vector pointing in the direction normal to the air
    % gap half way between the Poles for the purpose of extracting the
    % air-gap closing forces
    [gvector(1), gvector(2)] = pol2cart(design.thetap,1);
    gvector = unit(gvector);
    
    % determine a unit vector pointing in the direction tangential to the
    % radius half way between the Poles for the purpose of extracting
    % the cogging forces
    %
    % The dot product of orthogonal vectors is zero. In two dimensions the
    % slopes of perpendicular lines are negative reciprocals. Switch the x
    % and y coefficients and change the sign of one of them.
    % 
    % For vector V1 = <a, b>, the reciprocal V2 = <b, -a>. Now make it a
    % unit vector by dividing by the magnitude.
    coggingvector = unit([gvector(2), -gvector(1)]);
    
    armirongroup = 2;
        
    for i = 1:numel(design.feapos)

        % Draw the sim, i.e by creating the FemmProblem structure
        [design.FemmProblem, design.coillabellocs, design.yokenodeids] = ...
                            slottedfemmprob_radial(design, ...
                                'ArmatureType', design.ArmatureType, ...
                                'NWindingLayers', design.CoilLayers, ...
                                'Position', (design.feapos(i)+1) * design.thetap + design.FirstSlotCenter, ...
                                'ArmatureBackIronGroup', armirongroup, ...
                                'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
                                'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
                                'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
                                'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
                                'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
                                'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
                                'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);

        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);
        % analyse the problem
        ansfilename = analyse_mfemm(femfilename, ...
                                    simoptions.usefemm, ...
                                    simoptions.quietfemm);
                            
        
        if (exist('fpproc_interface_mex', 'file')==3) && ~simoptions.usefemm
            
            solution = fpproc(ansfilename);
            % activate field smoothing so values are interpolated across
            % mesh elements
            solution.smoothon();
            
            if i == 1
                design = corelosssetup(design, design.feapos, solution);
            end
            
            % get the integral of the vector potential in a slot, if two
            % layers get both layers
            design = slotintAdata(design, simoptions, design.feapos(i), solution);
            
            design = slotintBdata(design, simoptions, design.feapos(i), solution);

            for ii = 1:numel(design.CoreLoss)
                
                p = solution.getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                             design.CoreLoss(ii).meshy(:) );

                design.CoreLoss(ii).By(:,i,:) = reshape(p(2,:)', size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Bx(:,i,:) = reshape(p(3,:)', size(design.CoreLoss(ii).meshx));

            end
            
            if i == 1
                % get the forces
                solution.clearblock();
                solution.groupselectblock(1)
                
                design.gforce = dot([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
                                     gvector);
                                
                design.gvar = design.g;
                
                % get the cross-sectional area of the armature iron for
                % calcuation of material masses later
                solution.clearblock();
                solution.groupselectblock(armirongroup);
                design.ArmatureIronAreaPerPole = solution.blockintegral(5)/2;
                
                % get the cross-sectional area of the coil winding bundle
                design.CoilArea = solution.blockintegral ( 5, design.coillabellocs(1,1), design.coillabellocs(1,2) );
                
            end
            
            % get the cogging forces
            solution.clearblock();
            solution.groupselectblock(1)
            
            design.coggingforce(i) = dot([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
                                         coggingvector);

            design.coggingforce(i) = design.coggingforce(i) * design.Poles(1);
            
            % explicitly call the delete method on the solution
            delete(solution);
            clear solution;
            
        else
            if i == 1
                design = corelosssetup(design, design.feapos, ansfilename);
            end
            
            % open the solution in FEMM
            opendocument(ansfilename);

            % get the integral of the vector potential in a slot, if two
            % layers get both layers
            design = slotintAdata(design, simoptions, design.feapos(i));

            design = slotintBdata(design, simoptions, design.feapos(i));

            for ii = 1:numel(design.CoreLoss)
                
                p = mo_getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                       design.CoreLoss(ii).meshy(:) );

                design.CoreLoss(ii).By(:,i,:) = reshape(p(:,2), size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Bx(:,i,:) = reshape(p(:,3), size(design.CoreLoss(ii).meshx));

            end
            
            if i == 1
                % get the forces
                mo_clearblock();
                mo_groupselectblock(1);
                design.gforce = dot([mo_blockintegral(18)/2, mo_blockintegral(19)/2], ...
                                    gvector);
                design.gvar = design.g;
                
                % get the cross-sectional area of the armature iron for
                % calcuation of material masses later
                mo_clearblock();
                mo_groupselectblock(armirongroup);
                design.ArmatureIronAreaPerPole = mo_blockintegral(5)/2;
                
            end

            mo_close;
            
        end
        
        % clean up the FEA files
        delete(femfilename);
        delete(ansfilename);
    
    end
    
    % now get more force points
    pos = linspace(design.FEMMTol - design.g, 0, simoptions.NForcePoints-1);
    pos(end) = 2*design.g;
    
    if simoptions.GetVariableGapForce
        design.gforce = [design.gforce, rotorforces_RADIAL_SLOTTED(design.FemmProblem, 2, pos)];
    else
        design.gforce = [design.gforce, repmat(design.gforce, 1, numel(pos))];
    end
    design.gvar = [design.gvar, design.g + pos];
    
    % perform an inductance sim
    Lcurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
    [InductanceFemmProblem, Lcoillabellocs] = slottedLfemmprob_radial(design, ...
        'ArmatureType', design.ArmatureType, ...
        'NWindingLayers', design.CoilLayers, ...
        'CoilCurrent', Lcurrent, ...
        'MagnetRegionMeshSize', simoptions.femmmeshoptions.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
        'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
        'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
        'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
        'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize);
    
    % write the fem file to disk
    writefemmfile(femfilename, InductanceFemmProblem);
    % analyse the problem
    ansfilename = analyse_mfemm(femfilename, ...
                                simoptions.usefemm, ...
                                simoptions.quietfemm);
    
    if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
        solution = fpproc(ansfilename);
        [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
        % get the mutual inductance by dividing the flux in a neighbouring
        % phase by the applied current in the first phase
        temp = solution.getcircuitprops('2');
        % explicitly call the delete method on the solution
        delete(solution);
        clear solution;
    else
        opendocument(ansfilename);
        % Now get the resistance and inductance of the machine coils
        [design.CoilResistance, design.CoilInductance] = RandLfromFEMMcircuit('1');
        % get the mutual inductance by dividing the flux in a neighbouring
        % phase by the applied current in the first phase
        temp = mo_getcircuitproperties('2');
        mo_close;
    end
    design.CoilInductance(2) = abs(temp(3)) / Lcurrent;
    
    % clean up the FEA files
    delete(femfilename);
    delete(ansfilename);
    
    % now recalculate coil resistance
    design.MTL = rectcoilmtl( design.ls, ...
                              design.yd * design.thetas * design.Rcm, ...
                              design.thetac * design.Rcm );

    design.CoilResistance = 1.7e-8 * design.CoilTurns * design.MTL ./ (pi * (design.Dc/2)^2);
    
end


function design = slotintBdata(design, simoptions, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)),2);
        slotxpos = design.coillabellocs(1:(size(design.coillabellocs,1)),1);
                
        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % x and y directed flux in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
            else
                % x and y directed flux in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),2), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,2)];
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,1)];
              
        for i = 1:size(slotypos,1)

            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i,1));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i,2));
                design.intBdata.slotIntB(end,3,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = solution.blockintegral(9);
                
            else
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i,1));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i,2));
                design.intBdata.slotIntB(end,3,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = mo_blockintegral(9);
            end
        end

    end
    
    % store the relative coil/slot positions, we use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol(slotxpos(:,1), slotypos(:,1));
    design.intBdata.slotPos = [design.intBdata.slotPos; 
                               (-thetaslotypos./design.thetap) + design.FirstSlotCenter + feapos];
    
    % sort the data in ascending position order
    [design.intBdata.slotPos, idx] = sort(design.intBdata.slotPos);
    design.intBdata.slotIntB = design.intBdata.slotIntB(idx,:,:);
    
end


function design = slotintAdata(design, simoptions, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)),2);
        slotxpos = design.coillabellocs(1:(size(design.coillabellocs,1)),1);

        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % vector potential in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
            else
                % vector potential in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),2), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,2)];
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)))+1,1)];
              
        for i = 1:size(slotypos,1)

            if exist('fpproc_interface_mex', 'file') == 3 && ~simoptions.usefemm
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i,1));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
                
                % vector potential flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i,2));
                design.intAdata.slotIntA(end,2,1) = solution.blockintegral(1);
                
            else
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i,1));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
                
                % vector potential in slot on left hand side inner
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i,2));
                design.intAdata.slotIntA(end,2,1) = mo_blockintegral(1);
                
            end
        end

    end
    
    % store the relative coil/slot positions. We use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol(slotxpos(:,1), slotypos(:,1));
    design.intAdata.slotPos = [design.intAdata.slotPos; 
                              (-thetaslotypos./design.thetap) + design.FirstSlotCenter + feapos];
    
    % sort the data in ascending position order
    [design.intAdata.slotPos, idx] = sort(design.intAdata.slotPos);
    design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);
    
end



function design = corelosssetup(design, feapos, solution)
    
    % get the number of positions
    npos = numel(feapos);
    
    if isa(solution, 'fpproc')
        % it's an xfemm fpproc object

        % get the volume of each element under consideration
        design.CoreLoss(1).dV = solution.getgroupareas (2) .* design.ls;
        
        % get the location we will use to estimate the flux in the element
        temp = solution.getgroupcentroids (2);
        design.CoreLoss(1).meshx = temp(:,1);
        design.CoreLoss(1).meshy = temp(:,2);
        
        % The values along the first dimension (down the columns) of the
        % arrays will contain values sampled in the x-direction in the fea
        % simulation. In this case values along the dimension hbi, or ht in
        % the case of the teeth. The values along the second dimension are
        % the values sampled at each value xRVWp. The values along the
        % third dimension of the arrays will contain values which are
        % sampled from the y direction of the fea simulation, in this case
        % along the dimension taupm, or tsb and tc in the case of the
        % teeth.
        design.CoreLoss(1).Bx = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).By = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).Bz = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);

        % add xstep which is the size of the steps in xRVWp denormalised. This
        % is later used to find the value of dB/dx at each position
        design.CoreLoss(1).xstep = (feapos(2) - feapos(1)) * design.thetap;
        
    else
        error ('FEMM not currently supported for this, use xfemm.')
    end

end

