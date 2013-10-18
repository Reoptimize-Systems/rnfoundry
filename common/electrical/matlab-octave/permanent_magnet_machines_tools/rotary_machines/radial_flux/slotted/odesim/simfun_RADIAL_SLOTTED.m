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

    design.Hc = design.tc;
    
    % \Tau_{cs} is the thickness of the winding, i.e. the pitch of a
    % winding slot
    design.Wc = design.thetac * design.Rcm;
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    if ~isfield(design, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
     
        design.CoreLoss(2) = design.CoreLoss(1);
    end
    
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
    
    femfilename = [tempname, '_simfun_TORUS_SLOTTED.fem'];
    
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
        [design.FemmProblem, design.outermagsep, design.coillabellocs, design.yokenodeids] = ...
                            slottedfemmprob_radial(design, ...
                                'StatorType', design.StatorType, ...
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

        if i == 1
            design = corelosssetup(design, design.feapos);
        end
        
        % write the fem file to disk
        writefemmfile(femfilename, design.FemmProblem);
        % analyse the problem
        ansfilename = analyse_mfemm(femfilename, ...
                                    simoptions.usefemm, ...
                                    simoptions.quietfemm);
        
        if (exist('fpproc_interface_mex', 'file')==3) && ~simoptions.usefemm
            
            solution = fpproc(ansfilename);
            solution.smoothon();
            
            % get the integral of the vector potential in a slot, if two
            % layers get both layers
            design = slotintAdata(design, simoptions, design.feapos(i), solution);
            
            design = slotintBdata(design, simoptions, design.feapos(i), solution);

            for ii = 1:numel(design.CoreLoss)
                
                p = solution.getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                             design.CoreLoss(ii).meshy(:) );

                design.CoreLoss(ii).By(:,i,:) = reshape(p(2,:)', size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Bx(:,i,:) = reshape(p(3,:)', size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Hy(:,i,:) = reshape(p(6,:)', size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Hx(:,i,:) = reshape(p(7,:)', size(design.CoreLoss(ii).meshx));

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
                design.CoreLoss(ii).Hy(:,i,:) = reshape(p(:,6), size(design.CoreLoss(ii).meshx));
                design.CoreLoss(ii).Hx(:,i,:) = reshape(p(:,7), size(design.CoreLoss(ii).meshx));

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
        'StatorType', design.StatorType, ...
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



function design = corelosssetup(design, feapos)

    % get the number of positions
    npos = numel(feapos);

    % obtain the coordinates of the four corners of the first stator
    coords = getnodecoords_mfemm(design.FemmProblem);
    coords = coords(design.yokenodeids(1,1:4)+1,:);
    
    % CoreLoss(1) will contain the data for the yoke while
    % CoreLoss(2) will contain the data for the teeth central section
    % without the shoes
    % CoreLoss(3) will contain the data for the teeth shoes
    
    % yoke
    design.CoreLoss(1).dy = design.ty / 10;
    design.CoreLoss(1).coreycoords = ((design.CoreLoss(1).dy/2):design.CoreLoss(1).dy:(design.ty-design.CoreLoss(1).dy/2))';
    
    if strcmp(design.StatorType, 'si')
        design.CoreLoss(1).coreycoords = ...
            design.CoreLoss(1).coreycoords + coords(2,1) + design.ty;
    elseif strcmp(design.StatorType, 'so')
        design.CoreLoss(1).coreycoords = ...
            design.CoreLoss(1).coreycoords + coords(1,1) - design.ty;
    end

    design.CoreLoss(1).dx = min([0.01 / 5, design.thetap / 10]);
    design.CoreLoss(1).corexcoords = ((design.CoreLoss(1).dx/2):design.CoreLoss(1).dx:(design.thetap-design.CoreLoss(1).dx/2))';
    [design.CoreLoss(1).meshx, design.CoreLoss(1).meshy] = meshgrid(design.CoreLoss(1).coreycoords,design.CoreLoss(1).corexcoords);
    Ri = design.CoreLoss(1).meshx - design.CoreLoss(1).dx / 2;
    Ro = design.CoreLoss(1).meshx + design.CoreLoss(1).dx / 2;
    theta = design.CoreLoss(1).dy;
    design.CoreLoss(1).dA = annularsecarea(Ri, Ro, theta);
    design.CoreLoss(1).dA = repmat(reshape(design.CoreLoss(1).dA, size(design.CoreLoss(1).dA,1), 1, size(design.CoreLoss(1).dA,2)), [1, npos, 1]);
    
    % convert mesh coordinates from polar to cart
    [design.CoreLoss(1).meshx, design.CoreLoss(1).meshy] = pol2cart( design.CoreLoss(1).meshy, design.CoreLoss(1).meshx );
    
    % The values along the first dimension (down the columns) of the arrays
    % will contain values sampled in the x-direction in the fea simulation.
    % In this case values along the dimension hbi, or ht in the case of the
    % teeth. The values along the second dimension are the values sampled
    % at each value xRVWp. The values along the third dimension of the
    % arrays will contain values which are sampled from the y direction of
    % the fea simulation, in this case along the dimension taupm, or tsb
    % and tc in the case of the teeth.
    design.CoreLoss(1).Bx = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    design.CoreLoss(1).By = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    design.CoreLoss(1).Bz = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    design.CoreLoss(1).Hy = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    design.CoreLoss(1).Hx = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    design.CoreLoss(1).Hz = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(1).dz = design.ls;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(1).xstep = (feapos(2) - feapos(1)) * design.thetap;
            
    % now set up the tooth body sampling positions
    design.CoreLoss(2).dy = design.tc / 10;
    design.CoreLoss(2).coreycoords =  ...
        ((design.CoreLoss(2).dy/2):design.CoreLoss(2).dy:(design.tc-design.CoreLoss(1).dy/2))';
    design.CoreLoss(2).coreycoords = design.CoreLoss(2).coreycoords + coords(1,1);
    
    design.thetat = design.thetas-design.thetac;
    design.CoreLoss(2).dx = min([0.01 / 5, design.thetat / 4]);
    % these will be the angular positions of the points
    halftoothcorex = ((design.CoreLoss(2).dx/2):design.CoreLoss(2).dx:(design.thetat/2-design.CoreLoss(2).dx/2))';
    fulltoothcorex = ((design.CoreLoss(2).dx/2):design.CoreLoss(2).dx:(design.thetat-design.CoreLoss(2).dx/2))';
    if design.ypd == 1
        
        design.CoreLoss(2).corexcoords = [halftoothcorex; ...
                                          zeros(numel(fulltoothcorex)*(design.ypn-1),1); ...
                                          halftoothcorex + design.thetap - (design.thetat/2) ];

        for i = 1:(design.ypn-1)
            design.CoreLoss(2).corexcoords( numel(halftoothcorex) ...
                                            + (1+((i-1)*numel(fulltoothcorex)):(i*numel(fulltoothcorex))),1 ) ...
                               = fulltoothcorex + i*design.thetas - design.thetat/2;
        end
                                  
    elseif design.ypd == 2
        
        design.CoreLoss(2).corexcoords = [halftoothcorex; ...
                                          zeros(numel(fulltoothcorex)*floor(design.ypn/2),1)];

        for i = 1:(floor(design.ypn/2))
            design.CoreLoss(2).corexcoords((i*numel(fulltoothcorex)):((i+1)*numel(fulltoothcorex)-1),1) = ...
                fulltoothcorex + i*design.thetas - design.thetat/2;
        end
        
    else
        error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end
    
    design.CoreLoss(2).corexcoords = [ design.CoreLoss(2).corexcoords; ];
                                   
    design.CoreLoss(2).coreycoords = [ design.CoreLoss(2).coreycoords; ];
    
    [design.CoreLoss(2).meshx, design.CoreLoss(2).meshy] = ... 
        meshgrid(design.CoreLoss(2).coreycoords,design.CoreLoss(2).corexcoords);
    
    % calculate the element areas
    Ri = design.CoreLoss(2).meshx - (design.CoreLoss(2).dy / 2);
    Ro = design.CoreLoss(2).meshx + (design.CoreLoss(2).dy / 2);
    theta = design.CoreLoss(2).dx;
    design.CoreLoss(2).dA = annularsecarea(Ri, Ro, theta);
    design.CoreLoss(2).dA = repmat(reshape(design.CoreLoss(2).dA, size(design.CoreLoss(2).dA,1), 1, size(design.CoreLoss(2).dA,2)), [1, npos, 1]);
    
    [design.CoreLoss(2).meshx, design.CoreLoss(2).meshy] = ...
        pol2cart(design.CoreLoss(2).meshy, design.CoreLoss(2).meshx);
    
    design.CoreLoss(2).Bx = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    design.CoreLoss(2).By = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    design.CoreLoss(2).Bz = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    design.CoreLoss(2).Hy = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    design.CoreLoss(2).Hx = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    design.CoreLoss(2).Hz = zeros([size(design.CoreLoss(2).meshx,1), npos, size(design.CoreLoss(2).meshx,2)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(2).dz = design.CoreLoss(1).dz;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(2).xstep = design.CoreLoss(1).xstep;
    
    % shoes
    if (design.thetac - design.thetasg) > (design.Rmm*design.FEMMTol)
        
        shoedx = min([0.01 / 5, (design.thetac - design.thetasg)/2 / 10]);
        
        shoedy = design.tsb / 10;
        
        design.thetacss = (design.thetac - design.thetasg) / 2;
        
        shoeangle = atan((design.tsb - design.tsg) / design.thetacss);
        
        xcoords = ((shoedx/2):shoedx:(design.thetacss-design.CoreLoss(1).dx/2));
        
        for i = 1:numel(xcoords)
            
            % copy over the information about the materials etc to the new
            % structure
            design.CoreLoss(2+i).fq = design.CoreLoss(2).fq;
            design.CoreLoss(2+i).Bq = design.CoreLoss(2).Bq;
            design.CoreLoss(2+i).Pq = design.CoreLoss(2).Pq;
            design.CoreLoss(2+i).kh = design.CoreLoss(2).kh;
            design.CoreLoss(2+i).kc = design.CoreLoss(2).kc;
            design.CoreLoss(2+i).ke = design.CoreLoss(2).ke;
            design.CoreLoss(2+i).beta = design.CoreLoss(2).beta;
            design.CoreLoss(2+i).dz = design.CoreLoss(1).dz;
            design.CoreLoss(2+i).xstep = design.CoreLoss(1).xstep;
            
            % choose y coords in a line across the shoe at the current x
            % position
            shoewid = design.tsg + xcoords(i)*tan(shoeangle);
            design.CoreLoss(2+i).dy = shoedy;
            shoeycoords = ((design.CoreLoss(2+i).dy/2):design.CoreLoss(2+i).dy:(shoewid-design.CoreLoss(2+i).dy/2))';
            if isempty(shoeycoords)
                if strcmp(design.StatorType, 'si')
                    design.CoreLoss(2+i).coreycoords = coords(1,1) + shoewid/2;
                elseif strcmp(design.StatorType, 'so')
                    design.CoreLoss(2+i).coreycoords = coords(2,1) - shoewid/2;
                end
            else
                if strcmp(design.StatorType, 'si')
                    design.CoreLoss(2+i).coreycoords = shoeycoords + coords(1,1);
                elseif strcmp(design.StatorType, 'so')
                    design.CoreLoss(2+i).coreycoords = coords(2,1) - fliplr(shoeycoords);
                end
            end
            
            % add a single x coordinate for the tooth shoe at the bottom of
            % the sim
            design.CoreLoss(2+i).dx = shoedx;
            design.CoreLoss(2+i).corexcoords = design.thetat/2 +  design.thetacss - xcoords(i);
            
            if design.ypd == 1
                
                for ii = 1:(design.ypn-1)
                    
                    design.CoreLoss(2+i).corexcoords = ...
                        [design.CoreLoss(2+i).corexcoords;
                         ii*design.thetas - design.thetat/2 - (design.thetacss) + xcoords(i); 
                         ii*design.thetas + design.thetat/2 + (design.thetacss) - xcoords(i) ];
                     
                end
                
                % add a single x coordinate for the tooth shoe at the top
                % of the pole
                design.CoreLoss(2+i).dx = shoedx;
                design.CoreLoss(2+i).corexcoords = [design.CoreLoss(2+i).corexcoords;
                                                    design.thetap - (design.thetat/2 +  design.thetacss - xcoords(i)) ];
                
            elseif design.ypd == 2
                
                for ii = 1:(floor(design.ypn/2))
                    
                    design.CoreLoss(2+i).corexcoords = ...
                        [design.CoreLoss(2+i).corexcoords;
                         ii*design.thetas - design.thetat/2 - (design.thetacss) + xcoords(i);
                         ii*design.thetas + design.thetat/2 + (design.thetacss) - xcoords(i) ];
                     
                end
                
            else
                error('denominator of slots per pole must be 1 or 2, other values not yet supported')
            end
            
            % generate the sample points
            [design.CoreLoss(2+i).meshx, design.CoreLoss(2+i).meshy] = ...
                meshgrid(design.CoreLoss(2+i).coreycoords,design.CoreLoss(2+i).corexcoords);

            % Calculate the element areas
            Ri = design.CoreLoss(2+i).meshx - design.CoreLoss(2+i).dx / 2;
            Ro = design.CoreLoss(2+i).meshx + design.CoreLoss(2+i).dx / 2;
            theta = design.CoreLoss(2+i).dy;
            design.CoreLoss(2+i).dA = annularsecarea(Ri, Ro, theta);
            design.CoreLoss(2+i).dA = repmat(reshape(design.CoreLoss(2+i).dA, size(design.CoreLoss(2+i).dA,1), 1, size(design.CoreLoss(2+i).dA,2)), [1, npos, 1]);

            [design.CoreLoss(2+i).meshx, design.CoreLoss(2+i).meshy] = ...
                pol2cart(design.CoreLoss(2+i).meshy, design.CoreLoss(2+i).meshx);
            
            design.CoreLoss(2+i).Bx = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            design.CoreLoss(2+i).By = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            design.CoreLoss(2+i).Bz = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            design.CoreLoss(2+i).Hy = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            design.CoreLoss(2+i).Hx = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            design.CoreLoss(2+i).Hz = zeros([size(design.CoreLoss(2+i).meshx,1), npos, size(design.CoreLoss(2+i).meshy,2)]);
            
        end
    end

end

