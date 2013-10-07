function [design, simoptions] = simfun_TORUS_SLOTTED(design, simoptions)
% generates simulation data for a slotted torus type machine in preparation
% for a time series ode simulation
%
% 
%
%

% Stator slot/coil/winding terminology is based on that presented in:
%
% J. J. Germishuizen and M. J. Kamper, “Classification of symmetrical
% non-overlapping three-phase windings,” in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.

    design.Hc = design.tc;
    
    % Single layer
    %   q = 2qc
    %   Qs = 2Qc
    % Double layer
    %   q = qc
    %   Qs = Qc
    %
    % yp - Average coil pitch as defined by (Qs/Poles)
    % yd - Actual coil pitch as defined by round(yp) +/- k
    % Qs – number of stator slots
    % Qc – number of winding coils
    % q – number of slots per pole and phase, (1)
    % qn – numerator of q
    % qd – denominator of q
    % qc – number of coils per pole and phase, (2)
    % qcn – numerator of qc
    % qcd – denominator of qc
    
    % number of slots per pole and phase
    if ~isfield(design, 'qc')
        design.qc = fr(design.Qs, design.Phases * design.Poles);
    else
        [design.Qs,junk] = rat(design.qc * design.Phases * design.Poles);
    end
    
    % number of slots per pole
    slotsperpole = design.Qs / design.Poles;
    
    % get the total slot width
    design.tausm = design.taupm / slotsperpole;
    
    % get the numerator and denominator of qc
    [design.qcn,design.qcd] = rat(design.qc);
    
    % Average coil pitch as defined by (Qs/Poles)
    design.yp = fr(design.Qs, design.Poles);
    
    % get the numerator and denominator of the coil pitch in slots
    [design.ypn,design.ypd] = rat(design.yp);
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    % \Tau_{cs} is the thickness of the winding, i.e. the pitch of a
    % winding slot
    design.Wc = design.taucs;
    
    if ~isfield(design, 'CoreLoss')
        % Coreloss will be the armature back iron data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
     
        design.CoreLoss(2) = design.CoreLoss(1);
    end
    
    [design, simoptions] = simfun_TORUS(design, simoptions);
    
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'ShoeGapRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'YokeRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'CoilRegionMeshSize', -1);
    
    % now get the flux linkage in the coils, we will do this by getting the
    % integral of the vector potential in a slot over a half pole of
    % displacement, giving a quarter of the wave period
    
    nfeapos = 10;
    design.FirstSlotCenter = design.taupm / slotsperpole / 2;
    design.feapos = linspace(0, 1, nfeapos);
    
%     design = corelosssetup(design, feapos);
    
    design.intAdata.slotPos = [];
    design.intAdata.slotIntA = [];
    design.intBdata.slotPos = [];
    design.intBdata.slotIntB = [];
    
    femfilename = 'temp_simfun_TORUS_SLOTTED.fem';
        
    for i = 1:numel(design.feapos)
        
        % Draw the sim 
        [design.FemmProblem, design.outermagsep, design.coillabellocs, design.yokenodeids] = ...
                            slottedfemmprob_torus(design, ...
                                'NStages', simoptions.ndrawnstages, ...
                                'NWindingLayers', design.CoilLayers, ...
                                'Position', design.feapos(i) * design.taupm + design.FirstSlotCenter, ...
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
        ansfilename = analyse_mfemm(femfilename);
        
        if (exist('fpproc_interface_mex', 'file')==3)
            
            solution = fpproc(ansfilename);
            
            % get the integral of the vector potential in a slot, if two
            % layers get both layers
            design = slotintAdata(design, design.feapos(i), solution);
            
            design = slotintBdata(design, design.feapos(i), solution);

            for ii = 1:numel(design.CoreLoss)
                
                p = solution.getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                             design.CoreLoss(ii).meshy(:) );

                design.CoreLoss(ii).By(:,i,:) = reshape(p(2,:)', size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Bx(:,i,:) = reshape(p(3,:)', size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Hy(:,i,:) = reshape(p(6,:)', size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Hx(:,i,:) = reshape(p(7,:)', size(design.CoreLoss(ii).meshx))';

            end
            
            if i == 1
                % get the forces
                solution.clearblock();
                solution.groupselectblock(simoptions.ndrawnstages+1)
                design.gforce = solution.blockintegral(18)/2;
                design.gvar = design.g;
            end
            
        else
            % open the solution in FEMM
            opendocument(ansfilename);

            % get the integral of the vector potential in a slot, if two
            % layers get both layers
            design = slotintAdata(design, design.feapos(i));

            design = slotintBdata(design, design.feapos(i));

            for ii = 1:numel(design.CoreLoss)
                
                p = mo_getpointvalues( design.CoreLoss(ii).meshx(:), ...
                                       design.CoreLoss(ii).meshy(:) );

                design.CoreLoss(ii).By(:,i,:) = reshape(p(:,2), size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Bx(:,i,:) = reshape(p(:,3), size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Hy(:,i,:) = reshape(p(:,6), size(design.CoreLoss(ii).meshx))';
                design.CoreLoss(ii).Hx(:,i,:) = reshape(p(:,7), size(design.CoreLoss(ii).meshx))';

            end
            
            if i == 1
                % get the forces
                mo_clearblock();
                mo_groupselectblock(simoptions.ndrawnstages+1)
                design.gforce = mo_blockintegral(18)/2;
                design.gvar = design.g;
            end

            mi_close;
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
        design.gforce = [design.gforce, rotorforces_TORUS_SLOTTED(design.FemmProblem, simoptions.ndrawnstages, 2, pos)];
    else
        design.gforce = [design.gforce, repmat(design.gforce, 1, numel(pos))];
    end
    design.gvar = [design.gvar, design.g + pos];
    
    % perform an inductance sim
    Lcurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
    [InductanceFemmProblem, Loutermagsep, Lcoillabellocs] = slottedLfemmprob_torus(design, ...
        'NStages', design.NStages, ...
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
    ansfilename = analyse_mfemm(femfilename);
    
    if exist('fpproc_interface_mex', 'file') == 3
        solution = fpproc(ansfilename);
        [design.CoilResistance, design.CoilInductance] = solution.circuitRL('1');
        % get the mutual inductance by dividing the flux in a neighbouring
        % phase by the applied current in the first phase
        temp = solution.getcircuitprops('2');
    else
        opendocument(ansfilename);
        % Now get the resistance and inductance of the machine coils
        [design.CoilResistance, design.CoilInductance] = RandLfromFEMMcircuit('1');
        % get the mutual inductance by dividing the flux in a neighbouring
        % phase by the applied current in the first phase
        temp = mo_getcircuitproperties('2');
        mo_close;
        mi_close;
    end
    design.CoilInductance(2) = abs(temp(3)) / Lcurrent;
    
    % clean up the FEA files
    delete(femfilename);
    delete(ansfilename);
    
    % now recalculate coil resistance
    design.tauci = (design.yd * design.tausm) - design.taucs;
    design.Rmm = (design.Rmo + design.Rmi) / 2;
    theta = design.tauci / design.Rmm;
    design.taucio = theta * (design.Rmo - design.Wc);
    design.taucii = theta * (design.Rmi + design.Wc);
    
    design.MTL = isotrapzcoilmtl(design.taucio, ...
                                 design.taucii, ...
                                 design.Rmo + design.Rmi, ...
                                 design.Wc);

    design.CoilResistance = 1.7e-8 * design.CoilTurns * design.MTL ./ (pi * (design.Dc/2)^2);
    
end


function design = slotintBdata(design, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)/2),2);
        slotxpos = [design.coillabellocs(1:(size(design.coillabellocs,1)/2),1), ...
                    design.coillabellocs((1:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2,1) ];
                
        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3
                % x and y directed flux in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
                
                % x and y directed flux in slot on right hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i));
                design.intBdata.slotIntB(end,1,2) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,2) = solution.blockintegral(9);
            else
                % x and y directed flux in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
                
                % x and y directed flux in slot on right hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i));
                design.intBdata.slotIntB(end,1,2) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,2) = mo_blockintegral(9);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = design.coillabellocs(1:2:(size(design.coillabellocs,1)/2),2);
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)/2),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2))+1,1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2+1,1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2,1)];
              
        for i = 1:numel(slotypos)

            if exist('fpproc_interface_mex', 'file') == 3
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = solution.blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i));
                design.intBdata.slotIntB(end,3,1) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = solution.blockintegral(9);
                
                % x and y directed flux in slot on right hand side inner
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,3), slotypos(i));
                design.intBdata.slotIntB(end,1,2) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,2,2) = solution.blockintegral(9);
                
                % x and y directed flux in slot on right hand side outer
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,4), slotypos(i));
                design.intBdata.slotIntB(end,3,2) = solution.blockintegral(8);
                design.intBdata.slotIntB(end,4,2) = solution.blockintegral(9);
                
            else
                % x and y directed flux in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i));
                design.intBdata.slotIntB(end+1,1,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,1) = mo_blockintegral(9);
                
                % x and y directed flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i));
                design.intBdata.slotIntB(end,3,1) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,4,1) = mo_blockintegral(9);
                
                % x and y directed flux in slot on right hand side inner
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,3), slotypos(i));
                design.intBdata.slotIntB(end,1,2) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,2,2) = mo_blockintegral(9);
                
                % x and y directed flux in slot on right hand side outer
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,4), slotypos(i));
                design.intBdata.slotIntB(end,3,2) = mo_blockintegral(8);
                design.intBdata.slotIntB(end,4,2) = mo_blockintegral(9);
            end
        end

    end
    
    % store the relative coil/slot positions, we use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    design.intBdata.slotPos = [design.intBdata.slotPos; 
                               (-slotypos./design.taupm) + design.FirstSlotCenter + feapos];
    
    % sort the data in ascending position order
    [design.intBdata.slotPos, idx] = sort(design.intBdata.slotPos);
    design.intBdata.slotIntB = design.intBdata.slotIntB(idx,:,:);
    
end


function design = slotintAdata(design, feapos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.coillabellocs(1:(size(design.coillabellocs,1)/2),2);
        slotxpos = [design.coillabellocs(1:(size(design.coillabellocs,1)/2),1), ...
                    design.coillabellocs((1:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2,1) ];
                
        for i = 1:numel(slotypos)
            if exist('fpproc_interface_mex', 'file') == 3
                % vector potential in slot on left hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
                
                % vector potential in slot on right hand side
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i));
                design.intAdata.slotIntA(end,1,2) = solution.blockintegral(1);
            else
                % vector potential in slot on left hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
                
                % vector potential in slot on right hand side
                mo_clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i));
                design.intAdata.slotIntA(end,1,2) = mo_blockintegral(1);
            end
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = design.coillabellocs(1:2:(size(design.coillabellocs,1)/2),2);
        slotxpos = [design.coillabellocs(1:2:(size(design.coillabellocs,1)/2),1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2))+1,1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2+1,1), ...
                    design.coillabellocs((1:2:(size(design.coillabellocs,1)/2)) + size(design.coillabellocs,1)/2,1)];
              
        for i = 1:numel(slotypos)

            if exist('fpproc_interface_mex', 'file') == 3
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = solution.blockintegral(1);
                
                % vector potential flux in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,2), slotypos(i));
                design.intAdata.slotIntA(end,2,1) = solution.blockintegral(1);
                
                % vector potential in slot on right hand side inner
                % (leftmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,3), slotypos(i));
                design.intAdata.slotIntA(end,1,2) = solution.blockintegral(1);
                
                % vector potential in slot on right hand side outer
                % (rightmost) layer
                solution.clearblock();
                solution.selectblock(slotxpos(i,4), slotypos(i));
                design.intAdata.slotIntA(end,2,2) = solution.blockintegral(1);
                
            else
                % vector potential in slot on left hand side outer
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,1), slotypos(i));
                design.intAdata.slotIntA(end+1,1,1) = mo_blockintegral(1);
                
                % vector potential in slot on left hand side inner
                % (rightmost) layer
                solution.clearblock();
                mo_selectblock(slotxpos(i,2), slotypos(i));
                design.intAdata.slotIntA(end,2,1) = mo_blockintegral(1);
                
                % vector potential in slot on right hand side inner
                % (leftmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,3), slotypos(i));
                design.intAdata.slotIntA(end,1,2) = mo_blockintegral(1);
                
                % vector potential in slot on right hand side outer
                % (rightmost) layer
                mo_clearblock();
                mo_selectblock(slotxpos(i,4), slotypos(i));
                design.intAdata.slotIntA(end,2,2) = mo_blockintegral(1);
            end
        end

    end
    
    % store the relative coil/slot positions. We use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    design.intAdata.slotPos = [design.intAdata.slotPos; 
                              (-slotypos./design.taupm) + design.FirstSlotCenter + feapos];
    
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
    design.CoreLoss(1).coreycoords = design.CoreLoss(1).coreycoords + coords(1,1) + design.tsb + design.tc;

    design.CoreLoss(1).dx = min([0.01 / 5, design.taupm / 10]);
    design.CoreLoss(1).corexcoords = ((design.CoreLoss(1).dx/2):design.CoreLoss(1).dx:(design.taupm-design.CoreLoss(1).dx/2))';
    [design.CoreLoss(1).meshx, design.CoreLoss(1).meshy] = meshgrid(design.CoreLoss(1).coreycoords,design.CoreLoss(1).corexcoords);
    % The values along the first dimension (down the columns) of the arrays
    % will contain values sampled in the x-direction in the fea simulation.
    % In this case values along the dimension hbi, or ht in the case of the
    % teeth. The values along the second dimension are the values sampled
    % at each value xRVWp. The values along the third dimension of the
    % arrays will contain values which are sampled from the y direction of
    % the fea simulation, in this case along the dimension taupm, or tsb
    % and tc in the case of the teeth.
    design.CoreLoss(1).Bx = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    design.CoreLoss(1).By = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    design.CoreLoss(1).Bz = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    design.CoreLoss(1).Hy = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    design.CoreLoss(1).Hx = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    design.CoreLoss(1).Hz = zeros([size(design.CoreLoss(1).meshx,2), npos, size(design.CoreLoss(1).meshx,1)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(1).dz = design.Rmo - design.Rmi;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(1).xstep = (feapos(2) - feapos(1)) * design.taupm;
            
    % now set up the tooth body sampling positions
    design.CoreLoss(2).dy = (design.tsb+design.tc) / 10;
    design.CoreLoss(2).coreycoords =  ((design.CoreLoss(1).dy/2):design.CoreLoss(1).dy:(design.tsb+design.tc-design.CoreLoss(1).dy/2))';
    design.CoreLoss(2).coreycoords = design.CoreLoss(2).coreycoords + coords(1,1);
    
    design.tautm = design.tausm-design.taucs;
    design.CoreLoss(2).dx = min([0.01 / 5, design.tautm / 4]);
    
    halftoothcorex = (design.CoreLoss(2).dx/2:design.CoreLoss(2).dx:design.tautm/2-design.CoreLoss(2).dx/2)';
    fulltoothcorex = (design.CoreLoss(2).dx/2:design.CoreLoss(2).dx:design.tautm-design.CoreLoss(2).dx/2)';
    if design.ypd == 1
        
        design.CoreLoss(2).corexcoords = [halftoothcorex; ...
                                          zeros(numel(halftoothcorex)*design.ypn-1,1); ...
                                          halftoothcorex + design.taupm - (design.tautm/2) ];
                                      
        for i = 1:(design.ypn-1)
            design.CoreLoss(2).corexcoords((i*numel(fulltoothcorex)):((i+1)*numel(fulltoothcorex)-1),1) = ...
                fulltoothcorex + i*design.tausm - design.tautm/2;
        end
                                  
    elseif design.ypd == 2
        
        design.CoreLoss(2).corexcoords = [halftoothcorex; ...
                                          zeros(numel(halftoothcorex)*floor(design.ypn/2),1)];

        for i = 1:(floor(design.ypn/2))
            design.CoreLoss(2).corexcoords((i*numel(fulltoothcorex)):((i+1)*numel(fulltoothcorex)-1),1) = ...
                fulltoothcorex + i*design.tausm - design.tautm/2;
        end
        
    else
        error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end
    
    design.CoreLoss(2).corexcoords = [ design.CoreLoss(2).corexcoords; 
                                       design.CoreLoss(2).corexcoords ];
                                   
    design.CoreLoss(2).coreycoords = [ design.CoreLoss(2).coreycoords; 
                                       design.CoreLoss(2).coreycoords + design.ty + design.tsb + design.tc + design.CoreLoss(2).dy/2 ];
    
    [design.CoreLoss(2).meshx, design.CoreLoss(2).meshy] = ... 
        meshgrid(design.CoreLoss(2).coreycoords,design.CoreLoss(2).corexcoords);
    
    design.CoreLoss(2).Bx = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    design.CoreLoss(2).By = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    design.CoreLoss(2).Bz = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    design.CoreLoss(2).Hy = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    design.CoreLoss(2).Hx = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    design.CoreLoss(2).Hz = zeros([size(design.CoreLoss(2).meshx,2), npos, size(design.CoreLoss(2).meshx,1)]);
    % make dz value equal to the depth of the simulation as it is used
    % to calculate the volume of material generating the loss (volume
    % will be dx X dy X dz)
    design.CoreLoss(2).dz = design.CoreLoss(1).dz;
    % add xstep which is the size of the steps in xRVWp denormalised. This
    % is later used to find the value of dB/dx at each position
    design.CoreLoss(2).xstep = design.CoreLoss(1).xstep;
    
    % shoes
    if (design.taucs - design.tausgm) > design.FEMMTol
        
        shoedx = min([0.01 / 5, (design.taucs - design.tausgm)/2 / 10]);
        
        shoedy = design.tsb / 10;
        
        design.taucssm = (design.taucs - design.tausgm) / 2;
        
        shoeangle = atan((design.tsb - design.tsg) / design.taucssm);
        
        xcoords = ((shoedx/2):shoedx:(design.taucssm-design.CoreLoss(1).dx/2));
        
        for i = 1:numel(xcoords)
            
            % copy over the information about the materials etc to the new
            % structure
            design.CoreLoss(2+i).fq = design.CoreLoss(2).fq;
            design.CoreLoss(2+i).Bq = design.CoreLoss(2).Bq;
            design.CoreLoss(2+i).Pq = design.CoreLoss(2).Pq;
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
                design.CoreLoss(2+i).coreycoords = [coords(1,1) + shoewid/2, coords(2,1) - shoewid/2];
            else
                design.CoreLoss(2+i).coreycoords = ((design.CoreLoss(2+i).dy/2):design.CoreLoss(2+i).dy:(shoewid-design.CoreLoss(2+i).dy/2))';
                design.CoreLoss(2+i).coreycoords = shoeycoords + coords(1,1);
                design.CoreLoss(2+i).coreycoords = [ design.CoreLoss(2+i).coreycoords; coords(2,1) - flipud(shoeycoords) ];
            end
            
            % add a single x coordinate for the tooth shoe at the bottom of
            % the sim
            design.CoreLoss(2+i).dx = shoedx;
            design.CoreLoss(2+i).corexcoords = design.tautm/2 +  design.taucssm - xcoords(i);
            
            if design.ypd == 1
                
                for ii = 1:(design.ypn-1)
                    
                    design.CoreLoss(2+i).corexcoords = ...
                        [design.CoreLoss(2+i).corexcoords;
                         ii*design.tausm - design.tautm/2 - (design.taucssm) + xcoords(i);
                         ii*design.tausm + design.tautm/2 + (design.taucssm) - xcoords(i) ];
                     
                end
                
            elseif design.ypd == 2
                
                for ii = 1:(floor(design.ypn/2))
                    
                    design.CoreLoss(2+i).corexcoords = ...
                        [design.CoreLoss(2+i).corexcoords;
                         ii*design.tausm - design.tautm/2 - (design.taucssm) + xcoords(i);
                         ii*design.tausm + design.tautm/2 + (design.taucssm) - xcoords(i) ];
                     
                end
                
            else
                error('denominator of slots per pole must be 1 or 2, other values not yet supported')
            end
            
            % generate the sample points
            [design.CoreLoss(2+i).meshx, design.CoreLoss(2+i).meshy] = meshgrid(design.CoreLoss(2+i).coreycoords,design.CoreLoss(2+i).corexcoords);

            design.CoreLoss(2+i).Bx = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            design.CoreLoss(2+i).By = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            design.CoreLoss(2+i).Bz = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            design.CoreLoss(2+i).Hy = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            design.CoreLoss(2+i).Hx = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            design.CoreLoss(2+i).Hz = zeros([size(design.CoreLoss(2+i).meshx,2), npos, size(design.CoreLoss(2+i).meshy,1)]);
            
        end
    end

end

