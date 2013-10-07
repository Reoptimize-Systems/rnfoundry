function [g, closingForces, BeamInfo, design] = feairgapclosure_FM(design, options, BeamInfo, OuterPoleWeight, beaminfofcn)
% determines the air gap closure using a finite element analyis s of a flat
% profile type machine structure
%
% Syntax
%
% [g, closingForces, BeamInfo, design] = feairgapclosure_FM(design, options, BeamInfo, OuterPoleWeight, beaminfofcn)
%
% 

    % determine if the design of the structure has changed since last time
    % we evaluated it, and if so generate the appropriate information
    if design.OuterStructureNBeams ~= BeamInfo.LastOuterStructureNBeams
        % store the previously generated previous meshes as they will be
        % wiped by the beaminfo function
        tempmeshes = BeamInfo.PreviousMeshes;
        tempnbeams = BeamInfo.PreviousMeshesNBeams;
        % generate information about the structure using the supplied
        % function
        BeamInfo = beaminfofcn(design, options);
        % Not the number of beams used in the structure, so later we can
        % check if it has changed
        BeamInfo.LastOuterStructureNBeams = design.OuterStructureNBeams;
        % strip any existing FE mesh, as the structure has changed and it
        % must be rebuilt
        design = rmiffield(design, 'fens');
        design = rmiffield(design, 'gcells');
        % copy the previous meshes back into the beam info stucture
        BeamInfo.PreviousMeshes = tempmeshes;
        BeamInfo.PreviousMeshesNBeams = tempnbeams;
    end
    
    % Get the initial machine airgaps for the translator Poles. These are
    % the air-gaps for every node, whereas the forces will be applied to
    % the sections between nodes
    g = repmat(design.g, [BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections+1, BeamInfo.StructSides]);

    OuterPoleWeight = -[OuterPoleWeight * cos(design.AngleFromHorizontal), OuterPoleWeight * sin(design.AngleFromHorizontal)];

    % calculate the proportion of the per-pole force on each beam
    forceRatio = forceratio_linear(design, span1(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod));
    
    % first do a quick analytical test of the structure to see if we should
    % bother with the detailed fe analysis
    Def = analyticalstructdef_FM(design, options, OuterPoleWeight(1) * design.Poles(2) / 2, BeamInfo, forceRatio);
                                
    if Def > (min(((1-options.gfactor) * 2), 1) * design.g)
        
        closingForces = [];
        
        g(:,:,1) = g(:,:,1) - Def;

        g(:,:,2) = g(:,:,2) - Def;
        
        return;
        
    end
    
    % create the mesh if not present in the design structure
    if ~(isfield(design, 'fens') && isfield(design, 'gcells'))
        
        if any(BeamInfo.PreviousMeshesNBeams == design.OuterStructureNBeams)
            % we have previously generated the same mesh, so reuse it
            meshind = find(BeamInfo.PreviousMeshesNBeams == design.OuterStructureNBeams);
            
            design.fens = BeamInfo.PreviousMeshes{meshind}{1};
            design.gcells = BeamInfo.PreviousMeshes{meshind}{2};
            design.OuterSupportCells = BeamInfo.PreviousMeshes{meshind}{3};
            BeamInfo.distn1 = BeamInfo.PreviousMeshes{meshind}{4};
            BeamInfo.distn2 = BeamInfo.PreviousMeshes{meshind}{5};
            BeamInfo.clampedn = BeamInfo.PreviousMeshes{meshind}{6};
            
        else

            % calculate some importent mesh dimensions
            [ x, y, z, n, ...
                znodes, ...
                zguidecon, ...
                zsupp, ...
                zframe, ...
                zextreme, ...
                guideyshift, ...
                tolerance ...
            ] = dimsfebeamdef_FM(BeamInfo);

            % create the finite element mesh for the design
            [design.fens, design.gcells, BeamInfo, design.OuterSupportCells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, ...
                zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo);

            % append the mesh to the list of previously generated meshes
            BeamInfo.PreviousMeshes{end+1} = {design.fens, ...
                                              design.gcells, ...
                                              design.OuterSupportCells, ...
                                              BeamInfo.distn1, ...
                                              BeamInfo.distn2, ...
                                              BeamInfo.clampedn};
                                          
            BeamInfo.PreviousMeshesNBeams(end+1) = design.OuterStructureNBeams;
        
        end
        
    end
    
    % initialise the deflections to zero
    Def = zeros([BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections+1, BeamInfo.StructSides+1]);
    oldDef = Def;
    maxDefDiff = 0;
    
    % initialise the forces to zero
    ForceSize = size(Def);
    ForceSize(2) = ForceSize(2) - 1;
    closingForces = zeros(ForceSize);
    
    defThresh = 0.01 * design.g;
    
    firstrun = true;
    interpgx = 0:(BeamInfo.width/BeamInfo.OuterPoleSupports.Sections):BeamInfo.width;
    iters = 0;
    
    while ((min(min(min(g))) > 0 && maxDefDiff > defThresh) || firstrun) && iters < 50
        
        % get new forces on sections for double sided machine

        for k = 1:BeamInfo.OuterPoleSupports.NoPerSide

            gxi = interpgx(1:end-1)+BeamInfo.width/BeamInfo.OuterPoleSupports.Sections/2;

            % The following assumes a polynomial has been fitted to
            % the force per unit area of translator surface versus
            % the air gap
            closingForces(k,:,1) = polyvaln(design.p_gforce, interp1(interpgx, g(k,:,1), gxi));
            % Get the force per unit length on a support beam
            closingForces(k,:,1) = closingForces(k,:,1) * design.PoleWidth * forceRatio;

            if size(g,3) > 1

                closingForces(k,:,3) = -polyvaln(design.p_gforce, interp1(interpgx, g(k,:,2), gxi));

                % Get the force per unit length on a support beam
                closingForces(k,:,3) = closingForces(k,:,3) * design.PoleWidth * forceRatio;

                closingForces(k,:,2) = closingForces(k,:,1) + closingForces(k,:,3);

            else

                closingForces(k,:,2) = -closingForces(k,:,1);

            end
        end

        % preallocate the cell array
        if size(g,3) > 1, thirddim = 3; else thirddim = 1; end
        
        ForceVecCells = cell([size(closingForces, 1), size(closingForces, 2), thirddim]);
        
        for ii = 1:size(closingForces, 1)

            for jj = 1:size(closingForces, 2)
                
                ForceVecCells{ii,jj,1} = [0, closingForces(ii,jj,1) + OuterPoleWeight(1), OuterPoleWeight(2)];
                
                if size(g,3) > 1
                    ForceVecCells{ii,jj,3} = [0, closingForces(ii,jj,3) + OuterPoleWeight(1), OuterPoleWeight(2)];
                    ForceVecCells{ii,jj,2} = 0;
                end

            end

        end

        % get the lateral deflections due to the forces
        
        if ~isfield(design, 'gcells')
            [Def, BeamInfo, design.fens, design.gcells, design.OuterSupportCells] = febeamdef_FM(BeamInfo, ForceVecCells, false);
        else
            [Def, BeamInfo] = febeamdef_FM(BeamInfo, ForceVecCells, false, design.fens, design.gcells, design.OuterSupportCells);
        end

        maxDefDiff = max(abs(Def(:)) - abs(oldDef(:)));
        
        oldDef = Def;

        % Reinitialise the machine airgaps
        g = repmat(design.g, [BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections+1, BeamInfo.StructSides]);
    
        % get the new air gaps
        g(:,:,1) = g(:,:,1) - Def(:,:,1) + Def(:,:,2);

        g(:,:,2) = g(:,:,2) + Def(:,:,3) - Def(:,:,2);

        firstrun = false;
        
        % count the number of iterations
        iters = iters + 1;

    end
    
end    
    