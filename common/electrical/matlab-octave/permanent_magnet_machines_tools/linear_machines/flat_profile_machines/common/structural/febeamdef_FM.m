function [airgapdisps, BeamInfo, fens, gcells] = febeamdef_FM(BeamInfo, MForce, doplot, fens, gcells, OuterSupportCells)
% febeamdef_FM: Determines the deflections in the structure of a non-tubular 
% machine using finite elements based on simple beam element types
%
%

    if nargin < 3
        doplot = false;
    elseif nargin > 6
        error('Incorrect number of arguments.')
    end

    [ x, y, z, n, ...
      znodes, ...
      zguidecon, ...
      zsupp, ...
      zframe, ...
      zextreme, ...
      guideyshift, ...
      tolerance ... 
     ] = dimsfebeamdef_FM(BeamInfo);
    
    %% Material Properties

    g = 9.81;

    if nargin > 3
        
        % the mesh has already been created and we merely need to update
        % the beam properties of the appropriate gcells
        
        % Get the moments of inertia and cross-sectional area of the support beams
        % from the BeamInfo structure
        BeamInfo.SupportBeams.label = 1;
        Parameters = BeamInfo.SupportBeams;
        Parameters.x1x2_vector = [1 0 0];
        gcells = set(gcells, 'beamprops', Parameters);
        
        
        % Get the moments of inertia and cross-sectional area of the support beams
        % from the BeamInfo structure
        BeamInfo.GuideRails.label = 4;
        Parameters = BeamInfo.GuideRails;
        Parameters.x1x2_vector = [1 0 0];
        gcells = set(gcells, 'beamprops', Parameters);
        
        
        BeamInfo.GuideBearings.label = 5;
        Parameters = BeamInfo.GuideBearings;
        Parameters.x1x2_vector = [0 0 1];
        gcells = set(gcells, 'beamprops', Parameters);
        
        
        BeamInfo.OuterWebs.label = 3;
        Parameters = BeamInfo.OuterWebs;
        Parameters.x1x2_vector = [0 0 1];
        gcells = set(gcells, 'beamprops', Parameters);
        
        
        BeamInfo.OuterPoleSupports.label = 2;
        Parameters = BeamInfo.OuterPoleSupports;
        Parameters.x1x2_vector = [0 0 1];
        gcells = set(gcells, 'beamprops', Parameters);

        if nargin < 6
            
            OuterSupportCells = cell(BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections, 3);

            for j = 1:BeamInfo.OuterPoleSupports.NoPerSide

                zoutersecs = -z + (j-1)*(2*z/n);

                outerelementlen = 2*x / BeamInfo.OuterPoleSupports.Sections;

                for i = 1:BeamInfo.OuterPoleSupports.Sections

                    selcell = subset(gcells,  gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), -y, zoutersecs])));

                    febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));
                    
                    OuterSupportCells{j,i,1} = febsel;
            
                    selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), y, zoutersecs])));
            
                    febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));
                    
                    OuterSupportCells{j,i,3} = febsel;

                end

            end
            
        end

        % As the beam parameters have changed

    else

%         % Create the mesh from scratch
%         [fens, gcells, BeamInfo, OuterSupportCells] = createmeshfebeamdef_FM(x, y, z, n, znodes, ...
%             zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo);

        % Create the mesh from scratch
        [fens, gcells, BeamInfo, OuterSupportCells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, ...
            zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo);

    end
    
    if ~isfield(BeamInfo, 'K')
        
        if ~isfield(BeamInfo.SupportBeams, 'mater')

            % the material used in the beams at the edges to which the field pole
            % support beams are fixed
            prop = property_linel_iso (struct('E', BeamInfo.SupportBeams.E, 'nu', BeamInfo.SupportBeams.nu));
            BeamInfo.SupportBeams.mater = mater_defor_ss_linel_triax (struct('property', prop));
            BeamInfo.SupportBeams.label = 1;

            % the material used in the field pole support beams
            prop = property_linel_iso (struct('E', BeamInfo.OuterPoleSupports.E, 'nu', BeamInfo.OuterPoleSupports.nu));
            BeamInfo.OuterPoleSupports.mater = mater_defor_ss_linel_triax (struct('property', prop));
            BeamInfo.OuterPoleSupports.label = 2;

            % the material used in the members which hold the two sides of the
            % translator apart
            prop = property_linel_iso (struct('E', BeamInfo.OuterWebs.E, 'nu', BeamInfo.OuterWebs.nu));
            BeamInfo.OuterWebs.mater = mater_defor_ss_linel_triax (struct('property', prop));
            BeamInfo.OuterWebs.label = 3;

            % the material of which the guide rails are constructed
            prop = property_linel_iso (struct('E', BeamInfo.GuideRails.E, 'nu', BeamInfo.GuideRails.nu));
            BeamInfo.GuideRails.mater = mater_defor_ss_linel_triax (struct('property', prop));
            BeamInfo.GuideRails.label = 4;

            % the very stiff material used in the connections between the guide rails
            % and the support beams
            prop = property_linel_iso (struct('E', BeamInfo.GuideBearings.E, 'nu', BeamInfo.GuideBearings.nu));
            BeamInfo.GuideBearings.mater = mater_defor_ss_linel_triax (struct('property', prop));
            BeamInfo.GuideBearings.label = 5;

        end
        
        %% Geometry and displacement fields

        BeamInfo.geom = field(struct ('name', 'geom', ...
                             'dim', 3, ...
                             'fens', fens));

        BeamInfo.ur  = field(struct ('name', 'ur', ...
                            'dim', 6, ...
                            'data', zeros(get(BeamInfo.geom,'nfens'),6)));

        %% Apply EBC's

        % The clamped end nodes are selected. These will not actually be fully
        % clamped but merely acting as if on simple supports
        if ~isfield(BeamInfo, 'clampedn')
            BeamInfo.clampedn = [ fenode_select(fens, struct('box', [x, -x, y+guideyshift, -y-guideyshift, -zextreme, -zextreme], 'inflate', 0.01)), ...
                                  fenode_select(fens, struct('box', [x, -x, y+guideyshift, -y-guideyshift,  zextreme,  zextreme], 'inflate', 0.01)) ];
        end
        
        ebc_fenids = BeamInfo.clampedn;

        ebc_prescribed = [1];

        % Movement of the nodes in their x, y and z coordinates (1, 2 and 3) is
        % prevented, as is rotation about their 1st axis (4)
        ebc_comp = [1,2,3];

        % ebc_comp = [];

        for i = 1:numel(ebc_fenids)

            ebc_val = 0; % ebc_fenids * 0;

            BeamInfo.ur = set_ebc(BeamInfo.ur, ebc_fenids(i), ebc_prescribed, ebc_comp, ebc_val);

            BeamInfo.ur = apply_ebc (BeamInfo.ur);

        end

        %% Number equations

        BeamInfo.ur = numbereqns (BeamInfo.ur);
        BeamInfo.neqns = get(BeamInfo.ur, 'neqns');

        %% The FE Blocks
    
        % the material used in the outer pole support beams
        BeamInfo.OuterPoleSupports.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.OuterPoleSupports.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.OuterPoleSupports.label)))));

        % the material used in the beams at the edges to which the outer pole
        % support beams are fixed
        BeamInfo.SupportBeams.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.SupportBeams.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.SupportBeams.label)))));

        % the material used in the members which hold the two sides of the
        % translator apart
        BeamInfo.OuterWebs.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.OuterWebs.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.OuterWebs.label)))));

        % the material of which the guide rails are constructed
        BeamInfo.GuideRails.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.GuideRails.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.GuideRails.label)))));
        % BeamInfo.GuideRails.mater

        % the very stiff material used in the connections between the guide rails
        % and the support beams
        % finite element block of only the guide connections
        BeamInfo.GuideBearings.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.GuideBearings.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.GuideBearings.label)))));

        %%  Assemble the stiffness matrix

        BeamInfo.K = start (sparse_sysmat, BeamInfo.neqns);

        BeamInfo.K = assemble (BeamInfo.K, stiffness(BeamInfo.OuterPoleSupports.feb, BeamInfo.geom, BeamInfo.ur));

        BeamInfo.K = assemble (BeamInfo.K, stiffness(BeamInfo.SupportBeams.feb, BeamInfo.geom, BeamInfo.ur));

        BeamInfo.K = assemble (BeamInfo.K, stiffness(BeamInfo.OuterWebs.feb, BeamInfo.geom, BeamInfo.ur));

        BeamInfo.K = assemble (BeamInfo.K, stiffness(BeamInfo.GuideRails.feb, BeamInfo.geom, BeamInfo.ur));
        
        BeamInfo.K = assemble (BeamInfo.K, stiffness(BeamInfo.GuideBearings.feb, BeamInfo.geom, BeamInfo.ur));

        if ~isfield(BeamInfo, 'StructWeightLoads') || numel(BeamInfo.StructWeightLoads) ~= 4
            % Preallocate an array of elevecset objects to hold the weight
            % loads
            BeamInfo.StructWeightLoads(4) = elevecset;
        end
        
        % calculate the weight loads due to the structureal sections own
        % masses (these will stay the same for the same design)
        fi = force_intensity(struct('magn', [0; ...
            -cos(BeamInfo.OuterPoleSupports.AngleFromHorizontal) * BeamInfo.OuterPoleSupports.rho .* BeamInfo.OuterPoleSupports.A .* g;
            -sin(BeamInfo.OuterPoleSupports.AngleFromHorizontal) * BeamInfo.OuterPoleSupports.rho .* BeamInfo.OuterPoleSupports.A .* g]));

        BeamInfo.StructWeightLoads(1) = distrib_loads(BeamInfo.OuterPoleSupports.feb, BeamInfo.geom, BeamInfo.ur, fi);

        fi = set_magn(fi, [0; ...
            -cos(BeamInfo.SupportBeams.AngleFromHorizontal) * BeamInfo.SupportBeams.rho .* BeamInfo.SupportBeams.A .* g; ...
            -sin(BeamInfo.SupportBeams.AngleFromHorizontal) * BeamInfo.SupportBeams.rho .* BeamInfo.SupportBeams.A .* g]);

        BeamInfo.StructWeightLoads(2) = distrib_loads(BeamInfo.SupportBeams.feb, BeamInfo.geom, BeamInfo.ur, fi);

        fi = set_magn(fi, [0; ...
            -cos(BeamInfo.OuterWebs.AngleFromHorizontal) * BeamInfo.OuterWebs.rho .* BeamInfo.OuterWebs.A .* g; ...
            -sin(BeamInfo.OuterWebs.AngleFromHorizontal) * BeamInfo.OuterWebs.rho .* BeamInfo.OuterWebs.A .* g]);

        BeamInfo.StructWeightLoads(3) = distrib_loads(BeamInfo.OuterWebs.feb, BeamInfo.geom, BeamInfo.ur, fi);

        fi = set_magn(fi, [0; ...
            -cos(BeamInfo.GuideRails.AngleFromHorizontal) * BeamInfo.GuideRails.rho .* BeamInfo.GuideRails.A .* g; ...
            -sin(BeamInfo.GuideRails.AngleFromHorizontal) * BeamInfo.GuideRails.rho .* BeamInfo.GuideRails.A .* g]);

        BeamInfo.StructWeightLoads(4) = distrib_loads(BeamInfo.GuideRails.feb, BeamInfo.geom, BeamInfo.ur, fi);
    
        % Get the nodes of interest to the air gap (the nodes on the
        % support beams) if they have not already been extracted
        
        if ~(isfield(BeamInfo, 'distn1') || isfield(BeamInfo, 'distn2')) || ...
                numel(BeamInfo.distn1) ~= BeamInfo.OuterPoleSupports.NoPerSide * (BeamInfo.OuterPoleSupports.Sections + 1)
            
            % The nodes will be extracted from the top to the bottom
            BeamInfo.distn1 = [];
            BeamInfo.distn2 = [];
            %     tolerance = 2*eps;

            for j = 1:BeamInfo.OuterPoleSupports.NoPerSide

                zoutersecs = z - (j-1)*(2*z/n);

                BeamInfo.distn1 = [BeamInfo.distn1, fenode_select(fens, struct('box', [x, -x, -y, -y, zoutersecs, zoutersecs], 'inflate', tolerance))];
                BeamInfo.distn2 = [BeamInfo.distn2, fenode_select(fens, struct('box', [x, -x,  y,  y, zoutersecs, zoutersecs], 'inflate', tolerance))];

            end

            BeamInfo.distn1 = BeamInfo.distn1';
            BeamInfo.distn2 = BeamInfo.distn2';

        end
        
    end

    %% Loads

    F = sysvec;
    F = start (F, BeamInfo.neqns);

    % % Apply a force at the middle
    % tipn = fenode_select(fens, struct('box', update_box([], [-x,-y, 0]), 'inflate', tolerance));
    % 
    % evs = loads(nodal_load(struct ('id', tipn, 'dir', 2, 'magn', loadmagn)), BeamInfo.ur);
    % 
    % F = assemble (F, evs);
    
    if ~iscell(MForce)

        % Convert to cell array
        matMForce = MForce;

        if numel(matMForce) == 1
            
            MForce = {[0, matMForce, 0]};

        else

            for i = 1:size(matMForce, 1)
                for j = 1:size(matMForce, 2)
                    MForce{i,j,1} = [0, matMForce(i,j,1), 0];
                    MForce{i,j,2} = [0, 0, 0];
                    MForce{i,j,3} = [0, matMForce(i,j,3), 0];
                end
            end

        end

    end

    if numel(MForce) == 1;
        % MForce is the force per unit length on an element. if only one value
        % of MForce has been supplied, the same force is applied to all
        % sections of all outer beams
        MForceVec = MForce{1};

        MForce = cell(BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections, 3);

        for i = 1:size(MForce, 1)
            for j = 1:size(MForce, 2)
                MForce{i,j,1} = MForceVec;
                MForce{i,j,2} = [0, 0, 0];
                MForce{i,j,3} = [MForceVec(1), -MForceVec(2), MForceVec(3)];
            end
        end
    end
    
    % preallocate an aray of load objects
    loads(BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections, 2) = elevecset;
    
    % create a force intensity object which we will use repeatedly
    fi = force_intensity(struct('magn', MForce{1,1,1}'));
    
    % Apply the forces to the outer elements
    for j = 1:BeamInfo.OuterPoleSupports.NoPerSide

        for i = 1:BeamInfo.OuterPoleSupports.Sections

            fi = set_magn(fi, MForce{j,i,1}');
            
            loads(j,i,1) = distrib_loads(OuterSupportCells{j,i,1}, BeamInfo.geom, BeamInfo.ur, fi);

            fi = set_magn(fi, MForce{j,i,3}');
            
            loads(j,i,2) = distrib_loads(OuterSupportCells{j,i,3}, BeamInfo.geom, BeamInfo.ur, fi);

        end

    end
    
    loads = [reshape(loads,1,[],1), BeamInfo.StructWeightLoads];
        
    F = assemble(F, loads);

    % Solve the system storing the displacement (u) and rotation fields (r) in
    % ur
    BeamInfo.ur = scatter_sysvec(BeamInfo.ur, BeamInfo.K \ F);

    % Get and display the deflection of the node where the load is applied
    % load_displacement_mm = 1000 * somel(gather(BeamInfo.ur, tipn, 'values', 'noreshape'), 1:3)
    % load_displacement_mm = 1000 * distndisp(:,1:3)

    distndisp = gather(BeamInfo.ur, BeamInfo.distn1, 'values', 'noreshape');
    airgapdisps(:,:,1) = reshape( distndisp(:,2)', BeamInfo.OuterPoleSupports.Sections+1, [])';

    airgapdisps(:,:,2) = zeros(size(airgapdisps));
    
    distndisp = gather(BeamInfo.ur, BeamInfo.distn2, 'values', 'noreshape');
    airgapdisps(:,:,3) = reshape( distndisp(:,2)', BeamInfo.OuterPoleSupports.Sections+1, [])';

    if doplot

        % Finite element block containing all the elements (useful for drawing)
        feb = feblock_defor_ss_beam3 (struct ('gcells', gcells));

        % Graphics display
        gv = graphic_viewer;

        gv = reset (gv,[]);

        % The size of the deflections in the viewer will be increase by a given
        % scale
        scale = 10;

        % get the displacements
        u = slice(BeamInfo.ur, (1:3), 'u');

        % Get the rotations
        rot = slice(BeamInfo.ur, (4:6),'rot');

        % Now draw the results
        draw(feb, gv, struct ('x', BeamInfo.geom, 'u', 0*u, 'rot', 0*rot, 'facecolor', 'none', 'drawscale', 0.2));

        % % Some other plots that can be useful
        % selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x, 0, 0])) );
        % selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((1-1) * outerelementlen), -y, zoutersecs])) );
        % selcells = subset(gcells, gcell_select(fens, gcells, struct('box', [x, -x,  y,  -y, z, z], 'inflate', tolerance, 'anynode', true)));
        % febsel = feblock_defor_ss_beam3 (struct ('gcells', selcells));
        % draw(febsel, gv, struct ('x', BeamInfo.geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue') )
        % draw(febsel, gv, struct ('x', BeamInfo.geom, 'u', scale * u, 'rot', scale * rot, 'facecolor', 'blue') )
        % draw(febguide, gv, struct ('x', BeamInfo.geom, 'u', 0 * u, 'rot',
        % 0 * rot, 'facecolor', 'blue') )

        draw(feb, gv, struct ('x', BeamInfo.geom, 'u', scale * u, 'rot', scale * rot, 'facecolor', 'red', 'drawscale', 0.2));

        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        
        hold on

        % Draw the node numbers
%         draw(fens,gv, struct ('x', BeamInfo.geom, 'u', 0*u, 'color', 'blue', 'offset', 0.025, 'fontname', 'Ariel', 'fontsize', 8));

        % Draw the coordinate axis
        draw_axes(gv, struct([]));
        
        hold off

    end

end