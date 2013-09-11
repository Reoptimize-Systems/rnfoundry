function [airgapdisps, feb, geom, ur] = fielddef_FM(width, depth, height, BeamInfo, MForce, doplot)
% fielddef_FM: Determines the deflections in the structure of a non-tubular 
% machine
%
%

    if nargin < 6
        doplot = false;
    end
    
    %% Material Properties

    g = 9.81;

    % the material used in the beams at the edges to which the field pole
    % support beams are fixed
    prop = property_linel_iso (struct('E', BeamInfo.SupportBeams.E, 'nu', BeamInfo.SupportBeams.nu));
    BeamInfo.SupportBeams.mater = mater_defor_ss_linel_triax (struct('property', prop));
    BeamInfo.SupportBeams.label = 1;

    % the material used in the field pole support beams
    prop = property_linel_iso (struct('E', BeamInfo.FieldPoleSupports.E, 'nu', BeamInfo.FieldPoleSupports.nu));
    BeamInfo.FieldPoleSupports.mater = mater_defor_ss_linel_triax (struct('property', prop));
    BeamInfo.FieldPoleSupports.label = 2;

    % the material used in the members which hold the two sides of the
    % translator apart
    prop = property_linel_iso (struct('E', BeamInfo.FieldWebs.E, 'nu', BeamInfo.FieldWebs.nu));
    BeamInfo.FieldWebs.mater = mater_defor_ss_linel_triax (struct('property', prop));
    BeamInfo.FieldWebs.label = 3;

    % the material of which the guide rails are constructed
    prop = property_linel_iso (struct('E', BeamInfo.GuideRails.E, 'nu', BeamInfo.GuideRails.nu));
    BeamInfo.GuideRails.mater = mater_defor_ss_linel_triax (struct('property', prop));
    BeamInfo.GuideRails.label = 4;

    % the very stiff material used in the connections between the guide rails
    % and the support beams 
    prop = property_linel_iso (struct('E', BeamInfo.GuideBearings.E, 'nu', BeamInfo.GuideBearings.nu));
    BeamInfo.GuideBearings.mater = mater_defor_ss_linel_triax (struct('property', prop));
    BeamInfo.GuideBearings.label = 5;


    %% Parameters
    
%     tolerance = min(abs(znodes(1:end-1,1) - znodes(2:end,1))) / 1000;
    tolerance = 100 * eps;
    
    % The centre of the machine is located a (0,0,0) so it is convenient to
    % half the dimensions to find the appropriate coordinates for the model
    x = width / 2;
    y = depth / 2;
    z = height / 2 ;

    % calculate appropriate points for the members keeping the two field parts
    % apart
    if BeamInfo.FieldWebs.NoPerSide == 1
        % if there's only one field web, it has to go in the middle
        zsupp = 0;
    else
        % otherwize arrange the supports evenly spaced from the extremities
        % inward
        zsupp = (-z:2*z/(BeamInfo.FieldWebs.NoPerSide-1):z)';
    end

    % each support beam beam will be divided into n elements
    n = BeamInfo.FieldPoleSupports.NoPerSide - 1;

    % calculate appropriate node heights for the beams on the field back iron
    zframe = (-z:2*z/n:z)';

    % there will also be connections to the linear guide rails
    zguidecon = (-z+(2*z/(BeamInfo.GuideBearings.NoPerGuide+1)):2*z/(BeamInfo.GuideBearings.NoPerGuide+1):z-(2*z/(BeamInfo.GuideBearings.NoPerGuide+1)))';

    % get the unique members of all the nodes, put them into a column vector
    % and sort it in ascending order
%     znodes = sortrows(consolidator([zsupp; zframe; zguidecon],[],[],tolerance));  
%     znodes = sortrows([zsupp; zframe; zguidecon]);  
%     znodes = znodes(abs(znodes(2:end) - znodes(1:end-1)) > tolerance);
    znodes = sortrows(uniquetol([zsupp; zframe; zguidecon],tolerance));

    %% Draw the main support beams

    % these are the beams to which the beams spanning the field back iron stack
    % length will be fixed

    % Get the moments of inertia and cross-sectional area of the support beams
    % from the BeamInfo structure
    Parameters = BeamInfo.SupportBeams;

    % set the x1x2_vector of the field support beam elements so that the
    % 'strongest' direction is in the y direction (as the first axis will be
    % pointing in the z-direction, the same direction as the beam elements)
    Parameters.x1x2_vector = [1 0 0];
    Parameters.label = 1;

    nodepairs = [ repmat([-x,-y], size(znodes, 1)-1, 1), znodes(1:end-1,1), repmat([-x,-y], size(znodes, 1)-1, 1), znodes(2:end,1);
                  repmat([x,-y],  size(znodes, 1)-1, 1), znodes(1:end-1,1), repmat([x,-y],  size(znodes, 1)-1, 1), znodes(2:end,1);
                  repmat([-x,y],  size(znodes, 1)-1, 1), znodes(1:end-1,1), repmat([-x,y],  size(znodes, 1)-1, 1), znodes(2:end,1);
                  repmat([x,y],   size(znodes, 1)-1, 1), znodes(1:end-1,1), repmat([x,y],   size(znodes, 1)-1, 1), znodes(2:end,1) ];

    nels = 2;

    for j = 1:numel(znodes)-1:3*(numel(znodes))

        if j == 1
            [fens, gcells] = Beam_member([nodepairs(j,1:3); nodepairs(j,4:6)], nels, @gcellset_beam3, Parameters);
        else

            [fens2, gcells2] = Beam_member([nodepairs(j,1:3); nodepairs(j,4:6)], nels, @gcellset_beam3, Parameters);

            [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
            gcells = cat(gcells,gcells2);

        end

        for i = 1:numel(znodes)-2

            [fens2, gcells2] = Beam_member([nodepairs(j+i,1:3); nodepairs(j+i,4:6)], nels, @gcellset_beam3, Parameters);

            [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
            gcells = cat(gcells,gcells2);

        end

    end

    %% Add the long load bearing guide rails for the translator 

    % Get the moments of inertia and cross-sectional area of the support beams
    % from the BeamInfo structure
    Parameters = BeamInfo.GuideRails;

    % set the x1x2_vector of the field support beam elements so that the
    % 'strongest' direction is in the y direction (as the first axis will be
    % pointing in the z-direction, the same direction as the beam elements)
    Parameters.x1x2_vector = [1 0 0];

    nels = 7;

    % zextreme = (design.poles(1) .* design.Wp);
    zextreme = BeamInfo.GuideRails.length / 2;

    % we must shift the support beams in the positive and negative y directions
    % to accomodate the connections to the linear guides
    yshift = 50/1000;

    nodepairs = [-x, (-y - yshift), -zextreme,  -x, (-y - yshift),   zguidecon(1); % 1
                 repmat([-x, (-y - yshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([-x, (-y - yshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 1
                 -x, (-y - yshift),  zguidecon(end), -x, (-y - yshift),   zextreme; %1
                 -x, ( y + yshift), -zextreme,  -x, ( y + yshift),  zguidecon(1); % 2
                 repmat([-x, ( y + yshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([-x, ( y + yshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 2
                 -x, ( y + yshift),  zguidecon(end), -x, ( y + yshift),  zextreme; %2
                  x, (-y - yshift), -zextreme,   x, (-y - yshift),  zguidecon(1); % 3
                  repmat([x, (-y - yshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([x, (-y - yshift)], numel(zguidecon)-1, 1), zguidecon(2:end);  % 3
                  x, (-y - yshift),  zguidecon(end),  x, (-y - yshift),  zextreme; %3
                  x, ( y + yshift), -zextreme,   x, ( y + yshift),  zguidecon(1); % 4
                  repmat([x, ( y + yshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([x, ( y + yshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 4
                  x, ( y + yshift),  zguidecon(end),  x, ( y + yshift),  zextreme ];

    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], nels, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);

    end         

    %% Add connections between guide rails and translator frames

    % These are the connection points between the guide rails and the field
    % frame. The will be made much stiffer than the rest of the structural
    % material

    Parameters = BeamInfo.GuideBearings;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x, (-y - yshift)], numel(zguidecon), 1), zguidecon, repmat([-x, -y], numel(zguidecon), 1), zguidecon;
                 repmat([-x, ( y + yshift)], numel(zguidecon), 1), zguidecon, repmat([-x,  y], numel(zguidecon), 1), zguidecon;
                 repmat([ x, (-y - yshift)], numel(zguidecon), 1), zguidecon, repmat([ x, -y], numel(zguidecon), 1), zguidecon;
                 repmat([ x, ( y + yshift)], numel(zguidecon), 1), zguidecon, repmat([ x,  y], numel(zguidecon), 1), zguidecon;];


    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], 3, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);  

    end

    %% Now add the side supports separating the two field parts

    % These are the supports that hold the two parts of the field apart on
    % either side of the machine

    Parameters = BeamInfo.FieldWebs;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x, -y], numel(zsupp), 1), zsupp, repmat([-x, y], numel(zsupp), 1), zsupp;
                 repmat([ x, -y], numel(zsupp), 1), zsupp, repmat([ x, y], numel(zsupp), 1), zsupp];


    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], 3, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);  

    end

    %% Add the main field magnet support beams 

    % These are the beams to which the magnets are attached

    Parameters = BeamInfo.FieldPoleSupports;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x,-y], numel(zframe), 1), zframe, repmat([x,-y], numel(zframe), 1), zframe;
                 repmat([-x, y], numel(zframe), 1), zframe, repmat([x, y], numel(zframe), 1), zframe];


    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], BeamInfo.FieldPoleSupports.Sections, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);  

    end

    %% The FE Blocks

    % the material used in the field pole support beams
    BeamInfo.FieldPoleSupports.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.FieldPoleSupports.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.FieldPoleSupports.label)))));

    % the material used in the beams at the edges to which the field pole
    % support beams are fixed
    BeamInfo.SupportBeams.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.SupportBeams.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.SupportBeams.label)))));

    % the material used in the members which hold the two sides of the
    % translator apart
    BeamInfo.FieldWebs.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.FieldWebs.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.FieldWebs.label)))));

    % the material of which the guide rails are constructed
    BeamInfo.GuideRails.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.GuideRails.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.GuideRails.label)))));
    % BeamInfo.GuideRails.mater

    % the very stiff material used in the connections between the guide rails
    % and the support beams 
    % finite element block of only the guide connections
    BeamInfo.GuideBearings.feb = feblock_defor_ss_beam3 (struct ('mater', BeamInfo.GuideBearings.mater, ...
            'gcells', subset(gcells, gcell_select(fens, gcells, struct('label', BeamInfo.GuideBearings.label)))));


    % Finite element block containing all the elements (useful for drawing)
    feb = feblock_defor_ss_beam3 (struct ('gcells', gcells));


    %% Geometry and displacement fields

    geom = field(struct ('name', 'geom', ...
                         'dim', 3, ...
                         'fens', fens));
                     
    ur  = field(struct ('name', 'ur', ...
                        'dim', 6, ...
                        'data', zeros(get(geom,'nfens'),6)));

    %% Apply EBC's

    % The clamped end nodes are selected. These will not actually be fully
    % clamped but merely acting as if on simple supports
    clampedn = [ fenode_select(fens, struct('box', [x, -x, y+yshift, -y-yshift, -zextreme, -zextreme], 'inflate', 0.01)), ...
                 fenode_select(fens, struct('box', [x, -x, y+yshift, -y-yshift,  zextreme,  zextreme], 'inflate', 0.01)) ];

    ebc_fenids = clampedn;

    ebc_prescribed = [1];

    % Movement of the nodes in their x, y and z coordinates (1, 2 and 3) is
    % prevented, as is rotation about their 1st axis (4)
    ebc_comp = [1,2,3];

    % ebc_comp = [];

    for i = 1:numel(ebc_fenids)

        ebc_val = 0; % ebc_fenids * 0;

        ur = set_ebc(ur, ebc_fenids(i), ebc_prescribed, ebc_comp, ebc_val);

        ur = apply_ebc (ur);

    end

    %% Number equations

    ur = numbereqns (ur);
    neqns = get(ur, 'neqns');

    %%  Assemble the stiffness matrix

    K = start (sparse_sysmat, neqns);

    K = assemble (K, stiffness(BeamInfo.FieldPoleSupports.feb, geom, ur));

    K = assemble (K, stiffness(BeamInfo.SupportBeams.feb, geom, ur));

    K = assemble (K, stiffness(BeamInfo.FieldWebs.feb, geom, ur));

    K = assemble (K, stiffness(BeamInfo.GuideRails.feb, geom, ur));

    K = assemble (K, stiffness(BeamInfo.GuideBearings.feb, geom, ur));

    %% Loads

    F = sysvec;
    F = start (F, neqns);

    % % Apply a force at the middle
    % tipn = fenode_select(fens, struct('box', update_box([], [-x,-y, 0]), 'inflate', tolerance));
    % 
    % evs = loads(nodal_load(struct ('id', tipn, 'dir', 2, 'magn', loadmagn)), ur);
    % 
    % F = assemble (F, evs);
    
    % The weight of all the members, except the bearings which weigh nothing
    
    fi = force_intensity(struct('magn', [0; ...
                                        -cos(BeamInfo.FieldPoleSupports.AngleFromHorizontal) * BeamInfo.FieldPoleSupports.rho .* BeamInfo.FieldPoleSupports.A .* g;
                                        -sin(BeamInfo.FieldPoleSupports.AngleFromHorizontal) * BeamInfo.FieldPoleSupports.rho .* BeamInfo.FieldPoleSupports.A .* g]));
    F = assemble(F, distrib_loads(BeamInfo.FieldPoleSupports.feb, geom, ur, fi));

    fi = force_intensity(struct('magn', [0; ...
                                         -cos(BeamInfo.SupportBeams.AngleFromHorizontal) * BeamInfo.SupportBeams.rho .* BeamInfo.SupportBeams.A .* g; ...
                                         -sin(BeamInfo.SupportBeams.AngleFromHorizontal) * BeamInfo.SupportBeams.rho .* BeamInfo.SupportBeams.A .* g]));
    F = assemble(F, distrib_loads(BeamInfo.SupportBeams.feb, geom, ur, fi));

    fi = force_intensity(struct('magn', [0; ...
                                         -cos(BeamInfo.FieldWebs.AngleFromHorizontal) * BeamInfo.FieldWebs.rho .* BeamInfo.FieldWebs.A .* g; ...
                                         -sin(BeamInfo.FieldWebs.AngleFromHorizontal) * BeamInfo.FieldWebs.rho .* BeamInfo.FieldWebs.A .* g]));
    F = assemble(F, distrib_loads(BeamInfo.FieldWebs.feb, geom, ur, fi));

    fi = force_intensity(struct('magn', [0; ...
                                         -cos(BeamInfo.GuideRails.AngleFromHorizontal) * BeamInfo.GuideRails.rho .* BeamInfo.GuideRails.A .* g; ...
                                         -sin(BeamInfo.GuideRails.AngleFromHorizontal) * BeamInfo.GuideRails.rho .* BeamInfo.GuideRails.A .* g]));
    F = assemble(F, distrib_loads(BeamInfo.GuideRails.feb, geom, ur, fi));
    
%     % Graphics display
%     gv = graphic_viewer;
%     
%     gv = reset (gv,[]);
%     
%     % The size of the deflections in the viewer will be increase by a given
%     % scale
%     scale = 100;
%     
%     % get the displacements
%     u = slice(ur, (1:3), 'u');
%     
%     % Get the rotations
%     rot = slice(ur, (4:6),'rot');
%     
%     % Now draw the results
%     draw(feb, gv, struct ('x', geom, 'u', 0*u, 'rot', 0*rot, 'facecolor','none', 'drawscale', 0.1));
%     % 
    
    if iscell(MForce)
        
        if numel(MForce) == 1;
            % MForce is the force per unit length on am element. if only one value
            % of MForce has been supplied, the same force is applied to all
            % sections of all beams
            MForce = repmat(MForce, BeamInfo.FieldPoleSupports.NoPerSide, BeamInfo.FieldPoleSupports.Sections);
            
            for i = 1:size(MForce, 1)
                for j = 1:size(MForce, 2)
                    tempforce = MForce{i,j,2};
                    tempforce(2) = -tempforce(2);
                    MForce{i,j,2} = tempforce;
                end
            end
        end

        % Apply the forces to the field elements
        for j = 1:BeamInfo.FieldPoleSupports.NoPerSide

            zfieldsecs = -z + (j-1)*(2*z/n);

            fieldelementlen = 2*x / BeamInfo.FieldPoleSupports.Sections;

            for i = 1:BeamInfo.FieldPoleSupports.Sections

                selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + fieldelementlen/2 + ((i-1) * fieldelementlen), -y, zfieldsecs])) );

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                %         draw(febsel, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue', 'drawscale', 0.1) )

                fi = force_intensity(struct('magn', MForce{j,i,1}'));

                F = assemble(F, distrib_loads(febsel, geom, ur, fi));

                selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + fieldelementlen/2 + ((i-1) * fieldelementlen), y, zfieldsecs])) );

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                %         draw(febsel, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue', 'drawscale', 0.1) )

                fi = force_intensity(struct('magn', MForce{j,i,2}'));

                F = assemble(F, distrib_loads(febsel, geom, ur, fi));

            end

        end

        
    else

        if numel(MForce) == 1;
            % MForce is the force per unit length on am element. if only one value
            % of MForce has been supplied, the same force is applied to all
            % sections of all beams
            MForce = repmat(MForce, BeamInfo.FieldPoleSupports.NoPerSide, BeamInfo.FieldPoleSupports.Sections);
            MForce(:,:,2) = -MForce(:,:,1);
        end

        % Apply the forces to the field elements
        for j = 1:BeamInfo.FieldPoleSupports.NoPerSide

            zfieldsecs = -z + (j-1)*(2*z/n);

            fieldelementlen = 2*x / BeamInfo.FieldPoleSupports.Sections;

            for i = 1:BeamInfo.FieldPoleSupports.Sections

                selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + fieldelementlen/2 + ((i-1) * fieldelementlen), -y, zfieldsecs])) );

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                %         draw(febsel, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue', 'drawscale', 0.1) )

                fi = force_intensity(struct('magn', [0, MForce(j,i,1), 0]'));

                F = assemble(F, distrib_loads(febsel, geom, ur, fi));

                selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + fieldelementlen/2 + ((i-1) * fieldelementlen), y, zfieldsecs])) );

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                %         draw(febsel, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue', 'drawscale', 0.1) )

                fi = force_intensity(struct('magn', [0, MForce(j,i,2), 0]'));

                F = assemble(F, distrib_loads(febsel, geom, ur, fi));

            end

        end

    end

    % Solve the system storing the displacement (u) and rotation fields (r) in
    % ur
    ur = scatter_sysvec(ur, K \ F);

    % The nodes will be extracted from the top to the bottom
    % distn1 = fenode_select(fens, struct('box', [x-(2*x)/BeamInfo.FieldPoleSupports.Sections, -x+(2*x)/BeamInfo.FieldPoleSupports.Sections, -y, -y, z, -z], 'inflate', tolerance));
    % distn2 = fenode_select(fens, struct('box', [x-(2*x)/BeamInfo.FieldPoleSupports.Sections, -x+(2*x)/BeamInfo.FieldPoleSupports.Sections,  y,  y, z, -z], 'inflate', tolerance));

    % The nodes will be extracted from the top to the bottom
    distn1 = [];
    distn2 = [];
%     tolerance = 2*eps;
    
    for j = 1:BeamInfo.FieldPoleSupports.NoPerSide

        zfieldsecs = z - (j-1)*(2*z/n);

        distn1 = [distn1, fenode_select(fens, struct('box', [x, -x, -y, -y, zfieldsecs, zfieldsecs], 'inflate', tolerance))];
        distn2 = [distn2, fenode_select(fens, struct('box', [x, -x,  y,  y, zfieldsecs, zfieldsecs], 'inflate', tolerance))];

    end

    % Get and display the deflection of the node where the load is applied
    % load_displacement_mm = 1000 * somel(gather(ur, tipn, 'values', 'noreshape'), 1:3)
    % load_displacement_mm = 1000 * distndisp(:,1:3)

    distndisp = gather(ur, distn1', 'values', 'noreshape');
    airgapdisps(:,:,1) = reshape( distndisp(:,2)', BeamInfo.FieldPoleSupports.Sections+1, [])';

    airgapdisps(:,:,2) = zeros(size(airgapdisps));

    distndisp = gather(ur, distn2', 'values', 'noreshape');
    airgapdisps(:,:,3) = reshape( distndisp(:,2)', BeamInfo.FieldPoleSupports.Sections+1, [])';

% %% 
%     keyboard
    if doplot
        
        % Graphics display
        gv = graphic_viewer;

        gv = reset (gv,[]);

        % The size of the deflections in the viewer will be increase by a given
        % scale
        scale = 100;

        % get the displacements
        u = slice(ur, (1:3), 'u');

        % Get the rotations
        rot = slice(ur, (4:6),'rot');

        % Now draw the results
        draw(feb, gv, struct ('x', geom, 'u', 0*u, 'rot', 0*rot, 'facecolor', 'none', 'drawscale', 0.1));

        % % Some other plots that can be useful
        % selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x, 0, 0])) );
        % selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + fieldelementlen/2 + ((1-1) * fieldelementlen), -y, zfieldsecs])) );
        % selcells = subset(gcells, gcell_select(fens, gcells, struct('box', [x, -x,  y,  -y, z, z], 'inflate', tolerance, 'anynode', true)));
        % febsel = feblock_defor_ss_beam3 (struct ('gcells', selcells));
        % draw(febsel, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue') )
        % draw(febsel, gv, struct ('x', geom, 'u', scale * u, 'rot', scale * rot, 'facecolor', 'blue') )
        % draw(febguide, gv, struct ('x', geom, 'u', 0 * u, 'rot', 0 * rot, 'facecolor', 'blue') )



        draw(feb, gv, struct ('x', geom, 'u', scale * u, 'rot', scale * rot, 'facecolor', 'red', 'drawscale', 0.1));

        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        
        hold on

        % Draw the node numbers
%         draw(fens,gv, struct ('x', geom, 'u', 0*u, 'color', 'blue', 'offset', 0.025, 'fontname', 'Ariel', 'fontsize', 8));

        % Draw the coordinate axis
        draw_axes(gv, struct([]));
        
        hold off

    end

end