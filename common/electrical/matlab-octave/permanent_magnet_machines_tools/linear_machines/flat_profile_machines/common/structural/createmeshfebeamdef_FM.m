function [fens, gcells, BeamInfo, OuterSupportCells] = createmeshfebeamdef_FM(x, y, z, n, znodes, zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo)

    if ~isfield(BeamInfo.SupportBeams, 'label')

        BeamInfo.SupportBeams.label = 1;
        BeamInfo.OuterPoleSupports.label = 2;
        BeamInfo.OuterWebs.label = 3;
        BeamInfo.GuideRails.label = 4;
        BeamInfo.GuideBearings.label = 5;

    end
    
    %% Draw the main support beams

    % these are the beams to which the beams spanning the outer back iron stack
    % length will be fixed

    % Get the moments of inertia and cross-sectional area of the support beams
    % from the BeamInfo structure
    Parameters = BeamInfo.SupportBeams;

    % set the x1x2_vector of the outer support beam elements so that the
    % 'strongest' direction is in the y direction (as the first axis will be
    % pointing in the z-direction, the same direction as the beam elements)
    Parameters.x1x2_vector = [1 0 0];

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

    % set the x1x2_vector of the outer support beam elements so that the
    % 'strongest' direction is in the y direction (as the first axis will be
    % pointing in the z-direction, the same direction as the beam elements)
    Parameters.x1x2_vector = [1 0 0];

    nels = 7;

    nodepairs = [-x, (-y - guideyshift), -zextreme,  -x, (-y - guideyshift),   zguidecon(1); % 1
        repmat([-x, (-y - guideyshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([-x, (-y - guideyshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 1
        -x, (-y - guideyshift),  zguidecon(end), -x, (-y - guideyshift),   zextreme; %1
        -x, ( y + guideyshift), -zextreme,  -x, ( y + guideyshift),  zguidecon(1); % 2
        repmat([-x, ( y + guideyshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([-x, ( y + guideyshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 2
        -x, ( y + guideyshift),  zguidecon(end), -x, ( y + guideyshift),  zextreme; %2
        x, (-y - guideyshift), -zextreme,   x, (-y - guideyshift),  zguidecon(1); % 3
        repmat([x, (-y - guideyshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([x, (-y - guideyshift)], numel(zguidecon)-1, 1), zguidecon(2:end);  % 3
        x, (-y - guideyshift),  zguidecon(end),  x, (-y - guideyshift),  zextreme; %3
        x, ( y + guideyshift), -zextreme,   x, ( y + guideyshift),  zguidecon(1); % 4
        repmat([x, ( y + guideyshift)], numel(zguidecon)-1, 1), zguidecon(1:end-1), repmat([x, ( y + guideyshift)], numel(zguidecon)-1, 1), zguidecon(2:end); % 4
        x, ( y + guideyshift),  zguidecon(end),  x, ( y + guideyshift),  zextreme ];

    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], nels, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);

    end

    %% Add connections between guide rails and translator frames

    % These are the connection points between the guide rails and the outer
    % frame. The will be made much stiffer than the rest of the structural
    % material

    Parameters = BeamInfo.GuideBearings;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x, (-y - guideyshift)], numel(zguidecon), 1), zguidecon, repmat([-x, -y], numel(zguidecon), 1), zguidecon;
        repmat([-x, ( y + guideyshift)], numel(zguidecon), 1), zguidecon, repmat([-x,  y], numel(zguidecon), 1), zguidecon;
        repmat([ x, (-y - guideyshift)], numel(zguidecon), 1), zguidecon, repmat([ x, -y], numel(zguidecon), 1), zguidecon;
        repmat([ x, ( y + guideyshift)], numel(zguidecon), 1), zguidecon, repmat([ x,  y], numel(zguidecon), 1), zguidecon;];


    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], 3, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);

    end

    %% Now add the side supports separating the two outer parts

    % These are the supports that hold the two parts of the outer apart on
    % either side of the machine

    Parameters = BeamInfo.OuterWebs;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x, -y], numel(zsupp), 1), zsupp, repmat([-x, y], numel(zsupp), 1), zsupp;
        repmat([ x, -y], numel(zsupp), 1), zsupp, repmat([ x, y], numel(zsupp), 1), zsupp];


    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], 3, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);

    end

    %% Add the main outer support beams

    % These are the beams supporting the Poles on the outer part to
    % which the outer magnets/coils are attached

    Parameters = BeamInfo.OuterPoleSupports;

    Parameters.x1x2_vector = [0 0 1];

    nodepairs = [repmat([-x,-y], numel(zframe), 1), zframe, repmat([x,-y], numel(zframe), 1), zframe;
        repmat([-x, y], numel(zframe), 1), zframe, repmat([x, y], numel(zframe), 1), zframe];

    for i = 1:size(nodepairs,1)

        [fens2, gcells2] = Beam_member([nodepairs(i,1:3); nodepairs(i,4:6)], BeamInfo.OuterPoleSupports.Sections, @gcellset_beam3, Parameters);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, tolerance);
        gcells = cat(gcells,gcells2);

    end
    
    
    %% Get the beam segment finite element blocks

    OuterSupportCells = cell(BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections, 3);

    for j = 1:BeamInfo.OuterPoleSupports.NoPerSide

        zoutersecs = -z + (j-1)*(2*z/n);

        outerelementlen = 2*x / BeamInfo.OuterPoleSupports.Sections;

        for i = 1:BeamInfo.OuterPoleSupports.Sections

            %                     OuterSupportCells(j,i,1) = gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), -y, zoutersecs]));

            selcell = subset(gcells,  gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), -y, zoutersecs])));

            febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

            OuterSupportCells{j,i,1} = febsel;

            %             OuterSupportCells{j,i,1} = feblock_defor_ss_beam3 (struct ('gcells', selcell));

            %                     OuterSupportCells(j,i,3) = gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), y, zoutersecs]));


            selcell = subset(gcells, gcell_select(fens, gcells, struct('nearestto', [-x + outerelementlen/2 + ((i-1) * outerelementlen), y, zoutersecs])));

            febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

            OuterSupportCells{j,i,3} = febsel;

            %             OuterSupportCells{j,i,1} = feblock_defor_ss_beam3 (struct ('gcells', selcell));

        end

    end

end
        