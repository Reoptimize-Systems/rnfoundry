function [FemmProblem, nodes, nodeids, links] = rectmagaperiodic(FemmProblem, ypole, ymag, xmag, xoffset, pos, varargin)
% planarrectmagaperiodic: 

    Inputs.MagDirections = {0, 180};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpaceMaterial = 1;
    Inputs.SpaceGroup = 0;
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;

    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isnumeric(Inputs.MagDirections) && isvector(Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
    % get the number of existing nodes, segments, boundaries etc. if any
    elcount = elementcount_mfemm(FemmProblem);
    
    % add a periodic boundary for the top and bottom of the magnets
    % region
    [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, ...
                                                    'Rect Mags Periodic', 4);
    
    % construct the segments and label locations for the magnet regions
    [nodes, nodeids, links, rectcentres, spacecentres] = ...
        rectregionsyperiodic(xmag, ymag, (ypole-ymag), xoffset, pos, Inputs.Tol, elcount.NNodes);

    % add all the nodes to the problem
    for i = 1:size(nodes,1)
        FemmProblem = addnodes_mfemm(FemmProblem, nodes(i,1), nodes(i,2), 'InGroup', Inputs.MagnetGroup);
    end
    
    % Periodic boundary at bottom
    FemmProblem = addsegments_mfemm(FemmProblem, links(1,1), links(1,2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
        'InGroup', Inputs.MagnetGroup);
        
    % Add all the segments except the top and bottom which will have
    % boundary properties we must set manually
    FemmProblem = addsegments_mfemm(FemmProblem, links(2:end-1,1), links(2:end-1,2), ...
        'InGroup', Inputs.MagnetGroup);

    % Periodic boundary at top
    FemmProblem = addsegments_mfemm(FemmProblem, links(end,1), links(end,2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
        'InGroup', Inputs.MagnetGroup);
    
    % Now add the labels
    
    if ~isempty(Inputs.MagnetMaterial)
        
        % add the magnet labels
        for i = 1:size(rectcentres, 1)

            if rectcentres(i,3)

                FemmProblem = addblocklabel_mfemm(FemmProblem, rectcentres(i,1), rectcentres(i,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{2}, ...
                                                'InGroup', Inputs.MagnetGroup);

            else

                FemmProblem = addblocklabel_mfemm(FemmProblem, rectcentres(i,1), rectcentres(i,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{1}, ...
                                                'InGroup', Inputs.MagnetGroup);

            end

        end
    
    end
    
    if ~isempty(Inputs.SpaceMaterial)
        
        % Add the other lables
        for i = 1:size(spacecentres, 1)

            FemmProblem = addblocklabel_mfemm(FemmProblem, spacecentres(i,1), spacecentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.SpaceMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpaceGroup);

        end
    
    end

    
end