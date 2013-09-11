function [FemmProblem, nodes, nodeids, links, magblockinds] = annularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos, varargin)
% generates a periodic magnets containing region with an annular sector
% shape (i.e. a region bounded by two arcs and two straight lines, or a
% sector of an annulus) suitible for modelling radial flux machines.
%
% Syntax
%
% [FemmProblem, nodes, nodeids, links] = ...
%           annularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, roffset, pos)
% [...] = annularsecmagaperiodic(..., 'Paramter', Value)
%
% Input
%
% 

    Inputs.MagDirections = {'theta', 'theta+180'};
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
    
    % add a periodic boundary for the edges of the magnets region
    [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, ...
                                                    'Rect Mags Periodic', 4);
    
    % construct the segments and label locations for the magnet regions
    % first in a linear fashoin which we will manipulate into the real
    % arced shape by modifying the node locations
    [nodes, nodeids, links, rectcentres, spacecentres] = ...
        rectregionsyperiodic(rmag, thetamag, (thetapole-thetamag), roffset, pos, Inputs.Tol, elcount.NNodes);

    % get the vertical links by finding those links where the difference in
    % y coordinates of the link nodes is not zero, these links must be made
    % into arc segments
    vertlinks = links(abs(diff( [nodes(links(:,1)+1-elcount.NNodes,1), nodes(links(:,2)+1-elcount.NNodes,1)], 1, 2 )) < Inputs.Tol,:);
    angles = diff( [nodes(vertlinks(:,1)+1-elcount.NNodes,2), nodes(vertlinks(:,2)+1-elcount.NNodes,2)], 1, 2);
    
    % get the horizontal links, these will be segments
    horizlinks = links(abs(diff( [nodes(links(:,1)+1-elcount.NNodes,1), nodes(links(:,2)+1-elcount.NNodes,1)], 1, 2 )) >= Inputs.Tol,:);
    
    % transform the node locations to convert the rectangulr region to the
    % desired arced region 
    [nodes(:,1), nodes(:,2)] = pol2cart(nodes(:,2), nodes(:,1));
    [rectcentres(:,1), rectcentres(:,2)] = pol2cart(rectcentres(:,2),rectcentres(:,1));
    [spacecentres(:,1), spacecentres(:,2)] = pol2cart(spacecentres(:,2),spacecentres(:,1));
    
    for i = 1:size(vertlinks,1)
        if angles(i) < 0
            vertlinks(i,:) = fliplr(vertlinks(i,:));
            angles(i) = abs(angles(i));
        end
    end
    
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
    FemmProblem = addsegments_mfemm(FemmProblem, horizlinks(2:end-1,1), horizlinks(2:end-1,2), ...
        'InGroup', Inputs.MagnetGroup);
    
    % Add all the arc segments
    FemmProblem = addarcsegments_mfemm(FemmProblem, vertlinks(:,1), vertlinks(:,2), rad2deg(angles), ...
        'InGroup', Inputs.MagnetGroup);

    % Periodic boundary at top
    FemmProblem = addsegments_mfemm(FemmProblem, links(end,1), links(end,2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name, ...
        'InGroup', Inputs.MagnetGroup);
    
    % Now add the labels
    
    if ~isempty(Inputs.MagnetMaterial)
        
        magblockinds = [];
        
        % add the magnet labels
        for i = 1:size(rectcentres, 1)

            if rectcentres(i,3)

                [FemmProblem, magblockinds(end+1,1)] = addblocklabel_mfemm(FemmProblem, rectcentres(i,1), rectcentres(i,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{2}, ...
                                                'InGroup', Inputs.MagnetGroup);
                
            else

                [FemmProblem, magblockinds(end+1,1)] = addblocklabel_mfemm(FemmProblem, rectcentres(i,1), rectcentres(i,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                                'MaxArea', Inputs.MeshSize, ...
                                                'MagDir', Inputs.MagDirections{1}, ...
                                                'InGroup', Inputs.MagnetGroup);

            end
            
            magblockinds(end, 2) = rectcentres(i,3);

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