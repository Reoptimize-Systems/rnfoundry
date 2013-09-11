function FemmProblem = axialfluxinnerrotor2dfemmprob(ypole, ymag, xmag, xbackiron, magsep, varargin)

    Inputs.NInnerParts = 1;
    Inputs.MagArrangement = 'NS';
    Inputs.FemmProblem = newproblem_mfemm('planar');
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetMaterial = 1;
    Inputs.BackIronMaterial = 1;
    Inputs.OuterRegionsMaterial = 1;
    Inputs.MagnetSpaceMaterial = 1;
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.BackIronGroup = [];
    Inputs.MagnetRegionMeshSize = -1;
    Inputs.BackIronRegionMeshSize = -1;
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.MagnetGroup) && isempty(Inputs.MagnetSpaceGroup) && isempty(Inputs.BackIronGroup)
        Inputs.MagnetGroup = (1:Inputs.NInnerParts) + 1;
        Inputs.MagnetSpaceGroup = Inputs.MagnetGroup;
        Inputs.BackIronGroup = Inputs.MagnetGroup;
    end
    
    if isscalar(Inputs.MagnetGroup) 
        Inputs.MagnetGroup = repmat(Inputs.MagnetGroup, 1, Inputs.NInnerParts);
    end
    
    if isscalar(Inputs.MagnetSpaceGroup) 
        Inputs.MagnetSpaceGroup = repmat(Inputs.MagnetSpaceGroup, 1, Inputs.NInnerParts);
    end
    
    if isscalar(Inputs.BackIronGroup) 
        Inputs.BackIronGroup = repmat(Inputs.BackIronGroup, 1, Inputs.NInnerParts);
    end
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(ypole, Inputs.Position, Inputs.FractionalPolePosition, Inputs.RotorAnglePosition);

    if xbackiron == 0
        stagewidth = magsep + xmag;
        nstagefactor = 1;
    else
        stagewidth = magsep + 2*xmag + xbackiron;
        nstagefactor = 2;
    end
    
    switch Inputs.MagArrangement

        case 'NN'

            MagDirections = repmat([180, 0], Inputs.NInnerParts * nstagefactor, 1);
            MagDirections(2:2:end,:) = fliplr(MagDirections(2:2:end,:));

        case 'NS'

            MagDirections = repmat([0, 180], Inputs.NInnerParts * nstagefactor, 1);

        otherwise

            error('ROTARY:axfluxrotor:badtype', 'Unknown magnet arrangement specification, should be NN or NS.')

    end
    
    stagepositions = (0:Inputs.NInnerParts-1) .* stagewidth;
    stagepositions = stagepositions - max(stagepositions)/2;
    
    FemmProblem = Inputs.FemmProblem;
    
    for i = 1:numel(stagepositions)
        
        % the rotor parts will be separated by 
        xoffset = stagepositions(i);        
                                                                           
        if xbackiron == 0
            
            % draw the outer rotor parts
            FemmProblem = wrappedrectmagaperiodic(FemmProblem, ypole, ...
                                                   ymag, xmag, xoffset, ...
                                                   Inputs.Position, 0, ...
                                                   'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                   'MagDirections', MagDirections(i,:), ...
                                                   'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                   'SpaceGroup', Inputs.MagnetSpaceGroup(i), ...
                                                   'MagnetGroup', Inputs.MagnetGroup(i), ...
                                                   'WrapperGroup', Inputs.BackIronGroup(i), ...
                                                   'Tol', Inputs.Tol, ...
                                                   'MeshSize', Inputs.MagnetRegionMeshSize);
                                                                       
        else

            % we must draw two sets of magnets and create back iron btween
            % them
            FemmProblem = wrappedrectmagaperiodic(FemmProblem, ypole, ...
                                                  ymag, xmag, xoffset - xmag/2 - xbackiron/2, ...
                                                  Inputs.Position, 0, ...
                                                  'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                  'MagDirections', MagDirections((2*i)-1,:), ...
                                                  'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                  'SpaceGroup', Inputs.MagnetSpaceGroup(i), ...
                                                  'MagnetGroup', Inputs.MagnetGroup(i), ...
                                                  'WrapperGroup', Inputs.BackIronGroup(i), ...
                                                  'Tol', Inputs.Tol, ...
                                                  'MeshSize', Inputs.MagnetRegionMeshSize);
                                              
            
                                                                       
            FemmProblem = wrappedrectmagaperiodic(FemmProblem, ypole, ...
                                                  ymag, xmag, xoffset + xmag/2 + xbackiron/2, ...
                                                  Inputs.Position, 0, ...
                                                  'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                  'MagDirections', MagDirections((2*i),:), ...
                                                  'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                  'SpaceGroup', Inputs.MagnetSpaceGroup(i), ...
                                                  'MagnetGroup', Inputs.MagnetGroup(i), ...
                                                  'WrapperGroup', Inputs.BackIronGroup(i), ...
                                                  'Tol', Inputs.Tol, ...
                                                  'MeshSize', Inputs.MagnetRegionMeshSize);
                                              
            backironcorners = [xoffset - xbackiron/2, 2*ypole; 
                               xoffset + xbackiron/2, 2*ypole; 
                               xoffset - xbackiron/2, 0; 
                               xoffset + xbackiron/2, 0; ];
             
            % get the node ids of the back iron corners
            nodeids = findnode_mfemm(FemmProblem, backironcorners);
                                               
            % add a new periodic boundary for the top and bottom of the
            % region
            FemmProblem = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Rect Mags Back Iron Periodic', 4);
        
            FemmProblem = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2), ...
                                    'BoundaryMarker', FemmProblem.BoundaryProps(end).Name, ...
                                    'InGroup', Inputs.BackIronGroup(i));
                                
            FemmProblem = addsegments_mfemm(FemmProblem, nodeids(3), nodeids(4), ...
                                    'BoundaryMarker', FemmProblem.BoundaryProps(end).Name, ...
                                    'InGroup', Inputs.BackIronGroup(i));   
                                
            % Add block label for the back iron
            labelloc = rectcentre(backironcorners(1,:), backironcorners(4,:));
            
            FemmProblem = addblocklabel_mfemm(FemmProblem, labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.BackIronMaterial).Name, ...
                                            'MaxArea', Inputs.BackIronRegionMeshSize, ...
                                            'InGroup', Inputs.BackIronGroup(i));
            

        end
        
    end

end