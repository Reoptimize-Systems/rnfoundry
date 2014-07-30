function [FemmProblem, magcornernodeids, linktb] = radialfluxrotor2dfemmprob(thetapole, thetamag, rmag, rbackiron, drawnrotors, rrotor, varargin)
% adds the outermost rotor parts of a radial flux machine to a FemmProblem
% Structure

    Inputs.MagArrangement = 'NS';
    Inputs.PolarisationType = 'constant';
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
    Inputs.OuterRegionGroup = [];
    Inputs.MagnetRegionMeshSize = -1;
    Inputs.BackIronRegionMeshSize = -1;
    Inputs.OuterRegionsMeshSize = [-1, -1];
    Inputs.Tol = 1e-5;
    Inputs.DrawOuterRegions = true;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.MagnetGroup) && isempty(Inputs.MagnetSpaceGroup) ...
            && isempty(Inputs.BackIronGroup) && isempty(Inputs.OuterRegionGroup)
        
        if drawnrotors(1) && drawnrotors(2)
            Inputs.MagnetGroup = [1,2];
        else
            Inputs.MagnetGroup = [1,1];
        end
        Inputs.MagnetSpaceGroup = Inputs.MagnetGroup;
        Inputs.BackIronGroup = Inputs.MagnetGroup;
        Inputs.OuterRegionGroup = Inputs.MagnetGroup;
        
    end
    
    if isscalar(Inputs.MagnetGroup) 
        Inputs.MagnetGroup = repmat(Inputs.MagnetGroup, 1, 2);
    end
    
    if isscalar(Inputs.MagnetSpaceGroup) 
        Inputs.MagnetSpaceGroup = repmat(Inputs.MagnetSpaceGroup, 1, 2);
    elseif isempty(Inputs.OuterRegionGroup)
        Inputs.MagnetSpaceGroup = zeros(1, 2);
    end
    
    if isscalar(Inputs.BackIronGroup) 
        Inputs.BackIronGroup = repmat(Inputs.BackIronGroup, 1, 2);
    elseif isempty(Inputs.BackIronGroup)
        Inputs.BackIronGroup = Inputs.MagnetGroup;
    end
    
    if isscalar(Inputs.OuterRegionGroup) 
        Inputs.OuterRegionGroup = repmat(Inputs.OuterRegionGroup, 1, 2);
    elseif isempty(Inputs.OuterRegionGroup)
        Inputs.OuterRegionGroup = zeros(1, 2);
    end
    
    if numel(Inputs.OuterRegionsMeshSize) ~= 2
        error('AXIALFLUXOUTERROTOR2DFEMMPROB:badoutermeshspec', ...
            'OuterRegionsMeshSize must be a two element vector')
    end
    
    if strncmpi(Inputs.PolarisationType, 'constant', 1);
        
        switch Inputs.MagArrangement
            
            case 'NN'
                
                innerMagDirections = {180, 0};
                
            case 'NS'
                
                innerMagDirections = {0, 180};
                
            otherwise
                
                error('ROTARY:radiafluxrotor:badtype', ...
                    'Unknown magnet arrangement specification, should be NN or NS.')
                
        end
    
    elseif strncmpi(Inputs.PolarisationType, 'radial', 1);
        
        switch Inputs.MagArrangement

            case 'NN'

                innerMagDirections = {'theta', 'theta+180'};

            case 'NS'

                innerMagDirections = {'theta+180', 'theta'};

            otherwise

                error('ROTARY:radiafluxrotor:badtype', ...
                      'Unknown magnet arrangement specification, should be NN or NS.')

        end
    
    else
        error('Unknown magnet polarisation type.');
    end

    if drawnrotors(1) && ~drawnrotors(2)
        roffset = [rrotor + rmag/2, 0];
    elseif ~drawnrotors(1) && drawnrotors(2)
        roffset = [0, rrotor - rmag/2];
    elseif drawnrotors(1) && drawnrotors(2)
        roffset = [rrotor(1)+rmag/2, rrotor(2)-rmag/2];
    else
        error('ROTARY:radiafluxrotor:badtype', ...
              'At least one rotor part must be drawn.');
    end
  
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(thetapole, ...
                                     Inputs.Position, ...
                                     Inputs.FractionalPolePosition, ...
                                     Inputs.RotorAnglePosition);
    
	FemmProblem = Inputs.FemmProblem;
    
    % Add a proscribed A boundary we will use for the outer regions
    % boundary types
    [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Outer Boundary Pros A', 0, 'A0', 0);

    linktb1 = false;
    linktb2 = false;
    
    if drawnrotors(1)
        % an outer rotor
        
        % there will be three wrapper thicknesses on either side, the back
        % iron, and some air region
        outerwrapperthickness = [ 0, rbackiron; ];
        
        if Inputs.DrawOuterRegions
            outerwrapperthickness = [ outerwrapperthickness;
                                      0, 2*(rbackiron + rmag);
                                      0, 4*(rbackiron + rmag) ];
        end

        % draw the outer rotor parts
        [FemmProblem, ~, ~, ~, nodeids1, linktb1] = wrappedannularsecmagaperiodic(FemmProblem, thetapole, ...
                                                                       thetamag, rmag, roffset(1), ...
                                                                       Inputs.Position, outerwrapperthickness, ...
                                                                       'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                                       'MagDirections', fliplr(innerMagDirections), ...
                                                                       'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                                       'SpaceGroup', Inputs.MagnetSpaceGroup(1), ...
                                                                       'MagnetGroup', Inputs.MagnetGroup(1), ...
                                                                       'WrapperGroup', Inputs.BackIronGroup(1), ...
                                                                       'Tol', Inputs.Tol, ...
                                                                       'MeshSize', Inputs.MagnetRegionMeshSize);

        if Inputs.DrawOuterRegions && ~linktb1
            % the last arc segment added will be the outermost segment, so give it the
            % proscribed A boundary
            FemmProblem.ArcSegments(end).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
            FemmProblem.ArcSegments(end-1).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
        end
    end
    
    if drawnrotors(2)
        % an inner rotor
        
        % mirror the wrapper thicknesses on the other side
        innerwrapperthickness = [ rbackiron, 0; ];
        if Inputs.DrawOuterRegions
            innerwrapperthickness = [ innerwrapperthickness;
                                      (roffset(2) * 0.2), 0;
                                      (roffset(2)-(roffset(2) * 0.2)) * 0.5, 0 ];
        end
        
        [FemmProblem, ~, ~, ~, nodeids2, linktb2] = wrappedannularsecmagaperiodic(FemmProblem, thetapole, ...
                                                                       thetamag, rmag, roffset(2), ...
                                                                       Inputs.Position, innerwrapperthickness, ...
                                                                       'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                                       'MagDirections', innerMagDirections, ...
                                                                       'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                                       'SpaceGroup', Inputs.MagnetSpaceGroup(2), ...
                                                                       'MagnetGroup', Inputs.MagnetGroup(2), ...
                                                                       'WrapperGroup', Inputs.BackIronGroup(2), ...
                                                                       'Tol', Inputs.Tol, ...
                                                                       'MeshSize', Inputs.MagnetRegionMeshSize);

        if Inputs.DrawOuterRegions && ~linktb2
            % Again, the last arc segment added should be the outermost boundary
            % segment this time on the rhs, so give it the proscribed boundary condition
            FemmProblem.ArcSegments(end).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
            FemmProblem.ArcSegments(end-1).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
        end

    end
    
    % Add the back iron block labels
    if ~isempty(Inputs.BackIronMaterial)
        
        if drawnrotors(1)
            [x,y] = pol2cart(thetapole, roffset(1) + rmag/2 + rbackiron/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                                        'BlockType', FemmProblem.Materials(Inputs.BackIronMaterial).Name, ...
                                                        'MaxArea', Inputs.BackIronRegionMeshSize, ...
                                                        'InGroup', Inputs.BackIronGroup(1));
        end
                                                
        if drawnrotors(2)
            [x,y] = pol2cart(thetapole, roffset(2) - rmag/2 - rbackiron/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                                        'BlockType', FemmProblem.Materials(Inputs.BackIronMaterial).Name, ...
                                                        'MaxArea', Inputs.BackIronRegionMeshSize, ...
                                                        'InGroup', Inputs.BackIronGroup(2));
        end
        
    end
             
    if Inputs.DrawOuterRegions && ~isempty(Inputs.OuterRegionsMaterial)
        
        % Add the outer region block lables
        if drawnrotors(1)
            
            [x,y] = pol2cart( thetapole,  roffset(1) + rmag/2 + rbackiron + outerwrapperthickness(2,2)/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                              'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                              'MaxArea', Inputs.OuterRegionsMeshSize(1), ...
                                              'InGroup', Inputs.OuterRegionGroup(1));
                                                    
            [x,y] = pol2cart(thetapole, roffset(1) + rmag/2 + rbackiron + outerwrapperthickness(2,2) + outerwrapperthickness(3,2)/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                              'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                              'MaxArea', Inputs.OuterRegionsMeshSize(2), ...
                                              'InGroup', Inputs.OuterRegionGroup(1));    

        end
        
        if drawnrotors(2) 
            
            [x,y] = pol2cart(thetapole, roffset(2) - rmag/2 - rbackiron - innerwrapperthickness(2,1)/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                              'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                              'MaxArea', Inputs.OuterRegionsMeshSize(1), ...
                                              'InGroup', Inputs.OuterRegionGroup(2));

            [x,y] = pol2cart(thetapole, roffset(2) - rmag/2 - rbackiron - innerwrapperthickness(2,1) - innerwrapperthickness(3,1)/2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, ...
                                              'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                              'MaxArea', Inputs.OuterRegionsMeshSize(2), ...
                                              'InGroup', Inputs.OuterRegionGroup(2));
        end
        
    end

    if drawnrotors(1) && drawnrotors(2)
        magcornernodeids = [nodeids1([1,end-1]), nodeids2([2,end])];
    elseif drawnrotors(1)
        magcornernodeids = nodeids1([1,end-1]);
    elseif drawnrotors(2)
        magcornernodeids = nodeids2([2,end]);
    end
    
    linktb = linktb1 || linktb2;
    
end