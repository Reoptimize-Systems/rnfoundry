function FemmProblem = axialfluxouterrotor2dfemmprob(ypole, ymag, xmag, xbackiron, magsep, varargin)
% axialfluxouterrotor2dfemmprob: adds the outermost rotor parts of an axial
% flux machine to a FemmProblem Structure

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
    Inputs.OuterRegionGroup = [];
    Inputs.MagnetRegionMeshSize = -1;
    Inputs.BackIronRegionMeshSize = -1;
    Inputs.OuterRegionsMeshSize = [-1, -1];
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.MagnetGroup) && isempty(Inputs.MagnetSpaceGroup) && ...
       isempty(Inputs.BackIronGroup) && isempty(OuterRegionGroup)
   
        Inputs.MagnetGroup = [0,0];
        Inputs.MagnetSpaceGroup = Inputs.MagnetGroup;
        Inputs.BackIronGroup = Inputs.MagnetGroup;
        Inputs.OuterRegionGroup = Inputs.MagnetGroup;
        
    end
    
    if isscalar(Inputs.MagnetGroup) 
        Inputs.MagnetGroup = repmat(Inputs.MagnetGroup, 1, 2);
    end
    
    if isscalar(Inputs.MagnetSpaceGroup) 
        Inputs.MagnetSpaceGroup = repmat(Inputs.MagnetSpaceGroup, 1, 2);
    end
    
    if isscalar(Inputs.BackIronGroup) 
        Inputs.BackIronGroup = repmat(Inputs.BackIronGroup, 1, 2);
    end
    
    if isscalar(Inputs.OuterRegionGroup) 
        Inputs.OuterRegionGroup = repmat(Inputs.OuterRegionGroup, 1, 2);
    end
    
    if numel(Inputs.OuterRegionsMeshSize) ~= 2
        error('AXIALFLUXOUTERROTOR2DFEMMPROB:badoutermeshspec', ...
            'OuterRegionsMeshSize must be a two element vector')
    end
    
    switch Inputs.MagArrangement
        
        case 'NN'
            
            rightMagDirections = [180, 0];
            
        case 'NS'
            
            rightMagDirections = [0, 180];
            
        otherwise
            
            error('ROTARY:axfluxrotor:badtype', 'Unknown magnet arrangement specification, should be NN or NS.')
            
    end
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(ypole, Inputs.Position, Inputs.FractionalPolePosition, Inputs.RotorAnglePosition);
    
%     if magsep <= xmag;
%        error('ROTARY:axfluxrotor:badpos', 'Magnet separation is less than magnet width, magnets will overlap in centre'); 
%     end
     
    % the rotor parts will be separated by 
    xoffset = magsep/2 + xmag/2;
     
    % there will be three wrapper thicknesses on either side, the back
    % iron, and some air region
    wrapperthickness = [ xbackiron, 0;
                         (xbackiron + xmag + magsep/2), 0;
                         4*(xbackiron + xmag + magsep/2), 0 ];
    
	FemmProblem = Inputs.FemmProblem;
                     
    % Add a proscribed A boundary we will use for the outer regions
    % boundary types
    [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Outer Boundary Pros A', 0, 'A0', 0);
    
    % draw the outer rotor parts
    [FemmProblem, nodes, nodeids, links] = wrappedrectmagaperiodic(FemmProblem, ypole, ...
                                                                   ymag, xmag, -xoffset, ...
                                                                   Inputs.Position, wrapperthickness, ...
                                                                   'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                                   'MagDirections', [0, 180], ...
                                                                   'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                                   'SpaceGroup', Inputs.MagnetSpaceGroup(1), ...
                                                                   'MagnetGroup', Inputs.MagnetGroup(1), ...
                                                                   'WrapperGroup', Inputs.BackIronGroup(1), ...
                                                                   'Tol', Inputs.Tol, ...
                                                                   'MeshSize', Inputs.MagnetRegionMeshSize);
                                                               
    % the last segment added will be the outermost segment, so give it the
    % proscribed A boundary
    FemmProblem.Segments(end).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
    
    % mirror the wrapper thicknesses on the other side
    wrapperthickness = fliplr(wrapperthickness);
    
    [FemmProblem, nodes, nodeids, links] = wrappedrectmagaperiodic(FemmProblem, ypole, ...
                                                                   ymag, xmag, xoffset, ...
                                                                   Inputs.Position, wrapperthickness, ...
                                                                   'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                                   'MagDirections', rightMagDirections, ...
                                                                   'SpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                                   'SpaceGroup', Inputs.MagnetSpaceGroup(2), ...
                                                                   'MagnetGroup', Inputs.MagnetGroup(2), ...
                                                                   'WrapperGroup', Inputs.BackIronGroup(2), ...
                                                                   'Tol', Inputs.Tol, ...
                                                                   'MeshSize', Inputs.MagnetRegionMeshSize);
    
	% Again, the last segment added should be the outermost boundary
	% segment this time on the rhs, so give it the proscribed boundary condition
    FemmProblem.Segments(end).BoundaryMarker = FemmProblem.BoundaryProps(boundind).Name;
    
    
    % Add the back iron block labels
    
    if ~isempty(Inputs.BackIronMaterial)
        
        FemmProblem = addblocklabel_mfemm(FemmProblem, -(magsep/2 + xmag + wrapperthickness(1,2)/2), ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.BackIronMaterial).Name, ...
                                                    'MaxArea', Inputs.BackIronRegionMeshSize, ...
                                                    'InGroup', Inputs.BackIronGroup(1));
                                                
        FemmProblem = addblocklabel_mfemm(FemmProblem, magsep/2 + xmag + wrapperthickness(1,2)/2, ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.BackIronMaterial).Name, ...
                                                    'MaxArea', Inputs.BackIronRegionMeshSize, ...
                                                    'InGroup', Inputs.BackIronGroup(2));
    end
             
    if ~isempty(Inputs.OuterRegionsMaterial)
        
    	% Add the outer region block lables
        FemmProblem = addblocklabel_mfemm(FemmProblem, magsep/2 + xmag + wrapperthickness(1,2) + wrapperthickness(2,2)/2, ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                                    'MaxArea', Inputs.OuterRegionsMeshSize(1), ...
                                                    'InGroup', Inputs.OuterRegionGroup(2));                                                               

        FemmProblem = addblocklabel_mfemm(FemmProblem, -(magsep/2 + xmag + wrapperthickness(1,2) + wrapperthickness(2,2)/2), ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                                    'MaxArea', Inputs.OuterRegionsMeshSize(1), ...
                                                    'InGroup', Inputs.OuterRegionGroup(1));

        FemmProblem = addblocklabel_mfemm(FemmProblem, magsep/2 + xmag + wrapperthickness(1,2) + wrapperthickness(2,2) + wrapperthickness(3,2)/2, ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                                    'MaxArea', Inputs.OuterRegionsMeshSize(2), ...
                                                    'InGroup', Inputs.OuterRegionGroup(2));                                                               

        FemmProblem = addblocklabel_mfemm(FemmProblem, -(magsep/2 + xmag + wrapperthickness(1,2) + wrapperthickness(2,2) + wrapperthickness(3,2)/2), ypole, ...
                                                    'BlockType', FemmProblem.Materials(Inputs.OuterRegionsMaterial).Name, ...
                                                    'MaxArea', Inputs.OuterRegionsMeshSize(2), ...
                                                    'InGroup', Inputs.OuterRegionGroup(1));
                                            
    end

end