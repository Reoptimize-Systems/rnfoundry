function [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob(ypole, ymag, xmag, xbackiron, magsep, varargin)
% torusrotor2dfemmprob draws a torus type rotor of any number of stages and
% with several possible magnet arrangements
%
%

    Inputs.NStages = 1;
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
    
    % parse the inputs
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.MagnetGroup) && isempty(Inputs.MagnetSpaceGroup) && ...
            isempty(Inputs.BackIronGroup) && isempty(Inputs.OuterRegionGroup)
        Inputs.MagnetGroup = (1:Inputs.NStages + 1);
        Inputs.MagnetSpaceGroup = Inputs.MagnetGroup + Inputs.MagnetGroup(end);
        Inputs.BackIronGroup = Inputs.MagnetGroup;
        Inputs.OuterRegionGroup = [0,0];
    end
    
    if isempty(Inputs.MagnetGroup)
        Inputs.MagnetGroup = (1:Inputs.NStages + 1);
    end
    
    if isempty(Inputs.MagnetSpaceGroup)
        Inputs.MagnetSpaceGroup = Inputs.MagnetGroup + Inputs.MagnetGroup(end);
    end
    
    if isempty(Inputs.BackIronGroup)
        Inputs.BackIronGroup = Inputs.MagnetGroup;
    end
    
    if isempty(Inputs.OuterRegionGroup)
        Inputs.OuterRegionGroup = [0,0];
    end
    
    if isscalar(Inputs.MagnetGroup) 
        Inputs.MagnetGroup = repmat(Inputs.MagnetGroup, 1, Inputs.NStages+1);
    end
    
    if isscalar(Inputs.MagnetSpaceGroup) 
        Inputs.MagnetSpaceGroup = repmat(Inputs.MagnetSpaceGroup, 1, Inputs.NStages+1);
    end
    
    if isscalar(Inputs.BackIronGroup) 
        Inputs.BackIronGroup = repmat(Inputs.BackIronGroup, 1, Inputs.NStages+1);
    end
    
    if isscalar(Inputs.OuterRegionGroup) 
        Inputs.OuterRegionGroup = repmat(Inputs.OuterRegionGroup, 1, 2);
    end
    
    FemmProblem = Inputs.FemmProblem;
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(ypole, Inputs.Position, Inputs.FractionalPolePosition, Inputs.RotorAnglePosition);
    
    if Inputs.NStages == 1
        
        outermagsep = magsep(1);
        xmag = xmag(1);
        xbackiron = xbackiron(1);
        innerstagewidth = 0;
        
    else
        
        % do some error checking
        if isscalar(xmag)
            xmag = [xmag, xmag];
        end
        
        if isscalar(xbackiron)
            xbackiron = [xbackiron, xbackiron];
        end
        
        if isscalar(magsep)
            magsep = [magsep, magsep];
        end
        
        if numel(xmag) ~= 2
            error('xmag should be a scalar or 2 element vector.');
        end
        
        if numel(xbackiron) ~= 2
            error('xbackiron should be a scalar or 2 element vector.');
        end
        
        if numel(magsep) ~= 2
            error('magsep should be a scalar or 2 element vector.');
        end
        
        if anyalldims([xmag, xbackiron, magsep] < 0)
            error('All dimensions must be positive')
        end
        

        FemmProblem = axialfluxinnerrotor2dfemmprob(ypole, ymag, xmag(2), xbackiron(2), magsep(2), ...
                                                    'FemmProblem', FemmProblem, ...
                                                    'MagArrangement', Inputs.MagArrangement, ...
                                                    'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                    'BackIronMaterial', Inputs.BackIronMaterial, ...
                                                    'OuterRegionsMaterial', Inputs.OuterRegionsMaterial, ...
                                                    'MagnetSpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                    'MagnetGroup', Inputs.MagnetGroup(2:end-1), ...
                                                    'MagnetSpaceGroup', Inputs.MagnetSpaceGroup(2:end-1), ...
                                                    'BackIronGroup', Inputs.BackIronGroup(2:end-1), ...
                                                    'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
                                                    'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
                                                    'NInnerParts', Inputs.NStages - 1, ...
                                                    'Position', Inputs.Position, ...
                                                    'Tol', Inputs.Tol);
        
        outerstagewidth = magsep(1) + xmag(2) + xbackiron(2)/2;

        if xbackiron(2) > 0
            innerstagewidth = magsep(2) + 2*xmag(2) + xbackiron(2);
        else
            innerstagewidth = magsep(2) + xmag(2);
        end
        
        outermagsep = (Inputs.NStages - 2) * innerstagewidth + 2*outerstagewidth;
        
    end
    
    if xbackiron(1) == 0
        
        % set back iron material to be same as outer region material
        Inputs.BackIronMaterial = Inputs.OuterRegionsMaterial;
        
    end
    
    % draw the outer rotor parts, these will always be present
    FemmProblem = axialfluxouterrotor2dfemmprob(ypole, ymag, xmag(1), xbackiron(1), outermagsep, ...
                                                'FemmProblem', FemmProblem, ...
                                                'MagArrangement', Inputs.MagArrangement, ...
                                                'MagnetMaterial', Inputs.MagnetMaterial, ...
                                                'BackIronMaterial', Inputs.BackIronMaterial, ...
                                                'OuterRegionsMaterial', Inputs.OuterRegionsMaterial, ...
                                                'MagnetSpaceMaterial', Inputs.MagnetSpaceMaterial, ...
                                                'MagnetGroup', Inputs.MagnetGroup([1,end]), ...
                                                'MagnetSpaceGroup', Inputs.MagnetSpaceGroup([1,end]), ...
                                                'BackIronGroup', Inputs.BackIronGroup([1,end]), ...
                                                'OuterRegionGroup', Inputs.OuterRegionGroup, ...
                                                'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
                                                'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
                                                'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
                                                'Position', Inputs.Position, ...
                                                'Tol', Inputs.Tol);

end