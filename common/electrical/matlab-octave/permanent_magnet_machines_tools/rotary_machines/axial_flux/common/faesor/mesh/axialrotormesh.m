function [fens,gcells,maglabels,shaftlabels,outersep] = axialrotormesh(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tdiscsep, tausupp, nmodules, nmodulesupports, varargin)

    Options.ShaftAxialLayersPerM = 20;
    Options.DiscAxialLayersPerM = 50;
    Options.SupportAxialLayersPerM = 50;
    Options.BackIronRadialPointsPerM = 15;
    Options.MagnetRadialPointsPerM = 15;
    Options.CircumPoints = [];
    Options.CircumPointsPerM = 300;
    Options.MagnetLabels = [];
    Options.ShaftLabels = [];
    Options.NStages = 1;
    Options.Subset = nmodules;
    Options.MeshType = 'H20';
    
    Options = parseoptions(Options, varargin);

    if isempty(Options.MagnetLabels)
        Options.MagnetLabels = 1:Options.NStages+1;
    elseif isscalar(Options.MagnetLabels)
        Options.MagnetLabels = repmat(Options.MagnetLabels, 1, Options.NStages+1);
    elseif numel(Options.MagnetLabels) ~= Options.NStages+1
        error('MagnetLabels must be a sclalar or vector of values of length NStages+1, being one magnet label per rotor disc.')
    end
    
    maglabels = Options.MagnetLabels;
    
    if isempty(Options.ShaftLabels)
        Options.ShaftLabels = max(Options.MagnetLabels)+1:max(Options.MagnetLabels)+Options.NStages;
    elseif isscalar(Options.ShaftLabels)
        Options.ShaftLabels = repmat(Options.ShaftLabels, 1, Options.NStages);
    elseif numel(Options.ShaftLabels) ~= Options.NStages
        error('ShaftLabels must be a sclalar or vector of values of length NStages, being one shaft section label per rotor stage.')
    end
    
    shaftlabels = Options.ShaftLabels;
    
    if isscalar(tbi)
        % outer and inner discs are the same thickness
        tbi = [tbi, tbi];
    elseif numel(tbi) ~= 2
        error('The back iron thickness (tbi) must be a scalar or two-element vector.');
    end
    
    if any(tbi <= 0)
        error('tbi cannot be less than or equal to zero')
    end
    
    if ~samesize(Rs, Rbi, Rmi, Rmo, Rbo, tsuppb, tdiscsep, tausupp, nmodules, nmodulesupports)
        error('Rs, Rbi, Rmi, Rmo, Rbo, tsuppb, tdiscsep, tausupp, nmodules and nmodulesupports must all be scalars');
    else
        if isscalar(Rs)
            if any([Rs, Rbi, Rmi, Rmo, Rbo, tdiscsep, nmodules] <= 0)
                error('Rs, Rbi, Rmi, Rmo, Rbo, tdiscsep, and nmodules must all be greater than zero.')
            end
        else
            error('Rs, Rbi, Rmi, Rmo, Rbo, tsuppb, tdiscsep, tausupp, nmodules and nmodulesupports must all be scalar values')
        end
    end
    
    if any(diff([Rs, Rbi, Rmi, Rmo, Rbo]) <= 0)
        error('Rs, Rbi, Rmi, Rmo and Rbo must be greater than zero and each be bigger than the the variable earlier in the list i.e. Rs < Rbi < Rmi < Rmo < Rbo.')
    end
    
    if Options.Subset < 1
        Options.Subset = 1;
    else
        Options.Subset = round(Options.Subset);
    end
    
    outerdisclayers = axialrotormesh_layercalc(Options.DiscAxialLayersPerM, tbi(1), 1);
    innerdisclayers = axialrotormesh_layercalc(Options.DiscAxialLayersPerM, tbi(2), 1);
    supportlayers = axialrotormesh_layercalc(Options.SupportAxialLayersPerM, tsuppb, 1);
    shaftlayers = axialrotormesh_layercalc(Options.ShaftAxialLayersPerM, tdiscsep, 1);
    if isempty(Options.CircumPoints)
        Options.CircumPoints = ceil(Options.CircumPointsPerM * pi * 2 * Rbo);
    end
    
    % calculate the separation between the inner surfaces of the outer
    % discs
    outersep = ((Options.NStages-1) * tbi(2)) + (Options.NStages * tdiscsep);
    
    % Make the outer rotor which can have support fins
    [fens, gcells, modulexs, modulegap] = ...
        axialouterrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi(1), tsuppb, tausupp, ...
                                  nmodules, nmodulesupports, outersep, ...
                                  'disclayers', outerdisclayers, ...
                                  'supportlayers', supportlayers, ...
                                  'circ', Options.CircumPoints, ...
                                  'magl', Options.MagnetLabels([1, end]), ...
                                  'BackIronRadialPointsPerM', Options.BackIronRadialPointsPerM, ...
                                  'MagnetRadialPointsPerM', Options.MagnetRadialPointsPerM, ...
                                  'Subset', Options.Subset);

    if Options.NStages > 1
        
        % create the internal rotor stages
        [fens2, gcells2] = axialinnerrotorfinsmesh3d(...
            Rs, Rbi, Rmi, Rmo, Rbo, tbi(2), outersep, nmodules, ...
            Options.NStages-1, Options.MagnetLabels(2:end-1), innerdisclayers, modulexs,  ...
            'BackIronRadialPointsPerM', Options.BackIronRadialPointsPerM, ...
            'MagnetRadialPointsPerM', Options.MagnetRadialPointsPerM, ...
            'Subset', Options.Subset, ...
            'ModuleGap', modulegap);

        [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, 0);

        gcells = cat(gcells, gcells2);

    end

    % link the discs together
    [fens2,gcells2] = axialrotorfinlinkmesh3d(...
        Rs, Rbi, tbi(2), outersep, modulexs, nmodules, Options.NStages, ...
        shaftlayers, 'ShaftLabel', Options.ShaftLabels, ...
        'Subset', Options.Subset, ...
        'ModuleGap', modulegap);

    [fens, gcells, gcells2] = merge_meshes(fens, gcells, fens2, gcells2, 2*eps);

    gcells = cat(gcells, gcells2);
    
%     if ~isempty(Options.Subset)
%         
%         
%         
%         extent = 1.1 * (tbi(1) + tsuppb);
%         
%         [x,y] = pol2cart(linspace(0, Options.Subset*2*pi/nmodules, 10)', repmat(1.1*Rbo, 10,1));
% 
%         % selection region convex hull points
%         convhullpoints = [ 0, 0, -extent;
%                            x, y, repmat(-extent, size(x, 1), 1); 
%                            x, y, repmat(extent, size(x, 1), 1);  
%                            0, 0, extent; ];
% 
% %         C = convhull(convhullpoints);
% 
% %         trimesh(C, convhullpoints(:,1), convhullpoints(:,2), convhullpoints(:,3));
% 
%         gcells = subset(gcells, gcell_select(fens, gcells, struct('convhull', convhullpoints)));
% 
%         fens = subset(fens, fenode_select(fens, struct('convhull', convhullpoints)));
%         
%     end
           
    % change the type of mesh output if desired
    switch Options.MeshType
        
        case 'H8'
            
            % do nothing
            
        case 'H20'
            
            [fens,gcells] = H8_to_H20(fens,gcells);
            
        case 'H27'
            
            [fens,gcells] = H8_to_H27(fens,gcells);
            
        case 'H64'
            
            [fens,gcells] = H8_to_H64(fens,gcells);
            
        case 'T4'
            
            [fens,gcells] = H8_to_T4_random(fens,gcells);
            
        otherwise
            
            error('unrecognised mesh type: %s, expected H8, H20, H27, H64 or T4', Options.MeshType);
            
    end

end


function layers = axialrotormesh_layercalc(layersperm, thickness, minlayers)

    if nargin < 3
        minlayers = 1;
    end
    
     layers = max(minlayers, round(layersperm * thickness));

end