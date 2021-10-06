function [faeprob, innersurfacelabel] = outerstructuremesh_RF(faeprob, Rsii, Rsoi, Rsoo, ls, lst, discends, varargin)
% creates a structured FAESOR mesh of the inner structure of a rotary
% machine
%
% Syntax
% 
% [fens,gcells] = outerstructuremesh_RF(faeprob, Rsi, Rso, ls)
% [fens,gcells] = outerstructuremesh_RF(..., forcelabel)
%
% Input
%
% Output
%
% fens - FAESOR fenodeset object containing the nodes of the rotor fin mesh
%
% gcells - FAESOR gcellset_H20 oject containing the H20 type cells making
%   up the rotor fin mesh
%
%

%     options.BackIronRadialPointsPerM = 15;
%     options.MagnetRadialPointsPerM = 15;
%     options.Subset = nmodules;
%     if nmodules == 1
%         options.ModuleGap = 0;
%     else
%         options.ModuleGap = 2*pi/(100*nmodules);
%     end
    
%     options = parseoptions(options, varargin);
    
    Options.AxialLayersPerM = 3;
    Options.RadialLayersPerM = 20;
    Options.RadialDiscLayersPerM = 15;
    Options.CircumPointsPerM = 20;
    Options.InnerSurfaceLabel = 2;
    Options.MeshType = 'H8';
    
    Options = parseoptions(Options, varargin);
    
    innersurfacelabel = Options.InnerSurfaceLabel;
    layers = max(5, ceil(Options.AxialLayersPerM * ls));
    ncircumpoints = max(5, ceil(Options.CircumPointsPerM * (2 * pi * mean([Rsoo, Rsoi]))));
    nradpoints = max(5, ceil(Options.RadialLayersPerM * (Rsoo - Rsoi)));
    disclayers = max(5, ceil(Options.RadialDiscLayersPerM * (Rsoo - Rsii)));
    
    mergetol = 2*eps;
    
    % 1. make set of Q4 quads in a ladder shape of height 2 pi that can be
    %    morphed into a cylinder later

    % construct rungs of ladder (the spokes)

    % get the radial positions 
    xsp = linspace(0, 2*pi, ncircumpoints);

    ysp = linspace(Rsoi, Rsoo, nradpoints);

    % create a set of 2D cells we can extrude into the 3D mesh
    [faeprob.fens, faeprob.gcells] = Q4_blockx(xsp, ysp, []);
        
    % label parts of interest
    
    % The inner surface cells
    celllist = gcell_select(faeprob.fens, faeprob.gcells, struct ('box', [-2*pi, 2*pi, ysp(1), ysp(2)], 'inflate', mergetol));
    
    % set the cells to have the inner surface label
    faeprob.gcells = setlabels(faeprob.gcells, Options.InnerSurfaceLabel, celllist);
    
    % extrude the mesh
    [faeprob.fens, faeprob.gcells] = Q4_extrude(faeprob.fens, faeprob.gcells, layers, ls, true);
    
    % add the supporting disc
    
    if discends(1)
        
        % get the radial positions 
        xsp = linspace(0, 2*pi, ncircumpoints);

        ysp = linspace(Rsii, Rsoo, nradpoints);

        % create a set of 2D cells we can extrude into the 3D mesh
        [disc1fens, disc1gcells] = Q4_blockx(xsp, ysp, []);

        % extrude the mesh
        [disc1fens, disc1gcells] = Q4_extrude(disc1fens, disc1gcells, disclayers, -lst, true);

        % put into cell arrays
        disc1fens = {disc1fens};
        disc1gcells = {disc1gcells};
        
    else
        disc1fens = {};
        disc1gcells = {};
    end
    
    if discends(2)
        
        % get the radial positions 
        xsp = linspace(0, 2*pi, ncircumpoints);

        ysp = linspace(Rsii, Rsoo, nradpoints);

        % create a set of 2D cells we can extrude into the 3D mesh
        [disc2fens, disc2gcells] = Q4_blockx(xsp, ysp, []);

        % extrude the mesh
        [disc2fens, disc2gcells] = Q4_extrude(disc2fens, disc2gcells, disclayers, lst, true);
        
        % shift mesh to the end of the rotor
        xyz = get (disc2fens,'xyz');
        
        xyz(:,3) = xyz(:,3) + ls;
        
        disc2fens = set(disc2fens, 'xyz', xyz);

        % put into cell arrays
        disc2fens = {disc2fens};
        disc2gcells = {disc2gcells};
        
    else
        disc2fens = {};
        disc2gcells = {};
    end
    
    [faeprob.fens,gcellsa] = merge_n_meshes( [ {faeprob.fens}, disc1fens, disc2fens ], ...
                                             [ {faeprob.gcells}, disc1gcells, disc2gcells ], ...
                                             mergetol );

	faeprob.gcells = cat(gcellsa{:});
    
    % get the locations of the nodes in the mesh
    xyz = get (faeprob.fens,'xyz');

    [x,y] = pol2cart(-xyz(:,1),xyz(:,2));

    xyz = [x, y, xyz(:,3)];

%         % shift to match rotor disc positon which is located at the origin
%         xyz(:,3) = xyz(:,3) - (tbi/2);

    % replace the coordinates of the nodes
    faeprob.fens = set(faeprob.fens, 'xyz', xyz);

    % Merge together all the nodes within a tolerance
    [faeprob.fens, faeprob.gcells] = merge_nodes(faeprob.fens, faeprob.gcells, 10 * mergetol);
    
    % change the type of mesh output if desired
    switch Options.MeshType
        
        case 'H8'
            
            % do nothing
            
        case 'H20'
            
            [faeprob.fens,faeprob.gcells] = H8_to_H20(faeprob.fens,faeprob.gcells);
            
        case 'H27'
            
            [faeprob.fens,faeprob.gcells] = H8_to_H27(faeprob.fens,faeprob.gcells);
            
        case 'H64'
            
            [faeprob.fens,faeprob.gcells] = H8_to_H64(faeprob.fens,faeprob.gcells);
            
        case 'T4'
            
            [faeprob.fens,faeprob.gcells] = H8_to_T4_random(faeprob.fens,faeprob.gcells);
            
        otherwise
            
            error('unrecognised mesh type: %s, expected H8, H20, H27, H64 or T4', Options.MeshType);
            
    end

end
