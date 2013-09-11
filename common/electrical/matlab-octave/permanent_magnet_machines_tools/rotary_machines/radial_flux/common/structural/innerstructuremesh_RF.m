function [faeprob, labels] = innerstructuremesh_RF(faeprob, Rsii, Rsio, Rsoi, Rsoo, thetasp, ls, nspokes, varargin)
% creates a structured FAESOR mesh of the inner structure of a rotary
% machine
%
% Syntax
% 
% [fens,gcells] = innerstructuremesh_RF((Rsii, Rsio, Rsoi, Rsoo, nspokes, layers, circumpoints)
% [fens,gcells] = innerstructuremesh_RF((..., forcelabel)
%
% Input
%
% Rs - inner radius of rotor shaft 
% 
% Rbi - outer radius of rotor shaft, and inner radius of the rotor fin back
%   iron
% 
% Rmi - inner radius of magnet positions
% 
% Rmo - outer radius of magnet positions
% 
% Rbo - outer rotor fin radius
% 
% tbi - thickness of the fin
% 
% disclayers - number of disclayers the mesh will be made up of (before refinement)
%   in the same dimension as tbi, the rotor fin thickness.
% 
% circumpoints - number of circumferential points the mesh will be made up
%   of (before refinement)
%
% forcelabel - (optional) label for cells in the magnet region
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
    Options.RadialInnerShaftLayersPerM = 20;
    Options.RadialOuterShaftLayersPerM = 20;
    Options.RadialSpokeLayersPerM = 10;
    Options.CircumSpokeLayers = 3;
    Options.CircumPointsPerM = 20;
    Options.OuterSurfaceLabel = 1;
    Options.MeshType = 'H8';
    
    Options = parseoptions(Options, varargin);
    
    labels.outersurface = Options.OuterSurfaceLabel;
    layers = max(2, ceil(Options.AxialLayersPerM * ls));
    ncircumpoints = max(5, ceil(Options.CircumPointsPerM * (2 * pi * mean([Rsoo, Rsoi]))));
    nradsipoints = max(5, ceil(Options.RadialInnerShaftLayersPerM * (Rsoo - Rsoi)));
    nradsppoints = max(5, ceil(Options.RadialSpokeLayersPerM * (Rsoi - Rsio)));
    nradsopoints = max(5, ceil(Options.RadialOuterShaftLayersPerM * (Rsio - Rsii)));
    
    mergetol = 100*eps(Rsoo);
    
    if nspokes > 1
    % to make spoked mesh,
    
        % 1. make set of Q4 quads in a ladder shape of height 2 pi that can be
        %    morphed into a cylinder later

        % construct rungs of ladder (the spokes)

        % get the positions of the cells in a spoke and shift to the first
        % spoke position
        xsp = linspace(-thetasp/2, thetasp/2, Options.CircumSpokeLayers) + (2*pi/(2*nspokes));

        tempxsp = xsp;
        
        ysp = [ linspace(Rsii, Rsio, nradsipoints), ...
                linspace(Rsio, Rsoi, nradsppoints), ...
                linspace(Rsoi, Rsoo, nradsopoints) ];

        % remove the duplicate points from the radial positions
        ysp = unique(ysp);

        spcenters = linspace(0, 2*pi, nspokes+1);

        spcenters = spcenters(1:end-1) + spcenters(2)/2;
        
        % create a set of 2D cells we can extrude into the 3D mesh
        [faeprob.fens, faeprob.gcells] = Q4_blockx(tempxsp, ysp, []);

        % create the links between the start and end
        [ilinkfens,ilinkgcells] = Q4_blockx(linspace(tempxsp(1) - (2*pi/(2*nspokes)) + thetasp/2, tempxsp(1), ceil(ncircumpoints/nspokes)), ysp(1:nradsipoints), []);

        [faeprob.fens,faeprob.gcells,ilinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, ilinkfens, ilinkgcells, 0);
        faeprob.gcells = cat(faeprob.gcells, ilinkgcells);

        [olinkfens,olinkgcells] = Q4_blockx(linspace(tempxsp(1) - (2*pi/(2*nspokes)) + thetasp/2, tempxsp(1), ceil(ncircumpoints/nspokes)), ysp(end-nradsopoints+1:end), []);

        [faeprob.fens,faeprob.gcells,olinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, olinkfens, olinkgcells, 0);
        faeprob.gcells = cat(faeprob.gcells, olinkgcells);
        
        for ind = 2:nspokes
            
            xshift = 2*pi/nspokes;
            tempxsp = tempxsp + (2*pi/nspokes);
            
            % create a set of 2D cells we can extrude into the 3D mesh
            [newfens, newgcells] = Q4_blockx(tempxsp, ysp, []);
            
            % merge the meshes
            [faeprob.fens,faeprob.gcells,newgcells] = merge_meshes(faeprob.fens, faeprob.gcells, newfens, newgcells, 0);
            faeprob.gcells = cat(faeprob.gcells, newgcells);
            
            % create the links between this spoke and the previous spoke
            [ilinkfens,ilinkgcells] = Q4_blockx(linspace(tempxsp(end) - xshift, ...
                                                         tempxsp(1), ...
                                                         ceil(ncircumpoints/nspokes)), ...
                                                ysp(1:nradsipoints), []);
            
            [faeprob.fens,faeprob.gcells,ilinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, ilinkfens, ilinkgcells, 0);
            faeprob.gcells = cat(faeprob.gcells, ilinkgcells);
            
            [olinkfens,olinkgcells] = Q4_blockx(linspace(tempxsp(end) - xshift, ...
                                                         tempxsp(1), ...
                                                         ceil(ncircumpoints/nspokes)), ...
                                                ysp(end-nradsopoints+1:end), []);
            
            [faeprob.fens,faeprob.gcells,olinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, olinkfens, olinkgcells, 0);
            faeprob.gcells = cat(faeprob.gcells, olinkgcells);
            
        end
        
%         xshift = 2*pi/nspokes;
%         tempxsp = tempxsp + (2*pi/nspokes);
            
        % create the links between the start and end
        [ilinkfens,ilinkgcells] = Q4_blockx(linspace(tempxsp(end), ...
                                                     tempxsp(end) + (2*pi/(2*nspokes)) - thetasp/2, ...
                                                     ceil(ncircumpoints/nspokes)), ...
                                            ysp(1:nradsipoints), []);

        [faeprob.fens,faeprob.gcells,ilinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, ilinkfens, ilinkgcells, 0);
        faeprob.gcells = cat(faeprob.gcells, ilinkgcells);
        
        [olinkfens,olinkgcells] = Q4_blockx(linspace(tempxsp(end), ...
                                                     tempxsp(end) + (2*pi/(2*nspokes)) - thetasp/2, ...
                                                     ceil(ncircumpoints/nspokes)), ...
                                            ysp(end-nradsopoints+1:end), []);

        [faeprob.fens,faeprob.gcells,olinkgcells] = merge_meshes(faeprob.fens, faeprob.gcells, olinkfens, olinkgcells, 0);
        faeprob.gcells = cat(faeprob.gcells, olinkgcells);
        
        [faeprob.fens, faeprob.gcells] = merge_nodes(faeprob.fens, faeprob.gcells, mergetol);
        
        tempxsp = xsp - (2*pi/nspokes);
        
        for ind = 1:nspokes
            
            xshift = 2*pi/nspokes;
            tempxsp = tempxsp + xshift;
            
            % select the spoke nodes, we'll shift them to make the spokes
            % the same width along their length 
            spnodes = fenode_select(faeprob.fens, struct ('box', [tempxsp(end) tempxsp(1) Rsoo Rsio], 'inflate', mergetol));
            
            xyz = get (faeprob.fens,'xyz');
            
            spxyz = xyz(spnodes,:);
            
            % get the centerline of the spoke
            spcentre = (tempxsp(end) + tempxsp(1)) / 2;
            
            spxyz(spxyz(:,1) > (spcentre + mergetol),1) = spxyz(spxyz(:,1) > (spcentre + mergetol),1) ...
                                             - thetasp .* ( 1 - Rsio ./ spxyz(spxyz(:,1) > (spcentre + mergetol),2) )/2;
            
            spxyz(spxyz(:,1) < (spcentre - mergetol),1) = spxyz(spxyz(:,1) < (spcentre - mergetol),1) ...
                                             + thetasp .* ( 1 - Rsio ./ spxyz(spxyz(:,1) < (spcentre - mergetol),2) )/2;
                                         
            xyz(spnodes,:) = spxyz;
            
            % replace the coordinates of the nodes
            faeprob.fens = set(faeprob.fens, 'xyz', xyz);
            
        end

    end
    
    % label parts of interest
    
    % The outer surface cells
    celllist = gcell_select(faeprob.fens, faeprob.gcells, struct ('box', [-4*pi, 4*pi, ysp(end), ysp(end-1)], 'inflate', mergetol));
    
    % set the cells to have the outer surface label
    faeprob.gcells = setlabels(faeprob.gcells, Options.OuterSurfaceLabel, celllist);
    
    % extrude the mesh
    [faeprob.fens, faeprob.gcells] = Q4_extrude(faeprob.fens, faeprob.gcells, layers, ls, true);
    
    % get the locations of the nodes in the mesh
    xyz = get (faeprob.fens,'xyz');

    [x,y] = pol2cart(-xyz(:,1),xyz(:,2));

    xyz = [x, y, xyz(:,3)];

%         % shift to match rotor disc positon which is located at the origin
%         xyz(:,3) = xyz(:,3) - (tbi/2);

    % replace the coordinates of the nodes
    faeprob.fens = set(faeprob.fens, 'xyz', xyz);

    % Merge together all the nodes within a tolerance
    [faeprob.fens, faeprob.gcells] = merge_nodes(faeprob.fens, faeprob.gcells, mergetol);
    
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
