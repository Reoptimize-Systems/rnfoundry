function [faeprob, rimlabel] = outerstructuremesh_RF_2(faeprob, Rsii, Rsio, Rsoi, Rsoo, thetasp, ls, lst, lp, nspokes, fixedends, varargin)
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
    Options.RimLabel = 1;
    Options.MeshType = 'H8';
    
    Options = parseoptions(Options, varargin);
    
    rimlabel = Options.RimLabel;
    spokelabel = rimlabel + 1;
    weblabel = spokelabel + 1;
    axlelabel = weblabel + 1;
    
%     layers = meshdivcalc(2, 30, Options.AxialLayersPerM, ls);
    ncircumpoints = meshdivcalc(5, 50, Options.CircumPointsPerM, 2 * pi * mean([Rsoo, Rsoi]));
    nradsipoints = meshdivcalc(5, 20, Options.RadialInnerShaftLayersPerM , Rsoo - Rsoi);
    nradsppoints = meshdivcalc(5, 20, Options.RadialSpokeLayersPerM, Rsoi - Rsio);
    nradsopoints = meshdivcalc(5, 20, Options.RadialOuterShaftLayersPerM, Rsio - Rsii);
    
    rimlayers = meshdivcalc(2, 30, Options.AxialLayersPerM, ls);
    spokelayers = meshdivcalc(2, 5, Options.AxialLayersPerM, lst);
    weblayers = meshdivcalc(2, 3, Options.AxialLayersPerM, lp);
    axlelayers = spokelayers;
    
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
        
        faeprob.gcells = setlabels(faeprob.gcells, spokelabel);

        % create the links between the start and end
        [ilinkfens,ilinkgcells] = Q4_blockx(linspace(tempxsp(1) - (2*pi/(2*nspokes)) + thetasp/2, tempxsp(1), ceil(ncircumpoints/nspokes/2)), ysp(1:nradsipoints), []);
        
        ilinkgcells = setlabels(ilinkgcells, axlelabel);

        [olinkfens,olinkgcells] = Q4_blockx(linspace(tempxsp(1) - (2*pi/(2*nspokes)) + thetasp/2, tempxsp(1), ceil(ncircumpoints/nspokes/2)), ysp(end-nradsopoints+1:end), []);

        olinkgcells = setlabels(olinkgcells, rimlabel);
        
        [weblinkfens,weblinkgcells] = Q4_blockx(linspace(tempxsp(1) - (2*pi/(2*nspokes)) + thetasp/2, tempxsp(1), ceil(ncircumpoints/nspokes/2)), ysp(nradsopoints:end-nradsopoints+1), []);

        weblinkgcells = setlabels(weblinkgcells, weblabel);
        
        [faeprob.fens,gcells] = merge_n_meshes( {faeprob.fens, ilinkfens, olinkfens, weblinkfens}, ...
                                                {faeprob.gcells, ilinkgcells, olinkgcells, weblinkgcells}, 0);
        faeprob.gcells = cat(gcells{:});
        
        for ind = 2:nspokes
            
            xshift = 2*pi/nspokes;
            tempxsp = tempxsp + (2*pi/nspokes);
            
            % create a set of 2D cells we can extrude into the 3D mesh
            [spokefens, spokegcells] = Q4_blockx(tempxsp, ysp, []);
            
            spokegcells = setlabels(spokegcells, spokelabel);
            
            % merge the meshes
            [faeprob.fens,faeprob.gcells,spokegcells] = merge_meshes(faeprob.fens, faeprob.gcells, spokefens, spokegcells, 0);
            faeprob.gcells = cat(faeprob.gcells, spokegcells);
            
            % create the links between this spoke and the previous spoke
            [ilinkfens,ilinkgcells] = Q4_blockx( linspace(tempxsp(end) - xshift, ...
                                                          tempxsp(1), ...
                                                          ceil(ncircumpoints/nspokes)), ...
                                                 ysp(1:nradsipoints), []);
                                             
            ilinkgcells = setlabels(ilinkgcells, axlelabel);
            
            [olinkfens,olinkgcells] = Q4_blockx( linspace(tempxsp(end) - xshift, ...
                                                          tempxsp(1), ...
                                                          ceil(ncircumpoints/nspokes)), ...
                                                 ysp(end-nradsopoints+1:end), []);
                                             
            olinkgcells = setlabels(olinkgcells, rimlabel);
            
            [weblinkfens,weblinkgcells] = Q4_blockx( linspace(tempxsp(end) - xshift, ...
                                                              tempxsp(1), ...
                                                              ceil(ncircumpoints/nspokes)), ...
                                                     ysp(nradsopoints:end-nradsopoints+1), []);

            weblinkgcells = setlabels(weblinkgcells, weblabel);

            [faeprob.fens,gcells] = merge_n_meshes( {faeprob.fens, ilinkfens, olinkfens, weblinkfens}, ...
                                                    {faeprob.gcells, ilinkgcells, olinkgcells, weblinkgcells}, 0);
            faeprob.gcells = cat(gcells{:});
            
        end
        
%         xshift = 2*pi/nspokes;
%         tempxsp = tempxsp + (2*pi/nspokes);
            
        % create the links between the start and end
        [ilinkfens,ilinkgcells] = Q4_blockx(linspace(tempxsp(end), ...
                                                     tempxsp(end) + (2*pi/(2*nspokes)) - thetasp/2, ...
                                                     ceil(ncircumpoints/nspokes/2)), ...
                                            ysp(1:nradsipoints), []);

        ilinkgcells = setlabels(ilinkgcells, axlelabel);
        
        [olinkfens,olinkgcells] = Q4_blockx(linspace(tempxsp(end), ...
                                                     tempxsp(end) + (2*pi/(2*nspokes)) - thetasp/2, ...
                                                     ceil(ncircumpoints/nspokes/2)), ...
                                            ysp(end-nradsopoints+1:end), []);
                                        
        olinkgcells = setlabels(olinkgcells, rimlabel);

        [weblinkfens,weblinkgcells] = Q4_blockx( linspace(tempxsp(end), ...
                                                          tempxsp(end) + (2*pi/(2*nspokes)) - thetasp/2, ...
                                                          ceil(ncircumpoints/nspokes/2)), ...
                                                  ysp(nradsopoints:end-nradsopoints+1), []);

        weblinkgcells = setlabels(weblinkgcells, weblabel);

        [faeprob.fens,gcells] = merge_n_meshes( {faeprob.fens, ilinkfens, olinkfens, weblinkfens}, ...
                                                {faeprob.gcells, ilinkgcells, olinkgcells, weblinkgcells}, 0);
        faeprob.gcells = cat(gcells{:});
        
        [faeprob.fens, faeprob.gcells] = merge_nodes(faeprob.fens, faeprob.gcells, mergetol);
        
        tempxsp = xsp - (2*pi/nspokes);
        
        for ind = 1:nspokes
            
            xshift = 2*pi/nspokes;
            tempxsp = tempxsp + xshift;
            
            % select the spoke nodes, we'll shift them to make the spokes
            % the same width along their length 
            spnodes = fenode_select(faeprob.fens, struct ('box', [tempxsp(end) tempxsp(1) Rsoo Rsio], ...
                                                          'inflate', mergetol));
            
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
    
    % Identify the outer rim cells
    celllist = gcell_select(faeprob.fens, faeprob.gcells, struct ('box', [-4*pi, 4*pi, Rsoi, Rsoo], ...
                                                                  'inflate', mergetol));
    
    faeprob.gcells = setlabels(faeprob.gcells, rimlabel, celllist);
    
    % Identify the axle cells
    celllist = gcell_select(faeprob.fens, faeprob.gcells, struct ('box', [-4*pi, 4*pi, Rsii, Rsio], ...
                                                                  'inflate', mergetol));
    
    faeprob.gcells = setlabels(faeprob.gcells, axlelabel, celllist);

    % extrude the rim
    [rimfens1, rimgcells1] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', rimlabel))), ...
                                       spokelayers, -lst, true );
                                   
    rimfens1 = translate(rimfens1, [0 0 -lp]);
    
    [rimfens2, rimgcells2] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', rimlabel))), ...
                                       weblayers, -lp, true );
    
	[rimfens3, rimgcells3] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', rimlabel))), ...
                                       rimlayers, ls, true );
                                   
	
    [webfens, webgcells] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', weblabel))), ...
                                       weblayers, -lp, true );
                                   
    [spokefens, spokegcells] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', spokelabel))), ...
                                       weblayers, -lst, true );  
                                   
	spokefens = translate(spokefens, [0 0 -lp]);
    
    [axlefens1, axlegcells1] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', axlelabel))), ...
                                       spokelayers, -lst, true );
                                   
    axlefens1 = translate(axlefens1, [0 0 -lp]);
    
    [axlefens2, axlegcells2] = Q4_extrude( faeprob.fens, ...
                                       subset(faeprob.gcells,  gcell_select(faeprob.fens, faeprob.gcells, struct('label', axlelabel))), ...
                                       weblayers, -lp, true );
    
    [faeprob.fens,gcells] = merge_n_meshes( {rimfens1, rimfens2, rimfens3, webfens, spokefens, axlefens1, axlefens2}, ...
                                            {rimgcells1, rimgcells2, rimgcells3, webgcells, spokegcells, axlegcells1, axlegcells2}, ...
                                            0 );
    faeprob.gcells = cat(gcells{:});                      
    
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


function divs = meshdivcalc(mindivs, maxdivs, layersperlength, length)

    divs = min(maxdivs,  max(mindivs, ceil(layersperlength * length)));


end
