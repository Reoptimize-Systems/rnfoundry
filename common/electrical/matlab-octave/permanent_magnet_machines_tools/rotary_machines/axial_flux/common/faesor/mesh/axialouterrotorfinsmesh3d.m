function [fens, gcells, modulexs, modulegap] = axialouterrotorfinsmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, tsuppb, tausupp, nmodules, nmodulesupports, separation, varargin)
% creates a structured FAESOR mesh of one fin of a rotor of an axial-flux
% type machine
%
% Syntax
% 
% [fens,gcells] = axialrotorfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, layers, circumpoints)
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
% separation - this is the space between the inner sufaces of the outer
%   fins
%
% layers - number of layers the mesh will be made up of (before refinement)
%   in the same dimension as tbi, the rotor fin thickness.
% 
% circumpoints - number of circumferential points the mesh will be made up
%   of (before refinement)
%
%
% Output
%
% fens - FAESOR fenodeset object containing the nodes of the rotor fin mesh
%
% gcells - FAESOR gcellset_H20 oject containing the H20 type cells making
%   up the rotor fin mesh
%
%

    Inputs.disclayers = 2;
    Inputs.supportlayers = 2;
    Inputs.circumpoints = 25;
    Inputs.maglabels = [1, 2];
    Inputs.BackIronRadialPointsPerM = 15;
    Inputs.MagnetRadialPointsPerM = 15;
    Inputs.Subset = nmodules;

    Inputs = parseoptions(Inputs, varargin);
    
    if numel(Inputs.maglabels) ~= 2
       error('If supplied, magnet region labels must be two-element vectors, one for each fin') 
    end
    
    % create the bottom rotor part
    [fens1,gcells1,modulexs,modulegap] = axialrotorendfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, ...
                                tsuppb, tausupp, nmodules, nmodulesupports, 'bottom', ...
                                Inputs.disclayers, ...
                                Inputs.supportlayers, ...
                                Inputs.circumpoints, ...
                                Inputs.maglabels(1), ...
                                'BackIronRadialPointsPerM', Inputs.BackIronRadialPointsPerM, ...
                                'MagnetRadialPointsPerM', Inputs.MagnetRadialPointsPerM, ...
                                'Subset', Inputs.Subset);

    % move the mesh in the -ve x direction 
    xyz = get(fens1, 'xyz');
    xyz(:,3) = xyz(:,3) - (separation + tbi)/2;
    fens1 = set(fens1, 'xyz', xyz);

    % create the top rotor part
    [fens2,gcells2,modulexs] = axialrotorendfinmesh3d(Rs, Rbi, Rmi, Rmo, Rbo, tbi, ...
                                tsuppb, tausupp, nmodules, nmodulesupports, 'top', ...
                                Inputs.disclayers, ...
                                Inputs.supportlayers, ...
                                Inputs.circumpoints, ...
                                Inputs.maglabels(2), ...
                                'BackIronRadialPointsPerM', Inputs.BackIronRadialPointsPerM, ...
                                'MagnetRadialPointsPerM', Inputs.MagnetRadialPointsPerM, ...
                                'Subset', Inputs.Subset, ...
                                'ModuleGap', modulegap);

    % move the mesh in the -ve x direction 
    xyz = get(fens2, 'xyz');
    xyz(:,3) = xyz(:,3) + (separation + tbi)/2;
    fens2 = set(fens2, 'xyz', xyz);
    
    % merge the two meshes into one
    [fens,gcells1,gcells2] = merge_meshes(fens1, gcells1, fens2, gcells2, 0);

    % concatenate the cells into one for convenience
    gcells = cat(gcells1, gcells2);

end