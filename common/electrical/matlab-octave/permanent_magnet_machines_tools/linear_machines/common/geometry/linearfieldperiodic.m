function [FemmProblem, nodes, links, info] = linearfieldperiodic (FemmProblem, tmag, zmag, zs, toffset, tsvc, tsve, zsvi, zsvo, varargin)
% generates a linear magnets containing region with periodic boundaries at
% top and bottom
%
% Syntax
%
% [FemmProblem, nodes, links, info] = ...
%       linearfieldperiodic (FemmProblem, tmag, zmag, zs, toffset, tsvc, tsve, zsvi, zsvo)
% [...] = linearfieldperiodic(..., 'Parameter', Value)
%
%
% Inputs
%
%  FemmProblem - FemmProblem structure to which the geometry will be added
%
%  tmag - pole width in radians
%
%  zmag - magnet width in radians
%
%  zs - radial thickness of the magnets
%
%  toffset - radial displacement of the magnet centers from the center
%
%  tsvc - 
%
%  tsve - 
%
%  zsvi - 
%
%  zsvo - 
%
%                  :
%   toffset        :
%<---------------->:
%                  :
%         |         _______
%         |        /       |
%         |_______/        |
%         |                |
%         |________________|
%         |                | ^
%         |      tmag      | : zmag
%         |<-------------->| :
%         |________________| v
%         |<-----> tsve    | ^
%         |_______   zsvo  | :
%       ^ | tsvc  \  ^     | : 
%       : |<------>\_:_____| ;
%  zsvi : |        / :     | : zs
%       v |_______/  v     | :
%         |                | :
%         |________________| v
%         |                |
%         |                |         
%
%  In addition, a number of optional parameters can be specified as
%  parameter-value pairs. Possible parameter-value pairs are:
%
%  'MagDirections' - either a 2 element numeric vector, or a 2 element cell
%    array of strings. If numeric, these are the directions in degrees of
%    the magnet magnetisation. If a cell array of strings, these are
%    evaluated in the FEMM or xfemm lua interpreter to yield the magnet
%    direction in the magnet region elements. Variables that can be used in
%    these strings are:
%
%    'theta': angle in degrees of a line connecting the center of each
%             element with the origin 
%
%    'R'    : length of a line connecting the center of each element with the
%             origin
%
%    'x'    : x position of each element
%
%    'y'    : y position of each elements
% 
%    The default is {0, 180}, resulting in magnets pointing parallel to the
%    x axis
%
%  'MagnetMaterial' - 
%
%  'MagnetGroup' - 
%
%  'SpacerMaterial' - 
%
%  'SpaceGroup' - 
%
%  'Tol' - 
%
%  'MeshSize' - 
%
% Output
%
%  FemmProblem - 
%
%  nodes - 
%
%  links - 
%
%  info - 
%
%

    Inputs.MagDirections = {0, 180};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpacerMaterial = 1;
    Inputs.SpaceGroup = 0;
    Inputs.CavityMaterial = [];
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;
    Inputs.BoundName = '';
    Inputs.NPolePairs = 1;
    Inputs.Flip = false;

    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if numel(Inputs.MagDirections) ~= 2*Inputs.NPolePairs
        if numel (Inputs.MagDirections) == 2
            Inputs.MagDirections = repmat (Inputs.MagDirections, 1, Inputs.NPolePairs);
        else
            error ('The number of magnet directions supplied is not appropriate');
        end
    end
    
    if isnumeric (Inputs.MagDirections) && isvector (Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
    if isempty (Inputs.CavityMaterial)
        Inputs.CavityMaterial = Inputs.SpacerMaterial;
    end
    
    % get the number of existing nodes, segments, boundaries etc. if any
    elcount = elementcount_mfemm (FemmProblem);
    
    % add a periodic boundary for the edges of the magnets region
    if isempty (Inputs.BoundName)
        [FemmProblem, info.BoundaryInds] = addboundaryprop_mfemm (FemmProblem, ...
                                                'Rect Mags Periodic', 4);
    else
        BoundaryProp = newboundaryprop_mfemm (Inputs.BoundName, 4, false);
        info.BoundaryInds = elcount.NBoundaryProps + 1;
        FemmProblem.BoundaryProps(info.BoundaryInds) = BoundaryProp;
    end
    
    % construct the segments and label locations for the magnet regions
    npoles = Inputs.NPolePairs * 2;
    [nodes, links, info.OuterNodeIDs] = field_body (npoles, zmag, tmag, toffset, zs, tsvc, tsve, zsvi, zsvo, Inputs.Flip);
    zpole = zmag + zs;
    
    % add all the nodes to the problem
    [FemmProblem, ~, info.NodeIDs] = addnodes_mfemm (FemmProblem, ...
                            nodes(:,1), nodes(:,2), 'InGroup', Inputs.MagnetGroup);
                        
    % Add all the links as segments
    FemmProblem = addsegments_mfemm (FemmProblem, links(1:end,1), links(1:end,2), ...
        'InGroup', Inputs.MagnetGroup);
    
    % Periodic boundary at bottom
    [FemmProblem, seginds] = addsegments_mfemm (FemmProblem, info.OuterNodeIDs(4), info.OuterNodeIDs(3), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds).Name, ...
        'InGroup', Inputs.MagnetGroup);
    
    info.BottomSegInd = seginds;

    % Periodic boundary at top
    [FemmProblem, seginds] = addsegments_mfemm (FemmProblem, info.OuterNodeIDs(1), info.OuterNodeIDs(2), ...
        'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds).Name, ...
        'InGroup', Inputs.MagnetGroup);
    
    info.TopSegInd = seginds;
    
    rectcentres = [ toffset, zmag/4;
                    repmat(toffset, npoles-1, 1), zmag/2 + zs + zmag/2 + linspace(0, (npoles-2)*zpole, npoles-1)';
                    toffset, zpole * npoles - zmag/4;
                  ];
              
	rectcentres(:,3) = [ 1; repmat([2,1]', npoles/2, 1)];
    
    % Now add the labels
    info.MagnetBlockInds = [];
    if ~isempty(Inputs.MagnetMaterial)
        
        % add the magnet labels
        for i = 1:size(rectcentres, 1)

            [FemmProblem, info.MagnetBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, rectcentres(i,1), rectcentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'MagDir', Inputs.MagDirections{rectcentres(i,3)}, ...
                                            'InGroup', Inputs.MagnetGroup);
                
            info.MagnetBlockInds(end, 2) = rectcentres(i,3);

        end
    
    end
    
    if Inputs.Flip
        spacercentres = [ repmat(mean([toffset+tmag/2-tsvc, toffset-tmag/2]), npoles, 1), ...
                          zmag/2 + zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                        ];
    else
        spacercentres = [ repmat(mean([toffset-tmag/2+tsvc, toffset+tmag/2]), npoles, 1), ...
                          zmag/2 + zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                        ];
    end
    
    info.SpacerBlockInds = [];
    if ~isempty(Inputs.SpacerMaterial)
        % Add the other labels
        
        for i = 1:size(spacercentres, 1)

            [FemmProblem, info.SpacerBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, spacercentres(i,1), spacercentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.SpacerMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpaceGroup);

        end
    
    end
    
    if Inputs.Flip
        cavitycentres = [ repmat(mean([toffset+tmag/2, toffset+tmag/2-tsvc]), npoles, 1), ...
                          zmag/2 + zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                        ];
    else
        cavitycentres = [ repmat(mean([toffset-tmag/2+tsvc, toffset-tmag/2]), npoles, 1), ...
                          zmag/2 + zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                        ];
    end
    
    info.CavityBlockInds = [];
    if ~isempty(Inputs.SpacerMaterial)
        % Add the other labels
        
        for i = 1:size(cavitycentres, 1)

            [FemmProblem, info.CavityBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, cavitycentres(i,1), cavitycentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.CavityMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpaceGroup);

        end
    
    end

    
end

function [nodes, links, outernodeids] = field_body (npoles, zmag, tmag, toffset, zs, tsvc, tsve, zsvi, zsvo, flip)

    [nodes, links, cavouternodeids] = all_disc_cavities (npoles, zmag, tmag, zs, tsvc, tsve, zsvi, zsvo);
    
    zpole = zmag + zs;
    
    % shift everything up a bit so middle of first magnet half lies on x=0
    % axis
    nodes(:,2) = nodes(:,2) + (zs + zmag) / 2;
    
    % add the outer nodes going clockwise from top left
    outernodeids = size (nodes, 1):size (nodes, 1)+3;
    
    nodes = [ nodes;
              0, zpole * npoles;
              tmag, zpole * npoles;
              tmag, 0;
              0, 0 ];
          
    % add inner (left) nodes and links
    
    lnodestartid = size(nodes, 1);
    
    % TODO: this all needs fixed, bottom half mag not included
    
    % lower mag/steel boundary nodes
    nodes = [ nodes;
        zeros(npoles,1), zmag/2 + linspace(0, zpole * (npoles-1), npoles)' ];
    
    % mag/steel boundary nodes
    nodes = [ nodes;
        zeros(npoles,1), zmag/2 + zs + linspace(0, zpole * (npoles-1), npoles)' ];
    
    % link up left side
    links = [ links;
              outernodeids(4), lnodestartid; % bottom corner to bottom mag top left
              cavouternodeids(1:end,3), (lnodestartid:(lnodestartid+npoles-1))';
              cavouternodeids(1:end,1), (lnodestartid+npoles:(lnodestartid+2*npoles-1))';
              lnodestartid+2*npoles-1, outernodeids(1);
              (lnodestartid+1:(lnodestartid+npoles-1))', (lnodestartid+npoles:(lnodestartid+2*npoles-2))';
              cavouternodeids(1:end,1), cavouternodeids(1:end,3);
            ];
    
    % add outer (right) nodes and links
    % lower mag nodes
    rnodestartid = size(nodes, 1);
    nodes = [ nodes; repmat(tmag, npoles,1), zmag/2 + linspace(0, zpole * (npoles-1), npoles)' ];
    
    % upper mag nodes
    nodes = [ nodes; repmat(tmag, npoles,1), zmag/2 + zs + linspace(0, zpole * (npoles-1), npoles)' ];
    
    % link up left side
    links = [ links;
              outernodeids(3), rnodestartid; % bottom right corner to bottom mag top right
              (rnodestartid:(rnodestartid+npoles-1))', cavouternodeids(1:end,2);
              cavouternodeids(1:end,2), (rnodestartid+npoles:(rnodestartid+2*npoles-1))';
              (rnodestartid+1:(rnodestartid+npoles-1))', (rnodestartid+npoles:(rnodestartid+2*npoles-2))';
              rnodestartid+2*npoles-1, outernodeids(2);
              ...cavouternodeids(1:end,1), cavouternodeids(1:end,3);
            ];
        
    % make horizontal links
    links = [ links;
              (lnodestartid:(lnodestartid+2*npoles-1))', (rnodestartid:(rnodestartid+2*npoles-1))';
             ]
         
    if flip
        nodes = mfemmdeps.reflect2d (nodes, 'Axis', 'y');
        nodes(:,1) = nodes(:,1) + tmag;
    end
    
    % shift all nodes right to account for specified offset
    nodes(:,1) = nodes(:,1) - tmag/2 + toffset;

end


function [nodes, links, outernodeids] = all_disc_cavities (npoles, zmag, tmag, zs, tsvc, tsve, zsvi, zsvo)
% Draws multiple cavities like so:
%
%
%   12      15       
%    x-------x
%             \
%              \
%           16  x-------x 13
%              /
%             /
%    x-------x
%   14      17
%
%
%
%   6       9       
%    x-------x
%             \
%              \
%           10  x-------x 7
%              /
%             /
%    x-------x
%   8      11
%
%
%
%   0       3       
%    x-------x
%             \
%              \
%            4  x-------x 1
%              /
%             /
%    x-------x
%   2       5
%
%


    [nodes, links, outernodeids] = single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo);
    
    for ind = 2:npoles
        
        [newnodes, newlinks, newouternodeids] = single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo);
        
        links = [ links; newlinks + size(nodes,1) ];
        
        outernodeids = [ outernodeids;
                         newouternodeids + size(nodes,1) ];
                     
        nodes = [ nodes; newnodes(:,1), newnodes(:,2) + (ind-1) * (zmag + zs) ];
        
    end
end


function [nodes, links, outernodeids] = single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo)
% Draws the boundary of the cavity like so:
%
%   0       3       
%    x-------x
%             \
%              \
%            4  x-------x 1
%              /
%             /
%    x-------x
%   2       5
%
%

    nodes = [ 0, zsvi/2;
              tmag, 0;
              0, -zsvi/2;
              tsve, zsvo/2;
              tsvc, 0;
              tsve, -zsvo/2 ];
          
    links = [0, 3;
             3, 4;
             4, 5;
             5, 2;
             ...4, 1
             ];
         
	outernodeids = [0,1,2];
    
end


