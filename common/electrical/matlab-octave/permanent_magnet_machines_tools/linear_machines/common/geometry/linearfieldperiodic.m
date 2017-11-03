function [FemmProblem, nodes, links, info] = linearfieldperiodic (FemmProblem, type, vars, varargin)
% generates region containing linear magnets with periodic boundaries at
% top and bottom
%
% Syntax
%
% [FemmProblem, nodes, links, info] = ...
%       linearfieldperiodic (FemmProblem, type, vars)
% [...] = linearfieldperiodic(..., 'Parameter', Value)
%
%
% Inputs
%
%  FemmProblem - FemmProblem structure to which the geometry will be added
%
%  vars - structure containing the variables describing the linear field
%    design. This must contain appropriate fields depending on the they
%    type of geometry to be created. The appropriate fields are shown for
%    each geometry type below:
%
%    Geom Type: 'simple' 
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
%      tmag - thickness of the magnets
%
%      zmag - height of the magnets
%
%      zs - height of magnet spacer
%
%      toffset - displacement of the magnet centers from the origin
%
%      tsvc - (optional) thickness of cavity at spacer center. If not
%        present will be set to tmag/2.
%
%      tsve -(optional) thickness of cavity at spacer edge. If not
%        present will be set to tmag/2.
%
%      zsvi - (optional) height of cavity at spacer center end. If not
%        present will be set to zs/2.
%
%      zsvo - (optional) height of cavity at spacer inner end. If not
%        present will be set to zs/2.
%
%
%  In addition, a number of optional parameters can be specified as
%  parameter-value pairs. Possible parameter-value pairs are:
%
%  'NPolePairs' - number of pairs of magnets and spacers to be drawn.
%    Default is 1.
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
%  'MagnetMaterial' - Index of the material in the FemmProblem.Materials
%    structure to use as magnet material. If empty, no magnet region block
%    labels are drawn. Default is 1 if not supplied.
%
%  'MagnetGroup' - number to be assigned as the group number of the maget
%    regions. Default is 0 if not supplied. 
%
%  'SpacerMaterial' - Index of the material in the FemmProblem.Materials
%    structure to use as spacer material. If empty, no spacer region block
%    labels are drawn. Default is 1 if not supplied.
%
%  'SpacerGroup' - number to be assigned as the group number of the spacer
%    regions. Default is 0 if not supplied. Default is 1 if not supplied.
%
%  'MeshSize' - Mesh size to use in the entire region.
%
%  'AddPeriodicBoundaries' - logical (true/false) flag indicating whether
%    to add periodic boundaries to the top and bottom of the drawing.
%    Default is true.
%
% Output
%
%  FemmProblem - input femmproblem structure with new elements added.
%
%  nodes - (n x 2) matrix of coordinates of the new nodes added.
%
%  links - (k x 2) matrix of node links, indexed into the nodes output
%
%  info - structure containing other information about the drawing. It can
%    contain the following fields:
%   
%    'OuterNodeIDs' - array of four integers containing the ids of nodes at
%      the four corners of the field drawing. Nodes areindexed from the
%      ottom left corner, and then anti-clockwise around the drawing.
%
%    'NodeIDs' - array of integers containing the ids of all new nodes
%      added.
%
%    'BoundaryInds' - array of integers containing the indices in the
%      FemmProblem.Boundaries structure of the newly added boundaries. Not
%      present if AddPeriodicBoundaries is false.
%
%    'TopSegInd' - index in the FemmProblem structure of the top horizontal
%      segment in the drawing
%
%    'BottomSegInd' - index in the FemmProblem structure of the bottom
%      horizontal segment in the drawing
%
%    'MagnetBlockInds' - indices in the FemmProblem structure of the magnet
%      region block labels. Empty if the 'MagnetMaterial' optional argument
%      is empty (so no labels are draw).
%
%    'SpacerBlockInds' - indices in the FemmProblem structure of the spacer
%      region block labels. Empty if the 'MagnetMaterial' optional argument
%      is empty (so no labels are draw).
%
%    'CavityBlockInds' - indices in the FemmProblem structure of the cavity
%      region block labels. Empty if there are no cavities.
%

    Inputs.MagDirections = {0, 180};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpacerMaterial = 1;
    Inputs.SpacerGroup = 0;
    Inputs.CavityMaterial = [];
    Inputs.CavityGroup = 0;
%     Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;
    Inputs.BoundName = '';
    Inputs.NPolePairs = 1;
    Inputs.Flip = false;
    Inputs.AddPeriodicBoundaries = true;

    Inputs = parse_pv_pairs (Inputs, varargin);
    
    check.isScalarInteger (Inputs.MagnetMaterial, true, 'MagnetMaterial');
    check.isScalarInteger (Inputs.MagnetGroup, true, 'MagnetGroup');
    check.isScalarInteger (Inputs.SpacerMaterial, true, 'SpacerMaterial');
    check.isScalarInteger (Inputs.SpacerGroup, true, 'SpacerGroup');
    check.isNumericScalar (Inputs.MeshSize, true, 'MeshSize');
    check.isScalarInteger (Inputs.NPolePairs, true, 'NPolePairs');
    check.isLogicalScalar (Inputs.Flip, true, 'Flip');
    check.isLogicalScalar (Inputs.AddPeriodicBoundaries, true, 'AddPeriodicBoundaries');
    
    % get the number of existing nodes, segments, boundaries etc. if any
    elcount = elementcount_mfemm (FemmProblem);
    
    switch type
        
        case 'simple'
            
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
            
            assert (isfield (vars, 'tmag'), 'tmag field is not present in vars');
            assert (isfield (vars, 'zmag'), 'zmag field is not present in vars');
            assert (isfield (vars, 'zs'), 'zs field is not present in vars');
            assert (isfield (vars, 'toffset'), 'toffset field is not present in vars');
            check.isNumericScalar (vars.tmag, true, 'vars.tmag');
            check.isNumericScalar (vars.zmag, true, 'vars.zmag');
            check.isNumericScalar (vars.zs, true, 'vars.zs');
            check.isNumericScalar (vars.toffset, true, 'vars.toffset');
            vars = setfieldifabsent (vars, 'tsvc', vars.tmag/2);
            vars = setfieldifabsent (vars, 'tsve', vars.tmag/2);
            vars = setfieldifabsent (vars, 'zsvi', vars.zs/2);
            vars = setfieldifabsent (vars, 'zsvo', vars.zs/2);
            check.isNumericScalar (vars.tsvc, true, 'vars.tsvc');
            check.isNumericScalar (vars.tsve, true, 'vars.tsve');
            check.isNumericScalar (vars.zsvi, true, 'vars.zsvi');
            check.isNumericScalar (vars.zsvo, true, 'vars.zsvo');
            
            % add a periodic boundary for the edges of the magnets region
            if Inputs.AddPeriodicBoundaries
                if isempty (Inputs.BoundName)
                    [FemmProblem, info.BoundaryInds] = addboundaryprop_mfemm (FemmProblem, ...
                                                            'Rect Mags Periodic', 4);
                else
                    BoundaryProp = newboundaryprop_mfemm (Inputs.BoundName, 4, false);
                    info.BoundaryInds = elcount.NBoundaryProps + 1;
                    FemmProblem.BoundaryProps(info.BoundaryInds) = BoundaryProp;
                end
                
                boundname = FemmProblem.BoundaryProps(info.BoundaryInds).Name;
            else
                boundname = '';
                info.BoundaryInds = [];
            end

            % construct the segments and nodes for the chosen geometry
            npoles = Inputs.NPolePairs * 2;
            [nodes, links, info.OuterNodeIDs] = simple_field_body (npoles, vars.zmag, vars.tmag, vars.toffset, vars.zs, vars.tsvc, vars.tsve, vars.zsvi, vars.zsvo, Inputs.Flip);
            zpole = vars.zmag + vars.zs;

            % add all the nodes to the problem
            [FemmProblem, ~, info.NodeIDs] = addnodes_mfemm (FemmProblem, ...
                                    nodes(:,1), nodes(:,2), 'InGroup', Inputs.MagnetGroup);

            % Add all the links as segments
            FemmProblem = addsegments_mfemm (FemmProblem, links(1:end,1), links(1:end,2), ...
                'InGroup', Inputs.MagnetGroup);

            % Periodic boundary at bottom
            [FemmProblem, seginds] = addsegments_mfemm (FemmProblem, info.OuterNodeIDs(4), info.OuterNodeIDs(3), ...
                'BoundaryMarker', boundname, ...
                'InGroup', Inputs.MagnetGroup);

            info.BottomSegInd = seginds;
            
            % Periodic boundary at top
            [FemmProblem, seginds] = addsegments_mfemm (FemmProblem, info.OuterNodeIDs(1), info.OuterNodeIDs(2), ...
                'BoundaryMarker', boundname, ...
                'InGroup', Inputs.MagnetGroup);

            info.TopSegInd = seginds;

            magcentres = [ vars.toffset, vars.zmag/4;
                            repmat(vars.toffset, npoles-1, 1), vars.zmag/2 + vars.zs + vars.zmag/2 + linspace(0, (npoles-2)*zpole, npoles-1)';
                            vars.toffset, zpole * npoles - vars.zmag/4;
                          ];

            magcentres(:,3) = [ 1; repmat([2,1]', npoles/2, 1)];
            
            if Inputs.Flip
                spacercentres = [ repmat(mean([vars.toffset+vars.tmag/2-vars.tsvc, vars.toffset-vars.tmag/2]), npoles, 1), ...
                                  vars.zmag/2 + vars.zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                                ];
            else
                spacercentres = [ repmat(mean([vars.toffset-vars.tmag/2+vars.tsvc, vars.toffset+vars.tmag/2]), npoles, 1), ...
                                  vars.zmag/2 + vars.zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                                ];
            end
            
            if Inputs.Flip
                cavitycentres = [ repmat(mean([vars.toffset+vars.tmag/2, vars.toffset+vars.tmag/2-vars.tsvc]), npoles, 1), ...
                                  vars.zmag/2 + vars.zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                                ];
            else
                cavitycentres = [ repmat(mean([vars.toffset-vars.tmag/2+vars.tsvc, vars.toffset-vars.tmag/2]), npoles, 1), ...
                                  vars.zmag/2 + vars.zs/2 + linspace(0, (npoles-1)*zpole, npoles)';
                                ];
            end
            
        otherwise
            
            error ('RENEWNET:linearfieldperiodic:badtype', ...
                'Unrecognised linar field type specifcation: %s', Inputs.Type);
    
    end
    
    % Now add the labels
    info.MagnetBlockInds = [];
    if ~isempty(Inputs.MagnetMaterial)
        
        % add the magnet labels
        for i = 1:size(magcentres, 1)

            [FemmProblem, info.MagnetBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, magcentres(i,1), magcentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.MagnetMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'MagDir', Inputs.MagDirections{magcentres(i,3)}, ...
                                            'InGroup', Inputs.MagnetGroup);
                
            info.MagnetBlockInds(end, 2) = magcentres(i,3);

        end
    
    end
    
    info.SpacerBlockInds = [];
    if ~isempty(Inputs.SpacerMaterial)
        % Add the other labels
        
        for i = 1:size(spacercentres, 1)

            [FemmProblem, info.SpacerBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, spacercentres(i,1), spacercentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.SpacerMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.SpacerGroup);

        end
    
    end
    
    info.CavityBlockInds = [];
    if ~isempty(Inputs.SpacerMaterial)
        % Add the other labels
        
        for i = 1:size(cavitycentres, 1)

            [FemmProblem, info.CavityBlockInds(end+1,1)] = addblocklabel_mfemm (FemmProblem, cavitycentres(i,1), cavitycentres(i,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.CavityMaterial).Name, ...
                                            'MaxArea', Inputs.MeshSize, ...
                                            'InGroup', Inputs.CavityGroup);

        end
    
    end

    
end

function [nodes, links, outernodeids] = simple_field_body (npoles, zmag, tmag, toffset, zs, tsvc, tsve, zsvi, zsvo, flip)

    [nodes, links, cavouternodeids] = simple_all_disc_cavities (npoles, zmag, tmag, zs, tsvc, tsve, zsvi, zsvo);
    
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
             ];
         
    if flip
        nodes = mfemmdeps.reflect2d (nodes, 'Axis', 'y');
        nodes(:,1) = nodes(:,1) + tmag;
    end
    
    % shift all nodes right to account for specified offset
    nodes(:,1) = nodes(:,1) - tmag/2 + toffset;

end


function [nodes, links, outernodeids] = simple_all_disc_cavities (npoles, zmag, tmag, zs, tsvc, tsve, zsvi, zsvo)
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


    [nodes, links, outernodeids] = simple_single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo);
    
    for ind = 2:npoles
        
        [newnodes, newlinks, newouternodeids] = simple_single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo);
        
        links = [ links; newlinks + size(nodes,1) ];
        
        outernodeids = [ outernodeids;
                         newouternodeids + size(nodes,1) ];
                     
        nodes = [ nodes; newnodes(:,1), newnodes(:,2) + (ind-1) * (zmag + zs) ];
        
    end
end


function [nodes, links, outernodeids] = simple_single_disc_cavity (tmag, tsvc, tsve, zsvi, zsvo)
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


