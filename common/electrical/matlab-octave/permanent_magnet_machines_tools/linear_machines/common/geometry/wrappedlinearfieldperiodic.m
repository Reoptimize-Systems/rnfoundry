function [FemmProblem, wrapperthickness, info] = wrappedlinearfieldperiodic (FemmProblem, fieldtype, vars, wrapperthickness, varargin)
% creates a FemmProblem geometry of a radial section containing two magnets
% optionally wrapped with additional layers and with periodic edges
%
%
% Syntax
%
% [FemmProblem, wrapperthickness, info] = ...
%   wrappedlinearfieldperiodic(FemmProblem, fieldtype, vars, wrapperthickness)
% [...] = wrappedlinearfieldperiodic(..., 'param', value, ...)
%
% Description
%
% wrappedlinearfieldperiodic creates a periodic linear machine field
% geometry. In addition, any number of rectangular regions can be added
% either to inside or outside of the main region (like wrappers for the
% main region).
%
% Inputs
%
%  FemmProblem - FemmProblem structure to which the geometry will be added.
% 
%  fieldtype - string containing the geometry type to be drawn. See vars
%    below.
%
%  vars - structure containing the variables describing the linear field
%    design. This must contain appropriate fields depending on the they
%    type of geometry to be created (specified in the fieldtype input
%    described above). The appropriate fields are shown for each avaialble
%    geometry type below:
%
%    Geom Type: 'simple' 
%
%    A geometry type consisting of a series of rectangular regions with
%    alternating materials. Optionally a 'cavity' region may be added as
%    shown with a different material in this space. See also the
%    'CavityMaterial' and 'CavityGroup' Inputs below.
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
%  wrapperthickness - either an (n x 2) matrix or column vector of wrapper
%    thicknesses. If an (n x 2) matrix. the first column specifies the
%    thickness of any desired wrapper on the left of the magnets, and the
%    second column the thickness of wrappers on the right hand side. The
%    wrappers are added moving progressively further from the magnet
%    position (either to the left or right) down the rows of the matrix.
%    Wrappers with thicknesses less than a tolerance are not added. The
%    default tolerance is 1e-5, but this value can be changed using the
%    appropriate optional parameter value pair (see below). If
%    wrapperthickness is a column vector the same thicknesses are used on
%    both sides.
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
%  'CavityMaterial' - Relevant only of 'simple' geometry. Index of the
%    material in the FemmProblem.Materials structure to use as cavity
%    material. If empty, no magnet region block labels are drawn. Default
%    is 1 if not supplied.
%
%  'CavityGroup' - Relevant only of 'simple' geometry. Number to be
%    assigned as the group number of the cavity regions. Default is 0 if
%    not supplied. Default is 1 if not supplied.
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
%  wrapperthickness - either an (n x 2) matrix or column vector of wrapper
%    thicknesses. If an (n x 2) matrix. the first column specifies the
%    thickness of any wrappers on the left of the magnets, and the
%    second column the thickness of wrappers on the right hand side. The
%    wrappers are added moving progressively further from the magnet
%    position (either to the left or right) down the rows of the matrix.
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
%    'InnerCentres' - 
%
%    'OuterCentres' - 
%
%
% See also: linearfieldperiodic
%

    Inputs.MagDirections = {0, 180};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpacerMaterial = 1;
    Inputs.SpacerGroup = 0;
    Inputs.CavityMaterial = [];
    Inputs.CavityGroup = 0;
    Inputs.WrapperGroup = 0;
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;
    Inputs.NPolePairs = 1;
    Inputs.Flip = false;
    Inputs.AddPeriodicBoundaries = true;
    
    % parse the input arguments
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty (FemmProblem)
        FemmProblem = newproblem_mfemm ('axi');
    end
    
    check.allowedStringInputs (fieldtype, {'simple'}, true, 'fieldtype');
    
    check.isScalarInteger (Inputs.MagnetMaterial, true, 'MagnetMaterial');
    check.isScalarInteger (Inputs.MagnetGroup, true, 'MagnetGroup');
    check.isScalarInteger (Inputs.SpacerMaterial, true, 'SpacerMaterial');
    check.isScalarInteger (Inputs.SpacerGroup, true, 'SpacerGroup');
    check.isNumericScalar (Inputs.MeshSize, true, 'MeshSize');
    check.isScalarInteger (Inputs.NPolePairs, true, 'NPolePairs');
    check.isLogicalScalar (Inputs.Flip, true, 'Flip');
    check.isLogicalScalar (Inputs.AddPeriodicBoundaries, true, 'AddPeriodicBoundaries');
    
    if isscalar(Inputs.WrapperGroup)
        Inputs.WrapperGroup = repmat (Inputs.WrapperGroup, size(wrapperthickness));
    elseif ~samesize(wrapperthickness, Inputs.WrapperGroup)
        error ('RENEWNET:wrappedlinearfield:badwrapper', ...
            'Your must supply either a scalar wrapper group, or a vector the smame size as the number of wrappers.')
    end
    
    if isnumeric(Inputs.MagDirections) && isvector(Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
    [FemmProblem, ~, ~, info] = linearfieldperiodic (FemmProblem, fieldtype, vars, ...
                'MagDirections', Inputs.MagDirections, ...
                'MagnetMaterial', Inputs.MagnetMaterial, ...
                'MagnetGroup', Inputs.MagnetGroup, ...
                'SpacerMaterial', Inputs.SpacerMaterial, ...
                'SpacerGroup', Inputs.SpacerGroup, ...
                ... 'Tol', Inputs.Tol, ...
                'MeshSize', Inputs.MeshSize, ...
                'NPolePairs', Inputs.NPolePairs, ...
                'Flip', Inputs.Flip, ...
                'AddPeriodicBoundaries', Inputs.AddPeriodicBoundaries, ...
                'CavityMaterial', Inputs.CavityMaterial, ...
                'CavityGroup', Inputs.CavityGroup );
    
    zpole = vars.zmag + vars.zs;
    
    % now add the back iron nodes and links, depending on their thicknesses
    if size(wrapperthickness,2) < 2
        wrapperthickness = [wrapperthickness, wrapperthickness];
    elseif size(wrapperthickness,2) > 2
        error('RENEWNET:wrappedlinearfield:badwrapper', ...
            'wrapperthickness must be a scaler or a (1 x 2) vector or (n x 2) matrix')
    end
    
    elcount = elementcount_mfemm(FemmProblem);
    
    info.InnerCentres = zeros(size(wrapperthickness)) * NaN;
    
    innert = vars.toffset - vars.tmag/2;
    
    % wrapperthickness(1,1) is the first inner region thickness
    if wrapperthickness(1,1) > Inputs.Tol
        % add the nodes and segments for the inner side
        %
        
        if wrapperthickness(1,1) > vars.toffset
            error('RENEWNET:wrappedlinearfield:badwrapper', ...
                'wrapper thickness cannot be greater than magnet leftmost position.');
        end
        
        innert = innert - wrapperthickness(1,1);
        
        % First node is to left of first node in 'nodes' matrix. this is at
        % the bottom of the sim
        [FemmProblem, ~, botnodeid] = addnodes_mfemm (FemmProblem, ...
                            innert, ...
                            0, ...
                            'InGroup', Inputs.WrapperGroup(1,1));
        
        % Second node is to left of penultimate node in 'nodes' matrix
        [FemmProblem, ~, topnodeid] = addnodes_mfemm (FemmProblem, ...
                             innert, ...
                             Inputs.NPolePairs*2*zpole, ...
                            'InGroup', Inputs.WrapperGroup(1,1));

        if Inputs.AddPeriodicBoundaries
            % add a new periodic boundary for the top and bottom of the region
            [FemmProblem, info.BoundaryInds(end+1)] = ...
                addboundaryprop_mfemm (FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);
            
            tbboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;
        else
            tbboundarymarker = '';
        end

        % Seg with Periodic boundary at top
        [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                        info.OuterNodeIDs(1), ...
                                        topnodeid, ...
                                        'BoundaryMarker', tbboundarymarker, ...
                                        'InGroup', Inputs.WrapperGroup(1,1));

        info.TopSegInd = [info.TopSegInd, segind];

        lastbotnodeid = info.OuterNodeIDs(4);
        
        % Seg at bottom
        [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                        lastbotnodeid, ...
                                        botnodeid, ...
                                        'BoundaryMarker', tbboundarymarker, ...
                                        'InGroup', Inputs.WrapperGroup(1,1));

        info.BottomSegInd = [info.BottomSegInd, segind];
        
        % Add a node at the mid-point of the wrapper
        [FemmProblem, ~, midnodeid] = addnodes_mfemm (FemmProblem, ...
                                        innert, ...
                                        Inputs.NPolePairs*zpole, ...
                                        'InGroup', Inputs.WrapperGroup(1,1));
                        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                            botnodeid, ...
                                            midnodeid, ...
                                            'InGroup', Inputs.WrapperGroup(1,1) );
                                        
        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                            midnodeid, ...
                                            topnodeid, ...
                                            'InGroup', Inputs.WrapperGroup(1,1) );
                    
        info.InnerCentres(1,:) = [innert+wrapperthickness(1,1)/2, Inputs.NPolePairs*zpole];
        
    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,1) = 0;
        topnodeid = info.OuterNodeIDs(1);
        botnodeid = info.OuterNodeIDs(4);
    end
    
    % now add all subsequent inner wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,1) > Inputs.Tol
            
            innert = innert - wrapperthickness(i,1);
        
            if innert < Inputs.Tol
                error('Inner wrapper radii must all be greater than tolerance.')
            end
            
            lastbotnodeid = botnodeid;
            lasttopnodeid = topnodeid;
            
            % First node is to left of first node in 'nodes' matrix. this is at
            % the bottom of the sim
            [FemmProblem, ~, botnodeid] = addnodes_mfemm (FemmProblem, ...
                                innert, ...
                                0, ...
                                'InGroup', Inputs.WrapperGroup(i,1));
   
                
            % Second node is to left of penultimate node in 'nodes' matrix
            [FemmProblem, ~, topnodeid] = addnodes_mfemm (FemmProblem, ...
                                 innert, ...
                                 Inputs.NPolePairs*2*zpole, ...
                                'InGroup', Inputs.WrapperGroup(i,1));

            if Inputs.AddPeriodicBoundaries
                % add a new periodic boundary for the top and bottom of the region
                [FemmProblem, info.BoundaryInds(end+1)] = ...
                    addboundaryprop_mfemm(FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);

                tbboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;
            else
                tbboundarymarker = '';
            end
            % Seg with Periodic boundary at top
            [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                            lasttopnodeid, ...
                                            topnodeid, ...
                                            'BoundaryMarker', tbboundarymarker, ...
                                            'InGroup', Inputs.WrapperGroup(i,1));

            info.TopSegInd = [info.TopSegInd, segind];
            
            % Seg at bottom
            [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                            lastbotnodeid, ...
                                            botnodeid, ...
                                            'BoundaryMarker', tbboundarymarker, ...
                                            'InGroup', Inputs.WrapperGroup(i,1));

            info.BottomSegInd = [info.BottomSegInd, segind];
            
            % Add a node at the mid-point of the wrapper
            [FemmProblem, ~, midnodeid] = addnodes_mfemm (FemmProblem, ...
                                            innert, ...
                                            Inputs.NPolePairs*zpole, ...
                                            'InGroup', Inputs.WrapperGroup(i,1));
                                    
            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                             botnodeid, ...
                                             midnodeid, ...
                                             'InGroup', Inputs.WrapperGroup(i,1) );

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                             midnodeid, ...
                                             topnodeid, ...
                                             'InGroup', Inputs.WrapperGroup(i,1) );
            
            info.InnerCentres(i,:) = [innert+wrapperthickness(i,1)/2, Inputs.NPolePairs*zpole];
            
        else
            wrapperthickness(i,1) = 0;
        end
        
    end
    
    % Now we dow the wrappers on the other side
    info.OuterCentres = zeros(size(wrapperthickness)) * NaN;
    
    outert = vars.toffset + vars.tmag/2;
    
    % wrapperthickness(1,2) is the first outer region thickness
    if wrapperthickness(1,2) > Inputs.Tol
        
        outert = outert + wrapperthickness(1,2);

        % First node is to right of second node in 'nodes' matrix. this is at
        % the bottom of the sim
        [FemmProblem, ~, botnodeid] = addnodes_mfemm (FemmProblem, ...
                                     outert, ...
                                     0, ...
                                     'InGroup', Inputs.WrapperGroup(1,2));


        % Second node is to left of penultimate node in 'nodes' matrix
        [FemmProblem, ~, topnodeid] = addnodes_mfemm (FemmProblem, ...
                                     outert, ...
                                     Inputs.NPolePairs*2*zpole, ...
                                    'InGroup', Inputs.WrapperGroup(1,2));

        if Inputs.AddPeriodicBoundaries
            % add a new periodic boundary for the top and bottom of the
            % region
            [FemmProblem, info.BoundaryInds(end+1)] = ...
                addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

            tbboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;
        else
            tbboundarymarker = '';
        end

        % Seg with Periodic boundary at top
        [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                         info.OuterNodeIDs(2), ...
                                         topnodeid, ...
                                         'BoundaryMarker', tbboundarymarker, ...
                                         'InGroup', Inputs.WrapperGroup(1,2) );

        info.TopSegInd = [info.TopSegInd, segind];

        lastbotnodeid = info.OuterNodeIDs(3);

        % Seg at bottom
        [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                         lastbotnodeid, ...
                                         botnodeid, ...
                                         'BoundaryMarker', tbboundarymarker, ...
                                         'InGroup', Inputs.WrapperGroup(1,2));

        info.BottomSegInd = [info.BottomSegInd, segind];
        
        % Add a node at the mid-point of the wrapper
        [FemmProblem, ~, midnodeid] = addnodes_mfemm (FemmProblem, ...
                                            outert, ...
                                            Inputs.NPolePairs*zpole, ...
                                            'InGroup', Inputs.WrapperGroup(1,2));

        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                            botnodeid, ...
                                            midnodeid, ...
                                            'InGroup', Inputs.WrapperGroup(1,2) );

        % Seg joining top and bottom (the two most recently added nodes)
        FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                            midnodeid, ...
                                            topnodeid, ...
                                            'InGroup', Inputs.WrapperGroup(1,2) );

        info.OuterCentres(1,:) = [outert-wrapperthickness(1,2)/2, Inputs.NPolePairs*zpole ];
        
    else
        % Set the region thickness to be exactly zero so this can be tested
        % later
        wrapperthickness(1,2) = 0;
        topnodeid = info.OuterNodeIDs(2);
        botnodeid = info.OuterNodeIDs(3);
    end
    
    % now add all subsequent right hand wrappers
    for i = 2:size(wrapperthickness, 1)
        
        if wrapperthickness(i,2) > Inputs.Tol
            
            outert = outert + wrapperthickness(i,2);
            
            lastbotnodeid = botnodeid;
            lasttopnodeid = topnodeid;

            % First node is to right of second node in 'nodes' matrix. this is at
            % the bottom of the sim
            [FemmProblem, ~, botnodeid] = addnodes_mfemm (FemmProblem, ...
                                         outert, ...
                                         0, ...
                                         'InGroup', Inputs.WrapperGroup(i,2));


            % Second node is to left of penultimate node in 'nodes' matrix
            [FemmProblem, ~, topnodeid] = addnodes_mfemm (FemmProblem, ...
                                         outert, ...
                                         Inputs.NPolePairs*2*zpole, ...
                                        'InGroup', Inputs.WrapperGroup(i,2));

            if Inputs.AddPeriodicBoundaries
                % add a new periodic boundary for the top and bottom of the
                % region
                [FemmProblem, info.BoundaryInds(end+1)] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

                tbboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;
            else
                tbboundarymarker = '';
            end

            % Seg with Periodic boundary at top
            [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                             lasttopnodeid, ...
                                             topnodeid, ...
                                             'BoundaryMarker', tbboundarymarker, ...
                                             'InGroup', Inputs.WrapperGroup(i,2) );

            info.TopSegInd = [info.TopSegInd, segind];
            
            % Seg at bottom
            [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                             lastbotnodeid, ...
                                             botnodeid, ...
                                             'BoundaryMarker', tbboundarymarker, ...
                                             'InGroup', Inputs.WrapperGroup(i,2) );

            info.BottomSegInd = [info.BottomSegInd, segind];
            
            [FemmProblem, ~, midnodeid] = addnodes_mfemm (FemmProblem, ...
                                            outert, ...
                                            Inputs.NPolePairs*zpole, ...
                                            'InGroup', Inputs.WrapperGroup(i,2));

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                             botnodeid, ...
                                             midnodeid, ...
                                             'InGroup', Inputs.WrapperGroup(i,2) );

            % Seg joining top and bottom (the two most recently added nodes)
            FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                             midnodeid, ...
                                             topnodeid, ...
                                             'InGroup', Inputs.WrapperGroup(i,2) );
            
            info.OuterCentres(i,:) = [outert-wrapperthickness(i,2)/2, Inputs.NPolePairs*zpole];
            
        else
            wrapperthickness(i,2) = 0;
        end
        
    end
    

end


