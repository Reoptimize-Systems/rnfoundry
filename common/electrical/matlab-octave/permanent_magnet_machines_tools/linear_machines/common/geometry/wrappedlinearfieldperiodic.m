function [FemmProblem, wrapperthickness, info] = wrappedlinearfieldperiodic (FemmProblem, tmag, zmag, zs, toffset, tsvc, tsve, zsvi, zsvo, wrapperthickness, varargin)
% creates a FemmProblem geometry of a radial section containing two magnets
% optionally wrapped with additional layers and with periodic edges
%
%
% Syntax
%
% [FemmProblem, wrapperthickness, innercentres, outercentres, nodeids] = ...
%   wrappedannularsecmagaperiodic(FemmProblem, thetapole, thetamag, rmag, ...
%                                 roffset, pos, wrapperthickness)
%
% [...] = wrappedannularsecmagaperiodic(..., 'param', value, ...)
%
% Description
%
% wrappedannularsecmagaperiodic creates a periodic geometry of two magnets
% with spaces in between, with a base position shown in the figure below.
% In addition, any number of annular sectors can be added either to inside
% or outside of the main region (like wrappers for the main region).
%       
%
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
%
% Inputs
%
%  FemmProblem - FemmProblem structure to which the geometry will be added
%, ts, zmag, zs, toffset, tsvc, tsve, zsvi, zsvo, wrapperthickness
%
%  tmag - magnet width in radians
%
%  rmag - radial thickness of the magnets
%
%  toffset - radial displacement of the magnet centers from the center
%
%  pos - the angular position of the magnets
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
%  'NPolePairs' - 
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
%    The default is {'theta', 'theta+180'}, resulting in radially
%    magnetized magnets of opposite polarity.
%
%  'MagnetMaterial' - 
%
%  'MagnetGroup' - 
%
%  'SpaceMaterial' - 
%
%  'SpaceGroup' - 
%
%  'Tol' - 
%
%  'MeshSize' - 
%
%  'Flip' - 
%
% Output
%
%  FemmProblem - 
%
%  wrapperthickness - 
%
%  info - structure containing information about the problem drawing
%
%

    Inputs.MagDirections = {'theta', 'theta+180'};
    Inputs.MagnetMaterial = 1;
    Inputs.MagnetGroup = 0;
    Inputs.SpacerMaterial = 1;
    Inputs.SpaceGroup = 0; 
    Inputs.WrapperGroup = 0;
    Inputs.Tol = 1e-5;
    Inputs.MeshSize = -1;
    Inputs.NPolePairs = 1;
    Inputs.Flip = false;
    
    % parse the input arguments
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isscalar(Inputs.WrapperGroup)
        Inputs.WrapperGroup = repmat (Inputs.WrapperGroup, size(wrapperthickness));
    elseif ~samesize(wrapperthickness, Inputs.WrapperGroup)
        error ('RENEWNET:wrappedlinearfield:badwrapper', ...
            'Your must supply either a scalar wrapper group, or a vector the smame size as the number of wrappers.')
    end
    
    if isnumeric(Inputs.MagDirections) && isvector(Inputs.MagDirections)
        Inputs.MagDirections = {Inputs.MagDirections(1), Inputs.MagDirections(2)};
    end
    
%     % remove input fields specific to this function and pass the
%     % remaining fields to planarrectmagaperiodic as a series of p-v pairs
%     planarrectmagaperiodicInputs = rmfield(Inputs, {'WrapperGroup'});
%     
%     % Convert the inputs structure to a series of p-v pairs
%     planarrectmagaperiodicInputs = struct2pvpairs(planarrectmagaperiodicInputs);
    
    [FemmProblem, nodes, links, info] = linearfieldperiodic (FemmProblem, tmag, zmag, zs, toffset, tsvc, tsve, zsvi, zsvo, ...
                'MagDirections', Inputs.MagDirections, ...
                'MagnetMaterial', Inputs.MagnetMaterial, ...
                'MagnetGroup', Inputs.MagnetGroup, ...
                'SpacerMaterial', Inputs.SpacerMaterial, ...
                'SpaceGroup', Inputs.SpaceGroup, ...
                'Tol', Inputs.Tol, ...
                'MeshSize', Inputs.MeshSize, ...
                'NPolePairs', Inputs.NPolePairs, ...
                'Flip', Inputs.Flip);
    
    zpole = zmag + zs;
    
    % now add the back iron nodes and links, depending on their thicknesses
    if size(wrapperthickness,2) < 2
        wrapperthickness = [wrapperthickness, wrapperthickness];
    elseif size(wrapperthickness,2) > 2
        error('RENEWNET:wrappedlinearfield:badwrapper', ...
            'wrapperthickness must be a scaler or a (1 x 2) vector or (n x 2) matrix')
    end
    
    elcount = elementcount_mfemm(FemmProblem);
    
    info.InnerCentres = zeros(size(wrapperthickness)) * NaN;
    
    % wrapperthickness(1,1) is the first inner region thickness
    if wrapperthickness(1,1) > Inputs.Tol
        % add the nodes and segments for the inner side
        %
        
        if wrapperthickness(1,1) > toffset
            error('RENEWNET:wrappedlinearfield:badwrapper', ...
                'wrapper thickness cannot be greater than magnet leftmost position.');
        end
        
        innert = toffset - tmag/2 - wrapperthickness(1,1);
        
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

        % add a new periodic boundary for the top and bottom of the region
        [FemmProblem, info.BoundaryInds(end+1)] = addboundaryprop_mfemm (FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);

        botboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;

        % Seg with Periodic boundary at top
        [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                        info.OuterNodeIDs(1), ...
                                        topnodeid, ...
                                        'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name, ...
                                        'InGroup', Inputs.WrapperGroup(1,1));

        info.TopSegInd = [info.TopSegInd, segind];

        lastbotnodeid = info.OuterNodeIDs(4);
        
        % Seg at bottom
        [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                        lastbotnodeid, ...
                                        botnodeid, ...
                                        'BoundaryMarker', botboundarymarker, ...
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

            % add a new periodic boundary for the top and bottom of the region
            [FemmProblem, info.BoundaryInds(end+1)] = addboundaryprop_mfemm(FemmProblem, 'Left Wrap Annular Sec Mags Periodic', 4);

            botboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;

            % Seg with Periodic boundary at top
            [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                            lasttopnodeid, ...
                                            topnodeid, ...
                                            'BoundaryMarker', botboundarymarker, ...
                                            'InGroup', Inputs.WrapperGroup(i,1));

            info.TopSegInd = [info.TopSegInd, segind];
            
            % Seg at bottom
            [FemmProblem, segind] = addsegments_mfemm (FemmProblem, ...
                                            lastbotnodeid, ...
                                            botnodeid, ...
                                            'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name, ...
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
    
    % wrapperthickness(1,2) is the first outer region thickness
    if wrapperthickness(1,2) > Inputs.Tol
        
        outert = toffset + tmag/2 + wrapperthickness(1,2);

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

        % add a new periodic boundary for the top and bottom of the
        % region
        [FemmProblem, info.BoundaryInds(end+1)] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

        botboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;

        % Seg with Periodic boundary at top
        [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                         info.OuterNodeIDs(2), ...
                                         topnodeid, ...
                                         'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name, ...
                                         'InGroup', Inputs.WrapperGroup(1,2) );

        info.TopSegInd = [info.TopSegInd, segind];

        lastbotnodeid = info.OuterNodeIDs(3);

        % Seg at bottom
        [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                         lastbotnodeid, ...
                                         botnodeid, ...
                                         'BoundaryMarker', botboundarymarker, ...
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

            % add a new periodic boundary for the top and bottom of the
            % region
            [FemmProblem, info.BoundaryInds(end+1)] = addboundaryprop_mfemm(FemmProblem, 'Right Wrap Annular Sec Mags Periodic', 4);

            botboundarymarker = FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name;

            % Seg with Periodic boundary at top
            [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                             lasttopnodeid, ...
                                             topnodeid, ...
                                             'BoundaryMarker', FemmProblem.BoundaryProps(info.BoundaryInds(end)).Name, ...
                                             'InGroup', Inputs.WrapperGroup(i,2) );

            info.TopSegInd = [info.TopSegInd, segind];
            
            % Seg at bottom
            [FemmProblem, segind] = addsegments_mfemm ( FemmProblem, ...
                                             lastbotnodeid, ...
                                             botnodeid, ...
                                             'BoundaryMarker', botboundarymarker, ...
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


