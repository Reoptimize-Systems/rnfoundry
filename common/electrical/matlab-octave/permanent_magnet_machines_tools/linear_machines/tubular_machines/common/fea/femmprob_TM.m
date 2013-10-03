function FemmProblem = femmprob_TM(design, varargin)
% RunFEMMSimNew_ACTIAM: A function for generating simulations of a single
% pole of the Slotless Tubular Permanent Magnet machine in the FEMM finite
% element analysis program. 
%
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/design.Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of design.Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RsoVRm - scalar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   mode - Scalar value specifying what simulation type is to be performed.
%          Can have the values 0, 1, 2 or 3.
%          
%          0: Magnet in center with no steel removed
%
%          1: Magnet in centre with steel removed in accordance with the
%          ratios Rs2VHmag, Rs1VHmag, Ws2VhalfWs and Ws1VhalfWs.
%
%          2: Steel in centre with no steel removed
%
%          3: Steel in centre with steel removed in accordance with the
%          ratios Rs2VHmag, Rs1VHmag, Ws2VhalfWs and Ws1VhalfWs.
%
%          If mode (and subsequent arguments) are omitted the default is
%          zero. If the arguments following mode are omitted they default
%          to 0.5 in all cases.
%
%   Rs2VHmag - Ratio of the height of the air region from the inner edge of the steel
%              piece from the outer radius of the shaft to the total height
%              of the magnet, i.e. Rs1 / (Rm - Rso), see diagram.
%
%   Rs1VHmag - Ratio of the height of the air region from the centre of the steel
%              piece from the outer radius of the shaft to the total height
%              of the magnet, i.e. Rs2 / (Rm - Rso), see diagram.
% 
%   Ws2VhalfWs - Ratio of the width of the air region from the centre of
%                the steel piece to half of the total width of the steel
%                i.e. 2* Ws2 / Ws, see diagram.
% 
%   Ws1VhalfWs - Ratio of the width of the air region from the centre of
%                the steel piece to half of the total width of the steel
%                i.e. 2* Ws1 / Ws, see diagram.
%
%    |         _______
%    |        /       |
%    |_______/        |
%    |                |
%    |________________|
%    |                |
%    |      Hmag      |
%    |<-------------->|
%    |________________|
%    |<-----> Rs2     | ^
% Ws1|_______   Ws2   | :
% ^  |       \  ^     | : half Ws
% :  |        \_:_____| ;
%    <--------> Rs1
    
%     newdocument(0);
    
    Inputs.CurrentDensities = [0,0,0];
    Inputs.CoilCurrents = [];
    Inputs.TMType = 'ACTM';
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % The current density will be in A/m^2, femm accepts a value in
    % MA/m^2, so we must adjust to account for this
	Inputs.CurrentDensities = Inputs.CurrentDensities ./ 1e6;
    
    if ~all(isfield(design, {'Rs1VHmag', 'Ws2VhalfWs', 'Ws1VhalfWs', 'Ws1VhalfWs'}))
        
        design.Rs2VHmag = 0.5;
        design.Rs1VHmag = 0.5;
        design.Ws2VhalfWs = 0.25;
        design.Ws1VhalfWs = 0.25;
        
        if ~isfield(design, 'mode')
            design.mode = 0;
        end
        
        FemmProblem = newproblem_mfemm('axi', 'MinAngle', 30, 'LengthUnits','meters');
        
    else
        % use smaller angle constraint
        FemmProblem = newproblem_mfemm('axi', 'MinAngle', 5, 'LengthUnits','meters');
    end
    
    if ~isfield(design, 'HcMag')
        design.HcMag = mgoe2hc(40);
    end
    
    if design.mode == 0 || design.mode(1) == 1
        design.Wp = design.Rm * design.WpVRm;
        design.Ws = design.Wp * (1 - design.WmVWp);
        design.Wm = design.Wp * design.WmVWp;
        design.Ws1 = design.Ws1VhalfWs * (design.Ws/2);
        design.Ws2 = design.Ws2VhalfWs * (design.Ws/2);
        design.Wm1 = (1-design.Ws1VhalfWs) * (design.Wm/2) + (design.Ws/2);
        design.Wm2 = (1-design.Ws2VhalfWs) * (design.Wm/2) + (design.Ws/2);
        design.Rso = design.RsoVRm * design.Rm;
        design.Rs1 = design.Rso + design.Rs1VHmag * (design.Rm-design.Rso);
        design.Rs2 = design.Rso + design.Rs2VHmag * (design.Rm-design.Rso);
    else
        design.Wp = design.Rm * design.WpVRm;
        design.Ws = design.Wp * (design.WmVWp); % is actually design.Wm
        design.Wm = design.Wp * (1-design.WmVWp); % Is actually design.Ws
        design.Ws1 = design.Ws1VhalfWs * (design.Ws/2);
        design.Ws2 = design.Ws2VhalfWs * (design.Ws/2);
        design.Wm1 = (1-design.Ws1VhalfWs) * (design.Wm/2) + (design.Ws/2);
        design.Wm2 = (1-design.Ws2VhalfWs) * (design.Wm/2) + (design.Ws/2);
        design.Rso = design.RsoVRm * design.Rm;
        design.Rs1 = design.Rso + design.Rs1VHmag * (design.Rm-design.Rso);
        design.Rs2 = design.Rso + design.Rs2VHmag * (design.Rm-design.Rso);
    end
    
%     design.Ro = RoVRm * design.Rm;
%     design.Ra = RaVRo * design.Ro;
    
    %% Problem Setup

    % Add some boundary types
%     mi_addboundprop('Pros A', 0, 0, 0, 90, 0, 0, 0, 0, 0)
    [ FemmProblem, prosAboundind, prosAboundname ] = addboundaryprop_mfemm(FemmProblem, 'Pros A', 0, 'Phi', 90);
    
    for i = 1:7
%         mi_addAntiPeriodicBoundary(['antiPeriodic ' num2str(i)]);
        [ FemmProblem, APboundinds(i), APboundnames{i} ] = addboundaryprop_mfemm(FemmProblem, 'antiPeriodic', 5);
    end
    
    % Add some materials
    % Steel
%     steelname = mi_addsteel;
    steelname = '1117 Steel';
    [FemmProblem, steelind] = addmaterials_mfemm(FemmProblem, steelname);

    % Copper wire
    wirename = 'wire';
    [FemmProblem, wireind] = addmagnetwire_mfemm(FemmProblem, wirename, design.Dc);
    
    % Next magnets, 40 MGOe
    [FemmProblem, magname, magmatind] = addmagnet_mfemm(FemmProblem, design.HcMag);
    
    switch Inputs.TMType
        
        case 'ACTM'
            
            sheathmatname = 'Air';
            coilspacematname = 'Air';
            
        case 'STPMSM'
            
            sheathmatname = steelname;
            coilspacematname = 'Air';
            
        case 'TPMSM'
            
            sheathmatname = steelname;
            coilspacematname = steelname;
            
        otherwise
            
            error('Unknown machine type');
            
    end
    
    % Now draw some nodes
    
    % Bottom left node, corner of sim
	rzCoords(1,:) = [0,0];
    % Bottom node of outer shaft radius
	rzCoords(2,:) = [design.Rso, 0];
    % Bottom node of inner radius of steel
	rzCoords(3,:) = [design.Rs1, 0];
    % Bottom node of outer steel radius
	rzCoords(4,:) = [design.Rm, 0];
    % Bottom node of inner coil radius
    rzCoords(5,:) = [design.Ri, 0];
    % Bottom node of inner sheath Radius
    rzCoords(6,:) = [design.Ro, 0];
    % Bottom node of outer sheath radius
    rzCoords(7,:) = [design.Ra, 0];
    % Bottom right node of sim marking end of air
	rzCoords(8,:) = [4 * design.Ra, 0];
    % Node marking air-steel junction at the shaft
	rzCoords(9,:) = [design.Rso, design.Ws1];
    % Node marking air-steel junction closer to design.Rm
	rzCoords(10,:) = [design.Rs2, design.Ws2];
    % Node marking start of mag close to shaft
	rzCoords(11,:) = [design.Rso, design.Ws/2];
    % Node marking edge of mag at design.Rm
	rzCoords(12,:) = [design.Rm, design.Ws/2];
    
    rzCoords(13,:) = [design.Rso, design.Wm1];
    rzCoords(14,:) = [design.Rs2, design.Wm2];
    
    [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, rzCoords(:,1), rzCoords(:,2));
    
    % Now add some segments, those that can be mirrored
    segPropStructArray = struct('BoundaryMarker', APboundnames, 'InGroup', 1);
                    
    n0 = [1, 2, 3, 4, 5, 6, 7] - 1;
    n1 = [2, 3, 4, 5, 6, 7, 8] - 1;
	[FemmProblem, APseginds ] = addsegments_mfemm(FemmProblem, n0, n1, segPropStructArray);

	n0 = [9, 10, 11, 11, 12, 9, 11, 13] - 1;
    n1 = [10, 3,  9, 12,  4, 2, 13, 14] - 1;
	[FemmProblem, internalseginds ] = addsegments_mfemm(FemmProblem, n0, n1, 'InGroup', 1);

    FemmProblem = mirrorsegments_mfemm(FemmProblem, ...
                                       1:numel(FemmProblem.Segments), ...
                                       design.FEMMTol, ...
                                       'TwoPoints', [-1, design.Wp/2, 1, design.Wp/2]);

    % Now add missing segments and nodes
    rzCoords(15,:) = [design.Rs1, design.Wp/2];

    FemmProblem = addnodes_mfemm(FemmProblem, rzCoords(15,1), rzCoords(15,2));
    
    newnodeids = findnode_mfemm(FemmProblem, ...
                                [rzCoords(1,1), design.Wp-rzCoords(1,2);
                                 rzCoords(7,1), design.Wp-rzCoords(7,2);
                                 rzCoords(8,1), design.Wp-rzCoords(8,2);
                                 rzCoords(12,1), design.Wp-rzCoords(12,2)
                                 rzCoords(13,1), design.Wp-rzCoords(13,2);
                                 rzCoords(15,:);
                                 rzCoords(14,1), design.Wp-rzCoords(14,2); ]);
    
    segPropStructArray = struct('BoundaryMarker', {'', '', prosAboundname, '', '', '', ''}, ...
                                'InGroup', 1);
    
    n0 = [ [ 1, 7, 8, 12, 13, 14 ] - 1, newnodeids(6)];
	[FemmProblem, APseginds ] = addsegments_mfemm(FemmProblem, n0, newnodeids, segPropStructArray);

    % Add labels 
    meshSz = min([sqrt(design.Ra^2 + design.Wp^2) ./ 150, (design.Ro - design.Rm)/10]);
    
    if design.mode == 0 || design.mode == 1
        
        if design.mode == 0
            %  Magnet in center with no steel removed
            innerdiscmatname = steelname;
        elseif design.mode == 1
            % Magnet in centre with steel removed
            innerdiscmatname = 'Air';
        end
        
        % translator internal steel labels
        % bottom disc
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso + ((design.Rs2-design.Rso)/2), ...
                                          design.Ws2 / 2, ...
                                          'BlockType', innerdiscmatname, ...
                                          'MaxArea', meshSz);
        % top disc                 
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso + ((design.Rs2-design.Rso)/2), ...
                                          design.Wp - (design.Ws2 / 2), ...
                                          'BlockType', innerdiscmatname, ...
                                          'MaxArea', meshSz);
        
        % translator outer disc parts
        % bottom disc
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs2 + ((design.Rm-design.Rs2)/2), ...
                                          design.Ws2 + (((design.Ws/2)-design.Ws2)/2), ...
                                          'BlockType', steelname, ...
                                          'MaxArea', meshSz);
        % top disc
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs2 + ((design.Rm-design.Rs2)/2), ...
                                          design.Wp - (design.Ws2 + (((design.Ws/2)-design.Ws2)/2)), ...
                                          'BlockType', steelname, ...
                                          'MaxArea', meshSz);

        % Magnet labels

        % outer
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs1+((design.Rm-design.Rs1)/2), ...
                                          design.Wp/2, ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', 90); 
        % inner
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso+((design.Rs1-design.Rso)/2), ...
                                          design.Wp/2, ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', 90);  
        
    elseif design.mode == 2 || design.mode == 3
        % Steel in centre 
        
        if design.mode == 2
            %  Steel in center with no steel removed
            innerdiscmatname = steelname;
        elseif design.mode == 3
            % Steel in centre with steel removed
            innerdiscmatname = 'Air';
        end
        
        % translator mag label
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso+((design.Rs2-design.Rso)/2), ...
                                          design.Ws2 / 2, ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', 90, ...
                                          'InGroup', 1);
                                      
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso+((design.Rs2-design.Rso)/2), ...
                                          design.Wp - (design.Ws2 / 2), ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', -90, ...
                                          'InGroup', 1);
        
        % Steel Label
        % bottom disc
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs2 + ((design.Rm-design.Rs2)/2), ...
                                          design.Ws2 + (((design.Ws/2)-design.Ws2)/2), ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', 90, ...
                                          'InGroup', 1);
        % top disc
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs2 + ((design.Rm-design.Rs2)/2), ...
                                          design.Wp - (design.Ws2 + (((design.Ws/2)-design.Ws2)/2)), ...
                                          'BlockType', magname, ...
                                          'MaxArea', meshSz, ...
                                          'MagDir', -90, ...
                                          'InGroup', 1);
                                      
        % outer
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rs1+((design.Rm-design.Rs1)/2), ...
                                          design.Wp/2, ...
                                          'BlockType', steelname, ...
                                          'MaxArea', meshSz); 
        % inner
        FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                          design.Rso+((design.Rs1-design.Rso)/2), ...
                                          design.Wp/2, ...
                                          'BlockType', innerdiscmatname, ...
                                          'MaxArea', meshSz); 
        
    end  
    
    %%
    
    if ((design.Wp/3)-(design.WcVWp*design.Wp)) < design.FEMMTol
        % We will draw the coils with no space between them
        rzCCoords(1,:) = [design.Ri, design.Wp/3];
        rzCCoords(2,:) = [design.Ro, design.Wp/3];
        rzCCoords(3,:) = [design.Ri, 2*design.Wp/3];
        rzCCoords(4,:) = [design.Ro, 2*design.Wp/3];
        
        [ FemmProblem, nodeinds, cnodeids ] = addnodes_mfemm(FemmProblem, rzCCoords(:,1), rzCCoords(:,2));

        tnodeids = findnode_mfemm(FemmProblem, [design.Ri, design.Wp; ...
                                                design.Ro, design.Wp ] );
                                            
        n0 =         [[5, 6]-1, cnodeids([1, 1, 2, 3]), tnodeids([1, 2])];
        n1 = cnodeids([1, 2,              2, 3, 4, 4,             3, 4]);
        FemmProblem = addsegments_mfemm(FemmProblem, n0, n1, 'InGroup', 2);

    else
        % We will draw the coils and air space (could be iron to make slotted machine) 
        rzCCoords(1,:) = [design.Ri, ((design.Wp/3)-design.Wc)/2];
        rzCCoords(2,:) = [design.Ro, ((design.Wp/3)-design.Wc)/2];
        rzCCoords(3,:) = [design.Ri, design.Wc+((design.Wp/3)-design.Wc)/2];
        rzCCoords(4,:) = [design.Ro, design.Wc+((design.Wp/3)-design.Wc)/2];
        rzCCoords(5,:) = [design.Ri, design.Wc+3*((design.Wp/3)-design.Wc)/2];
        rzCCoords(6,:) = [design.Ro, design.Wc+3*((design.Wp/3)-design.Wc)/2];
        
        [ FemmProblem, odeinds, nodeids ] = addnodes_mfemm(FemmProblem, rzCCoords(:,1), rzCCoords(:,2));
        
        n0 = nodeids([1, 1, 2, 3, 3, 4, 5]);     
        n1 = nodeids([2, 3, 4, 4, 5, 6, 6]);
        [FemmProblem, coilseginds ] = addsegments_mfemm(FemmProblem, n0, n1, 'InGroup', 2);
        
        % Mirror the coil
        FemmProblem = mirrorsegments_mfemm(FemmProblem, ...
                                   coilseginds, ...
                                   design.FEMMTol, ...
                                   'TwoPoints', [0, design.Wp/2, 4*design.Rm, design.Wp/2]);
        
        % Join Coil Parts
        nodeids = findnode_mfemm(FemmProblem, [rzCCoords(5,1), rzCCoords(5,2); ...
                                               rzCCoords(6,1), rzCCoords(6,2); ...
                                               rzCCoords(5,1), design.Wp-rzCCoords(5,2); ...
                                               rzCCoords(6,1), design.Wp-rzCCoords(6,2) ]);
        
        n0 = nodeids([1, 2]); 
        n1 = nodeids([3, 4]);
        [FemmProblem, coilseginds ] = addsegments_mfemm(FemmProblem, n0, n1, 'InGroup', 2);
        
        % join top and bottoms to top and bottom
        nodeids = findnode_mfemm(FemmProblem, [rzCCoords(1,1), rzCCoords(1,2); ...
                                               rzCCoords(2,1), rzCCoords(2,2); ...
                                               rzCCoords(1,1), design.Wp-rzCCoords(1,2); ...
                                               rzCCoords(2,1), design.Wp-rzCCoords(2,2); ...
                                               design.Ri, 0;
                                               design.Ro, 0;
                                               design.Ri, design.Wp;
                                               design.Ro, design.Wp ]);
                                           
        n0 = nodeids([1, 2, 3, 4]); 
        n1 = nodeids([5, 6, 7, 8]);
        [FemmProblem, coilseginds ] = addsegments_mfemm(FemmProblem, n0, n1, 'InGroup', 2);
                    
        % labels for spaces between coils
        mesh = design.Rm / 200;

        [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                      design.Ri+(design.Ro-design.Ri)/2, ...
                                                      design.Wp/3, ...
                                                      'BlockType', coilspacematname, ...
                                                      'MaxArea', mesh);
        
        [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                      design.Ri+(design.Ro-design.Ri)/2, ...
                                                      2*design.Wp/3, ...
                                                      'BlockType', coilspacematname, ...
                                                      'MaxArea', mesh);
        
        [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                      design.Ri+(design.Ro-design.Ri)/2, ...
                                                      ((design.Wp/3)-design.Wc)/4, ...
                                                      'BlockType', coilspacematname, ...
                                                      'MaxArea', mesh);
        
        [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                      design.Ri+(design.Ro-design.Ri)/2, ...
                                                      design.Wp-((design.Wp/3)-design.Wc)/4, ...
                                                      'BlockType', coilspacematname, ...
                                                      'MaxArea', mesh);
        
    end
    
    % Coil labels
    
    if ~isempty(Inputs.CoilCurrents)

        FemmProblem = addcircuit_mfemm(FemmProblem, 'Coil C', 'TotalAmps_re', Inputs.CoilCurrents(1));

        blockPropStructArray = struct('BlockType', wirename, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 30, ...
                                      'InCircuit', 'Coil C', ...
                                      'Turns', design.CoilTurns);

    else
        
        % Add a solid copper region with the specified current density
        cuMaterial = newmaterial_mfemm('copper J C', ... 
                                       'Mu_x',1 , ...
                                       'Mu_y', 1, ...
                                       'H_c', 0, ...
                                       'J_re', Inputs.CurrentDensities(1), ...
                                       'Sigma', 58, ...
                                       'Density', 8600);
                                    
        [FemmProblem, cumatind] = addmaterials_mfemm(FemmProblem, cuMaterial);
        
        blockPropStructArray = struct('BlockType', FemmProblem.Materials(cumatind).Name, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 30);
    end
    
    [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                  design.Ri+(design.Ro-design.Ri)/2, ...
                                                  design.Wp/6, ...
                                                  blockPropStructArray);
    
    if ~isempty(Inputs.CoilCurrents)
        
        FemmProblem = addcircuit_mfemm(FemmProblem, 'Coil B', 'TotalAmps_re', Inputs.CoilCurrents(2));

        blockPropStructArray = struct('BlockType', wirename, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 20, ...
                                      'InCircuit', 'Coil B', ...
                                      'Turns', design.CoilTurns);
    else
        % Add a solid copper region with the specified current density        
        cuMaterial = newmaterial_mfemm('copper J B', ... 
                                       'Mu_x',1 , ...
                                       'Mu_y', 1, ...
                                       'H_c', 0, ...
                                       'J_re', Inputs.CurrentDensities(2), ...
                                       'Sigma', 58, ...
                                       'Density', 8600);
                                    
        [FemmProblem, cumatind] = addmaterials_mfemm(FemmProblem, cuMaterial);
        
        blockPropStructArray = struct('BlockType', FemmProblem.Materials(cumatind).Name, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 20);
    end
 
    [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                  design.Ri+(design.Ro-design.Ri)/2, ...
                                                  0.5*design.Wp, ...
                                                  blockPropStructArray);
    
    if ~isempty(Inputs.CoilCurrents)

        FemmProblem = addcircuit_mfemm(FemmProblem, 'Coil A', 'TotalAmps_re', Inputs.CoilCurrents(3));

        blockPropStructArray = struct('BlockType', wirename, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 10, ...
                                      'InCircuit', 'Coil A', ...
                                      'Turns', design.CoilTurns);
    else
        % Add a solid copper region with the specified current density
        cuMaterial = newmaterial_mfemm('copper J A', ... 
                                       'Mu_x',1 , ...
                                       'Mu_y', 1, ...
                                       'H_c', 0, ...
                                       'J_re', Inputs.CurrentDensities(3), ...
                                       'Sigma', 58, ...
                                       'Density', 8600);
                                    
        [FemmProblem, cumatind] = addmaterials_mfemm(FemmProblem, cuMaterial);
        
        blockPropStructArray = struct('BlockType', FemmProblem.Materials(cumatind).Name, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 10);
        
    end
    
    [FemmProblem, blockind] = addblocklabel_mfemm(FemmProblem, ...
                                                  design.Ri+(design.Ro-design.Ri)/2, ...
                                                  design.Wp-(design.Wp/6), ...
                                                  blockPropStructArray);
    
    
    % air label for shaft 
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      rzCoords(2,1)/2, ...
                                      design.Wp/2, ...
                                      'BlockType', 'Air', ...
                                      'MaxArea', meshSz * 2, ...
                                      'InGroup', 1);
    
    % air gap label
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      design.Rm+(design.Ri-design.Rm)/2, ...
                                      design.Wp/2, ...
                                      'BlockType', 'Air', ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 1);
    
    % outer air region label
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      2*design.Ra, ...
                                      design.Wp/2, ...
                                      'BlockType', 'Air', ...
                                      'MaxArea', meshSz * 10, ...
                                      'InGroup', 1);
    
    % Armature Sheath label
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      design.Ro+((design.Ra-design.Ro)/2), ...
                                      design.Wp/2, ...
                                      'BlockType', sheathmatname, ...
                                      'MaxArea', meshSz, ...
                                      'InGroup', 1);
	
end