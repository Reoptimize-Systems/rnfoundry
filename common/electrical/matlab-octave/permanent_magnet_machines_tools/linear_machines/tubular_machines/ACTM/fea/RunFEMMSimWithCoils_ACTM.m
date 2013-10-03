function [Dc, filename] = RunFEMMSimWithCoils_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RsoVRm, WcVWp, Rm, Ntot, kcufill, J, mode, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% RunFEMMSimNew_ACTM: A function for generating simulations of a single
% pole of the Slotless Tubular Permanent Magnet machine in the FEMM finite
% element analysis program. 
%
% Syntax
% 
% Dc = RunFEMMSimWithCoils_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RsoVRm, ... 
%                      WcVWp, Rm, Ntot, kcufill, J, mode, Rs2VHmag, ...
%                      Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RiVRm - scalar value of Ri/Rm, inner coil radius to Rm ratio
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   Ntot - Number of turns in the winding
%
%   kcufill - copper fill-factor of the winding
%
%   J - current density in coil, see information on mode for details
%
%   mode - Scalar or (1 x 2) vector specifying what simulation type is to be
%          performed. The first value, mode(1,1), determines the position
%          and construction of the translator/field and can take the
%          following values:
%
%          mode(1,1): Can have the values 0, 1, 2 or 3.
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
%          If present, the second value mode(1,2) determines whether the
%          magnets are present in the simulation. It may be desireable to
%          remove the magnets in order to determine the inductance of the
%          coils. mode(1,2) can be either 0 or 1. if zero, the magnets are
%          not present, if 1 they are present. The default is for the
%          magnets to be present. 
%
%          If present, the third value of mode determines whether the coils
%          are made up of actual circuits of blocks of copper with an
%          applied current density. If mode(3) == 0 blocks of copper are
%          used, the current density in J is applied to the whole block. If
%          mode(3) == 1 wire coils are used, the current densities in J are
%          used to calculate the current in the wire so that 
%          I = J * conductor area
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
% Ws1|_______   Ws2   | ¦
% ^  |       \  ^     | ¦ half Ws
% ¦  |        \_¦_____| ¦
%    <--------> Rs1

    newdocument(0);
    
    if nargin < 12
        
        Rs2VHmag = 0.5;
        Rs1VHmag = 0.5;
        Ws2VhalfWs = 0.25;
        Ws1VhalfWs = 0.25;
        if nargin == 10
            mode = [0, 1];
        end
        
        mi_probdef(0,'meters','axi',1e-8, 0, 30, 0);
        
    elseif (nargin > 11 && nargin < 15) || nargin < 10
        
        error('Invalid number of input arguments')
        
    else
        % use smaller angle constraint
        mi_probdef(0,'meters','axi',1e-8, 0, 5, 0);
        
    end
    
    if kcufill == -1
        kcufill = 0.6;
    end
    
    if size(mode,2) < 2
        % mode(2) = 1 means the magnets are present, this is the default
        % condition
        mode(2) = 1;
        mode(3) = 0;
    elseif size(mode,2) < 3
        mode(3) = 0;
    end
    
    if mode(1) == 0 
        Wp = Rm * WpVRm;
        Ws = Wp * (1 - WmVWp);
        Wm = Wp * WmVWp;
        Ws1 = Ws1VhalfWs * (Ws/2);
        Ws2 = Ws2VhalfWs * (Ws/2);
        Wm1 = (1-Ws1VhalfWs) * (Wm/2) + (Ws/2);
        Wm2 = (1-Ws2VhalfWs) * (Wm/2) + (Ws/2);
        Rso = RsoVRm * Rm;
        Rs1 = Rso + Rs1VHmag * (Rm-Rso);
        Rs2 = Rso + Rs2VHmag * (Rm-Rso);
    else
        Wp = Rm * WpVRm;
        Ws = Wp * (WmVWp); % is actually Wm
        Wm = Wp * (1-WmVWp); % Is actually Ws
        Ws1 = Ws1VhalfWs * (Ws/2);
        Ws2 = Ws2VhalfWs * (Ws/2);
        Wm1 = (1-Ws1VhalfWs) * (Wm/2) + (Ws/2);
        Wm2 = (1-Ws2VhalfWs) * (Wm/2) + (Ws/2);
        Rso = RsoVRm * Rm;
        Rs1 = Rso + Rs1VHmag * (Rm-Rso);
        Rs2 = Rso + Rs2VHmag * (Rm-Rso);
    end
    
    Ri = RiVRm * Rm;
    Ro = RoVRm * Rm;
    Ra = 2.5 * Ro;
    Router = 4 * Ro;
    Wc = WcVWp * Wp;
    Hc = Ro - Ri;
    Dc = ConductorDiameter(Wc*Hc, kcufill, Ntot);
    
    %% Problem Setup
    % Setup some problem variables
%     mi_setUpProblem([], 'meters');
    % Add some boundary types
    mi_addboundprop('Pros A', 0, 0, 0, 90, 0, 0, 0, 0, 0)
    
    for i = 1:7
        mi_addAntiPeriodicBoundary(['antiPeriodic ' num2str(i)]);
    end
    
    % Add some materials
    % First Steel
    steelname = mi_addsteel;
    % Next magnets
    magnetname = mi_addmagnet;
    % Then Copper wire
    wirename = 'wire';
    mi_addmagnetwire(wirename, Dc, 1.7e-8);
    % Finally Air
    mi_addair;
    
    % Add wire
    mi_addwire;
    
    % Now draw some nodes
    
    % Bottom left node, corner of sim
	rzCoords(1,:) = [0,0];
    % Bottom node of outer shaft radius
	rzCoords(2,:) = [Rso, 0];
    % Bottom node of inner radius of steel
	rzCoords(3,:) = [Rs1, 0];
    % Bottom node of outer steel radius
	rzCoords(4,:) = [Rm, 0];
    % Bottom node of inner coil radius
    rzCoords(5,:) = [Ri, 0];
    % Bottom node of inner sheath Radius
    rzCoords(6,:) = [Ro, 0];
    % Bottom node of outer sheath radius
    rzCoords(7,:) = [Ra, 0];
    % Bottom right node of sim marking end of air
	rzCoords(8,:) = [Router, 0];
    % Node marking air-steel junction at the shaft
	rzCoords(9,:) = [Rso, Ws1];
    % Node marking air-steel junction closer to Rm
	rzCoords(10,:) = [Rs2, Ws2];
    % Node marking start of mag close to shaft
	rzCoords(11,:) = [Rso, Ws/2];
    % Node marking edge of mag at Rm
	rzCoords(12,:) = [Rm, Ws/2];

    % Nodes marking the cut out region in the magnet
    rzCoords(13,:) = [Rso, Wm1];
    rzCoords(14,:) = [Rs2, Wm2];
    
    for n = 1:size(rzCoords,1)

        mi_addnode(rzCoords(n,1), rzCoords(n,2));

    end
    
    % Now add some segments, those that can be mirrored
    
    % First the segments along the bottom of the sim
    segPropStructArray(1) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 1');
    segPropStructArray(2) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 2');
    segPropStructArray(3) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 3');
    segPropStructArray(4) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 4');
    segPropStructArray(5) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 5');
    segPropStructArray(6) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 6');
    segPropStructArray(7) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 7');
    
    mi_addsegment2([rzCoords(1,1),...
                    rzCoords(2,1),...
                    rzCoords(3,1),...
                    rzCoords(4,1),...
                    rzCoords(5,1),...
                    rzCoords(6,1),...
                    rzCoords(7,1)],...
                   [rzCoords(1,2),...
                    rzCoords(2,2),...
                    rzCoords(3,2),...
                    rzCoords(4,2),...
                    rzCoords(5,2),...
                    rzCoords(6,2),...
                    rzCoords(7,2)],...
                   [rzCoords(2,1),...
                    rzCoords(3,1),...
                    rzCoords(4,1),...
                    rzCoords(5,1),...
                    rzCoords(6,1),...
                    rzCoords(7,1),...
                    rzCoords(8,1)],...
                   [rzCoords(2,2),...
                    rzCoords(3,2),...
                    rzCoords(4,2),...
                    rzCoords(5,2),...
                    rzCoords(6,2),...
                    rzCoords(7,2),...
                    rzCoords(8,2)], segPropStructArray);
    
    % Next some standard segments joining up the translator parts
    segPropStructArray = mi_segpropstruct(8, 1);
    
    mi_addsegment2([rzCoords(9,1),...
                    rzCoords(10,1),...
                    rzCoords(11,1),...
                    rzCoords(11,1),...
                    rzCoords(12,1),...
                    rzCoords(9,1),...
                    rzCoords(11,1),...
                    rzCoords(13,1)],...
                   [rzCoords(9,2),...
                    rzCoords(10,2),...
                    rzCoords(11,2),...
                    rzCoords(11,2),...
                    rzCoords(12,2),...
                    rzCoords(9,2),...
                    rzCoords(11,2),...
                    rzCoords(13,2)],...
                    [rzCoords(10,1),...
                    rzCoords(3,1),...
                    rzCoords(9,1),...
                    rzCoords(11,1),...
                    rzCoords(4,1),...
                    rzCoords(2,1),...
                    rzCoords(12,1),...
                    rzCoords(14,1)],...
                   [rzCoords(10,2),...
                    rzCoords(3,2),...
                    rzCoords(9,2),...
                    rzCoords(11,2),...
                    rzCoords(4,2),...
                    rzCoords(2,2),...
                    rzCoords(12,2),...
                    rzCoords(14,2)], segPropStructArray);
                
    % Add labels which can be mirrored
    meshSz = sqrt((2*Rm)^2 + Wp^2) ./ 150;
        
    if mode(1) == 0
 
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, meshSz, 1);
        mi_addblocklabel2(Rso+((Rs2-Rso)/2),rzCoords(10,2)/2, blockPropStructArray);
        % Steel Label
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, meshSz, 1);
        mi_addblocklabel2(rzCoords(10,1)+((rzCoords(12,1)-rzCoords(10,1))/2),rzCoords(10,2)+((rzCoords(12,2)-rzCoords(10,2))/2), blockPropStructArray);
        
    elseif mode(1) == 1
        
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 1);
        mi_addblocklabel2(Rso+((Rs2-Rso)/2),rzCoords(10,2)/2, blockPropStructArray);
        % Steel Label
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, meshSz, 1);
        mi_addblocklabel2(rzCoords(10,1)+((rzCoords(12,1)-rzCoords(10,1))/2),rzCoords(10,2)+((rzCoords(12,2)-rzCoords(10,2))/2), blockPropStructArray);
        
    elseif mode(1) == 2 || mode(1) == 3
        
        if mode(2) == 1
            blockPropStructArray = mi_blockpropstruct(1, magnetname, 0, meshSz, 1, '', 0, 90);
            mi_addblocklabel2(Rso+((Rs2-Rso)/2),rzCoords(10,2)/2, blockPropStructArray);
            % Steel Label
            blockPropStructArray = mi_blockpropstruct(1, magnetname, 0, meshSz, 1, '', 0, 90);
            mi_addblocklabel2(rzCoords(10,1)+((rzCoords(12,1)-rzCoords(10,1))/2),rzCoords(10,2)+((rzCoords(12,2)-rzCoords(10,2))/2), blockPropStructArray);
        else
            blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 1, '', 0, 0);
            mi_addblocklabel2(Rso+((Rs2-Rso)/2),rzCoords(10,2)/2, blockPropStructArray);
            % Steel Label
            blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 1, '', 0, 0);
            mi_addblocklabel2(rzCoords(10,1)+((rzCoords(12,1)-rzCoords(10,1))/2),rzCoords(10,2)+((rzCoords(12,2)-rzCoords(10,2))/2), blockPropStructArray);
        end
        
    end
  
    
    mi_clearselected;
    mi_selectgroup(1);
    mi_mirror(0, Wp/2, 4*Rm, Wp/2);
    
    % Now add missing segments and nodes
    rzCoords(15,:) = [Rs1, Wp/2];
    mi_addnode(rzCoords(15,1), rzCoords(15,2));
    
    segPropStructArray = mi_segpropstruct(8, 1);
    segPropStructArray(6) = mi_segpropstruct(1, 1, 1, 0, 'Pros A');
    
    mi_addsegment2([rzCoords(1,1),...
                    rzCoords(13,1),...
                    rzCoords(11,1),...
                    rzCoords(14,1),...
                    rzCoords(15,1),...
                    rzCoords(8,1),...
                    rzCoords(7,1),...
                    rzCoords(12,1)],...
                   [rzCoords(1,2),...
                    rzCoords(13,2),...
                    rzCoords(11,2),...
                    rzCoords(14,2),...
                    rzCoords(15,2),...
                    rzCoords(8,2),...
                    rzCoords(7,2),...
                    rzCoords(12,2)],...
                    [rzCoords(1,1),...
                    rzCoords(13,1),...
                    rzCoords(11,1),...
                    rzCoords(15,1),...
                    rzCoords(14,1),...
                    rzCoords(8,1),...
                    rzCoords(7,1),...
                    rzCoords(12,1)],...
                   [Wp-rzCoords(1,2),...
                    Wp-rzCoords(13,2),...
                    Wp-rzCoords(11,2),...
                    rzCoords(15,2),...
                    Wp-rzCoords(14,2),...
                    Wp-rzCoords(8,2),...
                    Wp-rzCoords(7,2),...
                    Wp-rzCoords(12,2)], segPropStructArray);
    
    rinner = rzCoords(13,1)+((rzCoords(15,1)-rzCoords(13,1))/2); 
    router = rzCoords(15,1)+((rzCoords(12,1)-rzCoords(15,1))/2);
    
    if mode(1) == 0 || mode(1) == 1
        if mode(2) == 1
            blockPropStructArray = mi_blockpropstruct(1, magnetname, 0, meshSz, 1, '', 0, 90);
            mi_addblocklabel2(router,Wp/2, blockPropStructArray);
            mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);
        else
            blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 1, '', 0, 0);
            mi_addblocklabel2(router,Wp/2, blockPropStructArray);
            mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);
        end
    elseif mode(1) == 2
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, meshSz, 1, '', 0, 0);
        mi_addblocklabel2(router,Wp/2, blockPropStructArray);  
        mi_addblocklabel2(rinner,Wp/2, blockPropStructArray); 
    elseif mode(1) == 3
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, meshSz, 1, '', 0, 0);
        mi_addblocklabel2(router,Wp/2, blockPropStructArray); 
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 1);
        mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);    
    end
    
    % Now we draw the rest of the coils, exactly how this is done will
    % depend on WcVWp.
    
    if WcVWp > 0.99 * 1/3
        % We will draw the coils with no space between them
        rzCCoords(1,:) = [Ri, Wp/3];
        rzCCoords(2,:) = [Ro, Wp/3];
        
        mi_addnode2(rzCCoords);

        segPropStructArray = mi_segpropstruct(3, 2);

        mi_addsegment2([rzCoords(5,1),...
                        rzCoords(6,1),...
                        rzCCoords(1,1)],...
                        [rzCoords(5,2),...
                        rzCoords(6,2),...
                        rzCCoords(1,2)],...
                        [rzCCoords(1,1),...
                        rzCCoords(2,1),...
                        rzCCoords(2,1)],...
                        [rzCCoords(1,2),...
                        rzCCoords(2,2),...
                        rzCCoords(2,2)], segPropStructArray);
        
        % Mirror the coil
        mi_clearselected;
        mi_selectgroup(2);
        mi_mirror(0, Wp/2, 4*Rm, Wp/2);
        
        % Join Coil Parts
        segPropStructArray = mi_segpropstruct(2, 0);
        
        mi_addsegment2([rzCCoords(1,1),...
                        rzCCoords(2,1)],...
                       [rzCCoords(1,2),...
                        rzCCoords(2,2)],...
                       [rzCCoords(1,1),...
                        rzCCoords(2,1)],...
                       [Wp-rzCCoords(1,2),...
                        Wp-rzCCoords(2,2)], segPropStructArray);
        
    else
        % We will draw the coils and air space (could be iron to make slotted machine) 
        rzCCoords(1,:) = [Ri, ((Wp/3)-Wc)/2];
        rzCCoords(2,:) = [Ro, ((Wp/3)-Wc)/2];
        rzCCoords(3,:) = [Ri, Wc+((Wp/3)-Wc)/2];
        rzCCoords(4,:) = [Ro, Wc+((Wp/3)-Wc)/2];
        rzCCoords(5,:) = [Ri, Wc+3*((Wp/3)-Wc)/2];
        rzCCoords(6,:) = [Ro, Wc+3*((Wp/3)-Wc)/2];
        
        mi_addnode2(rzCCoords);

        segPropStructArray = mi_segpropstruct(9, 2);

        mi_addsegment2([rzCoords(5,1),...
                        rzCoords(6,1),...
                        rzCCoords(1,1),...
                        rzCCoords(1,1),...
                        rzCCoords(2,1),...
                        rzCCoords(3,1),...
                        rzCCoords(3,1),...
                        rzCCoords(4,1),...
                        rzCCoords(5,1)],...
                        [rzCoords(5,2),...
                        rzCoords(6,2),...
                        rzCCoords(1,2),...
                        rzCCoords(1,2),...
                        rzCCoords(2,2),...
                        rzCCoords(3,2),...
                        rzCCoords(3,2),...
                        rzCCoords(4,2),...
                        rzCCoords(5,2)],...
                        [rzCCoords(1,1),...
                        rzCCoords(2,1),...
                        rzCCoords(2,1),...
                        rzCCoords(3,1),...
                        rzCCoords(4,1),...
                        rzCCoords(4,1),...
                        rzCCoords(5,1),...
                        rzCCoords(6,1),...
                        rzCCoords(6,1)],...
                        [rzCCoords(1,2),...
                        rzCCoords(2,2),...
                        rzCCoords(2,2),...
                        rzCCoords(3,2),...
                        rzCCoords(4,2),...
                        rzCCoords(4,2),...
                        rzCCoords(5,2),...
                        rzCCoords(6,2),...
                        rzCCoords(6,2)], segPropStructArray);
        
        % Mirror the coil
        mi_clearselected;
        mi_selectgroup(2);
        mi_mirror(0, Wp/2, 4*Rm, Wp/2);
        
        % Join Coil Parts
        segPropStructArray = mi_segpropstruct(2, 0);
        
        mi_addsegment2([rzCCoords(5,1),...
                        rzCCoords(6,1)],...
                       [rzCCoords(5,2),...
                        rzCCoords(6,2)],...
                       [rzCCoords(5,1),...
                        rzCCoords(6,1)],...
                       [Wp-rzCCoords(5,2),...
                        Wp-rzCCoords(6,2)], segPropStructArray);
                    
        % Air labels between coils

        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 0);
        mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp/3, blockPropStructArray);
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 0);
        mi_addblocklabel2(Ri+(Ro-Ri)/2,2*Wp/3, blockPropStructArray);
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 0);
        mi_addblocklabel2(Ri+(Ro-Ri)/2,((Wp/3)-Wc)/4, blockPropStructArray);
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, meshSz, 0);
        mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp-((Wp/3)-Wc)/4, blockPropStructArray);
        
    end
    
    % Coil labels
    
    if J(1) == 0 || mode(3) == 1
        mi_addcircprop('Coil C', J(1)*pi*(Dc/2).^2, 1);
        blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 30, 'Coil C', Ntot);
    else
        % The current density will be in A/m^2, femm accepts a value in
        % MA/m^2, so we must adjust to account for this
        J(1) = J(1) ./ 1e6;
        % Add a solid copper region with the specified current density
        mi_addmaterial('copper J C', 1, 1, 0, J(1), 58, 0, 0, 1, 0, 0, 0, 0, 0);
        blockPropStructArray = mi_blockpropstruct(1, 'copper J C', 0, meshSz, 30);
    end

    mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp/6, blockPropStructArray);

    if J(2) == 0 || mode(3) == 1
        mi_addcircprop('Coil B', J(2)*pi*(Dc/2).^2, 1);
        blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 20, 'Coil B', Ntot);
    else
        % The current density will be in A/m^2, femm accepts a value in
        % MA/m^2, so we must adjust to account for this
        J(2) = J(2) ./ 1e6;
        % Add a solid copper region with the specified current density
        mi_addmaterial('copper J B', 1, 1, 0, J(2), 58, 0, 0, 1, 0, 0, 0, 0, 0);
        blockPropStructArray = mi_blockpropstruct(1, 'copper J B', 0, meshSz, 20);
    end

    mi_addblocklabel2(Ri+(Ro-Ri)/2,0.5*Wp, blockPropStructArray);

    if J(3) == 0 || mode(3) == 1
        mi_addcircprop('Coil A', J(3)*pi*(Dc/2).^2, 1);
        blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 10, 'Coil A', Ntot);
    else
        % The current density will be in A/m^2, femm accepts a value in
        % MA/m^2, so we must adjust to account for this
        J(3) = J(3) ./ 1e6;
        % Add a solid copper region with the specified current density
        mi_addmaterial('copper J A', 1, 1, 0, J(3), 58, 0, 0, 1, 0, 0, 0, 0, 0);
        blockPropStructArray = mi_blockpropstruct(1, 'copper J A', 0, meshSz, 10);
    end

    mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp-(Wp/6), blockPropStructArray);

%     blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 1, 'Coil C', Ntot);
%     mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp/6, blockPropStructArray);
%     blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 1, 'Coil B', Ntot);
%     mi_addblocklabel2(Ri+(Ro-Ri)/2,0.5*Wp, blockPropStructArray);
%     blockPropStructArray = mi_blockpropstruct(1, wirename, 0, meshSz, 1, 'Coil A', Ntot);
%     mi_addblocklabel2(Ri+(Ro-Ri)/2,Wp-(Wp/6), blockPropStructArray);
        
    % air label for shaft 
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, 2*meshSz, 1);
    mi_addblocklabel2(rzCoords(2,1)/2,Wp/2, blockPropStructArray);
    
    % air gap label
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, 0.8*meshSz, 1);
    mi_addblocklabel2(Rm+(Ri-Rm)/2,Wp/2, blockPropStructArray);
    
    % Air region 1 label
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, 7.5*meshSz, 1);
    mi_addblocklabel2(Ro+((Ra-Ro)/2),Wp/2, blockPropStructArray);
    
    % outer air region label
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 1, 20*meshSz, 1);
    mi_addblocklabel2(Ra+(Router-Ra)/2,Wp/2, blockPropStructArray);
    
    filename = [tempname, int2str(round(rand()*100000)), '_ACTM.fem'];
    mi_saveas(filename);

	
end