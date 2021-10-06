function [filename] = RunFEMMSimNew_ACTM(WmVWp, WpVRm, RsoVRm, Rm, mode, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% RunFEMMSimNew_ACTM: A function for generating simulations of a single
% pole of the Air-Cored Tubular Permanent Magnet machine in the FEMM finite
% element analysis program. 
%
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   mode - Scalar or (1 x 2) vector specifying what simulation type is to
%          be performed. 
%
%          mode(1) have the values 0, 1, 2 or 3.
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
%          mode(2) determines the field direction of the magnets, if
%          length(mode) == 1, or mode is omitted the magnet direction is
%          set to 90 degrees, otherwise -90 degrees.
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
    
    if nargin < 6
        
        Rs2VHmag = 0.5;
        Rs1VHmag = 0.5;
        Ws2VhalfWs = 0.25;
        Ws1VhalfWs = 0.25;
        if nargin == 4
            mode = [0 0];
        end
        
        mi_probdef(0,'meters','axi',1e-8, 0, 30, 0);
        
    elseif (nargin > 5 && nargin < 9) || nargin < 4
        
        error('Invalid number of input arguments')
    
    else
        % use smaller angle constraint
        mi_probdef(0,'meters','axi',1e-8, 0, 5, 0);    
    end
    
    if length(mode) == 1
        mode = [mode 0];
    end
    
    if mode(1) == 0 || mode(1) == 1
        Wp = Rm * WpVRm;
        Ws = Wp * (1 - WmVWp);
        Wm = Wp * WmVWp;
        Ws1 = Ws1VhalfWs * (Ws/2);
        Ws2 = Ws2VhalfWs * (Ws/2);
        Wm1 = (1-Ws1VhalfWs) * (Wm/2) + (Ws/2);
        Wm2 = (1-Ws2VhalfWs) * (Wm/2) + (Ws/2);
        Rsho = RsoVRm * Rm;
        Rs1 = Rsho + Rs1VHmag * (Rm-Rsho);
        Rs2 = Rsho + Rs2VHmag * (Rm-Rsho);
    else
        Wp = Rm * WpVRm;
        Ws = Wp * (WmVWp); % is actually Wm
        Wm = Wp * (1-WmVWp); % Is actually Ws
        Ws1 = Ws1VhalfWs * (Ws/2);
        Ws2 = Ws2VhalfWs * (Ws/2);
        Wm1 = (1-Ws1VhalfWs) * (Wm/2) + (Ws/2);
        Wm2 = (1-Ws2VhalfWs) * (Wm/2) + (Ws/2);
        Rsho = RsoVRm * Rm;
        Rs1 = Rsho + Rs1VHmag * (Rm-Rsho);
        Rs2 = Rsho + Rs2VHmag * (Rm-Rsho);
    end
    
    %% Problem Setup
    % Setup some problem variables
%     mi_setUpProblem([], 'meters');
    % Add some boundary types
    mi_addboundprop('Pros A', 0, 0, 0, 90, 0, 0, 0, 0, 0)
    
    for i = 1:5
        mi_addAntiPeriodicBoundary(['antiPeriodic ' num2str(i)]);
    end
    
    % Add some materials
    % First Steel
    steelname = mi_addsteel;
    % Next magnets
    magname = mi_addmagnet( mgoe2hc(40) );
    % Finally Air
    mi_addair;
    
    % Add wire
    mi_addwire;
    
    % Now draw some nodes
    
    % Bottom left node, corner of sim
	rzCoords(1,:) = [0,0];
    % Bottom node of outer shaft radius
	rzCoords(2,:) = [Rsho, 0];
    % Bottom node of inner radius of steel
	rzCoords(3,:) = [Rs1, 0];
    % Bottom node of outer steel radius
	rzCoords(4,:) = [Rm, 0];
    % Bottom right of sim marking 3 Rm into the air gap, we will use a
    % denser mesh in this region
    rzCoords(5,:) = [3 * Rm, 0]; % 3
    % Bottom right node of sim marking end of air
	rzCoords(6,:) = [6 * Rm, 0]; % 6
    % Node marking air-steel junction at the shaft
	rzCoords(7,:) = [Rsho, Ws1];
    % Node marking air-steel junction closer to Rm
	rzCoords(8,:) = [Rs2, Ws2];
    % Node marking start of mag close to shaft
	rzCoords(9,:) = [Rsho, Ws/2];
    % Node marking edge of mag at Rm
	rzCoords(10,:) = [Rm, Ws/2];
    
    rzCoords(11,:) = [Rsho, Wm1];
    rzCoords(12,:) = [Rs2, Wm2];
    
    for n = 1:size(rzCoords,1)

        mi_addnode(rzCoords(n,1), rzCoords(n,2));

    end
    
    % Now add some segments, those that can be mirrored
    
    segPropStructArray(1) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 1');
    segPropStructArray(2) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 2');
    segPropStructArray(3) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 3');
    segPropStructArray(4) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 4');
    segPropStructArray(5) = mi_segpropstruct(1, 1, 1, 0, 'antiPeriodic 5');
    
    mi_addsegment2([rzCoords(1,1),...
                    rzCoords(2,1),...
                    rzCoords(3,1),...
                    rzCoords(4,1),...
                    rzCoords(5,1)],...
                   [rzCoords(1,2),...
                    rzCoords(2,2),...
                    rzCoords(3,2),...
                    rzCoords(4,2),...
                    rzCoords(5,2)],...
                   [rzCoords(2,1),...
                    rzCoords(3,1),...
                    rzCoords(4,1),...
                    rzCoords(5,1),...
                    rzCoords(6,1)],...
                   [rzCoords(2,2),...
                    rzCoords(3,2),...
                    rzCoords(4,2),...
                    rzCoords(5,2),...
                    rzCoords(6,2)], segPropStructArray);
                
    segPropStructArray = mi_segpropstruct(8, 1);
    
    mi_addsegment2([rzCoords(7,1),...
                    rzCoords(8,1),...
                    rzCoords(9,1),...
                    rzCoords(9,1),...
                    rzCoords(10,1),...
                    rzCoords(7,1),...
                    rzCoords(9,1),...
                    rzCoords(11,1)],...
                   [rzCoords(7,2),...
                    rzCoords(8,2),...
                    rzCoords(9,2),...
                    rzCoords(9,2),...
                    rzCoords(10,2),...
                    rzCoords(7,2),...
                    rzCoords(9,2),...
                    rzCoords(11,2)],...
                    [rzCoords(8,1),...
                    rzCoords(3,1),...
                    rzCoords(7,1),...
                    rzCoords(10,1),...
                    rzCoords(4,1),...
                    rzCoords(2,1),...
                    rzCoords(11,1),...
                    rzCoords(12,1)],...
                   [rzCoords(8,2),...
                    rzCoords(3,2),...
                    rzCoords(7,2),...
                    rzCoords(10,2),...
                    rzCoords(4,2),...
                    rzCoords(2,2),...
                    rzCoords(11,2),...
                    rzCoords(12,2)], segPropStructArray);
                
    % Add labels which can be mirrored
    if mode(2) == 0
        magangle = 90;
    else
        magangle = -90;
    end
        
    mesh = min(2*Rm/50, sqrt((2*Rm)^2 + Wp^2) ./ 250);
    
    if mode(1) == 0
        % translator steel label
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, mesh, 1);
        mi_addblocklabel2(Rsho+((Rs2-Rsho)/2),rzCoords(8,2)/2, blockPropStructArray);
        % Steel Label
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, mesh, 1);
        mi_addblocklabel2(rzCoords(8,1)+((rzCoords(10,1)-rzCoords(8,1))/2),rzCoords(8,2)+((rzCoords(10,2)-rzCoords(8,2))/2), blockPropStructArray);
    elseif mode(1) == 1
        % translator air labels
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, mesh, 1);
        mi_addblocklabel2(Rsho+((Rs2-Rsho)/2),rzCoords(8,2)/2, blockPropStructArray);
        % Steel Label
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, mesh, 1);
        mi_addblocklabel2(rzCoords(8,1)+((rzCoords(10,1)-rzCoords(8,1))/2),rzCoords(8,2)+((rzCoords(10,2)-rzCoords(8,2))/2), blockPropStructArray);
    elseif mode(1) == 2 || mode(1) == 3
        % translator mag labels
        blockPropStructArray = mi_blockpropstruct(1, magname, 0, mesh, 1, '', 0, magangle);
        mi_addblocklabel2(Rsho+((Rs2-Rsho)/2),rzCoords(8,2)/2, blockPropStructArray);
        % Steel Label
        blockPropStructArray = mi_blockpropstruct(1, magname, 0, mesh, 1, '', 0, magangle);
        mi_addblocklabel2(rzCoords(8,1)+((rzCoords(10,1)-rzCoords(8,1))/2),rzCoords(8,2)+((rzCoords(10,2)-rzCoords(8,2))/2), blockPropStructArray);
        
    end
  
    
    mi_clearselected;
    mi_selectgroup(1);
    mi_mirror(0, Wp/2, 4*Rm, Wp/2);
    
    % Now add missing segments and nodes
    rzCoords(13,:) = [Rs1, Wp/2];
    mi_addnode(rzCoords(13,1), rzCoords(13,2));
    
    segPropStructArray = mi_segpropstruct(7, 1);
    segPropStructArray(6) = mi_segpropstruct(1, 1, 1, 0, 'Pros A');
    
    mi_addsegment2([rzCoords(1,1),...
                    rzCoords(11,1),...
                    rzCoords(10,1),...
                    rzCoords(12,1),...
                    rzCoords(13,1),...
                    rzCoords(6,1),...
                    rzCoords(5,1)],...
                   [rzCoords(1,2),...
                    rzCoords(11,2),...
                    rzCoords(10,2),...
                    rzCoords(12,2),...
                    rzCoords(13,2),...
                    rzCoords(6,2),...
                    rzCoords(5,2)],...
                    [rzCoords(1,1),...
                    rzCoords(11,1),...
                    rzCoords(10,1),...
                    rzCoords(13,1),...
                    rzCoords(12,1),...
                    rzCoords(6,1),...
                    rzCoords(5,1)],...
                   [Wp-rzCoords(1,2),...
                    Wp-rzCoords(11,2),...
                    Wp-rzCoords(10,2),...
                    rzCoords(13,2),...
                    Wp-rzCoords(12,2),...
                    Wp-rzCoords(6,2),...
                    Wp-rzCoords(5,2)], segPropStructArray);
    
    router = rzCoords(13,1)+((rzCoords(10,1)-rzCoords(13,1))/2);
    rinner = rzCoords(11,1)+((rzCoords(13,1)-rzCoords(11,1))/2);
    
    if mode(1) == 0 || mode(1) == 1
        % Magnet labels
        blockPropStructArray = mi_blockpropstruct(1, magname, 0, mesh, 1, '', 0, magangle);
        mi_addblocklabel2(router,Wp/2, blockPropStructArray);  
        mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);  
    elseif mode(1) == 2
        % Steel labels
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, mesh, 1, '', 0, 0);
        mi_addblocklabel2(router,Wp/2, blockPropStructArray);  
        mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);
    elseif mode(1) == 3
        % Steel labels
        blockPropStructArray = mi_blockpropstruct(1, steelname, 0, mesh, 1, '', 0, 0);
        mi_addblocklabel2(router,Wp/2, blockPropStructArray); 
        blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, mesh, 1);
        mi_addblocklabel2(rinner,Wp/2, blockPropStructArray);    
    end
    
    % air label for shaft 
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, 2*mesh, 1);
    mi_addblocklabel2(rzCoords(2,1)/2,Wp/2, blockPropStructArray);
    
    % air label above translator
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 0, 7.5*mesh, 1);
    mi_addblocklabel2(2.999*Rm,Wp/2, blockPropStructArray);
    
    % second air label above translator
    blockPropStructArray = mi_blockpropstruct(1, 'Air', 1, 20*mesh, 1);
    mi_addblocklabel2(5.99*Rm,Wp/2, blockPropStructArray);
    
    filename = [tempname, int2str(round(rand()*100000)), '_ACTM.fem'];
    mi_saveas(filename);


end