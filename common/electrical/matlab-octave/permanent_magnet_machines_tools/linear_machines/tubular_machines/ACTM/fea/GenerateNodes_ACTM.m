function rzCoords = GenerateNodes_ACTM(WmVWp, WpVRm, RsVRm, RshoVRm, Rm)
%This function generates two tables containing the x and y coordinates
%respectively of the necessary nodes for the simulation. This does not include
%the property labels.
%Note that we are redefining the description of Ro in relation to Rm as this is 
%currently extremly unintuitive, it is now defined through the ratio Ro/Rm which
%is much clearer. The description of the armature thickness is through the ratio
%Ra/Ro where Ra is the outer radius of the armature.	
	Wp = Rm * WpVRm;
	%Wm = Wp * WmVWp;
	Ws = Wp * (1 - WmVWp);
	Rs1 = Rm * RsVRm;
    Rsho = RshoVRm * Rm;
    %Rs = Rso + ((Rm-Rso)/2); % steel half way point
    
    % Bottom left node, corner of sim
	rzCoords(1,:) = [0,0];
    % Bottom node of outer shaft radius
	rzCoords(2,:) = [Rsho, 0];
    % Bottom node of inner radius of steel
	rzCoords(3,:) = [Rs1, 0];
    % Bottom node of outer steel radius
	rzCoords(4,:) = [Rm, 0];
    % Bottom right node of sim marking end of air
	rzCoords(5,:) = [3* Rm, 0];
    % Node marking air-steel junction at the shaft
	rzCoords(6,:) = [Rsho, Ws1];
    % Node marking air-steel junction closer to Rm
	rzCoords(7,:) = [Rs2, Ws2];
    % Node marking start of mag close to shaft
	rzCoords(8,:) = [Rsho, Ws/2];
    % Node marking edge of mag at Rm
	rzCoords(9,:) = [Rm, Ws/2];
    
	for n = 1:size(rzCoords,2)
	
		mi_addnode(rzCoords(n,1), rzCoords(n,2));
		
	end
	
	
end