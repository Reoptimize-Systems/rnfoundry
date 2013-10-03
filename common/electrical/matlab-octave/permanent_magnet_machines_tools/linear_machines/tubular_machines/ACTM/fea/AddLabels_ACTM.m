function AddLabels_ACTM(XsegCoords, YsegCoords, AirMesh, MagMesh, TransSteelMesh, ShaftMesh, mode)
% If mode is zero we make a standard sim, i.e. only steel and magnet in stator
% otherwise non-magnetic parts are added
	mi_clearselected;

    % lower steel in translator
	n = 17;
    
	mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_setblockprop('1117 Steel', 0, TransSteelMesh, '', 0, 0, 0);

	mi_clearselected;

    % magnet in translator
	n = 11;

	mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_setblockprop('NdFeB 40 MGOe', 0, MagMesh, '', 90, 0, 0);

	mi_clearselected;

    % upper steel in translator
	n = 19;

	mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_setblockprop('1117 Steel', 0, TransSteelMesh, '', 0, 0, 0);

	mi_clearselected;

    % air region above translator
	n = 4;

	mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_setblockprop('Air', 0, AirMesh, '', 0, 0, 0);

	mi_clearselected;
	
	%shaft material
	n = 14;

	mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

	mi_setblockprop('Air', 0, ShaftMesh, '', 0, 0, 0);

	mi_clearselected;

    % inner translator steel/air
	if mode == 0 || mode == 1 
        % inner translator steel
		n = 18;

		mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_setblockprop('1117 Steel', 0, TransSteelMesh, '', 0, 0, 0);

		mi_clearselected;

        % inner translator steel
		n = 16;

		mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_setblockprop('1117 Steel', 0, TransSteelMesh, '', 0, 0, 0);

		mi_clearselected;
	else
        % inner translator air
		n = 18;

		mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_setblockprop('Air', 0, TransSteelMesh, '', 0, 0, 0);

		mi_clearselected;

        % inner translator air
		n = 16;

		mi_addblocklabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_selectlabel(XsegCoords(n)*0.999, YsegCoords(n));

		mi_setblockprop('Air', 0, TransSteelMesh, '', 0, 0, 0);

		mi_clearselected;

	end

end