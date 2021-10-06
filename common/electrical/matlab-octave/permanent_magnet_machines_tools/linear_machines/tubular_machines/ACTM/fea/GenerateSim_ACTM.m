function [Rm, Area] = GenerateSim_ACTM(filename, WmVWp, WpVRm, RsVRm, RsoVRm, RmorArea, mode)
%if mode = 0 we generate a standard sim, RmorArea is Area
%if mode = 1 we generate a standard sim, RmorArea is Rm
%if mode = 2 we use RsVRm, RmorArea is Area
%if mode = 3 we use RsVRm, RmorArea is Rm	

	if mode == 0 || mode == 1

		RsVRm = 0.5;

	end

    if mode == 0 || mode == 2

        Area = RmorArea;

        Rm = GetRm(WpVRm, Area);

        %Check Wp and Rm size is sensible and make sim bigger if not
        %this ensures the accuracy of the simulations, as the grid has
        %a resolution of only around 0.001 cm

        Wp = WpVRm * Rm;

        %If Wp is less than 1 cm we make the area bigger until
        %it is at least 1 cm
        while  Wp < 1

            Area = Area * 110/100;

            Rm = GetRm(WpVRm, Area);

            Wp = WpVRm * Rm;

        end

        Ro = RoVRm * Rm;

        %If The stator gap is less than 1 cm we make the area bigger until
        %it is at least 1 cm
        while  (Ro - Rm) < 1

            Area = Area * 110/100;

            Rm = GetRm(WpVRm, Area);

            Ro = RoVRm * Rm;

        end

        Ro = RoVRm * Rm;

        %If The steel spacer is less than 0.3 cm we make the area bigger until
        %it is at least 0.3 cm
        while  (Wp*(1-WmVWp)) < 0.6

            Area = Area * 110/100;

            Rm = GetRm(WpVRm, Area);

            Ro = RoVRm * Rm;

            Wp = WpVRm * Rm;

        end

        %If The magnet width is less than 0.3 cm we make the area bigger until
        %it is at least 0.3 cm
        while  (Wp*WmVWp) < 0.3

            Area = Area * 110/100;

            Rm = GetRm(WpVRm, Area);

            Ro = RoVRm * Rm;

            Wp = WpVRm * Rm;

        end

    else
        Rm = RmorArea;
        Wp = WpVRm * Rm;
    end

    %GapMesh = ChooseMesh(Wp, (Ro-Rm), 1/120)
	AirMesh = ChooseMesh(Wp, (2*Rm), (1/80)); %GapMesh % %ChooseGapMesh(Wp, (Ra-Ro))

    %if dealing with narrow regions, solver seems to be confused by large discontinuites
	%in mesh sizes between regions, we therefore ensure that the mesh sizes are similar when
	%WpVRm is small. We could just always choose the same size mesh, but this can result in
	%a failure due to the size of the problem in some cases. Reducing the mesh by an order 
        %of magnitude when we can reduces the likelyhood of this problem.
	if WpVRm < 0.15 

           TransSteelMesh = AirMesh * 2; %* 5 %ChooseMesh(Wp*(1-WmVWp), Rm, 0.2)
	       MagMesh = TransSteelMesh; %AirMesh %* 5 %ChooseMesh((Wp*WmVWp), Rm, 0.2)
           ArmSteelMesh = AirMesh * 2;

        else

            TransSteelMesh = AirMesh * 5; %ChooseMesh(Wp*(1-WmVWp), Rm, 0.2)
            MagMesh = AirMesh * 5;
	
	end
	
	ShaftMesh = MagMesh * 2; %* 5
	
	%Now proceed with sim
	SetUpProblem(filename);
	
	AddMaterialsAndBoundaries;

	xCoords, yCoords = GenerateNodes_ACTM(WmVWp, WpVRm, RsVRm, RsoVRm, Rm);

	XsegCoords, YsegCoords = GenerateSegments_ACTM(xCoords, yCoords);

	SetSegmentProperties_ACTM(XsegCoords, YsegCoords, Rm);

	AddLabels_ACTM(XsegCoords, YsegCoords, Rm, AirMesh, MagMesh, TransSteelMesh, ShaftMesh, mode);

	mi_zoomnatural;

	mi_saveas(filename);

end
