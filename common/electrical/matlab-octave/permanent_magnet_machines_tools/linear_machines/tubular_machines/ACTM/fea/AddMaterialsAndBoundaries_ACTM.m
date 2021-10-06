function AddMaterialsAndBoundaries_ACTM()

	AddAir;
	
	AddMagnet;
	
	AddSteel;
	
	AddWire;
	
	mi_addAntiPeriodicBoundary('antiperiodic1');
	
	mi_addAntiPeriodicBoundary('antiperiodic2');
	
	mi_addAntiPeriodicBoundary('antiperiodic3');
	
	mi_addAntiPeriodicBoundary('antiperiodic4');
	
	mi_addAntiPeriodicBoundary('proscribed A');
	
end