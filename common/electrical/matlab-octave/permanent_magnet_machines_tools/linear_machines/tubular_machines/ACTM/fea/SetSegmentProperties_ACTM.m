function SetSegmentProperties_ACTM(XsegCoords, YsegCoords, Rm)

	elementsize = 0.01 * Rm;
	
	%first check nothing's selected
	mi_clearselected;
	
	%First we will set the antiperiodic boundaries

	n = 3;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	n = 5;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('antiperiodic3', elementsize, 1, 0, 0);

	mi_clearselected;


	n = 7;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	n = 1;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('antiperiodic1', elementsize, 1, 0, 0);

	mi_clearselected;


	n = 2;
	
	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	n = 6;
	
	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('antiperiodic2', elementsize, 1, 0, 0);

	mi_clearselected;



	%Now set the outer segment property
	n = 4;
	
	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('proscribed A', elementsize, 0, 0, 0);

	mi_clearselected;



	%Now set the axis segment properties, solver will assume rotation around the 
	%x = 0 line if no boundary set, so leave this blank
	n = 8;
	
	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;


	n = 9;
	
	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 10;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 11;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 12;

	mi_selectsegment(XsegCoords(n), YsegCoords(n)+0.1);

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 13;

	mi_selectsegment(XsegCoords(n), YsegCoords(n)+0.1);

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 14;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 15;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 16;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 17;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 18;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;

	n = 19;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('', elementsize, 1, 0, 0);

	mi_clearselected;
	
	n = 21;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('antiperiodic4', elementsize, 1, 0, 0);

	mi_clearselected;
	
	n = 22;

	mi_selectsegment(XsegCoords(n), YsegCoords(n));

	mi_setsegmentprop('antiperiodic4', elementsize, 1, 0, 0);

	mi_clearselected;

end