function [XsegCoords, YsegCoords] = GenerateSegments_ACTM(xCoords, yCoords)

    YsegCoords = zeros(1,22);
    XsegCoords = YsegCoords;
    
	%First we generate segments 1 to 13
    for n = 1:13

        mi_addsegment(xCoords(n), yCoords(n), xCoords(n+1), yCoords(n+1));

        if xCoords(n) < xCoords(n+1)

            XsegCoords(n) = xCoords(n) + ((xCoords(n+1) - xCoords(n)) / 2);

        else

            XsegCoords(n) = xCoords(n+1) + ((xCoords(n) - xCoords(n+1)) / 2);

        end


        if yCoords(n) < yCoords(n+1)

            YsegCoords(n) = yCoords(n) + ((yCoords(n+1) - yCoords(n)) / 2);
            %print(YsegCoords(n))
        else

            YsegCoords(n) = yCoords(n+1) + ((yCoords(n) - yCoords(n+1)) / 2);
            %print(YsegCoords(n))
        end
    end

	% Add segments 14 to 19 manually
    
	mi_addsegment(xCoords(14), yCoords(14), xCoords(9), yCoords(9));
		
	XsegCoords(14) = xCoords(14);
		
	YsegCoords(14) = yCoords(14) + ((yCoords(14) - yCoords(1)) / 2);

	 
	mi_addsegment(xCoords(1), yCoords(1), xCoords(14), yCoords(14));
		
	XsegCoords(15) = xCoords(1);

	YsegCoords(15) = yCoords(1) + ((yCoords(14) - yCoords(1)) / 2);


	mi_addsegment(xCoords(2), yCoords(2), xCoords(13), yCoords(13));
		
	XsegCoords(16) = xCoords(2);
		
	YsegCoords(16) = yCoords(2) + ((yCoords(13) - yCoords(2)) / 2);


	mi_addsegment(xCoords(3), yCoords(3), xCoords(12), yCoords(12));
		
	XsegCoords(17) = xCoords(3);
		
	YsegCoords(17) = yCoords(3) + ((yCoords(12) - yCoords(3)) / 2);


	mi_addsegment(xCoords(10), yCoords(10), xCoords(7), yCoords(7));
		
	XsegCoords(18) = xCoords(10);
		
	YsegCoords(18) = yCoords(10) + ((yCoords(7) - yCoords(10)) / 2);


	mi_addsegment(xCoords(11), yCoords(11), xCoords(6), yCoords(6));
		
	XsegCoords(19) = xCoords(11);
		
	YsegCoords(19) = yCoords(11) + ((yCoords(6) - yCoords(11)) / 2);

	%back of sim
    mi_addsegment(xCoords(15), yCoords(15), xCoords(16), yCoords(16));
		
	XsegCoords(20) = xCoords(15);

	YsegCoords(20) = yCoords(16) / 2;
	
	%top of shaft
    mi_addsegment(xCoords(16), yCoords(16), xCoords(8), yCoords(8));
		
	XsegCoords(21) = xCoords(8) / 2;

	YsegCoords(21) = yCoords(16);
	
	%bottom of shaft
    mi_addsegment(xCoords(15), yCoords(15), xCoords(1), yCoords(1));
		
	XsegCoords(22) = xCoords(1) / 2;

	YsegCoords(22) = yCoords(15);

	
	mi_clearselected;

end
