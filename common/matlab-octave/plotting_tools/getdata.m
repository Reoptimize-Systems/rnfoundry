function getdata (option)
global A opt leg label

[fname,pname] = uigetfile('*.xls','Select File');
dbfile= strcat(pname,fname);
if length(dbfile) == 0 return; end

switch option
	case 1
		[A,leg]= readdata(dbfile);
	case 2
		[A,leg]= xlsread(dbfile);
	case 3
		[A,leg]= xlsread(dbfile,-1)
end

%Generar cell string array para titulos si no los hay
if isempty(leg) leg= cellstr( num2str( flipud(rot90([1:size(A,2)])) ) ); end

opt.xc= 1;
switch size(A,2)
	case 1
		opt.yc= 1; 
		opt.ec= 1;
		opt.zc= 1;
	case 2
		opt.yc= 2; 
		opt.ec= 2;
		opt.zc= 2;
	otherwise
		opt.yc= [2:size(A,2)-1];
		opt.ec= [fix(size(A,2)/2):size(A,2)];
		opt.zc= [3:size(A,2)];
end

return
