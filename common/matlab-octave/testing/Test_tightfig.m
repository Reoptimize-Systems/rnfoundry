x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(X)+cos(Y);

if exist ('htestfig', 'var') ...
        && isvalid (htestfig)
    
    close (htestfig);
    
end

htestfig = figure;

[~,hcntr] = contour(X,Y,Z);

xlabel ('x label');
ylabel ('y label');

hax = get (hcntr, 'Parent');

hcbar = colorbar ();

legend ('test', 'Location', 'southoutside')

title ('here is a title');
tightfig ()