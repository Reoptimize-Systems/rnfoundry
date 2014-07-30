% Test_internalslotnodelinks

%% shoe blunt edge with gap

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% shoe blunt edge with gap with insulation

ycoil = [0.5, 1.2];
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'InsulationThickness', 0.02);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off

%% shoe blunt edge with gap

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% shoe blunt edge with gap with insulation

ycoil = 0.5 + (0.5) * rand (1, 2);
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'InsulationThickness', 0.02);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off

%% shoe blunt edge with gap with insulation with high shoe control frac

ycoil = [ 0.2, 1.5 ];
yshoegap = 0.9 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.5 * xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks (ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                           'InsulationThickness', 0.02, ...
                           'ShoeCurveControlFrac', 0.95 );

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off

%% Shoe sharp point with gap

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% Shoe sharp point with gap with insulation

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0.3 * ycoil(1);

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'InsulationThickness', 0.02);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off


%% shoe blunt edge joined tips

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 1;
tol = 1e-5;
coillayers = 1;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
        'SplitX', true, 'CoilBaseFraction', 0.7);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% shoe blunt edge joined tips with insulation

ycoil =  0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'InsulationThickness', 0.02);
                          
plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' );
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off


%% Shoe sharp point joined tips

ycoil = 0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.5;
xshoegap = 0;
coillayers = 4;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, 'CoilBaseFraction', 0.7);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% Shoe sharp point joined tips with insulation

ycoil = 0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.5;
xshoegap = 0;
coillayers = 4;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'CoilBaseFraction', 0, 'InsulationThickness', 0.02);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' );
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off

%% No shoe with gap

ycoil = 0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0;
xshoegap = 0;
coillayers = 5;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% No shoe with gap with insulation

ycoil = 0.5 + (0.5) * rand (1, 2);
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0;
xshoegap = 0;
coillayers = 5;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc, vertlinkinds, toothlinkinds, inslabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, ...
                          'InsulationThickness', 0.02);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' );
plot ( inslabelloc(1), inslabelloc(2), 'xm');
hold off


%% Subfunctions

handles = internalslotnodelinks ();

ycoil = 1;
ycoilgap = ycoil;
ycoilbase = ycoil*5;
yshoegap = 0.95 * ycoil;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
xcoilbase = 0.1 * xcoil;
xcoilbody = xcoil - xcoilbase;

% shoecurvepoints
shoecontrolfrac = 0.5;
[x, y, Qx, Qy, Px, Py] = feval (handles{1}, shoecontrolfrac, xcore, xcoil, xshoebase, xshoegap, ycoilbase, ycoilgap, yshoegap);

plot (x, y, 'LineWidth', 2);

% inscurvepoints
tins = 0.01;
[insx, insy, insQx, insQy] = feval (handles{3}, Qx, Qy, Px, Py, tins);

hold all

plot (insx, insy, 'LineWidth', 2);

hold off

[x, y, Qx, Qy, m, c, Px, Py] = feval (handles{2}, xcore, xcoilbase, xcoilbody, ycoilbase, ycoilgap);
[insx, insy, insQx, insQy] = feval (handles{3}, Qx, Qy, Px, Py, tins);
plot (x, y, 'LineWidth', 2);
hold all
plot (insx, insy, 'LineWidth', 2);
hold off
axis equal