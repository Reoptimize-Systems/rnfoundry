% Test_internalslotnodelinks

%% shoe blunt edge with gap

ycoil = 1;
yshoegap = 0.3 * ycoil;

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

%% Shoe sharp point with gap

ycoil = 1;
yshoegap = 0.3 * ycoil;

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

%% shoe blunt edge joined tips

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, 'SplitX', true);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off

%% Shoe sharp point joined tips

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.5;
xshoegap = 0;
coillayers = 20;
tol = 1e-5;

[nodes, links, cornernodes, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol, 'CoilBaseFraction', 0.7);

plotnodelinks(nodes, links);
hold on
plot ( [ coillabelloc(:,1) ], ...
       [ coillabelloc(:,2) ], ...
       'xk' )
hold off
%%

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0;
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

% shoecurvepoints
shoecontrolfrac = 0.5;
[x, y] = feval (handles{1}, shoecontrolfrac, xcore, xcoil, xshoebase, xshoegap, ycoilbase, ycoilgap, yshoegap);

plot (x, y);
