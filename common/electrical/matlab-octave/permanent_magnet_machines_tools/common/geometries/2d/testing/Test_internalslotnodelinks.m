% Test_internalslotnodelinks

ycoil = 1;
yshoegap = 0.3 * ycoil;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 2;
tol = 1e-5;

[nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);

%%

ycoil = 1;
yshoegap = 0.3 * ycoil;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);

%%

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0.8 * xshoebase;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);

%%

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0.1;
xshoegap = 0;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);

%%

ycoil = 1;
yshoegap = 0;

xcore = 0.2;
xcoil = 1.0;
xshoebase = 0;
xshoegap = 0;
coillayers = 1;
tol = 1e-5;

[nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, coillabelloc] = ...
    internalslotnodelinks(ycoil, yshoegap, xcore, xcoil, xshoebase, xshoegap, coillayers, tol);

plotnodelinks(nodes, links);

