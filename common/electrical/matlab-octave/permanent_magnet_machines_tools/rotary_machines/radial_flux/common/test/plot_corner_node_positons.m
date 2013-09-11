
tempnodes = getnodecoords_mfemm(FemmProblem);
plot(nodes(cornernodes(1)+1,1),nodes(cornernodes(1)+1,2), 'bo')
hold on
plot(nodes(cornernodes(2)+1,1),nodes(cornernodes(2)+1,2), 'bx')

plot(nodes(cornernodes(3)+1,1),nodes(cornernodes(3)+1,2), 'rx')

plot(nodes(cornernodes(4)+1,1),nodes(cornernodes(4)+1,2), 'ro')

plot(nodes(:,1), nodes(:,2), 'g+')
plot(tempnodes(:,1), tempnodes(:,2), 'k+')
hold off

%%

tempnodes = getnodecoords_mfemm(FemmProblem);

plot(tempnodes(cornernodes(1)+elcount.NNodes+1,1),tempnodes(cornernodes(1)+elcount.NNodes+1,2), 'bo')
hold on
plot(tempnodes(cornernodes(2)+elcount.NNodes+1,1),tempnodes(cornernodes(2)+elcount.NNodes+1,2), 'bx')

plot(tempnodes(cornernodes(3)+elcount.NNodes+1,1),tempnodes(cornernodes(3)+elcount.NNodes+1,2), 'rx')

plot(tempnodes(cornernodes(4)+elcount.NNodes+1,1),tempnodes(cornernodes(4)+elcount.NNodes+1,2), 'ro')

plot(tempnodes(:,1), tempnodes(:,2), 'g+')
hold off

%%
tempnodes = getnodecoords_mfemm(FemmProblem);

plot(tempnodes(bottomnodes(1)+1,1),tempnodes(bottomnodes(1)+1,2), 'bo')
hold on
plot(tempnodes(bottomnodes(2)+1,1),tempnodes(bottomnodes(2)+1,2), 'bx')

plot(tempnodes(:,1), tempnodes(:,2), 'g+')
hold off

%%
tempnodes = getnodecoords_mfemm(FemmProblem);

plot(tempnodes(thisslotcornernodes(1)+1,1),tempnodes(thisslotcornernodes(1)+1,2), 'bo')
hold on
plot(tempnodes(thisslotcornernodes(2)+1,1),tempnodes(thisslotcornernodes(2)+1,2), 'bx')

plot(tempnodes(thisslotcornernodes(3)+1,1),tempnodes(thisslotcornernodes(3)+1,2), 'rx')

plot(tempnodes(thisslotcornernodes(4)+1,1),tempnodes(thisslotcornernodes(4)+1,2), 'ro')

plot(tempnodes(:,1), tempnodes(:,2), 'g+')

hold off

%%
tempnodes = getnodecoords_mfemm(FemmProblem);

plot(tempnodes(lastslotcornernodes(1)+1,1),tempnodes(lastslotcornernodes(1)+1,2), 'bo')
hold on
plot(tempnodes(lastslotcornernodes(2)+1,1),tempnodes(lastslotcornernodes(2)+1,2), 'bx')

plot(tempnodes(lastslotcornernodes(3)+1,1),tempnodes(lastslotcornernodes(3)+1,2), 'rx')

plot(tempnodes(lastslotcornernodes(4)+1,1),tempnodes(lastslotcornernodes(4)+1,2), 'ro')

plot(tempnodes(:,1), tempnodes(:,2), 'g+')

hold off

%%
tempnodes = getnodecoords_mfemm(FemmProblem);

plot(tempnodes(outernodes(1)+1,1),tempnodes(outernodes(1)+1,2), 'bo')
hold on
plot(tempnodes(outernodes(2)+1,1),tempnodes(outernodes(2)+1,2), 'bx')

plot(tempnodes(outernodes(3)+1,1),tempnodes(outernodes(3)+1,2), 'rx')

plot(tempnodes(outernodes(4)+1,1),tempnodes(outernodes(4)+1,2), 'ro')

plot(tempnodes(:,1), tempnodes(:,2), 'g+')
hold off

