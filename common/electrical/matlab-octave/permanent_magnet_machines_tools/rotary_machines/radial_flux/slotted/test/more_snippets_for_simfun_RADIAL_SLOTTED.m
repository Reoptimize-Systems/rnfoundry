
tempmaxflpos = maxflpos;

pos = linspace(-1,2,1000);

fl = fluxlinkagefrmintAslm(intAslm, ...
                          coilpitch, ...
                          pos, ...
                          design.CoilTurns, ...
                          design.ls, ...
                          design.CoilArea, ...
                          tempmaxflpos);

plot(pos, fl);
hold on
% x = linspace(0,1,7) + tempmaxflpos - (design.thetas/(2*design.thetap));
x = linspace(0,1,7) + tempmaxflpos + (design.thetas/(4*design.thetap));
x = [x;x];
y = repmat([max(fl); min(fl)], 1, size(x,2));
for i = 1:size(x,2)
    plot(x(:,i), y(:,i), ':r');
end

hold off

%%

plotslm(intAslm(1))
plotslm(intAslm(2))

%%

plot(design.psilookup(1,:), design.psilookup(2,:));