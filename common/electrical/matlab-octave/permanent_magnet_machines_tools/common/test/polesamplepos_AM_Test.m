
Wp = 0.25;

t = linspace(0, 4, 300);

x = sin(2*pi*(1/3)*t);

[pos,tpos] = polesamplepos_AM(x, t, Wp, 12);

plot(tpos, pos, 'DisplayName', 'pppos vs tttpos', 'XDataSource', 'tttpos', 'YDataSource', 'pppos'); figure(gcf)
