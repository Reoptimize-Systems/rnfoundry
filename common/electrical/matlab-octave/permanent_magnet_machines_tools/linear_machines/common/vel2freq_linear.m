function freq = vel2freq_linear (vel, polewidth)
% estimates electrical frequency at given velocity and machine pole width

    freq = vel / (2 * polewidth);

end