
t0 = 0;
tspan = 100;
fluxlinkage = sin (2 * pi * t0);

myodederiv = odederiv (t0, fluxlinkage);

odeopts = odeset ('OutputFcn', @(t,y,flag) test_odederiv_outputfcn (t, y, flag, myodederiv), ...
                  'MaxStep', 0.025 ...
                  );

[T, Y] = ode15s (@(t,y) test_num_deriv (t, y, myodederiv), [t0,t0+tspan],  fluxlinkage, odeopts);

myodederiv = odederiv (T(1), sin (T(1)));
flux_linkage_out = [];
dy = [];
for ind = 1:numel(T)

    [dy(ind), flux_linkage_out(ind)] = test_num_deriv (T(ind), Y(ind,:), myodederiv);
    
    myodederiv.update (T(ind), flux_linkage_out(ind))

end

plot (T, Y, T, flux_linkage_out, T, dy);

% figure; plot (linspace (0, 10, 1000), sin (2 * pi * linspace (0, 10, 1000)))