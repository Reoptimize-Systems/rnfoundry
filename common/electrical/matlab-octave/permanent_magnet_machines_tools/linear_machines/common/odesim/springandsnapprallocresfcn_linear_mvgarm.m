function results = springandsnapprallocresfcn_linear_mvgarm(design, simoptions, T, Y)

    % Now preallocate arrays of the correct sizes
    results.ydot = zeros(ceil(size(Y,1)/simoptions.skip), size(simoptions.IC, 2));
    results.dpsidxR = zeros(ceil(length(T)/simoptions.skip), design.Phases);
    results.EMF = zeros(ceil(length(T)/simoptions.skip), design.Phases);
    results.Fpto = zeros(ceil(length(T)/simoptions.skip), 1);
    results.xT = zeros(ceil(length(T)/simoptions.skip), 1);
    results.vT = zeros(ceil(length(T)/simoptions.skip), 1);
    results.excitation_force_heave = zeros(ceil(length(T)/simoptions.skip), 1);
    results.excitation_force_surge = zeros(ceil(length(T)/simoptions.skip), 1);
    results.radiation_force_heave = zeros(ceil(length(T)/simoptions.skip), 1);
    results.radiation_force_surge = zeros(ceil(length(T)/simoptions.skip), 1);
    results.buoyancy_force = zeros(ceil(length(T)/simoptions.skip), 1);
    results.FBDh = zeros(ceil(length(T)/simoptions.skip), 1);
    results.FBDs = zeros(ceil(length(T)/simoptions.skip), 1);
    results.Ffea_heave = zeros(ceil(length(T)/simoptions.skip), 1);
    results.Ffea_surge = zeros(ceil(length(T)/simoptions.skip), 1);
    results.Faddbuoy = zeros(ceil(length(T)/simoptions.skip), 2);
    results.Fa = zeros(ceil(length(T)/simoptions.skip), 2);

end
