function exportBodyHydroDataToCSV (body, waves, filename)

    T = body.hydroData.simulation_parameters.T';
    
    % get the real valued excitation forces
%     F_ex_real = squeeze (permute (body.hydroForce.fExt.re(1:3,:,:), [3, 1, 2]));
    
    cellstr2txtfile (filename, {'Wave Period (s), Wave Frequency (rad/s), F_ex_real_x (N), F_ex_real_y (N), F_ex_real_z (N), Tau_ex_real_x (N), Tau_ex_real_y (N), Tau_ex_real_z (N)'});
    
    dlmwrite (filename, [1./(waves.w(:)/(tau)), waves.w(:), body.hydroForce.fExt.re], '-append');

end