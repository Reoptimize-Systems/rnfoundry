function exportBodyHydroDataToCSV (body, filename)

    T = body.hydroData.simulation_parameters.T';
    
    % get the real valued excitation forces
    F_ex_real = squeeze (permute (body.hydroData.hydro_coeffs.excitation.re(1:3,:,:), [3, 1, 2]));
    
    cellstr2txtfile (filename, {'Wave Period (s), F_ex_real_x (N), F_ex_real_y (N), F_ex_real_z (N)'});
    
    dlmwrite (filename, [T, F_ex_real], '-append');

end