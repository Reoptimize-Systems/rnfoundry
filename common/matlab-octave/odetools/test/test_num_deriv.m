function [dy,flux_linkage]  = test_num_deriv (t, y, myodederiv)

    flux_linkage = sin (2 * pi * t);
    
    dy = myodederiv.derivative (t, flux_linkage);

end