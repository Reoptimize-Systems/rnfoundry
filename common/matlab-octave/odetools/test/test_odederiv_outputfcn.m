function status = test_odederiv_outputfcn (t, y, flag, myodederiv)

    status = 0;
    
    if isempty (flag)
        
%         [~,flux_linkage]  = test_num_deriv (t(end), y(end), myodederiv);
        
        myodederiv.update (t(end), y(end));
        
    elseif strcmp (flag, 'init')
        
    elseif strcmp (flag, 'done')
        
    end
    
end