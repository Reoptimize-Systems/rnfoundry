function [dx,y]  = test_state_space_odefcn (t, x, SS, K)

    u = -K * x; % 0.001;
    
    dx = SS.derivatives (u);
    
    y = SS.outputs (u);

end