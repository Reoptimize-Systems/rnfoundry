function [dx,y]  = test_state_space_odefcn (t, x, SS, K)

    u = -K * x;
    
    dx = SS.derivatives (u);
    
    y = SS.outputs (u);

end