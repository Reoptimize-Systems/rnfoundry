
% testing based on:
%
% http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlStateSpace

A = [ 0   1   0
     980  0  -2.8
      0   0  -100 ];

B = [ 0
      0
      100 ];

C = [ 1 0 0; eye(3)];

D = zeros ( size (C, 1), size (B, 2) );

x0 = [0.01 0 0];

SS = stateSpace (A, B, C, D, x0);

odeopts = odeset ('OutputFcn', @SS.outputfcn ...
                  ... , 'MaxStep', 0.0001 ...
                  , 'RelTol', 1e-6 ...
                  ... , 'AbsTol', [1e-5, 1e-5, 1e-5] ...
                  );

t0 = 0;
tspan = 2;

if exist ('ss', 'class')
    % compare to control systems state-space

    tss = t0:0.01:tspan;
    u = zeros(size(t));

    p1 = -10 + 10i;
    p2 = -10 - 10i;
    p3 = -50;

    K = place(A,B,[p1 p2 p3]);

    sys_cl = ss(A-B*K,B,C(1,:),0);

    [yss,tss,x] = lsim(sys_cl,u,tss,x0);
    
    figure;
    
    title ('method comparison');
    xlabel('Time (sec)')
    ylabel('Ball Position (m)')
    
    plot ( [t0, t0+tspan], [0, 0], 'Color', 'k', 'LineWidth', 0.1);
    hold on
    plot (tss, yss);
    hold off
    
else
    K = [-280.7143e+000 , -7.7857e+000,  -300.0000e-003];
end


[T, X] = ode45 (@(t,x) test_state_space_odefcn (t, x, SS, K), [t0,t0+tspan],  x0, odeopts);

% recalculate outputs
SS.reset ()

dx = [];
y = [];

for ind = 1:numel (T)
    
    [dx(:,ind),y(:,ind)]  = test_state_space_odefcn (T(ind), X(ind,:)', SS, K);
    
    SS.update (X(ind,:)');
    
end

%% plot the results

hold on
plot (T, y(1,:))
hold off
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')
% title ('ode stateSpace class method')

%% compare to simulink 

% load_system ('state_space_simulink_compare');

sim ('state_space_simulink_compare')

% close_system ('state_space_simulink_compare');

%% plot simulink results

% figure;
% plot ( [0, 2], [0, 0], 'Color', 'k', 'LineWidth', 0.1);
hold on
plot (simulink_state_space_y.time, simulink_state_space_y.signals.values);
hold off
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')

%%

legend ('ss lsim method', 'ode stateSpace classdef method', 'Simulink State Space Block Method');





              


