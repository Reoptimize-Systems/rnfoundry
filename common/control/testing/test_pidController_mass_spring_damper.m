% test_pidContoller_mass_spring_damper
%
% the following tests are based on the example provided at:
%
% http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID#9
%
%


m = 1; % kg
b = 10; % N s/m
k = 20; % N/m
F = 1; % N
setx = 0.05;
step = 1;

tstart = 0;
tend = 2;
maxstep = 0.01;
% t = 0:0.01:2;

initial_position = 0; 
initial_speed    = 0; 
 
x0 = [initial_position  initial_speed]; 

tspan = [ tstart, tend ];

%% open loop response

odeopts = odeset ('MaxStep', maxstep);

sol = ode45 ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, F, setx), ...
              tspan, ...
              x0, ...
              odeopts );

figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [setx, setx], 'k:');
hold off

title ('No control (ode45)');
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*setx]);

%% P control


Kp = 300;

C = pidController (Kp, 0, 0, 'InitialTime', tstart);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode45 ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('P control (ode45), Kp: %g', Kp));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PD control

Kp = 300;
Kd = 10;

C = pidController (Kp, 0, Kd, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode45 ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PD control (ode45), Kp: %g, Kd: %g', Kp, Kd));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PI control

Kp = 30;
Ki = 70;

C = pidController (Kp, Ki, 0, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode45 ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PI control (ode45), Kp: %g, Ki: %g', Kp, Ki));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PID control


Kp = 350;
Ki = 300;
Kd = 50;

C = pidController (Kp, Ki, Kd, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode45 ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PID control (ode45), Kp: %g, Ki: %g, Kd: %g', Kp, Ki, Kd));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% open loop response

odeopts = odeset ('MaxStep', maxstep);

sol = ode.rkfixed ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, F, setx), ...
              tspan, ...
              x0, ...
              odeopts );

figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [setx, setx], 'k:');
hold off

title ('No control (ode.rkfixed)');
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*setx]);

%% P control


Kp = 300;

C = pidController (Kp, 0, 0, 'InitialTime', tstart);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode.rkfixed ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('P control (ode.rkfixed), Kp: %g', Kp));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PD control

Kp = 300;
Kd = 10;

C = pidController (Kp, 0, Kd, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode.rkfixed ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PD control (ode.rkfixed), Kp: %g, Kd: %g', Kp, Kd));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PI control

Kp = 30;
Ki = 70;

C = pidController (Kp, Ki, 0, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode.rkfixed ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PI control (ode.rkfixed), Kp: %g, Ki: %g', Kp, Ki));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);

%% PID control

Kp = 350;
Ki = 300;
Kd = 50;

C = pidController (Kp, Ki, Kd, 'InitialTime', tstart, 'MaxOut', 1000, 'MinOut', -1000);

odeopts = odeset ( 'MaxStep', maxstep, ...
                   'OutputFcn', @(t,x,flag) ode.mass_spring_damper_outputFcn (t, x, flag, m, b, k, C, step ) );

sol = ode.rkfixed ( @(t,x) ode.mass_spring_damper (t, x, m, b, k, C, step), ...
              tspan, ...
              x0, ...
              odeopts );
          
figure;
plot (sol.x, sol.y(1,:));
hold on 
plot ([0,tend], [step, step], 'k:');
hold off

title ( sprintf ('PID control (ode.rkfixed), Kp: %g, Ki: %g, Kd: %g', Kp, Ki, Kd));
xlabel ('Time [s]');
ylabel ('x');
set (gca, 'YLim', [0, 2*step]);