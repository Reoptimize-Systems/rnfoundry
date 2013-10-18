% Test_igbtpwminverterlossses.m


I_CN = [ 15; 15; 15; 15; 75 ];
V_CEN = [ 2.5; 2.5; 2.5; 2.5; 2.5 ];
V_CE0 = [ 1; 1; 1; 1; 1 ];
V_FN = [ 1.8; 1.8; 1.8; 1.8; 2.2 ];
V_F0 = [ 0.7; 0.7; 0.7; 0.7; 0.7 ];
t_r = [ 200; 200; 200; 200; 200 ] .* 1e-9;
t_f = [ 200; 200; 200; 200; 300 ] .* 1e-9;
t_rr = [ 200; 200; 200; 200; 200 ] .* 1e-9;
Q_rr = [ 200; 200; 200; 200; 1100 ] .* 1e-9;

% 
I_out = [ 3.9; 5.2; 5; 2.85; 24 ]; 

% masimum collector current
I_CM = I_out * sqrt(2);

V_cc = [ 580; 580; 540; 580; 580 ];
F_s = [ 6000; 5700; 10800; 5400; 5700 ];
M = [ 100; 90; 90; 90; 95 ] ./ 100;
pf = [ 0.8; 0.8; 0.8; 0.8; 0.85 ];

t_rrN = t_rr ./ (0.8 + 0.2 .* (I_CM ./ I_CN));
t_fN = t_f ./ (2/3 + (1/3) .* (I_CM ./ I_CN));
t_rN = t_r ./ (I_CM ./ I_CN);
Q_rrN = Q_rr;

[ Ptot, Pi, Pd, Pon, Poff, Prr ] = ...
    igbtpwminverterlossses(I_CN, V_CEN, V_CE0, V_FN, V_F0, t_rN, t_fN, t_rrN, Q_rrN, I_CM, V_cc, F_s, M, pf)

