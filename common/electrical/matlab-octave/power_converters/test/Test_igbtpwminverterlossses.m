% Test_igbtpwminverterlossses.m


I_CN = [ 15; 15; 15; 15; 75 ];
V_CEN = [ 2.5; 2.5; 2.5; 2.5; 2.5 ];
V_CE0 = 
V_FN = [ 1.8; 1.8; 1.8; 1.8; 2.2 ];
V_F0 = [ 0.7; 0.7; 0.7; 0.7; 0.7 ];
t_rN = [ 200; 200; 200; 200; 200 ] .* 1e-9;
t_fN = [ 200; 200; 200; 200; 300 ] .* 1e-9;
t_rrN = [ 200; 200; 200; 200; 200 ] .* 1e-9;
Q_rrN = [ 200; 200; 200; 200; 1100 ] .* 1e-9;

I_CM,
V_cc = [ 580; 580; 540; 580; 580 ];
F_s = [ 6000; 5700; 10800; 5400; 5700 ];
M = [ 100; 90; 90; 90; 95 ] ./ 100;

[ Ptot, Pi, Pd, Pon, Poff, Prr ]= igbtpwminverterlossses(I_CN,V_CEN,V_CE0,V_FN,V_F0,t_rN,t_fN,t_rrN,Q_rrN,I_CM,V_cc,F_s,M)

